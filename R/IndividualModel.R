# ================================ INDIVIDUAL MODEL ============================== #
###: Get variability levels
###: 
###: Get a summary of the variability levels (inter-individual and/or intra-individual variability) present in the current project.
###: @return A collection of the variability levels present in the currently loaded project.
###: @examples
###: \dontrun{
###: getVariabilityLevels()
###: }
###: @export
getVariabilityLevels = function(){
  arguments = list()
  output = .processRequest("monolix", "getavailablelevels", arguments, "asynchronous")
  return(output)
}


###: Get individual parameter model
###: 
###: Get a summary of the information concerning the individual parameter model. The available informations are:
###: \itemize{
###: \item name: (\emph{string}) name of the individual parameter
###: \item distribution: (\emph{string}) distribution of the parameter values. The distribution type can be "normal", "logNormal", or "logitNormal".
###: \item formula: (\emph{string}) formula applied on individual parameters distribution
###: \item variability: a list giving, for each variability level, if individual parameters have variability or not
###: \item covariateModel: a list giving, for each individual parameter, if the related covariates are used or not.
###: If no covariate is used, this field is empty.
###: \item correlationBlocks : a list giving, for each variability level, the blocks of the correlation matrix of the random effects. 
###: A block is represented by a vector of individual parameter names. If there is no block, this field is empty.
###: }
###: @return A list of individual parameter model properties.
###: @examples
###: \dontrun{
###: indivModel = getIndividualParameterModel()
###: indivModel
###:  -> $name
###:       c("ka","V","Cl")
###:     $distribution
###:       c(ka = "logNormal", V = "normal", Cl = "logNormal")
###:     $formula
###:       "\\tlog(ka) = log(ka_pop) + eta_ka\\n\\n\\tlV = V_pop + eta_V\\n\\n\\tlog(Cl) = log(Cl_pop) + eta_Cl\\n\\n"
###:     $variability
###:       list( id = c(ka = TRUE, V = FALSE, Cl = TRUE) )
###:     $covariateModel
###:       list( ka = c(age = TRUE, sex = FALSE, wt = TRUE),
###:             V = c(age = FALSE, sex = FALSE, wt = FALSE),
###:             Cl = c(age = FALSE, sex = FALSE, wt = FALSE) )
###:     $correlationBlocks
###:       list( id = c("ka","V","Tlag") ) 
###: }
###: @seealso \code{\link{setIndividualParameterDistribution}} \code{\link{setIndividualParameterVariability}} \code{\link{setCovariateModel}}
###: @export
getIndividualParameterModel = function(){
  arguments = list()
  output = .processRequest("monolix", "getindividualparametermodels", arguments, "asynchronous")
  
  if (!is.null(output)){
    if (length(output$correlationBlocks) == 0)
      output$correlationBlocks <- NULL
    if (length(output$covariateModel) == 0)
      output$covariateModel <- NULL
  }
  
  return(output)
}

###: Set individual parameter distribution
###: 
###: Set the distribution of the estimated parameters.
###: Available distributions are "normal", "logNormal" and "logitNormal".\cr
###: Call \code{\link{getIndividualParameterModel}} to get a list of the available individual parameters within the current project.
###: @param ... A list of comma-separated pairs \{parameterName = (\emph{string})"distribution"\}.
###: @examples
###: \dontrun{
###: setIndividualParameterDistribution(V = "logNormal")
###: setIndividualParameterDistribution(Cl = "normal", V = "logNormal")
###: }
###: @seealso \code{\link{getIndividualParameterModel}}
###: @export
setIndividualParameterDistribution = function(...){
  arguments = list(...)
  if (is.null(names(arguments)) && length(arguments) == 1) 
    arguments <- arguments[[1]]
  parameterNames = names(arguments)
  
  if (length(parameterNames) == 0){
    if (length(arguments) > 0){
      .error("Unexpected type encountered. Please give a list of strings indexed by individual parameter names.")
      return(invisible(FALSE)) 
    } else {
      return(invisible(TRUE))
    }
  }
  else {
    for (i in 1:length(parameterNames)){
      if (parameterNames[i] == ""){
        .error(paste0("No individual parameter name found at position ",i,"."))
        return(invisible(FALSE))
      }
      else if (is.character(arguments[[parameterNames[i]]]) == FALSE){
        .error(paste0("Unexpected type encountered at position ",i,". Please give a string corresponding to an available individual parameter distribution name."))
        return(invisible(FALSE)) 
      }
    }
  }
  
  output = .processRequest("monolix", "setindividualdistribution", arguments, "synchronous", type = "STATUS")
  return(invisible(output))
}

###: Individual variability management
###: 
###: Add or remove inter-individual and/or intra-individual variability from some of the individual parameters present in the project.\cr
###: Call \code{\link{getIndividualParameterModel}} to get a list of the available parameters within the current project.
###: @param ... A list of comma-separated pairs \{variabilityLevel = \{individualParameterName = (\emph{bool})hasVariability\} \}.
###: @examples
###: \dontrun{
###: setIndividualParameterVariability(ka = TRUE, V = FALSE)
###: setIndividualParameterVariability(id = list(ka = TRUE), iov1 = list(ka = FALSE))
###: }
###: @seealso \code{\link{getIndividualParameterModel}}
###: @export
setIndividualParameterVariability = function(...){
  arguments = list(...)
  if (is.null(names(arguments)) && length(arguments) == 1){
    arguments <- arguments[[1]]
  }
  
  checkType <- function(input){
    indivParamNames = names(input)
    if (length(indivParamNames) == 0){
      if (length(input) > 0){
        .error("Unexpected type encountered. Please give a list of booleans indexed by individual parameter names.")
        return(invisible(FALSE)) 
      } else {
        return(invisible(TRUE))
      }
    }
    for (j in 1:length(indivParamNames)){
      if (indivParamNames[[j]] == "" || (is.logical(input[[indivParamNames[j]]]) == FALSE)){
        .error(paste0("Unexpected type encountered at position ",i,". Please give a list of comma-separated pairs {individualParameterName = (bool)hasVariability}."))
        return(invisible(FALSE)) 
      }
    }
    return(invisible(TRUE))
  }
  
  # one variability level case :
  if (length(arguments) > 0 && is.logical(arguments[[1]])){
    if (!checkType(arguments))
      return(invisible(FALSE))
  }
  
  # general case :
  else {
    levels = names(arguments)
    if (length(levels) == 0){
      if (length(arguments) > 0){
        .error("Unexpected type encountered. Please give a list indexed by variability level names.")
        return(invisible(FALSE)) 
      } else {
        return(invisible(TRUE))
      }
    }
    else {
      for (i in 1:length(levels)){
        if (levels[i] == ""){
          .error(paste0("No variability level found at position ",i,"."))
          return(invisible(FALSE))
        }
        else if (!checkType(arguments[[levels[i]]])){
          return(invisible(FALSE))
        }
      }
    }
  }
  
  output = .processRequest("monolix", "setindividualvariability", arguments, "synchronous", type = "STATUS")
  return(invisible(output))
}

###: Set covariate model
###: 
###: Set which are the covariates influencing individual parameters present in the project.
###: Call \code{\link{getIndividualParameterModel}} to get a list of the individual parameters present within the current project. 
###: and \code{\link{getCovariateInformation}} to know which are the available covariates for a given level of variability and a given individual parameter.
###: @param ... A list of comma-separated pairs \{parameterName = \{ covariateName = (\emph{bool})isInfluent, ...\} \}
###: @examples
###: \dontrun{
###: setCovariateModel( ka = c( Wt = FALSE, tWt = TRUE, lcat2 = TRUE),
###:                    Cl = c( SEX = TRUE )
###:                    )
###: }
###: @seealso \code{\link{getCovariateInformation}}
###: @export
setCovariateModel = function(...){
  arguments = list(...)
  if (length(arguments) == 1 && is.null(names(arguments)))
    arguments <- arguments[[1]]
  paramNames = names(arguments)
  
  if (length(paramNames) == 0){
    if (length(arguments) > 0){
      .error("Unexpected type encountered. Please give a list indexed by covariate names.")
      return(invisible(FALSE)) 
    } else {
      return(invisible(TRUE))
    }
  }
  else {
    for (i in 1:length(paramNames)){
      if (paramNames[i] == ""){
        .error(paste0("No individual parameter name found at position ",i,"."))
        return(invisible(FALSE))
      }
      else {
        covList = arguments[[paramNames[i]]]
        if (is.list(covList) == FALSE && is.vector(covList) == FALSE){
          .error(paste0("Unexpected type encountered at position ",i,". Please give a list, indexed by individual parameter names, of boolean lists."))
          return(invisible(FALSE)) 
        }
        else {
          covNames = names(covList)
          if (length(covNames) == 0){
            .error(paste0("Unexpected type encountered at position ",i,". Please give a list, indexed by individual parameter names, of boolean lists."))
            return(invisible(FALSE)) 
          }
          else {
            for (j in 1:length(covNames)){
              if (covNames[j] == ""){
                .error(paste0("Unexpected type encountered at position ",i,". Please give a list, indexed by individual parameter names, of boolean lists."))
                return(invisible(FALSE)) 
              }
              else if (length(covList[[covNames[j]]]) == 0){
                arguments[[paramNames[i]]][[covNames[j]]] <- NULL
              }
              else if (is.logical(covList[[covNames[j]]]) == FALSE){
                .error(paste0("Unexpected type encountered at position ",i,". Please give a list, indexed by individual parameter names, of boolean lists."))
                return(invisible(FALSE)) 
              }
            }
          }
        }
      }
    }
    
    output = .processRequest("monolix", "setcovariatemodel", arguments, "synchronous", type = "STATUS")
    return(invisible(output))
  }
}

###: Set correlation block structure
###: 
###: Define the correlation block structure associated to some of the variability levels of the current project.
###: Call \code{\link{getVariabilityLevels}} to get a list of the variability levels and \code{\link{getIndividualParameterModel}} to get a list of the available individual parameters within the current project.
###: @param ... A list of comma-separated pairs \{variabilityLevel = vector< (\emph{array<string>})parameterNames\} > \}.
###: @examples
###: \dontrun{
###: setCorrelationBlocks(id = list( c("ka","V","Tlag") ), iov1 = list( c("ka","Cl"), c("Tlag","V") ) )
###: }
###: @seealso \code{\link{getVariabilityLevels}} \code{\link{getIndividualParameterModel}}
###: @export
setCorrelationBlocks = function(...){
  arguments = list(...)
  if (is.null(names(arguments)) && length(arguments) == 1){
    arguments <- arguments[[1]]
  }
  
  checkType <- function(input){
    if (is.list(input) == FALSE){
      .error("Unexpected type encountered. Please give a list, indexed by variability level, of lists of individual parameters vectors.")
      return(invisible(FALSE)) 
    }
    else {
      if (length(input) > 0){
        for (j in 1:length(input)){
          if (is.null(input[[j]])){}
          else if (is.vector(input[[j]]) == FALSE){
            .error(paste0("Unexpected type encountered at position ",j,". Please give a list, indexed by variability level, of lists of individual parameters vectors."))
            return(invisible(FALSE)) 
          }
          else{
            group = input[[j]]
            for (k in 1:length(group)){
              if (is.character(group[[k]]) == FALSE){
                .error(paste0("Unexpected type encountered at position ",j,". Please give a list, indexed by variability level, of lists of individual parameters vectors."))
                return(invisible(FALSE)) 
              }
            } 
          }
        } 
      }
    }
    return(invisible(TRUE))
  }
  
  # one variability level case :
  if (length(arguments) > 0 && is.list(arguments[[1]]) == FALSE){
    if (!checkType(arguments))
      return (invisible(FALSE))
  }
  
  # general case :
  else {
    levels = names(arguments)
    if (length(levels) == 0){
      if (length(arguments) > 0){
        .error("Unexpected type encountered. Please give a list indexed by variability level names.")
        return(invisible(FALSE)) 
      } else {
        arguments = list(list(c()))
      }
    }
    else {
      for (i in 1:length(levels)){
        if (levels[i] == ""){
          .error(paste0("No variability level found at position ",i,"."))
          return(invisible(FALSE))
        }
        else if (!checkType(arguments[[levels[i]]])){
          return(invisible(FALSE))
        }
      }
    }
  }
  
  output = .processRequest("monolix", "setcorrelationgroup", arguments, "synchronous", type = "STATUS")
  return(invisible(output))
}
# ================================================================================ #