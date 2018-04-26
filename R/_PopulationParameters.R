# ============================= POPULATION PARAMETERS ============================ #
###: Get population parameters information
###: 
###: Get the name, the initial value, the estimation method and, if relevant, MAP parameters value of the population parameters present in the project. 
###: It is available for fixed effects, random effects, error model parameters, and latent covariates probabilities.
###: @return A data frame giving, for each population parameter, the corresponding :
###: \itemize{
###: \item initialValue : (\emph{double}) initial value
###: \item method : (\emph{string}) estimation method
###: \item priorValue : (\emph{double}) [MAP] typical value
###: \item priorSD : (\emph{double}) [MAP] standard deviation
###: }
###: @examples
###: \dontrun{
###: info = getPopulationParameterInformation()
###: info
###:     name      initialValue   method   typicalValue  stdDeviation
###:   ka_pop           1.0      MLE             NA            NA
###:   V_pop           10.0      MAP           10.0           0.5
###:   omega_ka         1.0     FIXED            NA            NA
###: }
###: @seealso \code{\link{setPopulationParameterInformation}}
###: @export
getPopulationParameterInformation = function(){
  arguments = list()
  info = .processRequest("monolix", "getpopulationparameterinformation", arguments, "asynchronous")
  
  if (!is.null(info)){
    output = data.frame(name = info$name, initialValue = as.numeric(info$initialValue), method = info$method, stringsAsFactors = FALSE)

    if (!is.null(info$priorValue)){
      nbParams = dim(output)[1]
      output$priorValue = rep(NA_real_,nbParams)
      output$priorSD = rep(NA_real_,nbParams)
      
      for (i in 1:length(info$priorValue)){
        index = which(output$name == names(info$priorValue)[[i]])
        if (index != 0){
          output[index,"priorValue"] <- as.numeric(info$priorValue[[i]])
          output[index,"priorSD"] <- as.numeric(info$priorSD[[i]]) 
        }
      }
    }
  }
  
  return(output)
}

###: Population parameters initialization and estimation method
###:
###: Set the initial value, the estimation method and, if relevant, the MAP parameters of one or several of the population parameters present within the current project (fixed effects + individual variances + error model parameters).
###: Available methods are:
###: \itemize{
###: \item "FIXED": Fixed
###: \item "MLE": Maximum Likelihood Estimation
###: \item "MAP": Maximum A Posteriori
###: }
###: Call \code{\link{getPopulationParameterInformation}} to get a list of the initializable population parameters present within the current project.
###: @param ... A list of comma-separated pairs \{paramName = list( initialValue = (\emph{double}), method = (\emph{string})"method"\}.
###: In case of "MAP" method, the user can specify the associated typical value and standard deviation values by using an additional list elements \{paramName = list( priorValue = (\emph{double})1, priorSD = (\emph{double})2 )\}.
###: By default, the prior value corresponds to the the population parameter and the prior standard deviation is set to 1.
###: @examples
###: \dontrun{
###: setPopulationParameterInformation(Cl_pop = list(initialValue = 0.5, method = "FIXED"), V_pop = list(intialValue = 1), ka_pop = list( method = "MAP", priorValue = 1.5, priorSD = 0.25 ) )
###: }
###: @seealso \code{\link{getPopulationParameterInformation}}
###: @export
setPopulationParameterInformation = function(...){
  arguments = list(...)
  
  if (is.null(names(arguments)) && length(arguments) == 1 && is.data.frame(arguments[[1]]) && "name" %in% names(arguments[[1]])){
    input = list()
    for (i in 1:nrow(arguments[[1]])){
      input[[arguments[[1]][i,"name"]]] <- as.list.data.frame(arguments[[1]][i,])
      input[[arguments[[1]][i,"name"]]]$name <- NULL
    }
    arguments <- input
  }

  parameterNames = names(arguments)
  if (length(parameterNames) == 0){
    if (length(arguments) > 0){
      .error("Unexpected type encountered. Please give a list of initialization properties indexed by population parameter names.")
      return(invisible(FALSE))
    } else {
      return(invisible(TRUE))
    }
  }
  else {
    for (i in 1:length(parameterNames)){
      if (parameterNames[i] == ""){
        .error(paste0("No population parameter name found at position ",i,"."))
        return(invisible(FALSE))
      }
      else if (is.list(arguments[[parameterNames[i]]])){
        fieldNames = names(arguments[[parameterNames[i]]])
        if (length(fieldNames) > 0){
          for (j in 1:length(fieldNames)){
            fieldValue = arguments[[parameterNames[i]]][[fieldNames[[j]]]]
            
            if (fieldNames[[j]] == "initialValue"){
              if (is.numeric(fieldValue) == FALSE){
                .error(paste0("Unexpected type encountered at position ",i," for field \"initialValue\". Please give a double."))
                return(invisible(FALSE))
              }
            } else if (fieldNames[[j]] == "method"){
              if (is.character(fieldValue) == FALSE){
                .error(paste0("Unexpected type encountered at position ",i," for field \"method\". Please give a string corresponding to a valid estimation method name."))
                return(invisible(FALSE))
              }
            } else if (fieldNames[[j]] == "priorValue"){
              if (is.numeric(fieldValue) == FALSE){
                .error(paste0("Unexpected type encountered at position ",i," for field \"priorValue\". Please give a double."))
                return(invisible(FALSE))
              }
            } else if (fieldNames[[j]] == "priorSD"){
              if (is.numeric(fieldValue) == FALSE){
                .error(paste0("Unexpected type encountered at position ",i," for field \"priorSD\". Please give a double."))
                return(invisible(FALSE))
              }
            } else {
              .error(paste0("\"",fieldNames[[j]],"\" is not a valid initialization property name."))
              return(invisible(FALSE))
            }
          }
        }
      } else {
        .error(paste0("Unexpected type encountered at position ",i,". Please give a list of initialization properties."))
        return(invisible(FALSE)) 
      }
    }
    
    output = .processRequest("monolix", "initializepopulationparameter", arguments, "synchronous", type = "STATUS")
    return(invisible(output))
  }
}

###: Initialize population parameters with the last estimated ones
###:
###: Set the initial value of all the population parameters present within the current project (fixed effects + individual variances + error model parameters) to the ones previously estimated.
###: These the values will be used in the population parameter estimation algorithm during the next scenario run.\cr
###: WARNING: If there is any set after a run, it will not be possible to set the initial values as the structure of the project has changed since last results.
###: @examples
###: \dontrun{
###: setInitialEstimatesToLastEstimates()
###: }
###: @seealso \code{\link{getEstimatedPopulationParameters}} \code{\link{getPopulationParameterInformation}} 
###: @export
setInitialEstimatesToLastEstimates = function(){
  arguments = list()
  output = .processRequest("monolix", "uselastestimates", arguments, "synchronous", type = "STATUS")
  return(invisible(output))
}