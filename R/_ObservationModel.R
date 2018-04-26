# =============================== OBSERVATION MODELS ============================= #
###: Get observations information
###: 
###: Get the name, the type and the values of the observations present in the project.
###: @return A list containing the name of the observations, their type and their values (id, time and observationName (and occasion if present in the data set)).
###: @examples
###: \dontrun{
###: info = getObservationInformation()
###: info
###:   -> $name
###:      c("concentration")
###:   -> $type
###:      c(concentration = "continuous")
###:   -> $concentration
###:        id   time concentration
###:         1    0.5     0.0
###:         .    .      .
###:         N    9.0    10.8
###: }
###: @export
getObservationInformation = function(){
  arguments = list()
  result = .processRequest("monolix", "getobservationsinformation", arguments, "asynchronous")
  output = NaN
  
  if (!is.null(result$value) && length(result$value) > 0){
    output = list(name = result$name, type = result$type)
    lvlNames = names(result$ids)
    
    for (iModel in 1:length(result$value)){
      name = result$name[[iModel]]
      time = c()
      value = c()
      ids = list()
      
      for (iSO in 1:length(result$value[[iModel]]$time)){
        time = append(time, result$value[[iModel]]$time[[iSO]])
        value = append(value, result$value[[iModel]]$value[[iSO]])
        
        for (iLvl in 1:length(lvlNames))
          ids[[lvlNames[[iLvl]]]] = append( ids[[lvlNames[[iLvl]]]],
                                            rep(result$ids[[lvlNames[[iLvl]]]][[iSO]], length(result$value[[iModel]]$time[[iSO]])) )
      }
      
      output[[name]] = list(time = time)
      output[[name]][[name]] <- value
      output[[name]] = as.data.frame( append(ids, output[[name]]) )
    }
  }
  
  return(output)
}

###: Get continuous observation models information
###: 
###: Get a summary of the information concerning the continuous observation models in the project. The following informations are provided.
###: \itemize{
###: \item prediction: (\emph{vector<string>}) name of the associated prediction
###: \item formula: (\emph{vector<string>}) formula applied on the observation
###: \item distribution: (\emph{vector<string>}) distribution of the observation in the Gaussian space. The distribution type can be "normal", "logNormal", or "logitNormal".
###: \item limits: (\emph{vector< pair<double,double> >}) lower and upper limits imposed to the observation. 
###: Used only if the distribution is logitNormal. If there is no logitNormal distribution, this field is empty.
###: \item errormodel: (\emph{vector<string>}) type of the associated error model
###: \item autocorrelation: (\emph{vector<bool>}) defines if there is auto correlation
###: }
###: Call \code{\link{getObservationInformation}} to get a list of the continuous observations present in the current project.
###: @return A list associating each continuous observation to its model properties.
###: @examples
###: \dontrun{
###: obsModels = getContinuousObservationModel()
###: obsModels
###:  -> $prediction
###:       c(Conc = "Cc")
###:     $formula
###:       c(Conc = "Conc = Cc + (a+b*Cc)*e")
###:     $distribution
###:       c(Conc = "logitNormal")
###:     $limits
###:       list(Conc = c(0,11.5))
###:     $errormodel
###:       c(Conc = "combined1")
###:     $autocorrelation
###:       c(Conc = TRUE)
###: }
###: @seealso \code{\link{getObservationInformation}} \code{\link{setObservationDistribution}} \code{\link{setObservationLimits}}
###: \code{\link{setErrorModel}} \code{\link{setAutocorrelation}}
###: @export
getContinuousObservationModel = function(){
 arguments = list()
 output = .processRequest("monolix", "getobservationmodels", arguments, "asynchronous")

  if (length(output) != 0){
   if (length(output$limits) == 0)
     output$limits <- NULL
   else
     output$limits = lapply(output$limits, as.numeric)
  }
 else
   output = NULL
 
  return(output)
}

###: Set observation model distribution limits
###: 
###: Set the minimum and the maximum values between which some of the observations can be found.
###: Used only if the distribution of the error model is "logitNormal", else wise it will not be taken into account
###: @param ... A list of comma-separated pairs \{observationModel = [(double)min,(double)max] \}
###: @examples
###: \dontrun{
###: setObservationLimits( Conc = c(-Inf,Inf), Effect = c(0,Inf) )
###: }
###: @seealso \code{\link{getContinuousObservationModel}} \code{\link{getObservationInformation}}
###: @export
setObservationLimits = function(...){
  arguments = list(...)
  if (length(arguments) == 0)
    return(invisible(TRUE))
  else if (length(arguments) == 1 && is.vector(arguments[[1]]) && length(names(arguments[[1]])) > 0)
    arguments = as.list(arguments[[1]])
  
  obsModelNames = names(arguments)
  if (length(obsModelNames) == 0){
    .error("No observation model names found.")
    return(invisible(FALSE))
  }
  else {
    for (i in 1:length(obsModelNames)){
      if (obsModelNames[i] == ""){
        .error(paste0("No observation model name found at position ",i,"."))
        return(invisible(FALSE))
      }
      else { 
        limits = arguments[[obsModelNames[i]]]
        if ((!is.list(limits) && !is.vector(limits)) || length(limits) != 2){
          .error(paste0("Unexpected type encountered at position ",i,
                         ". Please give a list of pairs of doubles [min,max]."))
          return(invisible(FALSE)) 
        }
        else if (!is.numeric(limits[[1]]) || !is.numeric(limits[[2]]) || limits[[1]] > limits[[2]]){
          .error(paste0("Unexpected type encountered at position ",i,
                         ". Please give a list of pairs of doubles [min,max]."))
          return(invisible(FALSE)) 
        }
      }
    }
    
    output = .processRequest("monolix", "setobservationmodellimits", arguments, "synchronous", type = "STATUS")
    return(invisible(output))
  }
}

###: Set observation model distribution
###: 
###: Set the distribution in the Gaussian space of some of the observation models.
###: Available distribution types are "normal", "logNormal", or "logitNormal".
###: Call \code{\link{getObservationInformation}} to get a list of the available observation models within the current project.
###: @param ... A list of comma-separated pairs \{observationModel = (\emph{string})"distribution"\}.
###: @examples 
###: \dontrun{
###: setObservationDistribution(Conc = "normal")
###: setObservationDistribution(Conc = "normal", Effect = "logNormal")
###: }
###: @seealso \code{\link{getContinuousObservationModel}}
###: @export
setObservationDistribution = function(...){
  arguments = list(...)
  if (length(arguments) == 0)
    return(invisible(TRUE))
  else if (length(arguments) == 1 && is.vector(arguments[[1]]) && length(names(arguments[[1]])) > 0)
    arguments = as.list(arguments[[1]])
  
  obsModelNames = names(arguments)
  if (length(obsModelNames) == 0){
    .error("No observation model names found.")
    return(invisible(FALSE))
  }
  else {
    for (i in 1:length(obsModelNames)){
      if (obsModelNames[i] == ""){
        .error(paste0("No observation model name found at position ",i,"."))
        return(invisible(FALSE))
      }
      else if (!is.character(arguments[[obsModelNames[i]]])){
        .error(paste0("Unexpected type encountered at position ",i,
                       ". Please give a string corresponding to an available observation distribution name."))
        return(invisible(FALSE)) 
      }
    }
    
    output = .processRequest("monolix", "setobservationdistribution", arguments, "synchronous", type = "STATUS")
    return(invisible(output))
  }
}

###: Set error model
###: 
###: Set the error model type to be used with some of the observation models. 
###: Call \code{\link{getObservationInformation}} to get a list of the observation models present in the current project.
###: 
###: Available error model types are :
###: \tabular{ll}{
###: "constant" \tab obs = pred + a*err\cr
###: "proportional" \tab obs = pred + (b*pred)*err\cr
###: "combined1" \tab obs = pred + (b*pred^c + a)*err\cr
###: "combined2" \tab obs = pred + sqrt(a^2 + (b^2)*pred^(2c))*err\cr
###: }\cr
###: Error model parameters will be initialized to 1 by default.
###: Call \code{\link{setPopulationParameterInformation}} to modify their initial value.\cr
###: The value of the exponent parameter is fixed by default when using the "combined1" and "combined2" models.\cr 
###: Use \code{\link{setPopulationParameterInformation}} to enable its estimation.\cr
###: @param ... A list of comma-separated pairs \{observationModel = (\emph{string})errorModelType\}.
###: @examples
###: \dontrun{
###: setErrorModel(Conc = "constant", Effect = "combined1")
###: }
###: @seealso \code{\link{getContinuousObservationModel}} \code{\link{setPopulationParameterInformation}}
###: @export
setErrorModel = function(...){
  arguments = list(...)
  if (length(arguments) == 0)
    return(invisible(TRUE))
  else if (length(arguments) == 1 && is.vector(arguments[[1]]) && length(names(arguments[[1]])) > 0)
    arguments = as.list(arguments[[1]])
  
  obsModelNames = names(arguments)
  if (length(obsModelNames) == 0){
    .error("No observation model names found.")
    return(invisible(FALSE))
  }
  else {
    for (i in 1:length(obsModelNames)){
      if (obsModelNames[i] == ""){
        .error(paste0("No observation model name found at position ",i,"."))
        return(invisible(FALSE))
      }
      else if (!is.character(arguments[[obsModelNames[i]]])){
        .error(paste0("Unexpected type encountered at position ",i,
                       ". Please give a string corresponding to a valid error model type."))
        return(invisible(FALSE))
      }
    }
    
    output = .processRequest("monolix", "seterrormodeltype", arguments, "synchronous", type = "STATUS")
    return(invisible(output))
  }
}

###: Set auto-correlation
###: 
###: Add or remove auto-correlation from the error model used on some of the observation models. \cr
###: Call \code{\link{getObservationInformation}} to get a list of the observation models present in the current project.
###: @param ... Sequence of comma-separated pairs \{(\emph{string})"observationModel",(\emph{boolean})hasAutoCorrelation\}.
###: @examples
###: \dontrun{
###: setAutocorrelation(Conc = TRUE)
###: setAutocorrelation(Conc = TRUE, Effect = FALSE)
###: }
###: @seealso \code{\link{getContinuousObservationModel}}
###: @export
setAutocorrelation = function(...){
  arguments = list(...)
  if (length(arguments) == 0)
    return(invisible(TRUE))
  else if (length(arguments) == 1 && is.vector(arguments[[1]]) && length(names(arguments[[1]])) > 0)
    arguments = as.list(arguments[[1]])
  
  obsModelNames = names(arguments)
  if (length(obsModelNames) == 0){
    .error("No observation model names found.")
    return(invisible(FALSE))
  }
  else {
    for (i in 1:length(obsModelNames)){
      if (obsModelNames[i] == ""){
        .error(paste0("No observation model name found at position ",i,"."))
        return(invisible(FALSE))
      }
      else if (!is.logical(arguments[[obsModelNames[i]]])){
        .error(paste0("Unexpected type encountered at position ",i,
                      ". Please give a boolean equaling TRUE to add autocorrelation to the error model used on the preceding observation model, FALSE to remove it."))
        return(invisible(FALSE))
      }
    }
    
    output = .processRequest("monolix", "setautocorrelation", arguments, "synchronous", type = "STATUS")
    return(invisible(output))
  }
}
# ================================================================================ #