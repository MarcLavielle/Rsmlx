# ===================================== RESULTS ================================== #
###: Get tasks with results
###: 
###: Get a list of the tasks which have results to provide. A task is the association of:
###: \itemize{
###: \item an algorithm (string)
###: \item a vector of methods (string) relative to this algorithm for the standardErrorEstimation and the loglikelihoodEstimation, TRUE or FALSE for the other one.
###: }
###: @return The list of tasks with results, indexed by algorithm names.
###: @examples
###: \dontrun{
###: tasks = getLaunchedTasks()
###: tasks
###:  -> $populationParameterEstimation = TRUE
###:     $conditionalModeEstimation = TRUE
###: 	   $standardErrorEstimation = "linearization"
###: }
###: @export
getLaunchedTasks = function(){
  arguments = list()
  output = .processRequest("monolix", "getlaunchedtasks", arguments, "asynchronous")
  return(output)
}

# Estimates ---------------------------------------------------------------------- #
###: Get last estimated population parameter value
###: 
###: Get the last estimated value of some of the population parameters present within the current project (fixed effects + individual variances + correlations + latent probabilities + error model parameters).\cr
###: WARNING: Estimated population parameters values cannot be accessible until the SAEM algorithm has been launched once.
###: @param ... [optional] (\emph{array<string>}) Names of the population parameters whose value must be displayed. Call \code{\link{getPopulationParameterInformation}} to get a list of the population parameters present within the current project.
###: If this field is not specified, the function will retrieve the values of all the available population parameters.
###: @return A named vector containing the last estimated value of each one of the population parameters passed in argument.
###: @examples
###: \dontrun{
###: getEstimatedPopulationParameters("V_pop") -> [V_pop = 0.5]
###: getEstimatedPopulationParameters("V_pop","Cl_pop") -> [V_pop = 0.5, Cl_pop = 0.25]
###: getEstimatedPopulationParameters() -> [V_pop = 0.5, Cl_pop = 0.25, ka_pop = 0.05]
###: }
###: @export
getEstimatedPopulationParameters = function(...){
  arguments = list(...)
  if ( length(arguments) > 0){
    for (i in 1:length(arguments)){
      if (is.character(arguments[[i]]) == FALSE){
        .error(paste0("Unexpected type encountered at position ",i,". Please give strings corresponding to a list of population parameter names."))
        return(invisible(FALSE))    
      }
    }
  }
  
  output = .processRequest("monolix", "getpopulationparameters", arguments, "asynchronous")
  return(output)
}

###: Get last estimated individual parameter values
###: 
###: Get the last estimated values for each subject of some of the individual parameters present within the current project.\cr
###: WARNING: Estimated individual parameters values cannot be accessible until the individual estimation algorithm has been launched once.\cr
###: NOTE: The user can choose to display only the individual parameter values estimated with a specific method.\cr
###: Existing individual estimation methods :
###: \tabular{ll}{
###: Conditional Mean SAEM \tab "saem" \cr
###: Conditional Mean \tab "conditionalMean" \cr
###: Conditional Mode \tab "conditionalMode"
###: }
###: WARNING: Only the methods which have been used during the last scenario run can provide estimation results.
###: @param ... (\emph{string}) Name of the individual parameters whose values must be displayed. Call \code{\link{getIndividualParameterModel}} to get a list of the individual parameters present within the current project. 
###: @param method [optional](\emph{string}) Individual parameter estimation method whose results should be displayed.
###: If there are latent covariate used in the model, the estimated modality is displayed too
###: If this field is not specified, the results provided by all the methods used during the last scenario run are displayed.
###: @return A data frame giving, for each wanted method, the last estimated values of the individual parameters of interest for each subject with the corresponding standard deviation values.
###: @examples
###: \dontrun{ 
###: indivParams = getEstimatedIndividualParameters() # retrieve the values of all the available individual parameters for all methods
###:   -> $saem
###:       id   Cl     V      ka
###:       1   0.28  7.71   0.29
###:       .   ...    ...    ...
###:       N   0.1047.62   1.51
###:       
###: indivParams = getEstimatedIndividualParameters("Cl", "V", method = "conditionalMean") # retrieve the values of the individual parameters "Cl" and "V" estimated by the conditional mode method
###: }
###: @seealso \code{\link{getEstimatedRandomEffects}} 
###: @export
getEstimatedIndividualParameters = function(..., method = ""){
  if (is.character(method) == FALSE){
    .error("Unexpected type encountered for \"method\" field. Please give a string corresponding to an individual estimation algorithm method.")
    method = ""   
  }
  
  parameterNames = list(...)
  if ( length(parameterNames) > 0){
    for (i in 1:length(parameterNames)){
      if (is.character(parameterNames[[i]]) == FALSE){
        .error(paste0("Unexpected type encountered at position ",i,". Please give strings corresponding to a list of individual parameter names."))
        return(invisible(FALSE))    
      }
    }
  }
  if (method != ""){
    arguments = list(list(method),parameterNames)
  } else if (length(parameterNames) != 0){
    arguments = list(list(),parameterNames)
  } else {
    arguments = list()
  }
  
  params = .processRequest("monolix", "getindividualparameters", arguments, "asynchronous")
  if (is.null(params))
    return(NULL)
  
  tryCatch(expr = {
    output = list()
    methodNames = names(params$parameters)
    
    for (iMeth in 1:length(methodNames)){
      method = methodNames[[iMeth]]
      
      if (method == "conditionalSD"){
        output[[method]] <- as.data.frame( append(params$ids,params$parameters[[method]]) )
      } else {
        tmp <- append(params$ids,params$parameters[[method]])
        
        if ("covariate" %in% names(params))
          tmp <- append(tmp,params$covariate)
        if ("bsmm" %in% names(params))
          tmp <- append(tmp,params$bsmm)
        
        output[[method]] <- as.data.frame(tmp)
      }
    }
  },
  error = function(err){
    output = list()
    .error("Bad output encountered.")
  })
  
  return(output)
}

###: Get estimated the random effects
###: 
###: Get the random effects for each subject of some of the individual parameters present within the current project.\cr
###: WARNING: Estimated random effects cannot be accessible until the individual estimation algorithm has been launched once.\cr
###: The user can choose to display only the random effects estimated with a specific method.\cr
###: NOTE: The random effects are defined in the gaussian referential, e.g. if ka is lognormally distributed around ka_pop, eta_i = log(ka_i)-log(ka_pop)
###: Existing individual estimation methods :
###: \tabular{ll}{
###: Conditional Mean SAEM \tab "saem" \cr
###: Conditional Mean \tab "conditionalMean" \cr
###: Conditional Mode \tab "conditionalMode"
###: }
###: WARNING: Only the methods which have been used during the last scenario run can provide estimation results. Please call \code{\link{getLaunchedTasks}} to get a list of the methods whose results are available.
###: @param ... (\emph{string}) Name of the individual parameters whose random effects must be displayed. Call \code{\link{getIndividualParameterModel}} to get a list of the individual parameters present within the current project. 
###: @param method [optional](\emph{string}) Individual parameter estimation method whose results should be displayed.
###: If this field is not specified, the results provided by all the methods used during the last scenario run are displayed.
###: @return A data frame giving, for each wanted method, the last estimated eta values of the individual parameters of interest for each subject with the corresponding standard deviation values.
###: @examples
###: \dontrun{ 
###: etaParams = getEstimatedRandomEffects() # retrieve the values of all the available random effects for all methods, without the associated standard deviations
###:   -> $saem
###:      id    Cl     V      ka
###:       1   0.28  7.71   0.29
###:       .   ...    ...    ...
###:       N   0.1047.62   1.51
###:       
###: etaParams = getEstimatedRandomEffects("Cl", "V", method = "conditionalMode") # retrieve the values of the individual parameters "Cl" and "V" estimated by the conditional mean from SAEM algorithm
###: }
###: @seealso \code{\link{getEstimatedIndividualParameters}}  
###: @export
getEstimatedRandomEffects = function(...,method = ""){
  if (is.character(method) == FALSE){
    .error("Unexpected type encountered for \"method\" field. Please give a string corresponding to an individual estimation algorithm method.")
    return(invisible(FALSE))    
  }
  
  parameterNames = list(...)
  if ( length(parameterNames) > 0){
    for (i in 1:length(parameterNames)){
      if (is.character(parameterNames[[i]]) == FALSE){
        .error(paste0("Unexpected type encountered at position ",i,". Please give strings corresponding to a list of individual parameter names."))
        return(invisible(FALSE))    
      }
    }
  }
  if (method != ""){
    arguments = list(list(method),parameterNames)
  } else if (length(parameterNames) != 0){
    arguments = list(list(),parameterNames)
  } else {
    arguments = list()
  }
  
  params = .processRequest("monolix", "getetaparameters", arguments, "asynchronous")
  if (is.null(params))
    return(NULL)
  
  tryCatch(expr = {
    toUser = function(x){
      out <- as.data.frame( append(params$ids,x) )
      return(out)
    }
    output <- lapply(params$parameters, function(x) toUser(x))
  },
  error = function(err){
    output = list()
    .error("Bad output encountered.")
  })
  
  return(output)
}

###: Get simulated individual parameters
###: 
###: Get the simulated values for each replicate of each subject of some of the individual parameters present within the current project.\cr
###: WARNING: Simulated individual parameters values cannot be accessible until the individual estimation with conditional mean algorithm has been launched once.\cr
###: @param ... (\emph{string}) Name of the individual parameters whose values must be displayed. Call \code{\link{getIndividualParameterModel}} to get a list of the individual parameters present within the current project. 
###: @return A list giving the last simulated values of the individual parameters of interest for each replicate of each subject.
###: @examples
###: \dontrun{ 
###: simParams = getSimulatedIndividualParameters() # retrieve the values of all the available individual parameters
###: simParams
###:      rep   id    Cl     V     ka
###:       1    1   0.022  0.37  1.79
###:       1    2   0.033  0.42  -0.92
###:       .    .    ...    ...  ...
###:       2    1   0.021  0.33  1.47
###:       .    .    ...    ...  ...
###: }
###: @seealso \code{\link{getSimulatedRandomEffects}}
###: @export
getSimulatedIndividualParameters = function(...){
  parameterNames = list(...)
  if ( length(parameterNames) > 0){
    for (i in 1:length(parameterNames)){
      if (is.character(parameterNames[[i]]) == FALSE){
        .error(paste0("Unexpected type encountered at position ",i,". Please give strings corresponding to a list of individual parameter names."))
        return(invisible(FALSE))    
      }
    }
  }

  results = .processRequest("monolix", "getsimulatedindividualparameters", parameterNames, "asynchronous")
  if (is.null(results)){
    return(NULL)
  }
  
  tryCatch(expr = {
    nbRep = length(results$parameters[[1]])
    nbSO = length(results$parameters[[1]][[1]])
    output = data.frame(matrix(unlist(results$parameters), nrow = nbRep*nbSO))
    names(output) = names(results$parameters)
    
    tmp = list(rep = rep((1:nbRep),each=nbSO))
    
    lvlNames = names(results$ids)
    for (i in length(lvlNames):1)
      tmp[[lvlNames[[i]]]] <- rep(results$ids[[i]],nbRep)
    
    output <- cbind(tmp,output)
  },
  
  error = function(err){
    output = list()
    .error("Bad output encountered.")
  })
  
  return(output)
}

###: Get simulated random effects
###: 
###: Get the simulated values for each replicate of each subject of some of the individual random effects present within the current project.\cr
###: WARNING: Simulated individual random effects values cannot be accessible until the individual estimation algorithm with conditional mean has been launched once.\cr
###: @param ... (\emph{string}) Name of the individual parameters whose values must be displayed. Call \code{\link{getIndividualParameterModel}} to get a list of the individual parameters present within the current project. 
###: @return A list giving the last simulated values of the individual random effects of interest for each replicate of each subject.
###: @examples
###: \dontrun{ 
###: simEtas = getSimulatedRandomEffects() # retrieve the values of all the available individual random effects
###: simEtas
###:      rep   id    Cl     V     ka
###:       1    1   0.022  0.37  1.79
###:       1    2   0.033  0.42  -0.92
###:       .    .    ...    ...  ...
###:       2    1   0.021  0.33  1.47
###:       .    .    ...    ...  ...
###: }
###: @seealso \code{\link{getIndividualParameterModel}}
###: @export
getSimulatedRandomEffects = function(...){
  parameterNames = list(...)
  if ( length(parameterNames) > 0){
    for (i in 1:length(parameterNames)){
      if (is.character(parameterNames[[i]]) == FALSE){
        .error(paste0("Unexpected type encountered at position ",i,". Please give strings corresponding to a list of individual parameter names."))
        return(invisible(FALSE))    
      }
    }
  }
  
  results = .processRequest("monolix", "getsimulatedetaparameters", parameterNames, "asynchronous")
  if (is.null(results)){
    return(NULL)
  }
  
  tryCatch(expr = {
    nbRep = length(results$parameters[[1]])
    nbSO = length(results$parameters[[1]][[1]])
    output = data.frame(matrix(unlist(results$parameters), nrow = nbRep*nbSO))
    names(output) = names(results$parameters)
    
    tmp = list(rep = rep((1:nbRep),each=nbSO))
    
    lvlNames = names(results$ids)
    for (i in length(lvlNames):1)
      tmp[[lvlNames[[i]]]] <- rep(results$ids[[i]],nbRep)
    
    output <- cbind(tmp,output)
  },
  
  error = function(err){
    output = list()
    .error("Bad output encountered.")
  })
  
  return(output)
}
# -------------------------------------------------------------------------------- #

# Algorithms Output -------------------------------------------------------------- #
###: Get SAEM algorithm iterations
###: 
###: Retrieve the successive values of some of the population parameters present within the current project (fixed effects + individual variances + correlations + latent probabilities + error model parameters) during the previous run of the SAEM algorithm.\cr
###: WARNING: Convergence history of population parameters values cannot be accessible until the SAEM algorithm has been launched once.
###: @param ... [optional] (\emph{array<string>}) Names of the population parameters whose convergence history must be displayed. Call \code{\link{getPopulationParameterInformation}} to get a list of the population parameters present within the current project.
###: If this field is not specified, the function will retrieve the values of all the available population parameters.
###: @return A list containing a pair composed by the number of exploratory and smoothing iterations and a data frame which associates each wanted population parameter to its successive values over SAEM algorithm iterations.
###: @examples
###: \dontrun{ 
###: report = getSAEMiterations()
###: report
###:   -> $iterationNumbers
###:        c(50,25)
###:      $estimates
###:           V    Cl
###:        0.25     0
###:         0.3   0.5
###:           .     .
###:        0.35  0.25
###: }
###: @export
getSAEMiterations = function(...){
  arguments = list(...)
  if ( length(arguments) > 0){
    for (i in 1:length(arguments)){
      if (is.character(arguments[[i]]) == FALSE){
        .error(paste0("Unexpected type encountered at position ",i,". Please give strings corresponding to a list of population parameter names."))
        return(invisible(FALSE))    
      }
    }
  }
  
  output = .processRequest("monolix", "getsaemconvergence", arguments, "asynchronous")
  
  if (!is.null(output))
    output$estimates <- as.data.frame(output$estimates, row.names = as.character(c(0:sum(output$iterationNumbers))))
  
  return(output)
}

###: Get the inverse of the Fisher Matrix
###: 
###: Get the inverse of the last estimated Fisher matrix computed either by all the Fisher methods used during the last scenario run or by the specific one passed in argument.\cr
###: WARNING: The Fisher matrix cannot be accessible until the Fisher algorithm has been launched once.\cr
###: The user can choose to display only the Fisher matrix estimated with a specific method.\cr
###: Existing Fisher methods :
###: \tabular{ll}{
###: Fisher by Linearization \tab "linearization" \cr
###: Fisher by Stochastic Approximation \tab "stochasticApproximation" \cr
###: }
###: WARNING: Only the methods which have been used during the last scenario run can provide results.
###: @param method [optional](\emph{string}) Fisher method whose results should be displayed.
###: If this field is not specified, the results provided by all the methods used during the last scenario run are displayed.
###: @return  A list whose each field contains the Fisher matrix computed by one of the available Fisher methods used during the ast scenario run.
###: A matrix is defined as a structure containing the following fields :
###: \tabular{ll}{
###: rownames \tab list of row names \cr
###: columnnames \tab list of column names \cr
###: rownumber \tab number of rows \cr
###: data \tab vector<...> containing matrix raw values (column major) \cr
###: }
###: @examples
###: \dontrun{
###: getCorrelationOfEstimates("linearization")
###:  -> list( linearization = list( data = c(1,0,0,0,1,-0.06,0,-0.06,1), rownumber = 3, rownames = c("Cl_pop","omega_Cl","a"), columnnames = c("Cl_pop","omega_Cl","a") )  )
###:  
###: getCorrelationOfEstimates() -> list( linearization = list(...), stochasticApproximation = list(...) )
###: }
###: @export
getCorrelationOfEstimates = function(method = ""){
  
  if (is.character(method) == FALSE){
    .error("Unexpected type encountered for \"method\" field. Please give a string corresponding to a Fisher algorithm method.")
    return(invisible(FALSE))    
  } else if (method != "") {
    arguments = list(method)
  } else {
    arguments = list()
  }
  
  output = .processRequest("monolix", "getfishermatrix", arguments, "asynchronous")
  
  if (!is.null(output)){
    for (i in 1:length(output))
      output[[i]] = .buildMatrix(output[[i]])
  }

  return(output)
}

###: Get standard errors of population parameters
###: 
###: Get the last estimated standard errors of population parameters computed either by all the Fisher methods used during the last scenario run or by the specific one passed in argument.\cr
###: WARNING: The standard errors cannot be accessible until the Fisher algorithm has been launched once.\cr
###: Existing Fisher methods :
###: \tabular{ll}{
###: Fisher by Linearization \tab "linearization" \cr
###: Fisher by Stochastic Approximation \tab "stochasticApproximation" \cr
###: }
###: WARNING: Only the methods which have been used during the last scenario run can provide results.
###: @param method [optional](\emph{string}) Fisher method whose results should be displayed.
###: If this field is not specified, the results provided by all the methods used during the last scenario run are retrieved
###: @return A list associating each retrieved Fisher algorithm method to the standard errors of population parameters computed during its last run.
###: @examples
###: \dontrun{
###: getEstimatedStandardErrors() -> list( linearization = [...], stochasticApproximation = [...] )
###: getEstimatedStandardErrors("linearization") -> list( linearization = [...] )
###: }
###: @export
getEstimatedStandardErrors = function(method = ""){
  if (is.character(method) == FALSE){
    .error("Unexpected type encountered for \"method\" field. Please give a string corresponding to a Fisher algorithm method.")
    return(invisible(FALSE))    
  } else if (method != "") {
    arguments = list(method)
  } else {
    arguments = list()
  }
  
  output = .processRequest("monolix", "getstandarderrors", arguments, "asynchronous")
  return(output)
}

###: Get Log-Likelihood values
###: 
###: Get the values computed by using a log-likelihood algorithm during the last scenario run, with or without a method-based filter.\cr
###: WARNING: The log-likelihood values cannot be accessible until the log-likelihood algorithm has been launched once.\cr
###: The user can choose to display only the log-likelihood values computed with a specific method.\cr
###: Existing log-likelihood methods :
###: \tabular{ll}{
###: Log-likelihood by Linearization \tab "linearization" \cr
###: Log-likelihood by Important Sampling \tab "importanceSampling" \cr
###: }
###: WARNING: Only the methods which have been used during the last scenario run can provide results.
###: @param method [optional](\emph{string}) Log-likelihood method whose results should be displayed.
###: If this field is not specified, the results provided by all the methods used during the last scenario run are retrieved.
###: @return  A list associating the name of each method passed in argument to the corresponding log-likelihood values computed by during the last scenario run.
###: @examples
###: \dontrun{
###: getEstimatedLogLikelihood()
###:  -> list(  linearization = [LL = -170.505, AIC = 350.280, BIC = 365.335] ,
###:            importanceSampling = [...] )
###:            
###: getEstimatedLogLikelihood("linearization")
###:  -> list(  linearization = [LL = -170.505, AIC = 350.280, BIC = 365.335] )
###:  }
###: @export
getEstimatedLogLikelihood = function(method = ""){
  if (is.character(method) == FALSE){
    .error("Unexpected type encountered for \"method\" field. Please give a string corresponding to a log-likelihood algorithm method.")
    return(invisible(FALSE))    
  } else if (method != "") {
    arguments = list(method)
  } else {
    arguments = list()
  }
  
  output = .processRequest("monolix", "getllvalues", arguments, "asynchronous")
  return(output)
}
# -------------------------------------------------------------------------------- #
# ================================================================================ #