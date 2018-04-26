# ============================== SCENARIO MANAGEMENT ============================= #
# Building Scenario -------------------------------------------------------------- #
###: Get current scenario
###:
###: Get the list of tasks that will be run at the next call to \code{\link{runScenario}}, the associated method (linearization true or false), and the associated list of plots. 
###: The list of tasks consist of the following tasks: populationParameterEstimation, conditionalDistributionSampling, conditionalModeEstimation, standardErrorEstimation, logLikelihoodEstimation, and plots.
###:
###: @return The list of tasks that corresponds to the current scenario, indexed by algorithm names.
###: @examples 
###: \dontrun{
###: scenario = getScenario()
###: scenario
###:  -> $tasks 
###:     populationParameterEstimation conditionalDistributionSampling conditionalModeEstimation standardErrorEstimation logLikelihoodEstimation plots
###:     TRUE                          TRUE							 TRUE                      FALSE                   FALSE                   FALSE 
###:     $linearization = T
###:     $plotList = "outputplot", "vpc"
###: }
###: @seealso \code{\link{setScenario}}
###: @export
getScenario = function(){
  arguments = list()
  output = .processRequest("monolix", "getscenario", arguments, "asynchronous")
  return(output)
}

###: Set scenario
###:
###: Clear the current scenario and build a new one from a given list of tasks, the linearization option and the list of plots. A task is the association of:
###: \itemize{
###: \item a task
###: \item a boolean
###: }
###: \subsection{NOTE}{by default the boolean is false, tus, the user can only state what will run during the scenario.}
###: 
###: \subsection{NOTE}{
###: Within a MONOLIX scenario, the order according which the different algorithms are run is fixed:
###: \tabular{ll}{
###: Algorithm \tab Algorithm Keyword \cr
###: Population Parameter Estimation \tab "populationParameterEstimation" \cr
###: Conditional Mode Estimation (EBEs) \tab "conditionalModeEstimation"  \cr
###: Sampling from the Conditional Distribution \tab "conditionalDistributionSampling" \cr
###: Standard Error and Fisher Information Matrix Estimation \tab "standardErrorEstimation" \cr
###: LogLikelihood Estimation \tab "logLikelihoodEstimation" \cr
###: Plots \tab "plots" \cr
###: }
###: }
###: @examples
###: \dontrun{
###: 
###: scenario = getScenario()
###: scenario$tasks = c(populationParameterEstimation = T, conditionalModeEstimation = T, conditionalDistributionSampling = T)
###: setScenario(scenario)
###: }
###: @seealso \code{\link{getScenario}}
###: @export
setScenario = function(...){
  arguments = list(...)
  if (is.null(names(arguments)) && length(arguments) == 1 && is.list(arguments[[1]])) 
    arguments <- arguments[[1]]
  
  if (length(arguments) == 0){
      return(invisible(TRUE))
  }
  else {
    scenario = list()
    
    tasks = names(arguments)
    
    for (i in 1:length(arguments)){
      
      if (is.null(tasks[i]) || tasks[i] == ""){
        .error(paste0("Unexpected type encountered at position ", i, ". Please give a string corresponding to a valid section name."))
        return(invisible(FALSE))
      }
      
      for(j in 1:length(arguments[[i]])){
        if (is.null(arguments[[i]][j]) || arguments[[i]][j] == ""){
          .error(paste0("Unexpected type encountered at position ", j, ". Please give a string corresponding to a valid section name."))
          return(invisible(FALSE))
        }
      }
      
      if (length(arguments[[i]]) != 0 && is.vector(arguments[[i]]) == FALSE){
        .error(paste0("Unexpected type encountered at position ",i,". Please give a vector of strings corresponding to valid section method names."))
        return(invisible(FALSE))
      } 
      scenario[[tasks[i]]] = arguments[[i]]
    }
  }

  output = .processRequest("monolix", "setScenario", scenario, "synchronous", type = "STATUS")
  return(invisible(output))
}
# -------------------------------------------------------------------------------- #

# Running Scenario --------------------------------------------------------------- #
###: Run current scenario
###:
###: Run the current scenario. By default, this task is processed sequentially.
###: Call 
###: \enumerate{
###: \item \code{\link{isRunning}} to check if the scenario is still running and get information about the current task,
###: \item \code{\link{abort}} to stop the execution.
###: }
###: To launch the function in the background, so that functions which do not modify the project ("get" functions for example) remains available, set the input argument "wait" to FALSE.
###: @param wait (\emph{bool}) Should R wait for run completion before giving back the hand to the user. Equals TRUE by default.
###: @examples
###: \dontrun{
###: runScenario() # sequential run
###: runScenario(wait = TRUE) # background run
###: }
###: @seealso \code{\link{setScenario}} \code{\link{getScenario}} \code{\link{abort}} \code{\link{isRunning}}
###: @export
runScenario = function(wait = TRUE){
  # on.exit(expr = self$abort())
  arguments = list()
  output = .processRequest("monolix", "runscenario", arguments, "synchronous", wait = wait, type = "STATUS")
  return(invisible(output))
}

###: Stop the current task run
###:
###: Stop the current task run.
###: @examples
###: \dontrun{
###: abort()
###: }
###: @seealso \code{\link{runScenario}}
###: @export
abort = function(){
  arguments = list()
  output = .processRequest("monolix", "abort", arguments, "asynchronous", type = "STATUS")
  return(invisible(output))
}

###: Get current scenario state
###:
###: Check if a scenario is currently running. If yes, information about the current running task are displayed.
###: @param verbose (\emph{bool}) Should information about the current running task be displayed in the console or not. Equals FALSE by default.
###: @return A boolean which equals TRUE if a scenario is currently running.
###: @examples
###: \dontrun{
###: isRunning()
###: }
###: @seealso \code{\link{runScenario}} \code{\link{abort}}
###: @export
isRunning = function(verbose = FALSE){
  arguments = list()
  output = .processRequest("monolix", "getworkflowstate", arguments, "asynchronous")
  
  tryCatch(expr = {
    if (verbose == TRUE){
      if (output[[1]] == TRUE){
        .info(paste0("Currently running task :\n",output[[2]]))
      } else {
        .info("No task currently running.")
      }
    }
    return(invisible(output[[1]]))
  },
  
  error = function(err){
    output = list()
    .error("Bad output encountered.")
  })
}

###: Get last run status
###:
###: Return an execution report about the last run with a summary of the error which could have occurred.
###: @return A structure containing 
###: \enumerate{
###: \item a boolean which equals TRUE if the last run has successfully completed, 
###: \item a summary of the errors which could have occurred.
###:}
###: @examples
###: \dontrun{
###: lastRunInfo = getLastRunStatus()
###: lastRunInfo$status
###:  -> TRUE
###: lastRunInfo$report
###:  -> ""
###:  }
###: @seealso \code{\link{runScenario}} \code{\link{abort}} \code{\link{isRunning}}
###: @export
getLastRunStatus = function(){
  arguments = list()
  output = .processRequest("monolix", "getlastrunstatus", arguments, "asynchronous")
  return(output)
}

###: Population parameter estimation
###:
###: Estimate the population parameters with the SAEM method. The associated method keyword is "saem".
###: By default, this task is not processed in the background of the R session. 
###: Notice that it does not impact the current scenario. Call 
###: \enumerate{
###: \item \code{\link{isRunning}} to check if the scenario is still running and get information about the current task,
###: \item \code{\link{abort}} to stop the execution.
###: }
###: To launch the function in the background, so that functions which do not modify the project ("get" functions for example) remains available, set the input argument "wait" to FALSE.
###: The initial values of the population parameters can be accessed by calling \code{\link{getPopulationParameterInformation}} and customized with \code{\link{setPopulationParameterInformation}}.\cr
###: The estimated population parameters are available using \code{\link{getEstimatedPopulationParameters}} function.
###: @param wait (\emph{bool}) Should R wait for run completion before giving back the hand to the user. Equals TRUE by default. 
###: @examples
###: \dontrun{
###: runPopulationParameterEstimation()
###: }
###: @seealso \code{\link{isRunning}} \code{\link{abort}}
###: @export
runPopulationParameterEstimation = function(wait = TRUE){
  if (is.logical(wait) == FALSE){
    .error("Unexpected type encountered for \"wait\" field. Please give a boolean.")
return(invisible(FALSE))
}
  arguments = list()
  output = .processRequest("monolix", "runsaem", arguments, "synchronous", wait = wait, type = "STATUS")
  return(invisible(output))
}

###: Estimation of the conditional modes (EBEs)
###:
###: Estimate the individual parameters using the conditional mode estimation algorithm (EBEs). The associated method keyword is "conditionalMode".
###: By default, this task is not processed in the background of the R session. 
###: Notice that it does not impact the current scenario. Call 
###: \enumerate{
###: \item \code{\link{isRunning}} to check if the scenario is still running and get information about the current task,
###: \item \code{\link{abort}} to stop the execution.
###: }
###: To launch the function in the background, so that functions which do not modify the project ("get" functions for example) remains available, set the input argument "wait" to FALSE.\cr
###: @param wait (\emph{bool}) Should R wait for run completion before giving back the hand to the user. Equals TRUE by default. 
###: @examples
###: \dontrun{
###: runConditionalModeEstimation()
###: }
###: @seealso \code{\link{isRunning}} \code{\link{abort}}
###: @export
runConditionalModeEstimation = function(wait = TRUE){
  if (is.logical(wait) == FALSE){
    .error("Unexpected type encountered for \"wait\" field. Please give a boolean.")
    return(invisible(FALSE))
  }
  
  arguments = list("conditionalMode")
  output = .processRequest("monolix", "runindivestim", arguments, "synchronous", wait = wait, type = "STATUS")
  return(invisible(output))
}

###: Sampling from the conditional distribution
###:
###: Estimate the individual parameters using conditional distribution sampling algorithm. The associated method keyword is "conditionalMean".
###: By default, this task is not processed in the background of the R session. 
###: Notice that it does not impact the current scenario. Call 
###: \enumerate{
###: \item \code{\link{isRunning}} to check if the scenario is still running and get information about the current task,
###: \item \code{\link{abort}} to stop the execution.
###: }
###: To launch the function in the background, so that functions which do not modify the project ("get" functions for example) remains available, set the input argument "wait" to FALSE.\cr
###: @param wait (\emph{bool}) Should R wait for run completion before giving back the hand to the user. Equals TRUE by default. 
###: @examples
###: \dontrun{
###: runConditionalDistributionSampling()
###: }
###: @seealso \code{\link{isRunning}} \code{\link{abort}}
###: @export
runConditionalDistributionSampling = function(wait = TRUE){
  if (is.logical(wait) == FALSE){
    .error("Unexpected type encountered for \"wait\" field. Please give a boolean.")
    return(invisible(FALSE))
  }
  
  arguments = list("conditionalMean")
  output = .processRequest("monolix", "runindivestim", arguments, "synchronous", wait = wait, type = "STATUS")
  return(invisible(output))
}

###: Standard error estimation
###:
###: Estimate the Fisher Information Matrix and the standard errors of the population parameters. By default, this task is not processed in the background of the R session. 
###: Notice that it does not impact the current scenario. Call 
###: \enumerate{
###: \item \code{\link{isRunning}} to check if the scenario is still running and get information about the current task,
###: \item \code{\link{abort}} to stop the execution.
###: }
###: To launch the function in the background, so that functions which do not modify the project ("get" functions for example) remains available, set the input argument "wait" to FALSE.\cr
###:
###: Existing methods:
###: \tabular{lll}{
###: \emph{Method} \tab \emph{Identifier} \cr
###: Estimate the FIM by Stochastic Approximation \tab linearization = F (default) \cr
###: Estimate the FIM by Linearization \tab linearization = T \cr
###: } 
###: The Fisher Information Matrix is available using \code{\link{getCorrelationOfEstimates}} function, while the standard errors are avalaible using \code{\link{getEstimatedStandardErrors}} function.
###: @param linearization option (\emph{boolean})[optional] method to be used. When no method is given, the stochastic approximation is used by default.
###: @param wait (\emph{bool}) Should R wait for run completion before giving back the hand to the user. Equals TRUE by default.
###: @examples
###: \dontrun{
###: runStandardErrorEstimation(linearization = T)
###: }
###: @seealso \code{\link{isRunning}} \code{\link{abort}}
###: @export
runStandardErrorEstimation = function(linearization=FALSE, wait = TRUE){
  if (is.logical(wait) == FALSE){
    .error("Unexpected type encountered for \"wait\" field. Please give a boolean.")
return(invisible(FALSE))
}
  # if (is.character(method) == FALSE){
  #   .error("Unexpected type encountered for \"method\" field. Please give a string corresponding to a fisher algorithm method.")
  #   return(invisible(FALSE))  
  # }
  # 
  # if (method != "")
  #   arguments = list(method)
  # else
  #   arguments = list()
  
  arguments <- ifelse(linearization==TRUE,list("linearization"),list("stochasticApproximation"))
  output = .processRequest("monolix", "runfisher", arguments, "synchronous", wait = wait, type = "STATUS")
  return(invisible(output))
}

###: Log-Likelihood estimation 
###:
###: Run the log-Likelihood estimation algorithm. By default, this task is not processed in the background of the R session. 
###: Notice that it does not impact the current scenario. Call 
###: \enumerate{
###: \item \code{\link{isRunning}} to check if the scenario is still running and get information about the current task,
###: \item \code{\link{abort}} to stop the execution.
###: }
###: To launch the function in the background, so that functions which do not modify the project ("get" functions for example) remains available, set the input argument "wait" to FALSE.\cr
###: Existing methods:
###: \tabular{lll}{
###: \emph{Method} \tab \emph{Identifier} \cr
###: Log-Likelihood estimation by linearization \tab linearization = T \cr
###: Log-Likelihood estimation by Importance Sampling (default) \tab linearization = F \cr
###: } 
###: The Log-likelihood outputs(-2LL, AIC, BIC) are available using \code{\link{getEstimatedLogLikelihood}} function
###: @param linearization option (\emph{boolean})[optional] method to be used. When no method is given, the importance sampling is used by default.
###: @param wait (\emph{bool}) Should R wait for run completion before giving back the hand to the user. Equals TRUE by default.
###: @examples
###: \dontrun{
###: runLogLikelihoodEstimation(linearization = T)
###: }
###: @seealso \code{\link{isRunning}} \code{\link{abort}}
###: @export
runLogLikelihoodEstimation = function(linearization=FALSE, wait = TRUE){
  # if (is.character(method) == FALSE){
  #   .error("Unexpected type encountered for \"method\" field. Please give a string corresponding to a log-likelihood algorithm method.")
  #   return(invisible(FALSE))  
  # }
  # 
  # arguments = list(method)
  arguments <- ifelse(linearization==TRUE,list("linearization"),list("importanceSampling"))
  output = .processRequest("monolix", "runll", arguments, "synchronous", wait = wait, type = "STATUS")
  return(invisible(output))
}
# -------------------------------------------------------------------------------- #
# ================================================================================ #