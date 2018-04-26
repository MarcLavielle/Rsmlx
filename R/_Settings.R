# ============================== SETTINGS MANAGEMENT ============================= #
.evalNumericalSettings = function(settingList){
  return( lapply(settingList, function(x) .evalNumerics(x)) )
}

# Project Settings --------------------------------------------------------------- #
###: Set project settings 
###: 
###: Set the value of one or several of the settings of the project. Associated settings are:
###: \tabular{lll}{
###: "directory" \tab (\emph{string}) \tab Path to the folder where simulation results will be saved. It should be a  writable directory.\cr
###: "exportResults" \tab (\emph{bool}) \tab Should results be exported.\cr
###: "seed" \tab (\emph{0< int <2147483647}) Seed used by random generators.\cr
###: "grid" \tab (\emph{int}) Number of points for the continuous simulation grid.\cr
###: "nbSimulations" \tab (\emph{int}) Simulation number.\cr
###: "dataAndModelNextToProject" \tab (\emph{bool}) Should data and model files be saved next to project.\cr
###: }
###: @param ... A collection of comma-separated pairs \{settingName = settingValue\}.
###: @examples
###: \dontrun{
###: setProjectSettings(directory = "/path/to/export/directory", seed = 12345)
###: }
###: @seealso \code{\link{getProjectSettings}}
###: @export
setProjectSettings = function(...){
  arguments = list(...)
  if (is.null(names(arguments)) && length(arguments) == 1)
    arguments <- arguments[[1]]
  settingNames = names(arguments)
  
  if ( length(settingNames) != length(arguments) ){
    .error("Unvalid input format. Please give a list of comma-separated pairs {settingName = settingValue}.")
    return(invisible(FALSE))
  }
  else if (length(settingNames) > 0){
    for (i in 1:length(settingNames)){
      if (.checkProjectSettingType(settingNames[i],arguments[[settingNames[i]]]) == FALSE){
        return(invisible(FALSE))
      }
    }
    output = .processRequest("monolix", "setprojectsettings", arguments, "synchronous", type = "STATUS")
    return(invisible(output))
  }
  else {
    return(invisible(TRUE))
  }
}

###: Get project settings
###: 
###: Get a summary of the project settings. Associated settings are:
###: \tabular{lll}{
###: "directory" \tab (\emph{string}) \tab Path to the folder where simulation results will be saved. It should be a writable directory.\cr
###: "exportResults" \tab (\emph{bool}) \tab Should results be exported.\cr
###: "seed" \tab (\emph{0< int <2147483647}) Seed used by random generators.\cr
###: "grid" \tab (\emph{int}) Number of points for the continuous simulation grid.\cr
###: "nbSimulations" \tab (\emph{int}) Simulation number.\cr
###: "dataAndModelNextToProject" \tab (\emph{bool}) Should data and model files be saved next to project.\cr
###: }
###: @param ... [optional] (string) Name of the settings whose value should be displayed. If no argument is provided, all the settings are returned.
###: @return An array which associates each setting name to its current value.
###: @examples
###: \dontrun{
###: getProjectSettings() # retrieve a list of all the project settings
###: getProjectSettings("directory","seed") # retrieve a list containing only the value of the settings whose name has been passed in argument
###: }
###: @seealso \code{\link{getProjectSettings}}
###: @export
getProjectSettings = function(...){
  arguments = list(...)
  
  if (length(arguments) > 0){
    for (i in 1:length(arguments)){
      if (is.character(arguments[[i]]) == FALSE){
        .error(paste0("Unexpected type encountered at position ",i,". Please give strings corresponding to valid settings names."))
        return(invisible(FALSE))    
      }
    }
  }
  
  output = .processRequest("monolix", "getprojectsettings", arguments, "asynchronous")
  output <- as.list(output)
  return(.evalNumericalSettings(output))
}

.checkProjectSettingType = function(settingName, settingValue){
  isValid = TRUE
  settingName = tolower(settingName)
  
  if (settingName == "directory"){
    if (is.character(settingValue) == FALSE){
      .error("Unexpected type encountered. Please give a string corresponding to the path to the wanted export directory.")
      isValid = FALSE
    }
  }
  else if (settingName == "exportresults"){
    if (is.logical(settingValue) == FALSE){
      .error("Unexpected type encountered. Please give a boolean equaling TRUE if results data should be exported and FALSE if not.")
      isValid = FALSE
    }
  }
  else if (settingName == "seed"){
    if (.isInteger(settingValue) == FALSE){
      .error("Unexpected type encountered. Please give a strictly positive integer lower or equal to 2147483646.")
      isValid = FALSE
    }
  }
  else if (settingName == "grid"){
    if (.isInteger(settingValue) == FALSE){
      .error("Unexpected type encountered. Please give a strictly positive integer.")
      isValid = FALSE
    }
  }
  else if (settingName == "nbsimulations"){
    if (.isInteger(settingValue) == FALSE){
      .error("Unexpected type encountered. Please give a strictly positive integer.")
      isValid = FALSE
    }
  }
  else if (settingName == "dataandmodelnexttoproject"){
    if (is.logical(settingValue) == FALSE){
      .error("Unexpected type encountered. Please give a boolean.")
      isValid = FALSE
    }
  }
  else {
    .error(paste0("\"",settingName,"\" is not a valid setting name."))
    isValid = FALSE
  }
  
  return(invisible(isValid))
}

# -------------------------------------------------------------------------------- #

# General Settings --------------------------------------------------------------- #
###: Set common settings for algorithms 
###: 
###: Set the value of one or several of the common settings for Monolix algorithms. Associated settings are:
###: \tabular{lll}{
###: "autoChains" \tab (\emph{bool}) \tab Automatically adjusted the number of chains to have at least a minimum number of subjects.\cr
###: "nbChains" \tab (\emph{int >0}) \tab Number of chains to be used if "autoChains" is set to FALSE.\cr
###: "minIndivForChains" \tab (\emph{int >0}) \tab Minimum number of individuals by chain.\cr
###: }
###: @param ... A collection of comma-separated pairs \{settingName = settingValue\}.
###: @examples
###: \dontrun{
###: setGeneralSettings(autoChains = FALSE, nbchains = 10)
###: }
###: @seealso \code{\link{getGeneralSettings}}
###: @export
setGeneralSettings = function(...){
  arguments = list(...)
  if (is.null(names(arguments)) && length(arguments) == 1)
    arguments <- arguments[[1]]
  settingNames = names(arguments)
  
  if ( length(settingNames) != length(arguments) ){
    .error("Unvalid input format. Please give a list of comma-separated pairs {settingName = settingValue}.")
    return(invisible(FALSE))
  }
  else if (length(settingNames) > 0){
    for (i in 1:length(settingNames)){
      if (.checkGeneralSettingType(settingNames[i],arguments[[settingNames[i]]]) == FALSE){
        return(invisible(FALSE))
      }
    }
    output = .processRequest("monolix", "setgeneralsettings", arguments, "synchronous", type = "STATUS")
    return(invisible(output))
  }
  else {
    return(invisible(TRUE))
  }
}

###: Get project general settings
###: 
###: Get a summary of the common settings for Monolix algorithms. Associated settings are:
###: \tabular{lll}{
###: "autoChains" \tab (\emph{bool}) \tab Automatically adjusted the number of chains to have at least a minimum number of subjects.\cr
###: "nbChains" \tab (\emph{int >0}) \tab Number of chains. \emph{Used only if "autoChains" is set to FALSE}.\cr
###: "minIndivForChains" \tab (\emph{int >0}) \tab Minimum number of individuals by chain.\cr
###: }
###: @param ... [optional] (string) Name of the settings whose value should be displayed. If no argument is provided, all the settings are returned.
###: @return An array which associates each setting name to its current value.
###: @examples
###: \dontrun{
###: getGeneralSettings() # retrieve a list of all the general settings
###: getGeneralSettings("nbChains","autoChains") # retrieve a list containing only the value of the settings whose name has been passed in argument
###: }
###: @seealso \code{\link{setGeneralSettings}}
###: @export
getGeneralSettings = function(...){
  arguments = list(...)
  
  if (length(arguments) > 0){
    for (i in 1:length(arguments)){
      if (is.character(arguments[[i]]) == FALSE){
        .error(paste0("Unexpected type encountered at position ",i,". Please give strings corresponding to valid settings names."))
        return(invisible(FALSE))    
      }
    }
  }
  
  output = .processRequest("monolix", "getgeneralsettings", arguments, "asynchronous")
  output <- as.list(output)
  return(.evalNumericalSettings(output))
}

.checkGeneralSettingType = function(settingName, settingValue){
  isValid = TRUE
  settingName = tolower(settingName)
  
  if (settingName == "autochains"){
    if (is.logical(settingValue) == FALSE){
      .error("Unexpected type encountered. \"autoChains\" must be a boolean.")
      isValid = FALSE
    }
  }
  else if (settingName == "nbchains" ){
    if (.isInteger(settingValue) == FALSE){
      .error("Unexpected type encountered. \"nbChains\" must be a strictly positive integer.")
      isValid = FALSE
    }
  }
  else if (settingName == "minindivforchains"){
    if (.isInteger(settingValue) == FALSE){
      .error("Unexpected type encountered. \"minIndivForChains\" must be a strictly positive integer.")
      isValid = FALSE
    }
  }
  else {
    .error(paste0("\"",settingName,"\" is not a valid setting name."))
    isValid = FALSE
  }
  
  return(invisible(isValid))
}
# -------------------------------------------------------------------------------- #

# MCMC Settings ------------------------------------------------------------------ #
###: Set settings associated to the MCMC algorithm
###: 
###: Set the value of one or several of the MCMC algorithm specific settings of the current project. Associated settings are:
###: \tabular{lll}{
###: "strategy" \tab (\emph{vector<int>[3]}) \tab Number of calls for each one of the three MCMC kernels.\cr
###: "acceptanceRatio" \tab (\emph{double}) \tab Target acceptance ratio.\cr
###: }
###: @param ... A collection of comma-separated pairs \{settingName = settingValue\}.
###: @examples
###: \dontrun{
###: setMCMCSettings(strategy = c(2,1,2))
###: }
###: @seealso \code{\link{getMCMCSettings}}
###: @export
setMCMCSettings = function(...){
  arguments = list(...)
  if (is.null(names(arguments)) && length(arguments) == 1)
    arguments <- arguments[[1]]
  settingNames = names(arguments)
  
  if ( length(settingNames) != length(arguments) ){
    .error("Unvalid input format. Please give a list of comma-separated pairs {settingName = settingValue}.")
    return(invisible(FALSE))
  }
  else if (length(settingNames) > 0){
    for (i in 1:length(settingNames)){
      if (.checkMCMCSettingType(settingNames[i],arguments[[settingNames[i]]]) == FALSE){
        return(invisible(FALSE))
      }
    }
    output = .processRequest("monolix", "setmcmcsettings", arguments, "synchronous", type = "STATUS")
    return(invisible(output))
  }
  else {
    return(invisible(TRUE))
  }
}

###: Get MCMC algorithm settings
###: 
###: Get the MCMC algorithm settings of the current project. Associated settings are:
###: \tabular{lll}{
###: "strategy" \tab (\emph{vector<int>[3]}) \tab Number of calls for each one of the three MCMC kernels.\cr
###: "acceptanceRatio" \tab (\emph{double}) \tab Target acceptance ratio.\cr
###: }
###: @param ... [optional] (string) Names of the settings whose value should be displayed. If no argument is provided, all the settings are returned.
###: @return An array which associates each setting name to its current value.
###: @examples
###: \dontrun{
###: getMCMCSettings() # retrieve a list of all the MCMC settings
###: getMCMCSettings("strategy") # retrieve a list containing only the value of the settings whose name has been passed in argument (here, the strategy)
###: }
###: @seealso \code{\link{setMCMCSettings}}
###: @export
getMCMCSettings = function(...){
  arguments = list(...)
  
  if (length(arguments) > 0){
    for (i in 1:length(arguments)){
      if (is.character(arguments[[i]]) == FALSE){
        .error(paste0("Unexpected type encountered at position ",i,". Please give strings corresponding to valid settings names."))
        return(invisible(FALSE))    
      }
    }
  }
  
  output = .processRequest("monolix", "getmcmcsettings", arguments, "asynchronous")
  output <- as.list(output)
  return(.evalNumericalSettings(output))
}

.checkMCMCSettingType = function(settingName, settingValue){
  settingName = tolower(settingName)
  
  if (settingName == "strategy"){
    if (length(settingValue) != 3){
      .error(paste0("Unexpected type encountered. MCMC strategy must be a vector of three integers corresponding to the number of calls for each one of the three MCMC kernels."))
      return(invisible(FALSE))
    }
    for (i in 1:length(settingValue)){
      if (.isInteger(settingValue[[i]]) == FALSE){
        .error(paste0("Unexpected type encountered. MCMC strategy must be a vector of three integers corresponding to the number of calls for each one of the three MCMC kernels."))
        return(invisible(FALSE))
      }
    }
  }
  else if (settingName == "acceptanceratio"){
    if (is.numeric(settingValue) == FALSE){
      .error(paste0("Unexpected type encountered. Target acceptance ratio must be a double between 0 and 1, exclusive."))
      return(invisible(FALSE))
    }
  }
  else {
    .error(paste0("\"",settingName,"\" is not a valid setting name."))
    isValid = FALSE
  }
  
  return(invisible(TRUE))
}
# -------------------------------------------------------------------------------- #

# Population Parameter Estimation Settings ------------------------------------------------------------------ #
###: Set population parameter estimation settings
###: 
###: Set the value of one or several of the population parameter estimation settings. Associated settings are:
###: \tabular{lll}{
###: "nbBurningIterations" \tab (\emph{int >=0}) \tab Number of iterations in the burn-in phase.\cr
###: "nbExploratoryIterations" \tab (\emph{int >=0}) \tab If "exploratoryAutoStop" is set to FALSE, it is the number of iterations in the exploratory phase. Else wise, if "exploratoryAutoStop" is set to TRUE, it is the maximum of iterations in the exploratory phase.\cr
###: "exploratoryAutoStop" \tab (\emph{bool}) \tab  Should the exploratory step automatically stop.\cr
###: "exploratoryInterval" \tab (\emph{int >0}) \tab Minimum number of interation in the exploratory phase. \emph{Used only if "exploratoryAutoStop" is TRUE}\cr
###: "exploratoryAlpha" \tab (\emph{0<= double <=1}) \tab Convergence memory in the exploratory phase. \emph{Used only if "exploratoryAutoStop" is TRUE}\cr
###: "nbSmoothingIterations" \tab (\emph{int >=0}) \tab   If "smoothingAutoStop" is set to FALSE, it is the number of iterations in the smoothing phase. Else wise, if "smoothingAutoStop" is set to TRUE, it is the maximum of iterations in the smoothing phase.\cr
###: "smoothingAutoStop" \tab (\emph{bool}) \tab   Should the smoothing step automatically stop.\cr
###: "smoothingInterval" \tab (\emph{int >0}) \tab Minimum number of interation in the smoothing phase. \emph{Used only if "smoothingAutoStop" is TRUE}.\cr
###: "smoothingAlpha" \tab (\emph{0.5< double <=1}) \tab Convergence memory in the smoothing phase. \emph{Used only if "smoothingAutoStop" is TRUE}.\cr
###: "smoothingRatio" \tab (\emph{0< double <1}) \tab Width of the confidence interval. \emph{Used only if "smoothingAutoStop" is TRUE}.\cr
###: "simulatedAnnealing" \tab (\emph{bool}) \tab  Should annealing be simulated.\cr
###: "tauOmega" \tab (\emph{double >0}) \tab Proportional rate on variance. \emph{Used only if "simulatedAnnealing" is TRUE}.\cr
###: "tauErrorModel" \tab (\emph{double >0}) \tab Proportional rate on error model. \emph{Used only if "simulatedAnnealing" is TRUE}. \cr
###: "variability" \tab (\emph{string}) \tab  Estimation method for parameters without variability: "firstStage" | "decreasing" | "none". \emph{Used only if arameters without variability are used in the project}.\cr
###: "nbOptimizationIterations" \tab (\emph{int >=1}) \tab Number of optimization iterations.\cr
###: "optimizationTolerance" \tab (\emph{double >0}) \tab Tolerance for optimization.\cr
###: }
###: @param ... A collection of comma-separated pairs \{settingName = SettingValue\}.
###: @examples
###: \dontrun{
###: setPopulationParameterEstimationSettings(exploratoryAutoStop = TRUE, tauOmega = 0.95)
###: }
###: @seealso \code{\link{getPopulationParameterEstimationSettings}}
###: @export
setPopulationParameterEstimationSettings = function(...){
  arguments = list(...)
  if (is.null(names(arguments)) && length(arguments) == 1)
    arguments <- arguments[[1]]
  settingNames = names(arguments)
  
  if ( length(settingNames) != length(arguments) ){
    .error("Unvalid input format. Please give a list of comma-separated pairs {settingName = settingValue}.")
    return(invisible(FALSE))
  }
  else if (length(settingNames) > 0){
    for (i in 1:length(settingNames)){
      if (.checkPopulationParameterEstimationSettingType(settingNames[i],arguments[[settingNames[i]]]) == FALSE){
        return(invisible(FALSE))
      }
    }
    output = .processRequest("monolix", "setpopulationparameterestimationsettings", arguments, "synchronous", type = "STATUS")
    return(invisible(output))
  }
  else {
    return(invisible(TRUE))
  }
}

###: Get population parameter estimation settings
###: 
###: Get the population parameter estimation settings. Associated settings are:
###: \tabular{lll}{
###: "nbBurningIterations" \tab (\emph{int >=0}) \tab Number of iterations in the burn-in phase.\cr
###: "nbExploratoryIterations" \tab (\emph{int >=0}) \tab If "exploratoryAutoStop" is set to FALSE, it is the number of iterations in the exploratory phase. Else wise, if "exploratoryAutoStop" is set to TRUE, it is the maximum of iterations in the exploratory phase.\cr
###: "exploratoryAutoStop" \tab (\emph{bool}) \tab  Should the exploratory step automatically stop.\cr
###: "exploratoryInterval" \tab (\emph{int >0}) \tab Minimum number of interation in the exploratory phase. \emph{Used only if "exploratoryAutoStop" is TRUE}\cr
###: "exploratoryAlpha" \tab (\emph{0<= double <=1}) \tab Convergence memory in the exploratory phase. \emph{Used only if "exploratoryAutoStop" is TRUE}\cr
###: "nbSmoothingIterations" \tab (\emph{int >=0}) \tab   If "smoothingAutoStop" is set to FALSE, it is the number of iterations in the smoothing phase. Else wise, if "smoothingAutoStop" is set to TRUE, it is the maximum of iterations in the smoothing phase.\cr
###: "smoothingAutoStop" \tab (\emph{bool}) \tab   Should the smoothing step automatically stop.\cr
###: "smoothingInterval" \tab (\emph{int >0}) \tab inimum number of interation in the smoothing phase. \emph{Used only if "smoothingAutoStop" is TRUE}.\cr
###: "smoothingAlpha" \tab (\emph{0.5< double <=1}) \tab Convergence memory in the smoothing phase. \emph{Used only if "smoothingAutoStop" is TRUE}.\cr
###: "smoothingRatio" \tab (\emph{0< double <1}) \tab Width of the confidence interval. \emph{Used only if "smoothingAutoStop" is TRUE}.\cr
###: "simulatedAnnealing" \tab (\emph{bool}) \tab  Should annealing be simulated.\cr
###: "tauOmega" \tab (\emph{double >0}) \tab Proportional rate on variance. \emph{Used only if "simulatedAnnealing" is TRUE}.\cr
###: "tauErrorModel" \tab (\emph{double >0}) \tab Proportional rate on error model. \emph{Used only if "simulatedAnnealing" is TRUE}. \cr
###: "variability" \tab (\emph{string}) \tab  Estimation method for parameters without variability: "firstStage" | "decreasing" | "none". \emph{Used only if arameters without variability are used in the project}.\cr
###: "nbOptimizationIterations" \tab (\emph{int >=1}) \tab Number of optimization iterations.\cr
###: "optimizationTolerance" \tab (\emph{double >0}) \tab Tolerance for optimization.\cr
###: }
###: @param ... [optional] (string) Name of the settings whose value should be displayed. If no argument is provided, all the settings are returned.
###: @return An array which associates each setting name to its current value.
###: @examples
###: \dontrun{
###: getPopulationParameterEstimationSettings() # retrieve a list of all the population parameter estimation settings
###: getPopulationParameterEstimationSettings("nbBurningIterations","smoothingInterval") # retrieve a list containing only the value of the settings whose name has been passed in argument (here, the number of burning iterations and the smoothing interval)
###: }
###: @seealso \code{\link{setPopulationParameterEstimationSettings}}
###: @export
getPopulationParameterEstimationSettings = function(...){
  arguments = list(...)
  
  if (length(arguments) > 0){
    for (i in 1:length(arguments)){
      if (is.character(arguments[[i]]) == FALSE){
        .error(paste0("Unexpected type encountered at position ",i,". Please give strings corresponding to valid settings names."))
        return(invisible(FALSE))    
      }
    }
  }
  
  output = .processRequest("monolix", "getpopulationparameterestimationsettings", arguments, "asynchronous")
  output <- as.list(output)
  return(.evalNumericalSettings(output))
}

.checkPopulationParameterEstimationSettingType = function(settingName, settingValue){
  isValid = TRUE
  settingName = tolower(settingName)
  
  if (settingName == "nbburningiterations"){
    if(.isInteger(settingValue) == FALSE){
      .error("Unexpected type encountered. The number of iterations in the burn-in phase must be an integer.")
      isValid = FALSE
    }
  }
  else if (settingName == "nbexploratoryiterations"){
    if(.isInteger(settingValue) == FALSE){
      .error("Unexpected type encountered. The number of iterations in the exploratory phase must be an integer.")
      isValid = FALSE
    }
  }
  else if (settingName == "exploratoryautostop"){
    if (is.logical(settingValue) == FALSE){
      .error("Unexpected type encountered. exploratoryautostop must be a boolean.")
      isValid = FALSE
    }
  }
  else if (settingName == "exploratoryinterval"){
    if (.isInteger(settingValue) == FALSE){
      .error("Unexpected type encountered. exploratoryinterval must be an integer.")
      isValid = FALSE
    }
  }
  else if (settingName == "exploratoryalpha"){
    if (is.numeric(settingValue) == FALSE){
      .error("Unexpected type encountered. exploratoryalpha must be a positive double.")
      isValid = FALSE
    }
  }
  else if (settingName == "nbsmoothingiterations"){
    if (.isInteger(settingValue) == FALSE){
      .error("Unexpected type encountered. The number of iterations in the smoothing phase must be an integer.")
      isValid = FALSE
    }
  }
  else if (settingName == "smoothingautostop"){
    if (is.logical(settingValue) == FALSE){
      .error("Unexpected type encountered. smoothingautostop must be a boolean.")
      isValid = FALSE
    }
  }
  else if (settingName == "smoothinginterval"){
    if (.isInteger(settingValue) == FALSE){
      .error("Unexpected type encountered. smoothinginterval must be an integer.")
      isValid = FALSE
    }
  }
  else if (settingName == "smoothingalpha"){
    if (is.numeric(settingValue) == FALSE){
      .error("Unexpected type encountered. smoothingalpha must be a double.")
      isValid = FALSE
    }
  }
  else if (settingName == "smoothingratio"){
    if (is.numeric(settingValue) == FALSE){
      .error("Unexpected type encountered. smoothingratio must be double between 0 and 1 exclusive.")
      isValid = FALSE
    }
  }
  else if (settingName == "simulatedannealing"){
    if (is.logical(settingValue) == FALSE){
      .error("Unexpected type encountered. simulatedannealing must be a boolean.")
      isValid = FALSE
    }
  }
  else if (settingName == "tauomega"){
    if (is.numeric(settingValue) == FALSE){
      .error("Unexpected type encountered. tauomega must be a strictly positive double.")
      isValid = FALSE
    }
  }
  else if (settingName == "tauerrormodel"){
    if (is.numeric(settingValue) == FALSE){
      .error("Unexpected type encountered. tauerrormodel must be a strictly positivedouble.")
      isValid = FALSE
    }
  }
  else if (settingName == "variability"){
    if (is.character(settingValue) == FALSE){
      .error("Unexpected type encountered. variability must be string corresponding to a valid variability strategy.")
      isValid = FALSE
    }
  }
  else if (settingName == "nboptimizationiterations"){
    if (.isInteger(settingValue) == FALSE){
      .error("Unexpected type encountered. nboptimizationiterations must be a strictly positive integer.")
      isValid = FALSE
    }
  }
  else if (settingName == "optimizationtolerance"){
    if (is.numeric(settingValue) == FALSE){
      .error("Unexpected type encountered. optimizationtolerance must be a strictly positive double.")
      isValid = FALSE
    }
  }
  else {
    .error(paste0("\"",settingName,"\" is not a valid setting name."))
    isValid = FALSE
  }
  
  return(invisible(isValid))
}
# --------------------------------------------------------------------------------------- #

# Conditional Mode Estimation Settings -------------------------------------------------- #
###: Set conditional mode estimation settings
###: 
###: Set the value of one or several of the conditional mode estimation settings. Associated settings are:
###: \tabular{lll}{
###: "nbOptimizationIterationsMode" \tab (\emph{int >=1}) \tab Maximum number of iterations. \cr
###: "optimizationToleranceMode" \tab (\emph{double >0}) \tab Optimization tolerance. \cr
###: }
###: @param ... A collection of comma-separated pairs \{settingName = settingValue\}.
###: @examples
###: \dontrun{
###: setConditionalModeEstimationSettings(nbOptimizationIterationsMode = 20, optimizationToleranceMode = 0.1)
###: }
###: @seealso \code{\link{getConditionalModeEstimationSettings}}
###: @export
setConditionalModeEstimationSettings = function(...){
  arguments = list(...)
  if (is.null(names(arguments)) && length(arguments) == 1)
    arguments <- arguments[[1]]
  settingNames = names(arguments)
  
  if ( length(settingNames) != length(arguments) ){
    .error("Unvalid input format. Please give a list of comma-separated pairs {settingName = settingValue}.")
    return(invisible(FALSE))
  }
  else if (length(settingNames) > 0){
    for (i in 1:length(settingNames)){
      if (.checkConditionalModeEstimationSettingType(settingNames[i],arguments[[settingNames[i]]]) == FALSE){
        return(invisible(FALSE))
      }
    }
    output = .processRequest("monolix", "setindivestimmodesettings", arguments, "synchronous", type = "STATUS")
    return(invisible(output))
  }
  else {
    return(invisible(TRUE))
  }
}

###: Get conditional mode estimation settings
###: 
###: Get the conditional mode estimation settings. Associated settings are:
###: \tabular{lll}{
###: "nbOptimizationIterationsMode" \tab (\emph{int >=1}) \tab Maximum number of iterations. \cr
###: "optimizationToleranceMode" \tab (\emph{double >0}) \tab Optimization tolerance. \cr
###: }
###: @param ... [optional] (string) Name of the settings whose value should be displayed. If no argument is provided, all the settings are returned.
###: @return An array which associates each setting name to its current value.
###: @examples
###: \dontrun{
###: getConditionalModeEstimationSettings() # retrieve a list of all the conditional mode estimation settings
###: getConditionalModeEstimationSettings("nbOptimizationIterationsMode") # retrieve a list containing only the value of the settings whose name has been passed in argument
###: }
###: @seealso \code{\link{setConditionalModeEstimationSettings}}
###: @export
getConditionalModeEstimationSettings = function(...){
  arguments = list(...)
  
  if (length(arguments) > 0){
    for (i in 1:length(arguments)){
      if (is.character(arguments[[i]]) == FALSE){
        .error(paste0("Unexpected type encountered at position ",i,". Please give strings corresponding to valid settings names."))
        return(invisible(FALSE))    
      }
    }
  }
  
  output = .processRequest("monolix", "getindivestimmodesettings", arguments, "asynchronous")
  output <- as.list(output)
  return(.evalNumericalSettings(output))
}


.checkConditionalModeEstimationSettingType = function(settingName, settingValue){
  isValid = TRUE
  settingName = tolower(settingName)
  
  if (settingName == "nboptimizationiterationsmode"){
    if (.isInteger(settingValue) == FALSE){
      .error(".Unexpected type encountered. nboptimizationiterationsmode must be an integer.")
      isValid = FALSE
    }
  }
  else if (settingName == "optimizationtolerancemode"){
    if (is.numeric(settingValue) == FALSE){
      .error("Unexpected type encountered. optimizationtolerancemode must be a double.")
      isValid = FALSE
    }
  }
  else {
    .error(paste0("\"",settingName,"\" is not a valid setting name."))
    isValid = FALSE
  }
  
  return(invisible(isValid))
}
# --------------------------------------------------------------------------------------- #

# Conditional Distribution Sampling Settings -------------------------------------------- #
###: Set conditional distribution sampling settings
###: 
###: Set the value of one or several of the conditional distribution sampling settings. Associated settings are:
###: \tabular{lll}{
###: "ratio" \tab (\emph{0< double <1}) \tab Width of the confidence interval.\cr
###: "nbMinIterations" \tab (\emph{int >=1}) \tab Minimum number of iterations.\cr
###: "nbSimulatedParameters" \tab (\emph{int >=1}) \tab Number of replicates.\cr
###: }
###: @param ... A collection of comma-separated pairs \{settingName = settingValue\}.
###: @examples
###: \dontrun{
###: setConditionalDistributionSamplingSettings(ratio = 0.05, nbMinIterations = 50)
###: }
###: @seealso \code{\link{getConditionalDistributionSamplingSettings}}
###: @export
setConditionalDistributionSamplingSettings = function(...){
  arguments = list(...)
  if (is.null(names(arguments)) && length(arguments) == 1)
    arguments <- arguments[[1]]
  settingNames = names(arguments)
  
  if ( length(settingNames) != length(arguments) ){
    .error("Unvalid input format. Please give a list of comma-separated pairs {settingName = settingValue}.")
    return(invisible(FALSE))
  }
  else if (length(settingNames) > 0){
    for (i in 1:length(settingNames)){
      if (.checkConditionalDistributionSamplingSettingType(settingNames[i],arguments[[settingNames[i]]]) == FALSE){
        return(invisible(FALSE))
      }
    }
    output = .processRequest("monolix", "setindivestimmeansettings", arguments, "synchronous", type = "STATUS")
    return(invisible(output))
  }
  else {
    return(invisible(TRUE))
  }
}

###: Get conditional distribution sampling settings
###: 
###: Get the conditional distribution sampling settings. Associated settings are:
###: \tabular{lll}{
###: "ratio" \tab (\emph{0< double <1}) \tab Width of the confidence interval.\cr
###: "nbMinIterations" \tab (\emph{int >=1}) \tab Minimum number of iterations.\cr
###: "nbSimulatedParameters" \tab (\emph{int >=1}) \tab Number of replicates.\cr
###: }
###: @param ... [optional] (string) Name of the settings whose value should be displayed. If no argument is provided, all the settings are returned.
###: @return An array which associates each setting name to its current value.
###: @examples
###: \dontrun{
###: getConditionalDistributionSamplingSettings() # retrieve a list of all the conditional distribution sampling settings
###: getConditionalDistributionSamplingSettings("ratio","nbMinIterations") # retrieve a list containing only the value of the settings whose name has been passed in argument
###: }
###: @seealso \code{\link{setConditionalDistributionSamplingSettings}}
###: @export
getConditionalDistributionSamplingSettings = function(...){
  arguments = list(...)
  
  if (length(arguments) > 0){
    for (i in 1:length(arguments)){
      if (is.character(arguments[[i]]) == FALSE){
        .error(paste0("Unexpected type encountered at position ",i,". Please give strings corresponding to valid settings names."))
        return(invisible(FALSE))    
      }
    }
  }
  
  output = .processRequest("monolix", "getindivestimmeansettings", arguments, "asynchronous")
  output <- as.list(output)
  return(.evalNumericalSettings(output))
}


.checkConditionalDistributionSamplingSettingType = function(settingName, settingValue){
  isValid = TRUE
  settingName = tolower(settingName)
  
  if (settingName == "ratio"){
    if (is.numeric(settingValue) == FALSE){
      .error("Unexpected type encountered. The ratio must be a double between 0 and 1, exclusive.")
      isValid = FALSE
    }
  }
  else if (settingName == "nbminiterations"){
    if (.isInteger(settingValue) == FALSE){
      .error(".Unexpected type encountered. The minimum number of iterations must be an integer.")
      isValid = FALSE
    }
  }
  else if (settingName == "nbsimulatedparameters"){
    if (.isInteger(settingValue) == FALSE){
      .error(".Unexpected type encountered. The number of simulated parameters must be an integer.")
      isValid = FALSE
    }
  }
  else {
    .error(paste0("\"",settingName,"\" is not a valid setting name."))
    isValid = FALSE
  }
  
  return(invisible(isValid))
}
# -------------------------------------------------------------------------------- #

# Standard Error Estimation Settings --------------------------------------------- #
###: Set standard error estimation settings
###: 
###: Set the value of one or several of the standard error estimation settings. Associated settings are:
###: \tabular{lll}{
###: "minIterations" \tab (\emph{int >=1}) \tab Minimum number of iterations.\cr
###: "maxIterations" \tab (\emph{int >=1}) \tab Maximum number of iterations.\cr
###: }
###: @param ... A collection of comma-separated pairs \{settingName = settingValue\}.
###: @examples
###: \dontrun{
###: setStandardErrorEstimationSettings(minIterations = 20, maxIterations = 250)
###: }
###: @seealso \code{\link{getStandardErrorEstimationSettings}}
###: @export
setStandardErrorEstimationSettings = function(...){
  arguments = list(...)
  if (is.null(names(arguments)) && length(arguments) == 1)
    arguments <- arguments[[1]]
  settingNames = names(arguments)
  
  if ( length(settingNames) != length(arguments) ){
    .error("Unvalid input format. Please give a list of comma-separated pairs {settingName = settingValue}.")
    return(invisible(FALSE))
  }
  else if (length(settingNames) > 0){
    for (i in 1:length(settingNames)){
      if (.checkStandardErrorEstimationSettingType(settingNames[i],arguments[[settingNames[i]]]) == FALSE){
        return(invisible(FALSE))
      }
    }
    output = .processRequest("monolix", "setstandarderrorestimationsettings", arguments, "synchronous", type = "STATUS")
    return(invisible(output))
  }
  else {
    return(invisible(TRUE))
  }
}

###: Get standard error estimation settings
###: 
###: Get the standard error estimation settings. Associated settings are:
###: \tabular{lll}{
###: "minIterations" \tab (\emph{int >=1}) \tab Minimum number of iterations. \cr
###: "maxIterations" \tab (\emph{int >=1}) \tab Maximum number of iterations. \cr
###: }
###: @param ... [optional] (string) Name of the settings whose value should be displayed. If no argument is provided, all the settings are returned.
###: @return An array which associates each setting name to its current value.
###: @examples
###: \dontrun{
###: getStandardErrorEstimationSettings() # retrieve a list of all the standard error estimation settings
###: getStandardErrorEstimationSettings("minIterations","maxIterations") # retrieve a list containing only the value of the settings whose name has been passed in argument
###: }
###: @seealso \code{\link{setStandardErrorEstimationSettings}}
###: @export
getStandardErrorEstimationSettings = function(...){
  arguments = list(...)
  
  if (length(arguments) > 0){
    for (i in 1:length(arguments)){
      if (is.character(arguments[[i]]) == FALSE){
        .error(paste0("Unexpected type encountered at position ",i,". Please give strings corresponding to valid settings names."))
        return(invisible(FALSE))    
      }
    }
  }
  
  output = .processRequest("monolix", "getstandarderrorestimationsettings", arguments, "asynchronous")
  output <- as.list(output)
  return(.evalNumericalSettings(output))
}


.checkStandardErrorEstimationSettingType = function(settingName, settingValue){
  isValid = TRUE
  settingName = tolower(settingName)
  
  if (settingName == "miniterations"){
    if (.isInteger(settingValue) == FALSE){
      .error("Unexpected type encountered. miniterations must be an integer.")
      isValid = FALSE
    }
  }
  else if (settingName == "maxiterations"){
    if (.isInteger(settingValue) == FALSE){
      .error("Unexpected type encountered. maxiterations must be an integer.")
      isValid = FALSE
    }
  }
  else {
    .error(paste0("\"",settingName,"\" is not a valid setting name."))
    isValid = FALSE
  }
  
  return(invisible(isValid))
}
# -------------------------------------------------------------------------------- #

# LogLikelihood Estimation Settings ---------------------------------------------- #
###: Set loglikelihood estimation settings
###: 
###: Set the value of the loglikelihood estimation settings. Associated settings are:
###: \tabular{lll}{
###: "nbFixedIterations" \tab (\emph{int >0}) \tab Monte Carlo size for the loglikelihood evaluation.\cr
###: "samplingMethod" \tab (\emph{string}) \tab Should the loglikelihood estimation use a given number of freedom degrees ("fixed") or test a sequence of degrees of freedom numbers before choosing the best one ("optimized").\cr
###: "nbFreedomDegrees" \tab (\emph{int >0}) \tab Degree of freedom of the Student t-distribution. \emph{Used only if "samplingMethod" is "fixed"}.\cr
###: "freedomDegreesSampling" \tab (\emph{vector<int(>0)>}) \tab Sequence of freedom degrees to be tested. \emph{Used only if "samplingMethod" is "optimized"}.\cr
###: }
###: @param ... A collection of comma-separated pairs \{settingName = settingValue\}.
###: @examples
###: \dontrun{
###: setLogLikelihoodEstimationSettings(nbFixedIterations = 20000)
###: }
###: @seealso \code{\link{getLogLikelihoodEstimationSettings}}
###: @export
setLogLikelihoodEstimationSettings = function(...){
  arguments = list(...)
  if (is.null(names(arguments)) && length(arguments) == 1)
    arguments <- arguments[[1]]
  settingNames = names(arguments)
  
  if ( length(settingNames) != length(arguments) ){
    .error("Unvalid input format. Please give a list of comma-separated pairs {settingName = settingValue}.")
    return(invisible(FALSE))
  }
  else if (length(settingNames) > 0){
    for (i in 1:length(settingNames)){
      if (.checkLogLikelihoodEstimationSettingType(settingNames[i],arguments[[settingNames[i]]]) == FALSE){
        return(invisible(FALSE))
      }
    }
    output = .processRequest("monolix", "setloglikelihoodestimationsettings", arguments, "synchronous", type = "STATUS")
    return(invisible(output))
  }
  else {
    return(invisible(TRUE))
  }
}

###: Get LogLikelihood algorithm settings
###: 
###: Get the loglikelihood estimation settings. Associated settings are:
###: \tabular{lll}{
###: "nbFixedIterations" \tab (\emph{int >0}) \tab Monte Carlo size for the loglikelihood evaluation.\cr
###: "samplingMethod" \tab (\emph{string}) \tab Should the loglikelihood estimation use a given number of freedom degrees ("fixed") or test a sequence of degrees of freedom numbers before choosing the best one ("optimized").\cr
###: "nbFreedomDegrees" \tab (\emph{int >0}) \tab Degree of freedom of the Student t-distribution. \emph{Used only if "samplingMethod" is "fixed"}.\cr
###: "freedomDegreesSampling" \tab (\emph{vector<int(>0)>}) \tab Sequence of freedom degrees to be tested. \emph{Used only if "samplingMethod" is "optimized"}.\cr
###: }
###: @param ... [optional] (string) Name of the settings whose value should be displayed. If no argument is provided, all the settings are returned.
###: @return An array which associates each setting name to its current value.
###: @examples
###: \dontrun{
###: getLogLikelihoodEstimationSettings() # retrieve a list of all the loglikelihood estimation settings
###: getLogLikelihoodEstimationSettings("nbFixedIterations","samplingMethod") # retrieve a list containing only the value of the settings whose name has been passed in argument (here, the number of fixed iterations and the method)
###: }
###: @seealso \code{\link{setLogLikelihoodEstimationSettings}}
###: @export
getLogLikelihoodEstimationSettings = function(...){
  arguments = list(...)
  
  if (length(arguments) > 0){
    for (i in 1:length(arguments)){
      if (is.character(arguments[[i]]) == FALSE){
        .error(paste0("Unexpected type encountered at position ",i,". Please give strings corresponding to valid settings names."))
        return(invisible(FALSE))    
      }
    }
  }
  
  output = .processRequest("monolix", "getloglikelihoodestimationsettings", arguments, "asynchronous")
  output <- as.list(output)
  return(.evalNumericalSettings(output))
}

.checkLogLikelihoodEstimationSettingType = function(settingName, settingValue){
  isValid = TRUE
  settingName = tolower(settingName)
  
  if (settingName == "nbfixediterations"){
    if (.isInteger(settingValue) == FALSE){
      .error("Unexpected type encountered. nbfixediterations must be an integer.")
      isValid = FALSE
    }
  }
  else if (settingName == "samplingmethod"){
    if (is.character(settingValue) == FALSE){
      .error("Unexpected type encountered. samplingmethod must be a string corresponding to a valid method for freedom degrees sampling (fixed or optimized).")
      isValid = FALSE
    }
  }
  else if (settingName == "nbfreedomdegrees"){
    if (.isInteger(settingValue) == FALSE){
      .error("Unexpected type encountered. nbfreedomdegrees must be an integer.")
      isValid = FALSE
    }
  }
  else if (settingName == "freedomdegreessampling"){
    for (i in 1:length(settingValue)){
      if (.isInteger(settingValue[[i]]) == FALSE){
        .error(paste0("Unexpected type encountered. LL freedom degrees sampling must be a vector of integers."))
        isValid = FALSE
        break  
      }
    }
  }
  else {
    .error(paste0("\"",settingName,"\" is not a valid setting name."))
    isValid = FALSE
  }
  
  return(invisible(isValid))
}
# -------------------------------------------------------------------------------- #

# Preferences -------------------------------------------------------------------- #
###: Set preferences
###: 
###: Set the value of one or several of the project preferences. Prefenreces are:
###: \tabular{lll}{
###: "relativePath" \tab (\emph{bool}) \tab Use relative path for save/load operations.\cr
###: "threads" \tab (\emph{int >0}) \tab Number of threads.\cr
###: "timeStamping" \tab (\emph{bool}) Create an archive containing result files after each run.\cr
###: "dpi" \tab (\emph{bool}) Apply high density pixel correction.\cr
###: "imageFormat" \tab (\emph{string}) Image format used to save monolix graphics.\cr
###: "delimiter" \tab (\emph{string}) Character use as delimiter in exported result files.\cr
###: "exportGraphics" \tab (\emph{bool}) Should graphics images be exported.\cr
###: "exportGraphicsData" \tab (\emph{bool}) Should graphics data be exported.\cr
###: }
###: @param ... A collection of comma-separated pairs \{preferenceName = settingValue\}.
###: @examples
###: \dontrun{
###: setPreferences(exportGraphics = FALSE, delimiter = ",")
###: }
###: @seealso \code{\link{getPreferences}}
###: @export
setPreferences = function(...){
  arguments = list(...)
  if (is.null(names(arguments)) && length(arguments) == 1)
    arguments <- arguments[[1]]
  settingNames = names(arguments)
  
  if ( length(settingNames) != length(arguments) ){
    .error("Unvalid input format. Please give a list of comma-separated pairs {preferenceName = preferenceValue}.")
    return(invisible(FALSE))
  }
  else if (length(settingNames) > 0){
    for (i in 1:length(settingNames)){
      if (.checkPreferenceType(settingNames[i],arguments[[settingNames[i]]]) == FALSE){
        return(invisible(FALSE))
      }
    }
    output = .processRequest("monolix", "setpreferences", arguments, "synchronous", type = "STATUS")
    return(invisible(output))
  }
  else {
    return(invisible(TRUE))
  }
}

###: Get project preferences
###: 
###: Get a summary of the project preferences. Preferences are:
###: \tabular{lll}{
###: "relativePath" \tab (\emph{bool}) \tab Use relative path for save/load operations.\cr
###: "threads" \tab (\emph{int >0}) \tab Number of threads.\cr
###: "timeStamping" \tab (\emph{bool}) Create an archive containing result files after each run.\cr
###: "dpi" \tab (\emph{bool}) Apply high density pixel correction.\cr
###: "imageFormat" \tab (\emph{string}) Image format used to save monolix graphics.\cr
###: "delimiter" \tab (\emph{string}) Character use as delimiter in exported result files.\cr
###: "exportGraphics" \tab (\emph{bool}) Should graphics images be exported.\cr
###: "exportGraphicsData" \tab (\emph{bool}) Should graphics data be exported.\cr
###: }
###: @param ... [optional] (string) Name of the preference whose value should be displayed. If no argument is provided, all the preferences are returned.
###: @return An array which associates each preference name to its current value.
###: @examples
###: \dontrun{
###: getPreferences() # retrieve a list of all the general settings
###: getPreferences("imageFormat","exportGraphics") # retrieve a list containing only the value of the preferences whose name has been passed in argument
###: }
###: @seealso \code{\link{setGeneralSettings}}
###: @export
getPreferences = function(...){
  arguments = list(...)
  
  if (length(arguments) > 0){
    for (i in 1:length(arguments)){
      if (is.character(arguments[[i]]) == FALSE){
        .error(paste0("Unexpected type encountered at position ",i,". Please give strings corresponding to valid preference names."))
        return(invisible(FALSE))    
      }
    }
  }
  
  output = .processRequest("monolix", "getpreferences", arguments, "asynchronous")
  output <- as.list(output)
  return(.evalNumericalSettings(output))
}

.checkPreferenceType = function(settingName, settingValue){
  isValid = TRUE
  settingName = tolower(settingName)
  
  if (settingName == "relativepath"){
    if (is.logical(settingValue) == FALSE){
      .error("Unexpected type encountered. \"relativePath\" must be a boolean.")
      isValid = FALSE
    }
  }
  else if (settingName == "threads" ){
    if (.isInteger(settingValue) == FALSE){
      .error("Unexpected type encountered. \"threads\" must be a strictly positive integer.")
      isValid = FALSE
    }
  }
  else if (settingName == "timestamping"){
    if (is.logical(settingValue) == FALSE){
      .error("Unexpected type encountered. \"timestamping\" must be a boolean.")
      isValid = FALSE
    }
  }
  else if (settingName == "dpi"){
    if (is.logical(settingValue) == FALSE){
      .error("Unexpected type encountered. \"dpi\" must be a boolean.")
      isValid = FALSE
    }
  }
  else if (settingName == "imageformat"){
    if (is.character(settingValue) == FALSE){
      .error("Unexpected type encountered. \"imageformat\" must be a string.")
      isValid = FALSE
    }
  }
  else if (settingName == "delimiter"){
    if (is.character(settingValue) == FALSE){
      .error("Unexpected type encountered. \"delimiter\" must be a string.")
      isValid = FALSE
    }
  }
  else if (settingName == "exportgraphics"){
    if (is.logical(settingValue) == FALSE){
      .error("Unexpected type encountered. \"exportgraphics\" must be a boolean.")
      isValid = FALSE
    }
  }
  else if (settingName == "exportgraphicsdata"){
    if (is.logical(settingValue) == FALSE){
      .error("Unexpected type encountered. \"exportgraphicsdata\" must be a boolean.")
      isValid = FALSE
    }
  }else {
    .error(paste0("\"",settingName,"\" is not a valid preference name."))
    isValid = FALSE
  }
  
  return(invisible(isValid))
}
# -------------------------------------------------------------------------------- #
# ================================================================================ #