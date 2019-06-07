.hiddenCall <- function(command){
  eval.parent(parse(text = command))
}

mlx.getLixoftConnectorsState <- function(quietly = TRUE) {
  r <- NULL
  .hiddenCall(paste0('r <- lixoftConnectors::getLixoftConnectorsState(quietly = ',quietly,')'))
  return(r)
}

mlx.initializeLixoftConnectors <- function() {
  r <- NULL
  .hiddenCall(paste0('r <- lixoftConnectors::initializeLixoftConnectors(software = "monolix", force = TRUE)'))
  return(r)
}

mlx.setStructuralModel <- function(modelFile = NULL) {
  .hiddenCall(paste0('r <- lixoftConnectors::setStructuralModel(modelFile = ',modelFile,')'))
}

mlx.getPopulationParameterInformation <- function() {
  r <- NULL
  .hiddenCall(paste0('r <- lixoftConnectors::getPopulationParameterInformation()'))
  return(r)
}

mlx.getScenario <- function() {
  r <- NULL
  .hiddenCall(paste0('r <- lixoftConnectors::getScenario()'))
  return(r)
}


mlx.getObservationInformation <- function() {
  r <- NULL
  .hiddenCall(paste0('r <- lixoftConnectors::getObservationInformation()'))
  return(r)
}

mlx.getLaunchedTasks <- function() {
  r <- NULL
  .hiddenCall(paste0('r <- lixoftConnectors::getLaunchedTasks()'))
  return(r)
}

mlx.getEstimatedLogLikelihood <- function() {
  r <- NULL
  .hiddenCall(paste0('r <- lixoftConnectors::getEstimatedLogLikelihood()'))
  return(r)
}

mlx.getEstimatedPopulationParameters <- function() {
  r <- NULL
  .hiddenCall(paste0('r <- lixoftConnectors::getEstimatedPopulationParameters()'))
  return(r)
}

mlx.getConditionalDistributionSamplingSettings <- function() {
  r <- NULL
  .hiddenCall(paste0('r <- lixoftConnectors::getConditionalDistributionSamplingSettings()'))
  return(r)
}

mlx.getConditionalModeEstimationSettings <- function() {
  r <- NULL
  .hiddenCall(paste0('r <- lixoftConnectors::getConditionalModeEstimationSettings()'))
  return(r)
}

mlx.getContinuousObservationModel <- function() {
  r <- NULL
  .hiddenCall(paste0('r <- lixoftConnectors::getContinuousObservationModel()'))
  return(r)
}

mlx.getCovariateInformation <- function() {
  r <- NULL
  .hiddenCall(paste0('r <- lixoftConnectors::getCovariateInformation()'))
  return(r)
}

mlx.getData <- function() {
  r <- NULL
  .hiddenCall(paste0('r <- lixoftConnectors::getData()'))
  return(r)
}
mlx.getEstimatedIndividualParameters <- function() {
  r <- NULL
  .hiddenCall(paste0('r <- lixoftConnectors::getEstimatedIndividualParameters()'))
  return(r)
}
mlx.getEstimatedRandomEffects <- function() {
  r <- NULL
  .hiddenCall(paste0('r <- lixoftConnectors::getEstimatedRandomEffects()'))
  return(r)
}
mlx.getEstimatedStandardErrors <- function() {
  r <- NULL
  .hiddenCall(paste0('r <- lixoftConnectors::getEstimatedStandardErrors()'))
  return(r)
}
mlx.getGeneralSettings <- function() {
  r <- NULL
  .hiddenCall(paste0('r <- lixoftConnectors::getGeneralSettings()'))
  return(r)
}
mlx.getIndividualParameterModel <- function() {
  r <- NULL
  .hiddenCall(paste0('r <- lixoftConnectors::getIndividualParameterModel()'))
  return(r)
}
mlx.getLogLikelihoodEstimationSettings <- function() {
  r <- NULL
  .hiddenCall(paste0('r <- lixoftConnectors::getLogLikelihoodEstimationSettings()'))
  return(r)
}
mlx.getPopulationParameterEstimationSettings <- function() {
  r <- NULL
  .hiddenCall(paste0('r <- lixoftConnectors::getPopulationParameterEstimationSettings()'))
  return(r)
}
mlx.getProjectSettings <- function() {
  r <- NULL
  .hiddenCall(paste0('r <- lixoftConnectors::getProjectSettings()'))
  return(r)
}
mlx.getSimulatedIndividualParameters <- function() {
  r <- NULL
  .hiddenCall(paste0('r <- lixoftConnectors::getSimulatedIndividualParameters()'))
  return(r)
}
mlx.getSimulatedRandomEffects <- function() {
  r <- NULL
  .hiddenCall(paste0('r <- lixoftConnectors::getSimulatedRandomEffects()'))
  return(r)
}
mlx.getStandardErrorEstimationSettings <- function() {
  r <- NULL
  .hiddenCall(paste0('r <- lixoftConnectors::getStandardErrorEstimationSettings()'))
  return(r)
}
mlx.runConditionalDistributionSampling <- function() {
  .hiddenCall(paste0('r <- lixoftConnectors::runConditionalDistributionSampling()'))
}
mlx.runConditionalModeEstimation <- function() {
  .hiddenCall(paste0('r <- lixoftConnectors::runConditionalModeEstimation()'))
}
mlx.runStandardErrorEstimation <- function(linearization=NULL) {
  .hiddenCall(paste0('r <- lixoftConnectors::runStandardErrorEstimation(linearization = ',linearization,')'))
}
mlx.runScenario <- function(wait=TRUE) {
  .hiddenCall(paste0('r <- lixoftConnectors::runScenario(wait = ',wait,')'))
}
mlx.setInitialEstimatesToLastEstimates <- function() {
  .hiddenCall(paste0('r <- lixoftConnectors::setInitialEstimatesToLastEstimates()'))
}

mlx.setPopulationParameterInformation <- function(a) {
  .hiddenCall(paste0('r <- lixoftConnectors::setPopulationParameterInformation(a)'))
}

mlx.loadProject <- function(projectFile=NULL) {
  .hiddenCall(paste0('r <- lixoftConnectors::loadProject(projectFile = "',projectFile,'")'))
}

mlx.setScenario <- function(a) {
  .hiddenCall(paste0('r <- lixoftConnectors::setScenario(a)'))
}
mlx.setConditionalDistributionSamplingSettings <- function(a) {
  .hiddenCall(paste0('r <- lixoftConnectors::setConditionalDistributionSamplingSettings(a)'))
}
mlx.setConditionalModeEstimationSettings <- function(a) {
  .hiddenCall(paste0('r <- lixoftConnectors::setConditionalModeEstimationSettings(a)'))
}
mlx.computePredictions <- function(a) {
  .hiddenCall(paste0('r <- lixoftConnectors::computePredictions(a)'))
}
mlx.setCovariateModel <- function(a) {
  .hiddenCall(paste0('r <- lixoftConnectors::setCovariateModel(a)'))
}


mlx.setCorrelationBlocks <- function(a) {
  .hiddenCall(paste0('r <- lixoftConnectors::setCorrelationBlocks(a)'))
}
mlx.setData <- function(a) {
  .hiddenCall(paste0('r <- lixoftConnectors::setData(a)'))
}
mlx.setErrorModel <- function(a) {
  .hiddenCall(paste0('r <- lixoftConnectors::setErrorModel(a)'))
}
mlx.setGeneralSettings <- function(g) {
  .hiddenCall(paste0('r <- lixoftConnectors::setGeneralSettings(g)'))
}
mlx.setLogLikelihoodEstimationSettings <- function(a) {
  .hiddenCall(paste0('r <- lixoftConnectors::setLogLikelihoodEstimationSettings(a)'))
}
mlx.setPopulationParameterEstimationSettings <- function(g) {
  .hiddenCall(paste0('r <- lixoftConnectors::setPopulationParameterEstimationSettings(g)'))
}
mlx.newProject <- function(data = NULL, modelFile = NULL) {
  .hiddenCall(paste0('r <- lixoftConnectors::newProject(data = data, modelFile = modelFile)'))
}

mlx.setProjectSettings <- function(directory=NULL) {
  .hiddenCall(paste0('r <- lixoftConnectors::setProjectSettings(directory = directory)'))
}
mlx.setStandardErrorEstimationSettings  <- function(a) {
  .hiddenCall(paste0('r <- lixoftConnectors::setStandardErrorEstimationSettings (a)'))
}
mlx.setObservationDistribution  <- function(a) {
  .hiddenCall(paste0('r <- lixoftConnectors::setObservationDistribution (a)'))
}


mlx.saveProject <- function(projectFile=NULL) {
  .hiddenCall(paste0('r <- lixoftConnectors::saveProject(projectFile = projectFile)'))
}

mlx.runPopulationParameterEstimation <- function() {
  .hiddenCall(paste0('r <- lixoftConnectors::runPopulationParameterEstimation()'))
}

mlx.runLogLikelihoodEstimation <- function(linearization = NULL) {
  .hiddenCall(paste0('r <- lixoftConnectors::runLogLikelihoodEstimation(linearization = linearization)'))
}




  