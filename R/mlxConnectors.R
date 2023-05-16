.hiddenCall <- function(command){
  eval.parent(parse(text = command))
}

mlx.getLixoftConnectorsState <- function(quietly = TRUE) {
  r <- NULL
  .hiddenCall(paste0('r <- lixoftConnectors::getLixoftConnectorsState(quietly = ',quietly,')'))
  return(r)
}

mlx.initializeLixoftConnectors <- function(software = "monolix", path="", force = TRUE) {
  r <- NULL
  .hiddenCall(paste0('r <- lixoftConnectors::initializeLixoftConnectors(software = "',software ,'", 
                     path = "',path,'", force=',force,')'))
  return(invisible(r))
}

mlx.setStructuralModel <- function(modelFile = NULL) {
  .hiddenCall(paste0('r <- lixoftConnectors::setStructuralModel(modelFile = "',modelFile,'")'))
}

mlx.getStructuralModel <- function() {
  r <- NULL
  .hiddenCall(paste0('r <- lixoftConnectors::getStructuralModel()'))
  return(r)
}

mlx.getTests <- function() {
  r <- NULL
  .hiddenCall(paste0('r <- lixoftConnectors::getTests()'))
  return(r)
}

mlx.getSAEMiterations <- function() {
  r <- NULL
  .hiddenCall(paste0('r <- lixoftConnectors::getSAEMiterations()'))
  return(r)
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
  for (k in 1:length(r)) {
    if (is.list(r[[k]]))
      r[[k]] <- unlist(r[[k]])
    if (!is.null(r[[k]]['-2LL'])) 
      names(r[[k]]) <- gsub("-2LL", "OFV", names(r[[k]]))
    i0 <- which(names(r[[k]])=='chosenDegree')
    if (length(i0)>0)
      r[[k]] <- r[[k]][-i0]
  }
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
  sn <- setdiff(r$name,names(r$covariate))
  if (length(sn)>0) {
    d <- mlx.getProjectSettings()$directory
    pind <- read.csv(file.path(d,"IndividualParameters/estimatedIndividualParameters.txt"))
    if (all(sn %in% names(pind)))
      r$covariate <- merge(r$covariate,pind[,c("id",sn)],by="id")
  }
  j.strat <- grep("stratification",r$type)
  if (length(j.strat) > 0) {
    strat.cov <- r$name[j.strat]
    r$covariate <- r$covariate %>% select(-strat.cov)
    r$type <- r$type[-j.strat] 
    r$name <- r$name[-j.strat] 
  }
  
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
  # for (k in 1:length(r)) {
  #   if (is.list(r[[k]]))
  #     r[[k]] <- unlist(r[[k]])
  # }
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
  if (is.factor(r$rep))  r$rep <- as.numeric(as.character(r$rep))
  return(r)
}
mlx.getSimulatedRandomEffects <- function() {
  r <- NULL
  .hiddenCall(paste0('r <- lixoftConnectors::getSimulatedRandomEffects()'))
  if (is.factor(r$rep))  r$rep <- as.numeric(as.character(r$rep))
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
mlx.runScenario <- function() {
  .hiddenCall('r <- lixoftConnectors::runScenario()')
}
mlx.setInitialEstimatesToLastEstimates <- function(fixedEffectsOnly = F) {
  .hiddenCall(paste0('r <- lixoftConnectors::setInitialEstimatesToLastEstimates(fixedEffectsOnly=fixedEffectsOnly)'))
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
mlx.setIndividualParameterModel <- function(a) {
  .hiddenCall(paste0('lixoftConnectors::setIndividualParameterModel(a)'))
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

mlx.setProjectSettings <- function(directory = NULL, dataandmodelnexttoproject = NULL) {
  if (!is.null(directory)) {
    .hiddenCall(paste0('r <- lixoftConnectors::setProjectSettings(directory = directory)'))
  } else {
    if (!is.null(dataandmodelnexttoproject)) {
      .hiddenCall(paste0('r <- lixoftConnectors::setProjectSettings(dataandmodelnexttoproject = dataandmodelnexttoproject)'))
    }
  }
}
mlx.setStandardErrorEstimationSettings  <- function(a) {
  .hiddenCall(paste0('r <- lixoftConnectors::setStandardErrorEstimationSettings (a)'))
}
mlx.setObservationDistribution  <- function(a) {
  .hiddenCall(paste0('r <- lixoftConnectors::setObservationDistribution (a)'))
}


mlx.saveProject <- function(projectFile=NULL) {
  if (is.null(projectFile)) {
    .hiddenCall(paste0('r <- lixoftConnectors::saveProject()'))
  } else {
    .hiddenCall(paste0('r <- lixoftConnectors::saveProject(projectFile = projectFile)'))
  }
}

mlx.runPopulationParameterEstimation <- function(parameters=NULL) {
  r <- NULL
  if (as.numeric(substr(packageVersion("lixoftConnectors"),1,4))>=2021 & !is.null(parameters))
    .hiddenCall(paste0('r <- lixoftConnectors::runPopulationParameterEstimation(parameters=parameters)'))
  else
    .hiddenCall(paste0('r <- lixoftConnectors::runPopulationParameterEstimation()'))
  .hiddenCall(paste0('r0 <- lixoftConnectors::runConditionalModeEstimation()'))
  return(r)
}

mlx.runLogLikelihoodEstimation <- function(linearization = FALSE) {
  .hiddenCall(paste0('r <- lixoftConnectors::runLogLikelihoodEstimation(linearization = linearization)'))
}

mlx.getLibraryModelName <- function(library) {
  .hiddenCall(paste0('r <- lixoftConnectors::getLibraryModelName(library)'))
}

mlx.saveFormattedFile <- function(path) {
  .hiddenCall(paste0('r <- lixoftConnectors::getFormatting()'))
  .hiddenCall(paste0('r$formattedFile <- path'))
  .hiddenCall(paste0('do.call(formatData, r)'))
}

mlx.getFormatting <- function() {
  r <- NULL
  .hiddenCall(paste0('r <- lixoftConnectors::getFormatting()'))
  return(r)
}

smlx.exportSimulatedData <- function(path) {
  .hiddenCall(paste0('r <- lixoftConnectors::exportSimulatedData(path)'))
}


smlx.importMonolixProject <- function(project) {
  .hiddenCall(paste0('r <- lixoftConnectors::importProject(project)'))
}

smlx.setNbReplicates <- function(nrep) {
  .hiddenCall(paste0('r <- lixoftConnectors::setNbReplicates(nrep)'))
}

smlx.getNbReplicates <- function() {
  .hiddenCall(paste0('r <- lixoftConnectors::getNbReplicates()'))
}

smlx.runSimulation <- function() {
  .hiddenCall(paste0('r <- lixoftConnectors::runSimulation()'))
}

smlx.getProjectSettings <- function() {
  .hiddenCall(paste0('r <- lixoftConnectors::getProjectSettings()'))
}

smlx.getSimulationResults <- function() {
  .hiddenCall(paste0('r <- lixoftConnectors::getSimulationResults()'))
}

smlx.getGroups <- function() {
  .hiddenCall(paste0('r <- lixoftConnectors::getGroups()'))
}

smlx.getCovariateElements <- function() {
  .hiddenCall(paste0('r <- lixoftConnectors::getCovariateElements()'))
}

smlx.getOccasionElements <- function() {
  .hiddenCall(paste0('r <- lixoftConnectors::getOccasionElements()'))
}

smlx.setProjectSettings <- function(...) {
  .hiddenCall(paste0('r <- lixoftConnectors::setProjectSettings(...)'))
}

smlx.getOutputElements <- function(...) {
  .hiddenCall(paste0('r <- lixoftConnectors::getOutputElements(...)'))
}

smlx.getTreatmentElements <- function(...) {
  .hiddenCall(paste0('r <- lixoftConnectors::getTreatmentElements(...)'))
}
