#' Covariate search functions
#' 
#' Search the associated covariate model to best fit the data set. Several type of covariate search is proposed
#' \itemize{
#' \item COSSAC: ranked SCM
#' \item SCM: SCM method
#' }
#' @param project a Monolix project
#' @param covariateToSearch [optional] defines the covariate to search. By default, all covariates are searched.
#' @param parameterToTest [optional] defines the parameters on which the covariates are applied. By default, all parameters are used.
#' @param settings [optional] defines the settings search. This is a list of all settings for the search. 
#' The settings are pInclusion, pElimination, useLinearization, rankedSCM, criteria, updateInit. By default, the values are .1, .05, FALSE, TRUE, 'LRT' and TRUE.
#' @export
covariateSearch <- function(project, covariateToSearch = NULL, parameterToTest = NULL, settings = NULL){
  ###################################################################################
  # Initialization
  ###################################################################################
  if(!file.exists(project)){
    message(paste0("ERROR: project : ", project,' does not exists'));
    return(invisible(FALSE))}
  
  # Check and initialize the settings 
  if(!is.null(settings)){
    if(!.checkCovariateSearchInput(inputName = "settings", inputValue = settings)){return(invisible(FALSE))}
  }
  if(is.null(settings$pInclusion)){settings$pInclusion <- 0.1 }
  if(is.null(settings$pElimination)){settings$pElimination <- 0.05 }
  if(is.null(settings$useLinearization)){settings$useLinearization <- FALSE }
  if(is.null(settings$rankedSCM)){settings$rankedSCM <- TRUE }
  if(is.null(settings$criteria)){settings$criteria <- 'LRT' }
  if(is.null(settings$updateInit)){settings$updateInit <- TRUE }
  
  loadProject(project)
  # Check if linearization is possible
  for(indexObservationModel in 1:length(getObservationInformation()$name)){settings$useLinearization <- settings$useLinearization & (getObservationInformation()$type[[indexObservationModel]]=="continuous")}
  
  # check the parameterToTest
  projectParameters <- getIndividualParameterModel()$name
  if(is.null(parameterToTest)){
    validParameters <- projectParameters
  }else{
    if(!.checkCovariateSearchInput(inputName = "parameterToTest", inputValue = parameterToTest)){return(invisible(FALSE))}
    validParameters <- intersect(projectParameters, parameterToTest)
  }
  indivParam = validParameters
  
  # Check the covariateToSearch
  projectCovariates <- getCovariateInformation()$name
  if(length(projectCovariates)==0){
    message("There is no covariate to search")
    return(invisible(FALSE))
  }else{
    if(is.null(covariateToSearch)){
      validCovariates <- projectCovariates
    }else{
      if(!.checkCovariateSearchInput(inputName = "covariateToSearch", inputValue = covariateToSearch)){return(invisible(FALSE))}
      validCovariates <- intersect(projectCovariates, covariateToSearch)
    }
  }
  covariate = validCovariates
  
  if(settings$rankedSCM){
    method <- 'COSSAC'
  }else{
    method <- 'SCM'
  }
  projectToSaveName <- toString(sub(pattern=".mlxtran", replacement=paste0('_covSearch_', method, '.mlxtran'), project))
  saveProject(projectToSaveName);loadProject(projectToSaveName);
 
  # Define the scenario associated to the type of test and the method
  .defineScenario(settings$useLinearization, settings$rankedSCM)
  summary.file = toString(sub(pattern=".mlxtran", replacement="_summary.txt", projectToSaveName))
  
  ######################################################################################################################
  # Initialization on a first run
  ######################################################################################################################
  summary <- c(date(),'\n'); nbRun <- 0; t_strat <- proc.time(); referenceOFV <- NULL;
  if(settings$criteria == 'BIC'){
    criteriaToDisplay <- 'BIC'
    if(settings$useLinearization){
      indexLL <- 3
    }else{
      indexLL <- 4
    }
  }else{
    indexLL <- 1
    criteriaToDisplay <- '-2LL'
  }
  
  # Make a first run if needed
  if((as.list(getLaunchedTasks())$logLikelihoodEstimation==FALSE)){
    bScenario <- runScenario(TRUE); nbRun = nbRun+1;
  }else{
    bScenario <- T
  }
  referenceOFV <- getEstimatedLogLikelihood()[[1]][indexLL]
  
  #############################################################################################################################
  # Forward inclusion step
  #############################################################################################################################
  lineDisplay <- paste0( " ========================================================\n Forward inclusion step (reference ",criteriaToDisplay," = ",format(referenceOFV, nsmall = 1),")\n")
  summary <- c(summary, lineDisplay); cat(lineDisplay)
  # Get all the possibilities
  remainingCovariateStructure <- .getRemainingCovariateStructure(covariate, indivParam)
  bFoward <- TRUE
  while(bFoward&(length(remainingCovariateStructure$indivParam)>0)){
    
    initialEstimates <- getPopulationParameterInformation();
    if(settings$rankedSCM){
      # We test the remaining most probable
      nbConfig <- 1;
    }else{
      # We test all the possibilities
      nbConfig <- length(remainingCovariateStructure$indivParam); 
    }
    OFvalues <- array(dim=nbConfig); estimatedPopParam <- list()
    
    # Check if it is interesting to do the run
    if(remainingCovariateStructure$pValue[1]<.5){
      for(indexConfig in 1:nbConfig){
        # get the parameter - covariate relationship
        evaluatedParameter <- remainingCovariateStructure$indivParam[indexConfig]; 
        evaluatedCovariate <- remainingCovariateStructure$covariate[indexConfig];
        # Initialize the population parameters estimates
        setPopulationParameterInformation(initialEstimates);
        # Run the scenario with te additional parameter - covariate relationship
        OFvalues[indexConfig] <- .setCovAndRun(evaluatedParameter, evaluatedCovariate, indexLL, forwardSearch = T); nbRun = nbRun+1;
        estimatedPopParam[[indexConfig]] <- getEstimatedPopulationParameters() 
        eval(parse(text=paste0('setCovariateModel(',evaluatedParameter,' = list(', evaluatedCovariate, ' = FALSE))')))
        
        lineDisplay <- paste0("Evaluation of covariate ", evaluatedCovariate, ' with parameter ',evaluatedParameter,', (', criteriaToDisplay,' = ', format(OFvalues[indexConfig], nsmall = 1),') \n');cat(lineDisplay)
        summary <- c(summary, lineDisplay);  
      }
      
      # Get the best model and its associated pValue
      indexMin <- which.min(OFvalues) 
      if(settings$criteria == "LRT"){
        pValue <- .getPvalueLRT(referenceOFV, OFV = OFvalues[indexMin], dof = 1)
      }else{
        if(OFvalues[indexMin] < referenceOFV){
          pValue <- 0
        }else{
          pValue <- 1
        }
      }
      if(pValue <= settings$pInclusion){
        evaluatedParameter <- remainingCovariateStructure$indivParam[indexMin]; evaluatedCovariate <- remainingCovariateStructure$covariate[indexMin]
        eval(parse(text=paste0('setCovariateModel(',evaluatedParameter,' = list(', evaluatedCovariate, ' = TRUE))')))
        referenceOFV <- OFvalues[indexMin]
        lineDisplay <- paste0(" +++++++++++++++++++++++ \n covariate ",evaluatedCovariate, ' was included with parameter ', evaluatedParameter,' (pVal ',format(pValue, nsmall = 1),', ',criteriaToDisplay ,' = ',format(referenceOFV, nsmall = 1),') \n'," +++++++++++++++++++++++ \n");
        summary <- c(summary, lineDisplay); cat(lineDisplay)
        if(settings$updateInit){
          newInitialConditions = estimatedPopParam[[indexMin]]
          for(indexParam in 1:length(newInitialConditions)){
            eval(parse(text=paste0('setPopulationParameterInformation(',names(newInitialConditions)[indexParam],' = list(initialValue = ',newInitialConditions[indexParam],'))')))
          }
        }
        # Update the pValue of the remaining pairs
        remainingCovariateStructure <- remainingCovariateStructure[-indexMin,]
        possibleCovariateStructure <- .getCovariateStructure(covariate, indivParam)
        if(length(remainingCovariateStructure$indivParam)>0){
          for(indexRemaining in 1:length(remainingCovariateStructure$indivParam)){
            indexLine <- which((as.character(possibleCovariateStructure$indivParam)==as.character(remainingCovariateStructure$indivParam[indexRemaining]))&(as.character(possibleCovariateStructure$covariate)==as.character(remainingCovariateStructure$covariate[indexRemaining])))
            remainingCovariateStructure$pValue[indexRemaining] <- possibleCovariateStructure$pValue[indexLine]
          }
          # Sort it
          indexSort <- sort(x=remainingCovariateStructure$pValue, decreasing = F, index.return = TRUE)$ix
          remainingCovariateStructure <- remainingCovariateStructure[indexSort,]
        }
      }else{
        setPopulationParameterInformation(initialEstimates);
        if(settings$rankedSCM){# keep searching in the list
          remainingCovariateStructure <- remainingCovariateStructure[-indexMin,]
        }else{ 
          lineDisplay <- paste0("+++++++++++++++++++++++\n No additional covariate \n +++++++++++++++++++++++\n ")
          summary <- c(summary, lineDisplay); cat(lineDisplay)
          bFoward <- FALSE
        }
      }
    }else{
      bFoward <- FALSE
    }
  }
  
  #############################################################################################################################
  # Backward elimination step
  #############################################################################################################################
  # Search the elimination possibilities
  lineDisplay <- paste0( " ========================================================\n Backward elimination step (reference ",criteriaToDisplay," = ",format(referenceOFV, nsmall = 1),")\n")
  summary <- c(summary, lineDisplay); cat(lineDisplay)
  bBackward <- TRUE
  remainingCovariateStructure <- .getCurrentCovariateStructure(covariate, indivParam)
  
  while(bBackward&(length(remainingCovariateStructure$indivParam)>0)){
    
    initialEstimates <- getPopulationParameterInformation();
    if(settings$rankedSCM){
      # We test the remaining most probable
      nbConfig <- 1;
    }else{
      # We test all the possibilities
      nbConfig <- length(remainingCovariateStructure$indivParam); 
    }
    OFvalues <- array(dim=nbConfig); estimatedPopParam <- list()
    
    # Check if it is interesting to do the run
    if(remainingCovariateStructure$pValue[1]>.01){
      for(indexConfig in 1:nbConfig){
        # get the parameter - covariate relationship
        evaluatedParameter <- remainingCovariateStructure$indivParam[indexConfig]; 
        evaluatedCovariate <- remainingCovariateStructure$covariate[indexConfig];
        # Initialize the population parameters estimates
        setPopulationParameterInformation(initialEstimates);
        # Run the scenario with te additional parameter - covariate relationship
        OFvalues[indexConfig] <- .setCovAndRun(evaluatedParameter, evaluatedCovariate, indexLL, forwardSearch = F); nbRun = nbRun+1;
        estimatedPopParam[[indexConfig]] <- getEstimatedPopulationParameters() 
        eval(parse(text=paste0('setCovariateModel(',evaluatedParameter,' = list(', evaluatedCovariate, ' = TRUE))')))
        
        lineDisplay <- paste0("Evaluation of covariate ", evaluatedCovariate, ' not on parameter ',evaluatedParameter,', (', criteriaToDisplay,' = ', format(OFvalues[indexConfig], nsmall = 1),') \n');cat(lineDisplay)
        summary <- c(summary, lineDisplay);  
      }
      
      # Get the best model and its associated pValue
      indexMin <- which.min(OFvalues) 
      if(settings$criteria == "LRT"){
        pValue <- .getPvalueLRT(referenceOFV, OFV = OFvalues[indexMin], dof = -1)
      }else{
        if(OFvalues[indexMin] > referenceOFV){
          pValue <- 0
        }else{
          pValue <- 1
        }
      }
      if(pValue >= settings$pElimination){
        evaluatedParameter <- remainingCovariateStructure$indivParam[indexMin]; evaluatedCovariate <- remainingCovariateStructure$covariate[indexMin]
        eval(parse(text=paste0('setCovariateModel(',evaluatedParameter,' = list(', evaluatedCovariate, ' = FALSE))')))
        referenceOFV <- OFvalues[indexMin]
        lineDisplay <- paste0(" +++++++++++++++++++++++ \n covariate ",evaluatedCovariate, ' was removed from parameter ', evaluatedParameter,' (pVal ',format(pValue, nsmall = 1),', ',criteriaToDisplay ,' = ',format(referenceOFV, nsmall = 1),') \n'," +++++++++++++++++++++++ \n");
        summary <- c(summary, lineDisplay); cat(lineDisplay)
        if(settings$updateInit){
          newInitialConditions = estimatedPopParam[[indexMin]]
          for(indexParam in 1:length(newInitialConditions)){
            eval(parse(text=paste0('setPopulationParameterInformation(',names(newInitialConditions)[indexParam],' = list(initialValue = ',newInitialConditions[indexParam],'))')))
          }
        }
        # Update the pValue of the remaining pairs
        remainingCovariateStructure <- remainingCovariateStructure[-indexMin,]
        possibleCovariateStructure <- .getCovariateStructure(covariate, indivParam)
        if(length(remainingCovariateStructure$indivParam)>0){
          for(indexRemaining in 1:length(remainingCovariateStructure$indivParam)){
            indexLine <- which((as.character(possibleCovariateStructure$indivParam)==as.character(remainingCovariateStructure$indivParam[indexRemaining]))&(as.character(possibleCovariateStructure$covariate)==as.character(remainingCovariateStructure$covariate[indexRemaining])))
            remainingCovariateStructure$pValue[indexRemaining] <- possibleCovariateStructure$pValue[indexLine]
          }
          # Sort it
          indexSort <- sort(x=remainingCovariateStructure$pValue, decreasing = T, index.return = TRUE)$ix
          remainingCovariateStructure <- remainingCovariateStructure[indexSort,]
        }
      }else{
        setPopulationParameterInformation(initialEstimates);
        if(settings$rankedSCM){
          # keep searching in the list
          remainingCovariateStructure <- remainingCovariateStructure[-indexMin,]
        }else{
          lineDisplay <- paste0("+++++++++++++++++++++++\n No covariate to remove\n +++++++++++++++++++++++\n ")
          summary <- c(summary, lineDisplay); cat(lineDisplay)
          bBackward <- FALSE
        }
      }
    }else{
      bBackward <- FALSE
    }
  }
  
  # Make the summary
  runScenario(TRUE); nbRun = nbRun+1;
  OFValue = getEstimatedLogLikelihood()[[1]][indexLL]
  
  summary.file = toString(sub(pattern=".mlxtran", replacement="_summary.txt", projectToSaveName))
  summary <- c(summary, paste0("========================================================\n FINAL MODEL (",criteriaToDisplay," = ",format(OFValue, nsmall = 1),") \n"))
  summary <- c(summary, paste0("-> target parameters: ",toString(indivParam)," \n"," -> searched covariates: ",toString(covariate)," \n\n"))
  backwardList <- .getCurrentCovariateStructure(covariate, indivParam)
  if(!is.null(backwardList$covariate)){# There are covariates 
    for(indexList in 1:length(backwardList[,1])){
      lineDisplay <- paste0('Covariate ',backwardList$covariate[indexList], ' on parameter ',backwardList$indivParam[indexList],'\n');summary <- c(summary, lineDisplay)
    }
  }
  summary <- c(summary, c(paste0("\n Done with ",toString(nbRun)," runs in ",toString(floor(proc.time()[3] - t_strat[3])),"s\n", date(),'\n'),"========================================================\n\n"))
  cat(summary); cat(summary, file = summary.file)
  saveProject(projectFile = projectToSaveName)
}  

###################################################################################
# Check the inputs 
###################################################################################
.checkCovariateSearchInput = function(inputName, inputValue){
  isValid = TRUE
  inputName = tolower(inputName)
  if(inputName == tolower("parameterToTest")){
    if(is.vector(inputValue) == FALSE){
      message("ERROR: Unexpected type encountered. parameterToTest must be a vector")
      isValid = FALSE
    }else if(length(intersect(getIndividualParameterModel()$name, inputValue))==0){
      message("ERROR: parameterToTest does not have valid parameter in its definition.")
      isValid = FALSE
    }
  }else if(inputName == tolower("covariateToSearch")){
    if(is.vector(inputValue) == FALSE){
      message("ERROR: Unexpected type encountered. covariateToSearch must be a vector")
      isValid = FALSE
    }else if(length(intersect(getCovariateInformation()$name, inputValue))==0){
      message("ERROR: covariateToSearch have no valid covariate in its definition.")
      isValid = FALSE
    }
  }else if(inputName == tolower("settings")){
    if(is.list(inputValue) == FALSE){
      message("ERROR: Unexpected type encountered. settings must be a list")
      isValid = FALSE
    }else {
      for (i in 1:length(inputValue)){
        if(!.checkCovariateSearchSettings(settingName = names(inputValue)[i], settingValue = inputValue[[i]])){
          isValid = FALSE
        }
      }
    }
  }
  return(invisible(isValid))
}

.checkCovariateSearchSettings = function(settingName, settingValue){
  isValid = TRUE
  settingName = tolower(settingName)
  if(settingName == tolower("pInclusion")){
    if(is.double(settingValue) == FALSE){
      message("ERROR: Unexpected type encountered. pInclusion must be a double.")
      isValid = FALSE
    }else{
      if(settingValue<0){
        message("ERROR: pInclusion must be strictly positive.")
        isValid = FALSE
      }
    }
  }else if(settingName == tolower("pElimination")){
    if(is.double(settingValue) == FALSE){
      message("ERROR: Unexpected type encountered. pElimination must be a double.")
      isValid = FALSE
    }else{
      if(settingValue<0){
        message("ERROR: pElimination must be strictly positive.")
        isValid = FALSE
      }
    }
  }else if(settingName == tolower("useLinearization")){
    if(is.logical(settingValue) == FALSE){
      message("ERROR: Unexpected type encountered. useLinearization must be a boolean.")
      isValid = FALSE
    }
  }else if(settingName == tolower("rankedSCM")){
    if(is.logical(settingValue) == FALSE){
      message("ERROR: Unexpected type encountered. rankedSCM must be a boolean.")
      isValid = FALSE
    }
  }else if(settingName == tolower("criteria")){
    if(is.character(settingValue) == FALSE){
      message("ERROR: Unexpected type encountered. criteria must be a string")
      isValid = FALSE
    }else if(length(intersect(tolower(settingValue), c(tolower('BIC'), tolower('LRT'))))==0){
      message("ERROR: method must be either 'BIC' or 'LRT'.")
      isValid = FALSE
    }
  }else if(settingName == tolower("updateInit")){
    if(is.logical(settingValue) == FALSE){
      message("ERROR: Unexpected type encountered. updateInit must be a boolean.")
      isValid = FALSE
    }
  }else{
    message("WARNING: ",settingName,' is not a valid setting')
  }
  return(isValid)
}

######################################################################################################################
# Get current covariate structure
# The output is a data.frame with indivParam, covariate and pValue (order by wald pValue is required)
######################################################################################################################
.getCurrentCovariateStructure <- function(covariate, indivParam){
  param <- cov <- pValue <-NULL
  for(indexParam in 1:length(indivParam)){
    hasCov <- eval(parse(text=paste0('which(getIndividualParameterModel()$covariateModel$',indivParam[indexParam],'==TRUE)')))
    if(length(hasCov)>0){
      covOnIndiv <- intersect(covariate, names(hasCov))
      if(length(covOnIndiv)>0){
        cov <- c(cov, covOnIndiv)
        param <- c(param, rep(indivParam[indexParam],length(covOnIndiv)))
        pValue <- c(pValue, rep(1,length(covOnIndiv)))
      }
    }
  }
  # Compute the pValue
  testValues <- .getTestsIndivParam()
  if(!is.null(testValues)){
    for(indexDF in 1:length(param)){
      indexTest <- which((as.character(testValues$parameter)==param[indexDF])&(as.character(testValues$covariate)==cov[indexDF]))
      if(length(indexTest)>0){
        pValue[indexDF] <- testValues[indexTest, 4]
      }
    }
  }
  
  if(is.null(cov)){
    currentCovariateStructure <- NULL
  }else{
    indexSort <- sort(x=pValue, decreasing = TRUE, index.return = TRUE)$ix
    currentCovariateStructure <- data.frame(indivParam = param[indexSort], covariate = cov[indexSort], pValue=pValue[indexSort])
  }
  
  return(currentCovariateStructure)
  
}

######################################################################################################################
# Get full covariate structure
# The output is a data.frame with indivParam, covariate and pValue (order by test pValue is required)
######################################################################################################################
.getCovariateStructure <- function(covariate, indivParam){
  dfOut <- data.frame(indivParam = rep(indivParam, times = 1, each = length(covariate)), 
                      covariate = rep(covariate, length(indivParam)), 
                      pValue = rep(.025, length(indivParam)*length(covariate)))
  
  indexFinal <- 1:length(dfOut[,3])
  
  testValues <- .getTestsRandomEffects() 
  if(!is.null(testValues)){
    for(indexDF in 1:length(dfOut[,1])){
      indexTest <- which((as.character(testValues$eta)==paste0('eta_',dfOut$indivParam[indexDF]))&(as.character(testValues$covariate)==as.character(dfOut$covariate[indexDF])))
      if(length(indexTest)>0){
        dfOut$pValue[indexDF] <- testValues[indexTest, 4]
      }
    }
    indexFinal <-  sort(dfOut$pValue, index.return= T)$ix
  }
  
  return(dfOut[indexFinal,])
}

######################################################################################################################
# Get the remaining covariate structure
# The output is a data.frame with indivParam, covariate and pValue (order by test pValue is required)
######################################################################################################################
.getRemainingCovariateStructure<- function(covariate, indivParam){
  currentCovariateStructure <- .getCurrentCovariateStructure(covariate, indivParam)
  possibleCovariateStructure <- .getCovariateStructure(covariate, indivParam)
  
  remainingCovariateStructure <- possibleCovariateStructure
  if(length(currentCovariateStructure[,1])>0){
    for(index in 1:length(currentCovariateStructure[,1])){
      indexLine <- which((as.character(remainingCovariateStructure$indivParam)==as.character(currentCovariateStructure$indivParam[index]))&(as.character(remainingCovariateStructure$covariate)==as.character(currentCovariateStructure$covariate[index])))
      remainingCovariateStructure <- remainingCovariateStructure[-indexLine,]
    }
  }
  return(remainingCovariateStructure)
}

######################################################################################################################
# Get the tests associated to the random effects w.r.t. the covariates
######################################################################################################################
.getTestsRandomEffects <- function(){
  testFile <- paste0(getProjectSettings()$directory,'/Tests/correlationRandomEffectsCovariates.txt')
  if(file.exists(testFile)){
    test <- read.table(file = testFile, header = T, sep = "")
    if(length(test[1,])==1){test <- read.table(file = testFile, header = T, sep = ",")}
    if(length(test[1,])==1){test <- read.table(file = testFile, header = T, sep = ";")}
    if(length(test[1,])==1){test <- read.table(file = testFile, header = T, sep = "\t")}
  }else{
    test <- NULL
  }
  return(test)
}

######################################################################################################################
# Get the tests associated to the individual parameters w.r.t. the covariates
######################################################################################################################
.getTestsIndivParam <- function(){
  testFile <- paste0(getProjectSettings()$directory,'/Tests/correlationIndividualParametersCovariates.txt')
  if(file.exists(testFile)){
    test <- read.table(file = testFile, header = T, sep = "")
    if(length(test[1,])==1){test <- read.table(file = testFile, header = T, sep = ",")}
    if(length(test[1,])==1){test <- read.table(file = testFile, header = T, sep = ";")}
    if(length(test[1,])==1){test <- read.table(file = testFile, header = T, sep = "\t")}
  }else{
    test <- NULL
  }
  return(test)
}

######################################################################################################################
# Compute the pValue based on a statistical test (LogLikelihood Ratio test) 
######################################################################################################################
#' @importFrom stats pchisq
.getPvalueLRT <- function(referenceOFV, OFV, dof){
  if(dof>0){
    pValue <- 1-pchisq(max(referenceOFV-OFV,.1), df=dof)
  }else{
    pValue <- 1-pchisq(max(-(referenceOFV-OFV),.1), df=-dof)
  }
  return(pValue)
}


#############################################################################################################################
# Get the scenario for the run
#############################################################################################################################
.defineScenario <- function(useLinearization, useTests){
  # Define the scenario associated to the type of test and the method
  scenario <- getScenario()
  if(useLinearization){
    scenario$linearization = T
  }else{
    scenario$linearization = F 
  }
  scenario$plotList <- "covariatemodeldiagnosis";
  
  if(useTests){# compute both conditional distribution and mode
    scenario$tasks <- c(populationParameterEstimation = T, conditionalDistributionSampling=T, conditionalModeEstimation = T, logLikelihoodEstimation = T, plots = T);
  }else{
    if(useLinearization){# compute only mode
      scenario$tasks <- c(populationParameterEstimation = T, conditionalModeEstimation = T, logLikelihoodEstimation = T)
    }else{# compute only mode
      scenario$tasks <- c(populationParameterEstimation = T, conditionalDistributionSampling=T, logLikelihoodEstimation = T)
    }
  }
  
  setScenario(scenario)
}

#############################################################################################################################
# Set the covariate, run the scenario and get the pValue
#############################################################################################################################
.setCovAndRun <- function(evaluatedParameter, evaluatedCovariate, indexLL, forwardSearch){
  if(forwardSearch){
    eval(parse(text=paste0('setCovariateModel(',evaluatedParameter,' = list(', evaluatedCovariate, ' = TRUE))')))
  }else{
    eval(parse(text=paste0('setCovariateModel(',evaluatedParameter,' = list(', evaluatedCovariate, ' = FALSE))')))
  }
  # Increase the associated omega to be able to explore the domain if it has variability
  if(length(intersect(x = evaluatedParameter, y = names(which(getIndividualParameterModel()$variability[[1]]==TRUE))))>0){
    eval(parse(text=paste0('setPopulationParameterInformation(omega_',evaluatedParameter,' = list(initialValue = 1))')))
  }
  bScenario <- runScenario(TRUE);
  
  if(bScenario){# The run is ok
    OFValue = getEstimatedLogLikelihood()[[1]][indexLL]
  }else{
    OFValue <- Inf
  }
  return(OFValue)
}