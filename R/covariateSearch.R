#' Covariate model building
#' 
#' Automatic search of the best covariate model. Automatic covariate model building is available directly in the lixoftConnectors package using the function
#' \code{runModelBuilding}. Please migrate, as this function will be deprecated in the future. \cr \cr
#
#' Two methods for covariate model building are proposed
#' \itemize{
#' \item SCM: stepwise covariate modeling method
#' In the forward selection, at each step, each of the remaining (i.e not yet included) parameter-covariate relationships are added to the model in an univariate model (one model per relationship), and run. Among all models, the model that improves some criteria (LRT, BIC or AIC) most is selected and taken forward to the next step.
#' During backward elimination, parameter-covariate relationships are removed in an univariate manner. 
#' \item COSSAC: COnditional Sampling for Stepwise Approach based on Correlation tests method 
#' COSSAC makes use of the information contained in the base model run to choose which covariate to try first (instead of trying all covariates "blindly" as in SCM). 
#' Indeed, the correlation between the individual parameters (or random effects) and the covariates hints at possibly relevant parameter-covariate relationships. If the EBEs (empirical Bayes estimates) are used, shrinkage may bias the result. 
#' COSSAC instead uses samples from the a posteriori conditional distribution (available as "conditional distribution" task in MonolixSuite2018) to calculate the correlation between the random effects and covariates. 
#' A p-value can be derived using the Pearson's correlation test for continuous covariate and ANOVA for categorical covariate. The p-values are used to sort all the random effect-covariate relationships. Relationships with the lowest p-value are added first, run and confirmed using a likelihood ratio test, AIC or BIC criteria. 
#' }
#' @param project a Monolix project
#' @param final.project [optional] string corresponding to the final Monolix project (default: 'runFinal.mlxtran' in covariate search output folder)
#' @param method [optional] string correspondig to the method. It can be 'COSSAC' or 'SCM'. By default, COSSAC' is used. 
#' @param covToTest [optional] vector of covariates to test. Cannot be used if testRelations is defined. By default, all covariates are tested. 
#' @param covToTransform [optional] vector of covariates to transform. The transformation consists in a log transform of the covariate with centering by the mean value (ex: WT is transformed into log(WT/mean) with mean the mean WT value over the individuals of the data set). Both the transformed and untransformed covariate are tested by the algorithm. By default, no covariate is transformed.
#' Note: 
#' adding a non-transformed covariate on a lognormally distributed parameter results in an exponential relationship: log(V) = log(Vpop) + beta*WT + eta <=> V = Vpop * exp(beta*WT) * exp(eta)
#' adding a log-transformed covariate on a lognormally distributed parameter results in a power law relationship: log(V) = log(Vpop) + beta*log(WT/70) + eta <=> V = Vpop * (WT/70)^beta * exp(eta)
#' @param paramToUse [optional] vector of parameters which may be function of covariates. Cannot be used if testRelations is defined. By default, all parameters are tested.
#' @param testRelations [optional] list of parameter-covariate relationships to test, ex: list(V=c("WT","SEX"),Cl=c("CRCL")). Cannot be used if covToTest or paramToUse is defined. By default, all parameter-covariate relationships are tested.
#' @param settings [optional] list of settings for the covariate search: 
#' \itemize{
#' \item \code{pInclusion}  [positive double] threshold on the LRT p-value to accept the model with the added parameter-covariate relationship during forward selection (default = .1). Only used if criteria="LRT".
#' \item \code{pElimination} [positive double] threshold on the LRT p-value to accept the model without the removed parameter-covariate relationship during the backward elimination (default = .05). Only used if criteria="LRT".
#' \item \code{criteriaThreshold} [positive double] the threshold on the AIC or BIC difference to accept the model with added/removed parameter-covariate relationship (default = 0). Only used if criteria="BIC" or "AIC.
#' \item \code{linearization} [boolean] whether the computation of the likelihood is based on a linearization of the model (default = FALSE).
#' \item \code{criteria} [string] criteria to optimize. It can be the "BIC", "AIC", or "LRT"  (default="LRT").
#' \item \code{direction} [string] method for covariate search. It can be "backward", "forward", or "both" (default = "both").
#' \item \code{updateInit} [boolean] whether to update or not the initial parameters using the estimates of the parent model (default = FALSE)
#' \item \code{saveRun} [boolean] whether to save or not each run (default = TRUE)
#' }
#' @seealso \code{getModelBuildingSettings} settings for model building with lixoftConnectors \cr
#' \code{runModelBuilding} run model building with lixoftConnectors \cr
#' \code{getModelBuildingResults} results for model building with lixoftConnectors
#' @examples
#' \dontrun{
#' # RsmlxDemo1.mlxtran is a Monolix project for modelling the pharmacokinetics (PK) of warfarin 
#' # using a PK model with parameters ka, V, Cl.
#' 
#' # In this example, three covariates (wt, age, sex) are available with the data
#' # covariatesearch will compute the best covariate model, in term of BIC, 
#' # for the three PK parameters using the three covariates. 
#' r1 <- covariateSearch(project="RsmlxDemo1.mlxtran")
#'   
#' # Instead of using the COSSAC method, we can use the SCM method:
#' r2 <- covariateSearch(project="RsmlxDemo1.mlxtran", method = 'SCM')
#' 
#' # Here, the covariate model is built using age and wt only, for V and Cl only:
#' r3 <- covariateSearch(project    = "RsmlxDemo1.mlxtran", 
#'                       paramToUse = c("V","Cl"), 
#'                       covToTest  = c("age","wt"))
#' }
#' 
#' # See http://monolix.lixoft.com/rsmlx/covariatesearch/ for detailed examples of covariatesearch
#' # Download the demo examples here: http://monolix.lixoft.com/rsmlx/installation
#'
#' @export
covariateSearch <- function(project, final.project=NULL, method = NULL, covToTest = NULL, covToTransform=NULL, paramToUse = NULL, testRelations = NULL, settings = NULL){
  warning("This function will be deprecated. Please migrate to lixoftConnectors::runModelBuilding().", 
          call. = FALSE, immediate. = TRUE)
  
  ###################################################################################
  # Initial check
  ###################################################################################
  # Check the validity of the project
  
  initRsmlx()
  
  r <- prcheck(project, f="cov", paramToUse=paramToUse, method=method)
  if (r$demo)
    return(r$res)
  project <- r$project
  
  # Check the validity of the final.project
  if(!is.null(final.project)){
    if(file.exists(final.project)){
      message(paste0("ERROR: project '", final.project, "' already exists"))
      return(invisible(FALSE))
    }else{
      if(length(grep(pattern = '.mlxtran', x = final.project))<1){
        message(paste0("ERROR: project '", project, "' should have a .mlxtran extension"))
        return(invisible(FALSE))
      }
    }
  }
  
  if(is.null(method)){
    method <- 'COSSAC'
  }
  if(!.checkCovariateSearchInput(inputName = "method", inputValue = method)){return(invisible(FALSE))}
  
  # check the paramToUse
  projectParameters <- mlx.getIndividualParameterModel()$name
  if(is.null(paramToUse)){
    validParameters <- projectParameters
  }else{
    if(!.checkCovariateSearchInput(inputName = "paramToUse", inputValue = paramToUse)){return(invisible(FALSE))}
    validParameters <- intersect(projectParameters, paramToUse)
  }
  indivParam = validParameters
  
  # check the covariate to transform
  meanCov <- NULL
  if(!is.null(covToTransform)){
    if(!.checkCovariateSearchInput(inputName = "covToTest", inputValue = covToTransform)){return(invisible(FALSE))}
    if(is.null(covToTest) && is.null(testRelations)){
      covToTest <- mlx.getCovariateInformation()$name
    }
    for(index in 1:length(covToTransform)){
      cov <- covToTransform[index]
      indexCov <- NULL
      eval(parse(text=paste0('indexCov <- which(names(mlx.getCovariateInformation()$type)=="',cov,'")')))
      if(length(intersect(mlx.getCovariateInformation()$type[indexCov],"continuous"))>0){
        eval(parse(text=paste0('meanCov <- mean(mlx.getCovariateInformation()$covariate$',cov,')')))
        newCov <- paste0("log_",cov)
        eval(parse(text=paste0('lixoftConnectors::addContinuousTransformedCovariate(',newCov,'="log(',cov,'/',toString(meanCov),')")')))
        cat(paste0(newCov," was added. \n" ))
        
        # Add it in the covariate to test
        if (!is.null(covToTest)) {
          covToTest = c(covToTest, newCov)
        }
      }else{
        warning(paste0("Covariate ",cov," can not be transformed as a continuous covariate"))        
      }
    }
  }
  
  # Check the covToTest
  projectCovariates <- mlx.getCovariateInformation()$name
  if(length(projectCovariates)==0){
    message("There is no covariate to test")
    return(invisible(FALSE))
  }else{
    if(is.null(covToTest)){
      validCovariates <- projectCovariates
    }else{
      if(!.checkCovariateSearchInput(inputName = "covToTest", inputValue = covToTest)){return(invisible(FALSE))}
      validCovariates <- intersect(projectCovariates, covToTest)
    }
  }
  covariate = validCovariates
  
  # Check the testRelations
  if(!is.null(testRelations)){
    if(!is.null(covToTest)||!is.null(paramToUse)){
      message(paste0("ERROR: testRelations can not be defined if either covToTest or paramToUse is not NULL"))
      return(invisible(FALSE))
    }
    if(!.checkCovariateSearchInput(inputName = "testRelations", inputValue = testRelations)){return(invisible(FALSE))}
  }
  
  # Check and initialize the settings 
  if(!is.null(settings)){
    if(!.checkCovariateSearchInput(inputName = "settings", inputValue = settings)){return(invisible(FALSE))}
  }
  if(is.null(settings$pInclusion)){settings$pInclusion <- 0.1 }
  if(is.null(settings$pElimination)){settings$pElimination <- 0.05 }
  if(is.null(settings$criteria)){settings$criteria <- 'LRT' }
  if(is.null(settings$criteriaThreshold)){settings$criteriaThreshold <- 0 }
  if(is.null(settings$linearization)){settings$linearization <- FALSE }
  if(is.null(settings$rankedSCM)){settings$rankedSCM <- TRUE }
  if(is.null(settings$criteria)){settings$criteria <- 'LRT' }
  if(is.null(settings$direction)){settings$direction <- 'both' }
  if(is.null(settings$updateInit)){settings$updateInit <- FALSE }
  if(is.null(settings$saveRun)){settings$saveRun <- TRUE }
  
  # Check if linearization is possible
  for(indexObservationModel in 1:length(mlx.getObservationInformation()$name)){settings$linearization <- settings$linearization & (mlx.getObservationInformation()$type[[indexObservationModel]]=="continuous")}
  
  
  ###################################################################################
  # Structure initialization
  ###################################################################################
  # Initialization of relationshipsToTest
  if(is.null(testRelations)){
    relationshipsToTest <- data.frame(indivParam = rep(indivParam, times = 1, each = length(covariate)), 
                                      covariate = rep(covariate, length(indivParam)))
  }else{
    relationshipsToTest <- data.frame(indivParam = names(testRelations)[1], covariate = testRelations[[1]])
    if(length(testRelations)>1){
      for(index in 2:length(testRelations)){
        relationshipsToTest <-  rbind(relationshipsToTest, data.frame(indivParam = names(testRelations)[index], covariate = testRelations[[index]]))
      }
    }
  }
  
  if(method=='COSSAC'){
    settings$rankedSCM <- T
  }else{
    settings$rankedSCM <- F
  }
  
  nbRun <- 0;
  outputDirectory <- mlx.getProjectSettings()$directory
  covariateSeachOutputFolder <- paste0(outputDirectory,'/covariateSearch_',method, '_',settings$criteria,'/')
  dir.create(path = covariateSeachOutputFolder, showWarnings = F, recursive = T)
  
  if(is.null(final.project)){
    final.project <- paste0(covariateSeachOutputFolder,'runFinal.mlxtran')
  }
  
  # Define the scenario associated to the type of test and the method
  .defineScenario(settings$linearization, settings$rankedSCM)
  
  summary.file <- paste0(covariateSeachOutputFolder,"covSearchSummary.txt")
  if(settings$rankedSCM){additionalDisplay <- ''}else{additionalDisplay <- "+++++++++++++++++++++++\n"}
  projectRun <- paste0(covariateSeachOutputFolder,'run_',toString(nbRun),'.mlxtran');mlx.saveProject(projectFile = projectRun);
  
  if(settings$criteria == 'BIC'){
    criteriaToDisplay <- 'BIC'
    indexLL <- 4 # Index in mlx.getEstimatedLogLikelihood()[[1]] function
    displayPvalue <- FALSE
  }else if(settings$criteria == 'AIC'){
    criteriaToDisplay <- 'AIC'
    indexLL <- 3
    displayPvalue <- FALSE
  }else{
    indexLL <- 1
    criteriaToDisplay <- 'OFV'
    displayPvalue <- TRUE
  }
  if(settings$linearization){
    indexLL = max(indexLL-1,1)
  }
  ######################################################################################################################
  # Initialization on a first run
  ######################################################################################################################
  summary <- c(date(),'\n');  t_strat <- proc.time(); referenceOFV <- NULL;
  
  # Make a first run
  bScenario <- mlx.runScenario(); nbRun = nbRun+1;
  referenceOFV <- mlx.getEstimatedLogLikelihood()[[1]][indexLL]
  
  #############################################################################################################################
  # Forward inclusion step
  #############################################################################################################################
  # Get all the possibilities
  bForward <- length(intersect(settings$direction,c("both","forward")))
  remainingCovariateStructure <- .getRemainingCovariateStructure(relationshipsToTest)
  if(bForward){
    lineDisplay <- paste0("========================================================\nForward inclusion step (reference objective function ",criteriaToDisplay," = ",format(referenceOFV, nsmall = 1),")\n")
    summary <- c(summary, lineDisplay); cat(lineDisplay);cat(summary, file = summary.file)
  }
  while(bForward&(length(remainingCovariateStructure$indivParam)>0)){
    initialEstimates <- mlx.getPopulationParameterInformation(); # Keep the initial estimates of the project
    if(settings$rankedSCM){
      idConfig <- which.min(x = remainingCovariateStructure$pValue); # We test the remaining most probable
    }else{
      idConfig <- 1:length(remainingCovariateStructure$indivParam); # We test all the possibilities
    }
    # Initialize the results of the current step
    OFvalues <- array(dim=c(length(idConfig),3)); estimatedPopParam <- list()
    
    # Check if it is interesting to do the run
    if(min(remainingCovariateStructure$pValue[idConfig])<.4){
      iterConfig <- 0  
      for(indexConfig in idConfig){
        iterConfig <- iterConfig+1
        # get the parameter - covariate relationship
        evaluatedParameter <- remainingCovariateStructure$indivParam[indexConfig];
        evaluatedCovariate <- remainingCovariateStructure$covariate[indexConfig];
        # Initialize the population parameters estimates and add the parameter - covariate relaionship
        mlx.setPopulationParameterInformation(initialEstimates);
        eval(parse(text=paste0('lixoftConnectors::setCovariateModel(',evaluatedParameter,' = list(', evaluatedCovariate, ' = TRUE))')))
        # Run the scenario with the additional parameter - covariate relationship
        if(settings$saveRun){projectRun <- paste0(covariateSeachOutputFolder,'run_',toString(nbRun),'.mlxtran');mlx.saveProject(projectRun);}
        OFvalues[iterConfig,1] <- .getOFV(indexLL); nbRun = nbRun+1;
        OFvalues[iterConfig,2] <- nbRun-1
        OFvalues[iterConfig,3] <- .getDof(covariate = evaluatedCovariate)
        estimatedPopParam[[iterConfig]] <- mlx.getEstimatedPopulationParameters()
        # Get back to the initial state before decision
        eval(parse(text=paste0('lixoftConnectors::setCovariateModel(',evaluatedParameter,' = list(', evaluatedCovariate, ' = FALSE))')))
        
        if(settings$rankedSCM == F){endofline='\n'}else{endofline=paste0('')}
        lineDisplay <- paste0("Run ",toString(OFvalues[iterConfig,2]),": Evaluation of adding covariate ", evaluatedCovariate, ' on parameter ',evaluatedParameter,' ==> ', criteriaToDisplay,' = ', format(as.double(OFvalues[iterConfig,1]), nsmall = 1),', ',criteriaToDisplay,'-difference = ',format(as.double(referenceOFV-OFvalues[iterConfig,1]), nsmall = 1),endofline);
        summary <- c(summary, lineDisplay);cat(lineDisplay);cat(summary, file = summary.file);  
      }
      
      # Get the best model and its associated pValue
      if(settings$criteria == "LRT"){
        pValueAll <- .getPvalueLRT(referenceOFV, OFV = as.double(OFvalues[,1]), dof = as.double(OFvalues[,3]))
        iterMin <- which.min(pValueAll)
        indexMin <- idConfig[iterMin]
        pValue <- pValueAll[iterMin]
      }else{
        iterMin <- which.min(as.double(OFvalues[,1]));
        indexMin <- idConfig[iterMin]
        if(as.double(OFvalues[iterMin,1])<referenceOFV-settings$criteriaThreshold){
          pValue <- 0
        }else{
          pValue <- 1
        }
      }
      
      if(pValue <= settings$pInclusion){ # If the pValue is valid wrt the inclusion criteria
        # Add the parameter and list
        evaluatedParameter <- remainingCovariateStructure$indivParam[indexMin]; evaluatedCovariate <- remainingCovariateStructure$covariate[indexMin]
        eval(parse(text=paste0('lixoftConnectors::setCovariateModel(',evaluatedParameter,' = list(', evaluatedCovariate, ' = TRUE))')))
        referenceOFV <- OFvalues[iterMin,1]
        if(settings$criteria == "LRT" & settings$rankedSCM==T){pValueDisplay <- paste0(', pVal = ',format(pValue, digits = 3),'')}else{pValueDisplay<-''}
        lineDisplay <- paste0(pValueDisplay,' \n',"=> covariate ", evaluatedCovariate, ' was included on parameter ', evaluatedParameter,' (new ref. ',criteriaToDisplay ,' = ',format(referenceOFV, nsmall = 1),') \n',additionalDisplay);
        summary <- c(summary, lineDisplay); cat(lineDisplay);cat(summary, file = summary.file)
        if(settings$updateInit){
          newInitialConditions = estimatedPopParam[[iterMin]]
          for(indexParam in 1:length(newInitialConditions)){
            eval(parse(text=paste0('lixoftConnectors::setPopulationParameterInformation(',names(newInitialConditions)[indexParam],' = list(initialValue = ',newInitialConditions[indexParam],'))')))
          }
        }
        # Update the pValue of the remaining pairs
        remainingCovariateStructure <- .updateRemainingPvalue(structure = remainingCovariateStructure[-indexMin,])
      }else{
        mlx.setPopulationParameterInformation(initialEstimates);
        if(settings$rankedSCM){# keep searching in the list
          if(settings$criteria == "LRT" & settings$rankedSCM==T){pValueDisplay <- paste0(', pVal = ',format(pValue, digits = 3),'')}else{pValueDisplay<-''}
          lineDisplay <- paste0(pValueDisplay,' \n',"=> covariate ",remainingCovariateStructure$covariate[indexMin], ' was not included on parameter ', remainingCovariateStructure$indivParam[indexMin],' (ref. ',criteriaToDisplay ,' = ',format(referenceOFV, nsmall = 1),') \n',additionalDisplay);
          summary <- c(summary, lineDisplay); cat(lineDisplay);cat(summary, file = summary.file)
          remainingCovariateStructure <- remainingCovariateStructure[-indexMin,]
        }else{ 
          lineDisplay <- paste0("\n+++++++++++++++++++++++\n No additional covariate \n+++++++++++++++++++++++\n")
          summary <- c(summary, lineDisplay); cat(lineDisplay);cat(summary, file = summary.file)
          bForward <- FALSE
        }
      }
    }else{
      bForward <- FALSE
    }
  }
  
  #############################################################################################################################
  # Backward elimination step
  #############################################################################################################################
  bBackward <- length(intersect(settings$direction,c("both","backward")))
  remainingCovariateStructure <- .getCurrentCovariateStructure(relationshipsToTest)
  
  if(bBackward){
    # Search the elimination possibilities
    lineDisplay <- paste0("========================================================\nBackward elimination step (reference ",criteriaToDisplay," = ",format(referenceOFV, nsmall = 1),")\n")
    summary <- c(summary, lineDisplay); cat(lineDisplay);cat(summary, file = summary.file)
  }
  
  while(bBackward&(length(remainingCovariateStructure$indivParam)>0)){
    initialEstimates <- mlx.getPopulationParameterInformation(); # Keep the initial estimates of the project
    if(settings$rankedSCM){
      idConfig <- which.max(x = remainingCovariateStructure$pValue); # We test the remaining most probable
    }else{
      idConfig <- 1:length(remainingCovariateStructure$indivParam); # We test all the possibilities
    }
    # Initialize the results of the current step
    OFvalues <- array(dim=c(length(idConfig),3)); estimatedPopParam <- list()
    
    # Check if it is interesting to do the run
    if(max(remainingCovariateStructure$pValue[idConfig])>.01){
      iterConfig <- 0  
      for(indexConfig in idConfig){
        iterConfig <- iterConfig+1
        # get the parameter - covariate relationship
        evaluatedParameter <- remainingCovariateStructure$indivParam[indexConfig];
        evaluatedCovariate <- remainingCovariateStructure$covariate[indexConfig];
        # Initialize the population parameters estimates and add the parameter - covariate relaionship
        mlx.setPopulationParameterInformation(initialEstimates);
        eval(parse(text=paste0('lixoftConnectors::setCovariateModel(',evaluatedParameter,' = list(', evaluatedCovariate, ' = FALSE))')))
        # Run the scenario with the additional parameter - covariate relationship
        if(settings$saveRun){projectRun <- paste0(covariateSeachOutputFolder,'run_',toString(nbRun),'.mlxtran');mlx.saveProject(projectRun);}
        OFvalues[iterConfig,1] <- .getOFV(indexLL); nbRun = nbRun+1;
        OFvalues[iterConfig,2] <- nbRun-1
        OFvalues[iterConfig,3] <- -.getDof(covariate = evaluatedCovariate)
        estimatedPopParam[[iterConfig]] <- mlx.getEstimatedPopulationParameters()
        # Get back to the initial state before decision
        eval(parse(text=paste0('lixoftConnectors::setCovariateModel(',evaluatedParameter,' = list(', evaluatedCovariate, ' = TRUE))')))
        if(settings$rankedSCM == F){endofline='\n'}else{endofline=''}
        lineDisplay <- paste0("Run ",toString(OFvalues[iterConfig,2]),": Evaluation of removing covariate ", evaluatedCovariate, ' from parameter ',evaluatedParameter,' ==> ', criteriaToDisplay,' = ', format(as.double(OFvalues[iterConfig,1]), nsmall = 1),', ', criteriaToDisplay,'-difference = ',format(as.double(referenceOFV-OFvalues[iterConfig,1]), nsmall = 1),endofline);
        summary <- c(summary, lineDisplay);cat(lineDisplay);cat(summary, file = summary.file);  
      }
      
      # Get the best model and its associated pValue
      if(settings$criteria == "LRT"){
        pValueAll <- .getPvalueLRT(referenceOFV, OFV = as.double(OFvalues[,1]), dof = as.double(OFvalues[,3]))
        iterMin <- which.min(pValueAll)
        indexMin <- idConfig[iterMin]
        pValue <- pValueAll[iterMin]
      }else{
        iterMin <- which.min(as.double(OFvalues[,1]));
        indexMin <- idConfig[iterMin]
        if(as.double(OFvalues[iterMin,1])<referenceOFV+settings$criteriaThreshold){
          pValue <- 1
        }else{
          pValue <- 0
        }
      }
      
      if(pValue >= settings$pElimination){
        # Add the parameter and list
        evaluatedParameter <- remainingCovariateStructure$indivParam[indexMin]; evaluatedCovariate <- remainingCovariateStructure$covariate[indexMin]
        eval(parse(text=paste0('lixoftConnectors::setCovariateModel(',evaluatedParameter,' = list(', evaluatedCovariate, ' = FALSE))')))
        referenceOFV <- OFvalues[iterMin,1]
        if(settings$criteria == "LRT" & settings$rankedSCM==T){pValueDisplay <- paste0(', pVal = ',format(pValue, digits = 3),'\n')}else{pValueDisplay<-''}
        lineDisplay <- paste0(pValueDisplay,"=> covariate ", evaluatedCovariate, ' was removed from parameter ', evaluatedParameter,' (new ref. ',criteriaToDisplay ,' = ',format(referenceOFV, nsmall = 1),') \n',additionalDisplay);
        summary <- c(summary, lineDisplay); cat(lineDisplay);cat(summary, file = summary.file)
        if(settings$updateInit){
          newInitialConditions = estimatedPopParam[[iterMin]]
          for(indexParam in 1:length(newInitialConditions)){
            eval(parse(text=paste0('lixoftConnectors::setPopulationParameterInformation(',names(newInitialConditions)[indexParam],' = list(initialValue = ',newInitialConditions[indexParam],'))')))
          }
        }
        # Update the pValue of the remaining pairs
        remainingCovariateStructure <- .updateRemainingPvalue(structure = remainingCovariateStructure[-indexMin,])
      }else{
        mlx.setPopulationParameterInformation(initialEstimates);
        if(settings$rankedSCM){
          if(settings$criteria == "LRT" & settings$rankedSCM==T){pValueDisplay <- paste0(', pVal = ',format(pValue, digits = 3),'\n')}else{pValueDisplay<-''}
          remainingCovariateStructure <- remainingCovariateStructure[-indexMin,]
          lineDisplay <- paste0(pValueDisplay,"=> covariate ",evaluatedCovariate, ' was kept on parameter ', evaluatedParameter,' (ref. ',criteriaToDisplay ,' = ',format(referenceOFV, nsmall = 1), ') \n',additionalDisplay);
          summary <- c(summary, lineDisplay); cat(lineDisplay);cat(summary, file = summary.file)
        }else{
          lineDisplay <- paste0("+++++++++++++++++++++++\n No covariate to remove\n +++++++++++++++++++++++\n")
          summary <- c(summary, lineDisplay); cat(lineDisplay);cat(summary, file = summary.file)
          bBackward <- FALSE
        }
      }
    }else{
      bBackward <- FALSE
    }
  }
  
  # Make the summary
  mlx.saveProject(projectFile = final.project)
  OFValue = .getOFV(indexLL); nbRun = nbRun+1;
  summary <- c(summary, paste0("========================================================\nFINAL MODEL (",criteriaToDisplay," = ",format(OFValue, nsmall = 1),")\n"))
  summary <- c(summary, paste0("--> target parameters: ",toString(unique(relationshipsToTest$indivParam)),"\n --> searched covariates: ",toString(unique(relationshipsToTest$covariate))," \n\n"))
  backwardList <- .getCurrentCovariateStructure(relationshipsToTest)
  if(!is.null(backwardList$covariate)){# There are covariates 
    for(indexList in 1:length(backwardList[,1])){
      lineDisplay <- paste0('Covariate ',backwardList$covariate[indexList], ' on parameter ',backwardList$indivParam[indexList],'\n');summary <- c(summary, lineDisplay)
    }
  }
  summary <- c(summary, c(paste0("\n => Done with ",toString(nbRun)," runs in ",toString(floor(proc.time()[3] - t_strat[3])),"s\n", date(),'\n'),"========================================================\n\n"))
  cat(summary); cat(summary, file = summary.file)
}  

###################################################################################
# Check the inputs 
###################################################################################
.checkCovariateSearchInput = function(inputName, inputValue){
  isValid = TRUE
  inputName = tolower(inputName)
  if(inputName == tolower("paramToUse")){
    if(is.vector(inputValue) == FALSE){
      message("ERROR: Unexpected type encountered. paramToUse must be a vector")
      isValid = FALSE
    }else if(prod(is.element(el = inputValue, set =mlx.getIndividualParameterModel()$name))==0){
      message("ERROR: paramToUse has at least one non-valid parameter name in its definition.")
      isValid = FALSE
    }
  }else if(inputName == tolower("covToTest")){
    if(is.vector(inputValue) == FALSE){
      message("ERROR: Unexpected type encountered. covToTest must be a vector")
      isValid = FALSE
    }else if(prod(is.element(el = inputValue, set = mlx.getCovariateInformation()$name))==0){
      message("ERROR: covToTest has at least one non-valid covariate name in its definition.")
      isValid = FALSE
    }
  }else if(inputName == tolower("testRelations")){
    if(is.list(inputValue) == FALSE){
      message("ERROR: Unexpected type encountered. testRelations must be a list")
      isValid = FALSE
    }else{
      for(indexList in 1:length(inputValue)){
        # Check the name
        if(!is.element(el = names(inputValue)[indexList],set = mlx.getIndividualParameterModel()$name)){
          message(paste0("ERROR: in testRelations, ", names(inputValue)[indexList], " is not a valid parameter name."))
          isValid = FALSE
          }
        # Check the values
          if(prod(is.element(el = inputValue[[indexList]],set = mlx.getCovariateInformation()$name))==0){
            message(paste0("ERROR: in testRelations, some elements of (",toString(inputValue[[indexList]]), ") are not valid covariates."))
            isValid = FALSE
          }
      }
    } 
  }else if(inputName == tolower("method")){
    if(is.character(inputValue) == FALSE){
      message("ERROR: Unexpected type encountered. method must be a string")
      isValid = FALSE
    }else if(prod(is.element(el = inputValue, set = c('COSSAC','SCM')))==0){
      message("ERROR: method should be either 'COSSAC' or 'SCM'.")
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
  }else if(settingName == tolower("criteriaThreshold")){
    if(is.double(settingValue) == FALSE){
      message("ERROR: Unexpected type encountered. criteriaThreshold must be a double.")
      isValid = FALSE
    }else{
      if(settingValue<0){
        message("ERROR: criteriaThreshold must be  positive.")
        isValid = FALSE
      }
    }
  }else  if(settingName == tolower("linearization")){
    if(is.logical(settingValue) == FALSE){
      message("ERROR: Unexpected type encountered. linearization must be a boolean.")
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
    }else if(length(intersect(tolower(settingValue), c(tolower('BIC'), tolower('LRT'), tolower('AIC'))))==0){
      message("ERROR: criteria must be either 'AIC', BIC' or 'LRT'.")
      isValid = FALSE
    }
  }else if(settingName == tolower("direction")){
    if(is.character(settingValue) == FALSE){
      message("ERROR: Unexpected type encountered. direction must be a string")
      isValid = FALSE
    }else if(length(intersect(tolower(settingValue), c(tolower('both'), tolower('forward'), tolower('backward'))))==0){
      message("ERROR: direction must be either 'both', 'forward', or 'backward'.")
      isValid = FALSE
    }
  }else if(settingName == tolower("updateInit")){
    if(is.logical(settingValue) == FALSE){
      message("ERROR: Unexpected type encountered. updateInit must be a boolean.")
      isValid = FALSE
    }
  }else if(settingName == tolower("saveRun")){
    if(is.logical(settingValue) == FALSE){
      message("ERROR: Unexpected type encountered. saveRun must be a boolean.")
      isValid = FALSE
    }
  }else{
    message("ERROR: ",settingName,' is not a valid setting')
    isValid = FALSE
  }
  return(isValid)
}

######################################################################################################################
# Get current covariate structure
# The output is a data.frame with indivParam, covariate and pValue
######################################################################################################################
.getCurrentCovariateStructure <- function(relationshipsToTest){
  param <- cov <- pValue <-NULL
  if(length(relationshipsToTest$indivParam)>0){
    for(indexParam in 1:length(relationshipsToTest$indivParam)){
      hasCov <- eval(parse(text=paste0('which(mlx.getIndividualParameterModel()$covariateModel$',relationshipsToTest$indivParam[indexParam],'==TRUE)')))
      if(length(hasCov)>0){
        covOnIndiv <- intersect(relationshipsToTest$covariate[indexParam], names(hasCov))
        if(length(covOnIndiv)>0){
          cov <- c(cov, covOnIndiv)
          param <- c(param, rep(as.character(relationshipsToTest$indivParam[indexParam]),length(covOnIndiv)))
          pValue <- c(pValue, rep(1,length(covOnIndiv)))
        }
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
    currentCovariateStructure <- data.frame(indivParam = param, covariate = cov, pValue=pValue)
  }
  
  return(currentCovariateStructure)
}

######################################################################################################################
# Get full covariate structure
# The output is a data.frame with indivParam, covariate and pValue 
######################################################################################################################
.getCovariateStructure <- function(relationshipsToTest){
  dfOut <- data.frame(indivParam = relationshipsToTest$indivParam, 
                      covariate = relationshipsToTest$covariate, 
                      pValue = rep(.025, length(relationshipsToTest$indivParam)))
  
  indexFinal <- 1:length(dfOut[,3])
  
  testValues <- .getTestsRandomEffects() 
  if(!is.null(testValues)){
    for(indexDF in 1:length(dfOut[,1])){
      indexTest <- which((as.character(testValues$eta)==paste0('eta_',dfOut$indivParam[indexDF]))&(as.character(testValues$covariate)==as.character(dfOut$covariate[indexDF])))
      if(length(indexTest)>0){
        if(testValues[indexTest, 4]=='<1e-16'){
          dfOut$pValue[indexDF] <- 1e-16
        }else{
          dfOut$pValue[indexDF] <- as.double(as.character(testValues[indexTest, 4]))
        }
      }
    }
  }
  
  return(dfOut[indexFinal,])
}

######################################################################################################################
# Get the remaining covariate structure
# The output is a data.frame with indivParam, covariate and pValue 
######################################################################################################################
.getRemainingCovariateStructure<- function(relationshipsToTest){
  currentCovariateStructure <- .getCurrentCovariateStructure(relationshipsToTest)
  
  remainingCovariateStructure <- .getCovariateStructure(relationshipsToTest)
  if(length(currentCovariateStructure[,1])>0){# We remove all the current covariate
    for(index in 1:length(currentCovariateStructure[,1])){
      indexLine <- which((as.character(remainingCovariateStructure$indivParam)==as.character(currentCovariateStructure$indivParam[index]))&(as.character(remainingCovariateStructure$covariate)==as.character(currentCovariateStructure$covariate[index])))
      remainingCovariateStructure <- remainingCovariateStructure[-indexLine,]
    }
  }
  return(remainingCovariateStructure)
}

######################################################################################################################
# Update the pValue of covariate structure
# The output is a data.frame with indivParam, covariate and pValue 
######################################################################################################################
.updateRemainingPvalue <- function(structure){
  structureOut <- structure
  
  testValues <- .getTestsRandomEffects() 
  if(!is.null(testValues)){
    for(indexDF in 1:length(structureOut[,1])){
      indexTest <- which((as.character(testValues$eta)==paste0('eta_',structureOut$indivParam[indexDF]))&(as.character(testValues$covariate)==as.character(structureOut$covariate[indexDF])))
      if(length(indexTest)>0){
        if(testValues[indexTest, 4]=='<1e-16'){
          structureOut$pValue[indexDF]  <- 1e-16
        }else{
          structureOut$pValue[indexDF]  <- as.double(as.character(testValues[indexTest, 4]))
        }
      }
    }
  }
  
  return(structureOut)
}

######################################################################################################################
# Get the tests associated to the random effects w.r.t. the covariates
######################################################################################################################
.getTestsRandomEffects <- function(){
  testFile <- paste0(mlx.getProjectSettings()$directory,'/Tests/correlationRandomEffectsCovariates.txt')
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
  testFile <- paste0(mlx.getProjectSettings()$directory,'/Tests/correlationIndividualParametersCovariates.txt')
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
  pValue = dof
  for(index in 1:length(dof)){
    if(dof[index]>0){
      pValue[index] <- 1-pchisq(max(referenceOFV-OFV[index],.1), df=dof[index])
    }else{
      pValue[index] <- 1-pchisq(max(-(referenceOFV-OFV[index]),.1), df=-dof[index])
    }
  }
  return(pValue)
}

#############################################################################################################################
# Get the scenario for the run
#############################################################################################################################
.defineScenario <- function(linearization, useTests){
  # Define the scenario associated to the type of test and the method
  scenario <- mlx.getScenario()
  if(linearization){
    scenario$linearization = T
  }else{
    scenario$linearization = F 
  }
  scenario$plotList <- "covariatemodeldiagnosis";
  
  if(useTests){# compute both conditional distribution and mode
    scenario$tasks <- c(populationParameterEstimation = T, conditionalDistributionSampling=T, conditionalModeEstimation = T, logLikelihoodEstimation = T, plots = T);
  }else{
    if(linearization){# compute only mode
      scenario$tasks <- c(populationParameterEstimation = T, conditionalModeEstimation = T, logLikelihoodEstimation = T)
    }else{# compute only mode
      scenario$tasks <- c(populationParameterEstimation = T, conditionalDistributionSampling = T, logLikelihoodEstimation = T)
    }
  }
  mlx.setScenario(scenario)
}

#############################################################################################################################
# Run the scenario and get the OFV
#############################################################################################################################
.getOFV <- function(indexLL){
  bScenario <- mlx.runScenario();
  
  if(bScenario){# The run is ok
    OFValue = mlx.getEstimatedLogLikelihood()[[1]][indexLL]
  }else{
    OFValue <- Inf
  }
  return(OFValue)
}

#############################################################################################################################
# Get the number of degree of freedom associated to a covariate
#############################################################################################################################
.getDof <- function(covariate){
  indexCov <- NULL
  eval(parse(text=paste0('indexCov <- which(names(mlx.getCovariateInformation()$type)=="',covariate,'")')))
  if(length(intersect(mlx.getCovariateInformation()$type[indexCov],c("categorical","categoricaltransformed")))>0){
    eval(parse(text=paste0('dof <- length(unique(mlx.getCovariateInformation()$covariate$',covariate,'))-1')))
  }else{
    dof = 1
  }
}
