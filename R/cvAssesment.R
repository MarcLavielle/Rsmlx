#' Convergence assesment function
#'
#' Allow to test differente initial conditions and different seeds.
#' @param project project name
#' @param type type of assessment. [Optional] It defines the scanario of the assesment. It can be 'FAST', 'LINEARIZED', or 'EXTENDED'. The default value is 'EXTENDED'
#' @param nbReplicates number of replicates. [Optional] The deualt value is 5
#' @param initialParameters data.frame defining the initial parameters for the assessment. [Optional] for the user that want to define it by himself.
#' @param initialParametersMinMax named list of initial parameters min and max for the assessment. [Optional]  for the user that want to define it by himself.
#' @examples
#' \dontrun{
#' cvAssesment(project) => run the assessment for the project with default values
#' cvAssesment(project, nbReplicates  = 10) => run the assessment for the project with 10 replicates
#' cvAssesment(project, type = 'FAST') => run the assessment with only SAEM, no standard error nor Log-likelihood will be computed
#' out <- cvAssesment(project, initialParameters = data.frame(ka_pop = seq(from=1, to = 5, by=1), V_pop_ = seq(from=1, to = 5, by=1)))
#' => run the assessment with only ka_pop and V_pop (if population parameters) as varying initial conditions. The initial conditions are defined in the data frame
#' out <- cvAssesment(project, initialParametersMinMax = list(ka_pop = c(min=1, max = 3), Cl_pop = c(min=.1, max = .5)) )
#' => run the assessment with only ka_pop and Cl_pop (if population parameters) as varying initial conditions. The initial conditions are uniformly drawn from the associated min and max
#' }
#' @export
cvAssesment<-function(project, type=NULL, nbReplicates = NULL, initialParameters = NULL, initialParametersMinMax = NULL){

  library(MlxConnectors)

  # Define the type of assessment, EXTENDED by default
  if(is.null(type)){
    type <- 'EXTENDED'
  }
  if(!.checkCvAssessmentInput(inputName = "type", inputValue = type)){return(invisible(FALSE))}

  # Define the number of replicates, 5 by default
  if(is.null(nbReplicates)){
    nbReplicates <- 5
  }
  if(!.checkCvAssessmentInput(inputName = "nbReplicates", inputValue = nbReplicates)){return(invisible(FALSE))}

  ###############################################################################################
  # Get the information of the project
  ###############################################################################################
  status <- loadProject(project)
  # Get the export directory
  expDirectory <- getExportSettings('Directory')
  # Define the scenario
  if((type=="LINEARIZATION")&(length(which(getAvailableTasks()$ll$methods=="Linearization"))==0)){type <- 'EXTENDED'}
  scenario <- switch (type,
                      "FAST" =  'saem = list( run = TRUE, methods = \"saem\" )',
                      "LINEARIZATION" = 'saem = list(run = TRUE,methods = \"saem\"), indivestim = list(run = TRUE, methods =\"indivestim_mode\"), fisher = list(run = TRUE, methods = \"fisher_lin\"), ll = list( run = TRUE, methods = \"ll_lin\")',
                      "EXTENDED"= 'saem = list(run = TRUE, methods = \"saem\" ), indivestim = list(run = TRUE, methods =\"indivestim_mean\"), fisher = list( run = TRUE, methods = \"fisher_sa\"), ll = list( run = TRUE, methods = \"ll_is\")'
  )

  ###############################################################################
  # Draw the replicates
  ###############################################################################
  # Get the population parameters
  populationParameter <- unlist(getAvailablePopulationParameters())[c(TRUE,FALSE)]
  if(!is.null(initialParameters)){
    if(!.checkCvAssessmentInput(inputName = "initialParameters", inputValue = initialParameters)){return(invisible(FALSE))}
    # The user defines all the initial parameters
    paramCvAss <- intersect(x=names(initialParameters), y = populationParameter)
    popParam_rep <- initialParameters[,match(paramCvAss, names(initialParameters))]
    if(is.vector(popParam_rep)){# COnvert it into a data frame
      popParam_rep <- data.frame(matrix(ncol = 1, nrow = length(popParam_rep), data = popParam_rep))
      colnames(popParam_rep) <- paramCvAss
    }
    nbReplicates <- length(popParam_rep[,1])
  }else{
    if(!is.null(initialParametersMinMax)){
      # The user defined the parameters along with the min and the max
      if(!.checkCvAssessmentInput(inputName = "initialParametersMinMax", inputValue = initialParametersMinMax)){return(invisible(FALSE))}
      indexAssessedParameters <- intersect(x = names(initialParametersMinMax), y = populationParameter)
      popParam_rep <- data.frame(matrix(ncol = length(indexAssessedParameters), nrow = nbReplicates))
      colnames(popParam_rep) <- indexAssessedParameters
      for(indexParam in 1:length(indexAssessedParameters)){
        paramBounds <- initialParametersMinMax[[match(names(initialParametersMinMax),indexAssessedParameters[indexParam])]]
        popParam_rep[,indexParam] <- runif(nbReplicates, min = min(paramBounds), max = max(paramBounds))
      }
    }else{
      isFixed <- unlist(getAvailablePopulationParameters())[ c(FALSE,TRUE) ]
      indexAssessedParameters <- intersect(which(isFixed==FALSE),grep(pattern="_pop", x = populationParameter))
      popParam_rep <- data.frame(matrix(ncol = length(indexAssessedParameters), nrow = nbReplicates))
      colnames(popParam_rep) <- populationParameter[indexAssessedParameters]
      for(indexParam in 1:length(indexAssessedParameters)){
        paramName <- populationParameter[indexAssessedParameters[indexParam]]
        offset <- runif(nbReplicates, min = -.5, max = .5)
        popParam_rep[,indexParam] <- getValueTransformedReferential(paramName, paramValue =  as.double(getPopulationParametersInitialValue(paramName)), offset)
      }
    }
  }

  ###############################################################################
  # Run the replicates
  ###############################################################################
  populationParameters <- names(getPopulationParametersInitialValue())
  # Removed the fixed ones
  isFixed <- which(unlist(getAvailablePopulationParameters())[ c(FALSE,TRUE) ]==TRUE)
  if(length(isFixed)>0){
    fixedPopulationParameters <- unlist(getAvailablePopulationParameters())[ c(TRUE,FALSE) ][isFixed]
    populationParameters <- populationParameters[populationParameters != fixedPopulationParameters]
  }
  nbPopParam <- length(populationParameters)
  assessmentPop <- data.frame(matrix(ncol = nbPopParam, nrow = nbReplicates)); colnames(assessmentPop)<- populationParameters
  assessmentStdErrors <- data.frame(matrix(ncol = nbPopParam, nrow = nbReplicates)); colnames(assessmentStdErrors)<- populationParameters
  assessmentLL <- data.frame(vector(length = nbReplicates)); colnames(assessmentLL)<- '-2LL'
  cvReport <- list()
  for(indexReplicate in 1:nbReplicates){
    cat(paste0("\nWorking on replicate n ",toString(indexReplicate),'\n'))
    exportDirectory_CvAss <- paste0(expDirectory,'cvAssesment_rep',toString(indexReplicate))
    # Set the scenario and the initial values
    eval(parse(text=paste0('setScenario(',scenario,')')))
    setRandomSeed()
    for(indexParam in 1:length(names(popParam_rep))){
      eval(parse(text=paste0('setPopulationParameterInitialValue(',names(popParam_rep)[indexParam],' = popParam_rep[indexReplicate,indexParam])')))
    }
    # Set the export directory
    setExportSettings(directory=exportDirectory_CvAss)
    runScenario(TRUE)
    # Get the result
    indexOrder <- match(populationParameters, names(getPopulationParametersValue()))
    assessmentPop[indexReplicate,] <- getPopulationParametersValue()[indexOrder]
    if(!(type=='FAST')){
      fisherValue <- switch (type,  "LINEARIZATION" = getStandardErrors()$fisher_lin, "EXTENDED" = getStandardErrors()$fisher_sa)
      fisherValue[sapply(fisherValue, is.null)] <- NA
      assessmentStdErrors[indexReplicate,] <- fisherValue
      assessmentLL[indexReplicate,1] <- switch (type,  "LINEARIZATION" = -2*getLogLikelihoodValues()$ll_lin[1], "EXTENDED" = -2*getLogLikelihoodValues()$ll_is[1])
    }
    cvReport[[indexReplicate]] <- getSAEMConvergenceReport()$convergencereport


    ###############################################################################
    # Prepare informations for the plot
    ###############################################################################
    # Define the type of plot
    if(type=="FAST"){displayLL <- FALSE; displaySA <- FALSE;}else{displayLL <- TRUE; displaySA <- TRUE;}
    #Define the commun x axis for the output graphics
    xVal <- 1:indexReplicate

    ###############################################################################
    # Plot the output figures
    # Where all the population parameters converge
    ###############################################################################
    nbFig <- nbPopParam+displayLL
    x_NbFig <- ceiling(max(sqrt(nbFig),1)); y_NbFig <- ceiling(nbFig/x_NbFig)
    par(mfrow = c(x_NbFig, y_NbFig), oma = c(0, 3, 1, 1), mar = c(3, 1, 0, 3), mgp = c(2, 1, 0), xpd = NA)
    for(indexFigure in 1:nbPopParam){
      paramName <- populationParameters[indexFigure]
      paramValue <- assessmentPop[xVal,indexFigure]
      if(displaySA){paramSA <- assessmentStdErrors[xVal,indexFigure]}else{paramSA<-rep(NaN, nbReplicates)}
      yVal_max <- getValueTransformedReferential(paramName, paramValue, offset = paramSA)
      yVal_min <- getValueTransformedReferential(paramName, paramValue, offset = -paramSA)
      marginDislplay <-  .001*mean(abs(paramValue))
      plot(xVal, paramValue, type="b", col="blue", ylab = paramName, ylim = c(min(yVal_min)-marginDislplay,max(yVal_max)+marginDislplay),xlab = "")
      if(!(type=='FAST')){
        lines(xVal, yVal_max, type="b", col="green", pch = 2); lines(xVal, yVal_min, type="b", col="green", pch = 2)
      }
    }
    if(displayLL){
      plot(xVal, assessmentLL[xVal,1], type="b", col="blue", ylab = "-2LL", ylim = c(min(assessmentLL[xVal,1])-1, max(assessmentLL[xVal,1])+1),xlab = "")
    }

    ####################################################################################
    # Plot the varying curves
    ####################################################################################
    nbFig <- length(cvReport[[1]])
    x_NbFig <- ceiling(max(sqrt(nbFig),1)); y_NbFig <- ceiling(nbFig/x_NbFig)
    par(mfrow = c(x_NbFig, y_NbFig),  oma = c(0, 3, 1, 1), mar = c(3, 1, 0, 3), mgp = c(2, 1, 0), xpd = NA)
    for(indexFigure in 1:nbFig){
      # get the bounds of the plot
      ymin <- 1e10; ymax <- -1e10; xlim <- 0; yF <- NULL
      for(idReplicate in 1:indexReplicate){
        xlim <- max(xlim,length(cvReport[[idReplicate]][1][[1]]))
        yVal <- cvReport[[idReplicate]][indexFigure][[1]]
        yF <- c(yF,yVal[length(yVal)])
        ymin <- min(ymin,min(yVal[which(!is.infinite(yVal))])); ymax <- max(ymax,max(yVal[which(!is.infinite(yVal))]))
      }
      ymin <- ymin-0.02*mean(abs(yF)); ymax <- ymax+0.02*mean(abs(yF))
      # Make the plot
      for(idReplicate in 1:indexReplicate){
        xVal <- seq(1,length(cvReport[[idReplicate]][1][[1]]))
        yVal <- cvReport[[idReplicate]][indexFigure][[1]]
        ylabel <- names(cvReport[[idReplicate]])[indexFigure]; ylabel <- gsub(ylabel, pattern = "-2LL", replacement = "Complete -2LL")
        couleur <- palette()[idReplicate+1]
        if(idReplicate==1){
          plot(xVal,yVal,type="b",col=couleur,ylab = ylabel,xlab = "",xlim = c(0,xlim), ylim = c(ymin,ymax))
        }else{
          lines(xVal,yVal,type="b",col=couleur)
        }
      }
      lines(1:xlim,mean(yF)*(1+0.05)+0*(1:xlim),type="l", pch=22, lty=1, col="black"); lines(1:xlim,mean(yF)*(1-0.05)+0*(1:xlim),type="l", pch=22, lty=1, col="black")
    }
  }
  return(list(assessmentPop,assessmentStdErrors,assessmentLL,cvReport))
}

.checkCvAssessmentInput= function(inputName, inputValue){
  isValid = TRUE
  inputName = tolower(inputName)
  if(inputName == "nbreplicates"){
    if((is.double(inputValue) == FALSE)&&(is.integer(inputValue) == FALSE)){
      message("ERROR: Unexpected type encountered. The number of replicates must be an integer.")
      isValid = FALSE
    }else{
      if(!(as.integer(inputValue) == inputValue)){
        message("ERROR: Unexpected type encountered. The number of replicates must be an integer.")
        isValid = FALSE
      }else if(inputValue<1){
        message("ERROR:the number of replicates must be a strictly positive integer.")
        isValid = FALSE
      }
    }
  }else if(inputName == "type"){
    if(is.character(inputValue) == FALSE){
      message("ERROR: Unexpected type encountered. The type must be a character.")
      isValid = FALSE
    }else{
      if(is.na(match(x = tolower(inputValue), table = c('fast','linearization','extended')))){
        message("ERROR: the type must be a character in the following list: 'FAST', 'LINEARIZATION', 'EXTENDED'.")
        isValid = FALSE
      }
    }
  }else if(inputName == "initialparameters"){
    if(is.data.frame(inputValue) == FALSE){
      message("ERROR: Unexpected type encountered. initialParameters must be a data frame.")
      isValid = FALSE
    }else if(length(intersect(x = names(getPopulationParametersInitialValue()), y=names(inputValue)))==0){
      message("ERROR: initialParameters have no population parameters in its definition.")
      isValid = FALSE
    }
  }else if(inputName == "initialparametersminmax"){
    if(is.list(inputValue) == FALSE){
      message("ERROR: Unexpected type encountered. initialParametersMinMax must be a list.")
      isValid = FALSE
    }else{
      for(indexList in 1:length(inputValue)){
        inputList <- inputValue[[indexList]]
        if(length(intersect(c("min","max"), names(inputList)))<2){
          message(paste0("ERROR: In initialParametersMinMix, min and max must be defined for list n?",toString(indexList)))
          isValid = FALSE
        }
      }
      if(isValid){
        if(length(intersect(names(getPopulationParametersInitialValue()), names(inputValue)))==0){
          message("ERROR: initialParametersMinMax have no population parameters in its definition.")
          isValid = FALSE
        }
      }

    }
  }
  return(invisible(isValid))
}

