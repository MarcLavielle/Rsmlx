#' Function that sample data set
#'
#' Sample data set from a project
#' @param project project where the data set will be sampled
#' @param dataFolder [optional] folder where to get the data set
#' @param [optional] settings of the sampling. These settings are
#' \itemize{
#' \item samples, the number of bootstrapped data set to generate (default value is 100)
#' \item sample_size, the number of subjects in each bootstrap data set (default value is the  number of individuals in the original data set).
#' \item regenerateData, boolean to regenerate or not the data sets (default value is F).
#' \item stratify_on, categorical covariate (or tranformed categorical covariate) of the project (default is NULL). It allows to respect the proportion of these covariate in the generated data set
#' }
#'
#' @export
generateBootstrap = function(project, settings, dataFolder = NULL){

  # define the scenario in order to only have SAEM
  setScenario(tasks =  c(populationParameterEstimation = T))

  # Prepare all the output folders
  exportDir <- getProjectSettings()$directory
  projectName <- substr(basename(project), 1, nchar(basename(project))-8)
  dir.create(file.path(exportDir, 'bootstrap/'), showWarnings = FALSE, recursive = T)

  # Get the data set information
  referenceDataset <- getData()
  cov <- getCovariateInformation()
  datasetFile <- referenceDataset$dataFile

  if(is.null(dataFolder)){
    dir.create(file.path(exportDir, 'bootstrap/data/'), showWarnings = FALSE)

    # Load the data set
    dataset <- read.table(file=datasetFile, header = T, sep = "", dec = ".")
    if(length(dataset[1,])==1){dataset <- read.table(file=datasetFile, header = T, sep = ",", dec = ".")}
    if(length(dataset[1,])==1){dataset <- read.table(file=datasetFile, header = T, sep = ";", dec = ".")}
    if(length(dataset[1,])==1){dataset <- read.table(file=datasetFile, header = T, sep = "\t", dec = ".")}

    indexID <- which(referenceDataset$headerTypes=="id")
    nameID <- unique(dataset[, indexID])
    nbIndiv <- length(nameID)
    if(is.na(settings$sample_size)){settings$sample_size = nbIndiv}

    validID <- list()
    if(is.na(settings$stratify_on)){
      nbCAT = 1
      indexPropCAT <- 1
      propCAT <- rep(settings$sample_size, nbCAT)
      validID[[indexPropCAT]] <- nameID
    }else{
      indexCAT <- which(cov$name == settings$stratify_on)
      catValues <- cov$covariate[,indexCAT+1]# +1 comes from the ID column
      nameCAT <- unique(catValues)
      nbCAT <- length(nameCAT)
      propCAT <- rep(settings$sample_size, nbCAT)
      validID <- list()
      for(indexPropCAT in 1:nbCAT){
        indexIdCat <- which(catValues==nameCAT[indexPropCAT])
        propCAT[indexPropCAT] <- max(1,floor(settings$sample_size*length(indexIdCat)/nbIndiv))
        validID[[indexPropCAT]] <- as.character(cov$covariate[indexIdCat,1])
      }
    }

    for(indexSample in 1:settings$samples){
      datasetFileName <- paste0(exportDir,'/bootstrap/data/dataset_',toString(indexSample),'.csv')
      if(!file.exists(datasetFileName)){
        ##################################################################################################################
        # Generate the data set
        ##################################################################################################################
        # Sample the IDs
        sampleIDs <- NULL
        for(indexValidID in 1:length(validID)){
          sampleIDs <- c(sampleIDs,  sample(x = validID[[indexValidID]], size = propCAT[indexValidID], replace = T) )
        }
        # get the datas
        data <- NULL
        for(indexSampleSize in 1:length(sampleIDs)){
          indexLine <- which(dataset[,indexID]==sampleIDs[indexSampleSize])
          dataToInsert <- dataset[indexLine,]
          dataToInsert[,indexID] <- indexSampleSize
          data <- rbind(data, dataToInsert)
        }

        write.table(x = data, file = datasetFileName,
                    eol = '\n', quote = FALSE, dec = '.',  row.names = FALSE, col.names = TRUE )
      }
      ##################################################################################################################
      # Generate the project file
      ##################################################################################################################
      # set the data file and the export directory
      bootData <- referenceDataset
      bootData$dataFile <- datasetFileName
      setData(bootData)
      saveProject(projectFile = paste0(exportDir,'/bootstrap/',projectName,'_bootstrap_',toString(indexSample),'.mlxtran'))
    }
  }else{
    dataFiles <- list.files(path = dataFolder, pattern = '*.txt|*.csv')
    for(indexSample in 1:length(dataFiles)){
      bootData <- referenceDataset
      bootData$dataFile <- paste0(dataFolder,dataFiles[indexSample])
      setData(bootData)
      saveProject(projectFile = paste0(exportDir,'/bootstrap/',projectName,'_bootstrap_',toString(indexSample),'.mlxtran'))
    }
  }

}

##################################################################################################################
# Clean the bootstrap folder
##################################################################################################################
cleanBootstrap <- function(project){
  # Prepare all the output folders
  exportDir <- getProjectSettings()$directory
  listProjectsToDelete <-list.files(path = paste0(exportDir,'/bootstrap/'),pattern = '*.mlxtran')

  if(length(listProjectsToDelete)>0){
    for(indexProject in 1:length(listProjectsToDelete)){
      projectBoot <-  paste0(exportDir,'/bootstrap/',listProjectsToDelete[indexProject])
      loadProject(projectBoot)
      exportDirToClean <- getProjectSettings()$directory
      unlink(exportDirToClean, recursive = TRUE)
      unlink(projectBoot, recursive = FALSE)
    }
    unlink(file.path(exportDir, '/bootstrap/data/'), recursive = FALSE)
  }
}

##################################################################################################################
# Run the bootstrap
##################################################################################################################
runBootstrap <- function(project, dataFolder = NULL, settings = NULL){

  loadProject(project)
  # Check and initialize the settings
  if(!is.null(settings)){
    if(!.checkBootstrapInput(inputName = "settings", inputValue = settings)){return(invisible(FALSE))}
  }
  if(is.null(settings$samples)){ settings$samples <- 100 }
  if(is.null(settings$sample_size)){ settings$sample_size <- NA}
  if(is.null(settings$regenerateData)){ settings$regenerateData <- F}
  if(is.null(settings$stratify_on)){ settings$stratify_on <- NA}

  if(!is.null(dataFolder)){
    dataFiles <- list.files(path = dataFolder, pattern = '*.txt|*.csv')
    if(length(dataFiles)>0){
      settings$samples <- length(dataFiles)
      settings$regenerateData <- F
    }else{
      message("WARNING: Folder ",dataFolder,' does not exist or does not contain any data set')
      return(invisible(FALSE))
    }
  }

  # Prepare all the output folders
  param <- getPopulationParameterInformation()$name[which(getPopulationParameterInformation()$method!="FIXED")]
  exportDir <- getProjectSettings()$directory
  projectName <- substr(basename(project), 1, nchar(basename(project))-8)

  if(settings$regenerateData){
    cleanBootstrap(project)
  }
  cat("Generating projects with bootstrap data sets\n")
  if(is.null(dataFolder)){
    cat("Generating data sets too\n")
    generateBootstrap(project, settings)
  }else{
    generateBootstrap(project, settings, dataFolder)
  }

  paramResults <- array(dim = c(settings$samples, length(param)))
  for(indexSample in 1:settings$samples){
    projectBoot <-  paste0(exportDir,'/bootstrap/',projectName,'_bootstrap_',toString(indexSample),'.mlxtran')
    loadProject(projectBoot)
    cat(paste0('Project ',toString(indexSample),'/',toString(settings$samples)))

    # Check if the run was done
    if(!file.exists(paste0(getProjectSettings()$directory,'/populationParameters.txt'))){
      cat(' => Running SAEM \n')
      runScenario()
    }else{
      cat(' => already computed \n')
    }
    paramResults[indexSample,] <-  getEstimatedPopulationParameters();
  }
  colnames(paramResults) <- names(getEstimatedPopulationParameters())

  # Plot the results
  nbFig <- length(param)
  x_NbFig <- ceiling(max(sqrt(nbFig),1)); y_NbFig <- ceiling(nbFig/x_NbFig)
  par(mfrow = c(x_NbFig, y_NbFig), oma = c(0, 3, 1, 1), mar = c(3, 1, 0, 3), mgp = c(1, 1, 0), xpd = NA)
  for(indexFigure in 1:nbFig){
    res <- paramResults[,indexFigure]
    resQ <- quantile(res,c(.05,.95))
    bxp <- boxplot(res, xlab = paste0(colnames(paramResults)[indexFigure],'\n 90% CI: [',toString(round(resQ[1],3)),', ',toString(round(resQ[2],3)),']'))
  }
  write.table(x = paramResults, file = paste0(exportDir,'/bootstrap/',projectName,'bootstrapResults.txt'),
              eol = "\n", sep = ",", col.names = T, quote = F, row.names = F)
  return(paramResults)
}

###################################################################################
# Check the inputs
###################################################################################
.checkBootstrapInput = function(inputName, inputValue){
  isValid = TRUE
  inputName = tolower(inputName)
  if(inputName == tolower("settings")){
    if(is.list(inputValue) == FALSE){
      message("ERROR: Unexpected type encountered. settings must be a list")
      isValid = FALSE
    }else {
      for (i in 1:length(inputValue)){
        if(!.checkBootstrapSettings(settingName = names(inputValue)[i], settingValue = inputValue[[i]])){
          isValid = FALSE
        }
      }
    }
  }
  return(invisible(isValid))
}

.checkBootstrapSettings = function(settingName, settingValue){
  isValid = TRUE
  settingName = tolower(settingName)
  if(settingName == tolower("samples")){
    if((is.double(settingValue) == FALSE)&&(is.integer(settingValue) == FALSE)){
      message("ERROR: Unexpected type encountered. samples must be an integer.")
      isValid = FALSE
    }else{
      if(!(as.integer(settingValue) == settingValue)){
        message("ERROR: Unexpected type encountered. samples must be an integer.")
        isValid = FALSE
      }else if(settingValue<1){
        message("ERROR: samples must be a strictly positive integer.")
        isValid = FALSE
      }
    }
  }else if(settingName == tolower("sample_size")){
    if((is.double(settingValue) == FALSE)&&(is.integer(settingValue) == FALSE)){
      message("ERROR: Unexpected type encountered. sample_size must be an integer.")
      isValid = FALSE
    }else{
      if(!(as.integer(settingValue) == settingValue)){
        message("ERROR: Unexpected type encountered. sample_size must be an integer.")
        isValid = FALSE
      }else if(settingValue<1){
        message("ERROR: samples must be a strictly positive integer.")
        isValid = FALSE
      }
    }
  }else if(settingName == tolower("regenerateData")){
    if((is.logical(settingValue) == FALSE)){
      message("ERROR: Unexpected type encountered. regenerateData must be an boolean")
      isValid = FALSE
    }
  }else if(settingName == tolower("stratify_on")){
    if(is.vector(inputValue) == FALSE){
      message("ERROR: Unexpected type encountered. covariateToSearch must be a vector")
      isValid = FALSE
    }else{
      if(length(intersect(getCovariateInformation()$name, inputValue))==0){
        message(paste0("ERROR: ",inputValue," is not a valid covariate of the project."))
        isValid = FALSE
      }else{
        indexCAT <- which(getCovariateInformation()$name==inputValue)
        catType <- getCovariateInformation()$type[indexCAT[1]]
        if(!((catType=="categorical")||(catType=="categoricaltransformed"))){
          message(paste0("ERROR: ",inputValue," is not a categorical covariate."))
          isValid = FALSE
        }
      }
    }
  }else{
    message("WARNING: ",settingName,' is not a valid setting')
  }
  return(isValid)
}
