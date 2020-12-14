#' Generate statistics for Visual Predictive Check
#' Creates data to further plot VPC from monolix project
#' 
#' @param project (\emph{string}) Monolix project
#' @param time (\emph{string}) (\emph{optional}) in {`time`, `timeSinceLastDose`} (default `time`).
#' `timeSinceLastDose` only possible when there is an `amount` column in the dataset
#' @param obsName (\emph{string}) (\emph{optional}) observation name. By default the first observation is considered.
#' 
#' \strong{parameters startification}
#' @param stratSplit (\emph{vector(string)}) (\emph{optional}) vector of covariate names used to split vpc plot.
#' (no split by default).
#' @param stratFilter (\emph{list}) (\emph{optional}) list of covariate names used to filter vpc plot
#' names are covariate names and values are category name in case of categorical covariate
#' and group id in case of continuous covariate.
#' (no filter by default).
#' @param stratScale (\emph{list}) (\emph{optional}) list of continuous covariate scaling
#' Names are continuous covariate names and values are vector of break doubles.
#' If not scaling defined, by default break is defined at the median.
#'
#' \strong{parameters for vpc computation}
#' @param level (\emph{int}) (\emph{optional}) level for prediction intervals computation (default 90).
#' @param higherPercentile (\emph{int}) (\emph{optional}) (\emph{continuous data}) (integer)
#' Higher percentile for empirical and predicted percentiles computation (default 90). For continuous data only.
#' @param useCorrpred (\emph{bool}) (\emph{optional}) (\emph{continuous data}) if TRUE, pcVPC are computed using Uppsala prediction correction (default FALSE).
#' For continuous data only.
#' @param useCensored (\emph{bool}) (\emph{optional}) (\emph{continuous data}) Choose to use BLQ data or to ignore it to compute the VPC (default TRUE).
#' For continuous data only.
#' @param censoring (\emph{string}) (\emph{optional}) (\emph{continuous data}) in {`simulated`, `loq`}, BLQ data can be simulated, or can be equal to the limit of quantification (LOQ) (default `simulated`). For continuous data only.
#' @param event.averageNumberEvents (\emph{bool}) (\emph{optional}) (\emph{event data}) If TRUE compute average number events in addition to survival curve (default FALSE)
#' If you do not need it, leave it to FALSE, it will speed up computations
#' 
#' \strong{parameters for x-axis  bins}
#' @param xBinsSettings (\emph{binSettingsClass}) (\emph{optional}) (\emph{continuous and discrete data}) a binSettingsClass object \link{getBinsSettings}.
#' a list of settings for time axis binning.
#' @param yBinsSettings (\emph{binSettingsClass}) (\emph{optional}) (\emph{countable discrete data}) a binsSettingsClass object \link{getBinsSettings}.
#' a list of settings for y axis binning.
#' @return A list:
#' \itemize{
#'   \item \code{vpcPercentiles}: a dataframe with vpc percentiles (bins, empirical and theorical vpc) 
#'   \item \code{observations}: a dataframe with observations
#'   \item \code{obsName}: Name of observation
#'   \item \code{timeName}: Name of time
#' }
#' @export
#' @examples
#' \dontrun{
#'   stratSplit <- c("sex", "wt")
#'   stratFilter = list("age" = 1)
#'   stratScale = list(age = c(65), wt = c(70))
#'   xBinsSettings <- getBinsSettings(is.fixedBins = TRUE, fixedBins = c(0, 24, 48, 72, 96, 120))
#'   vpcStats(project="RsmlxDemo1.mlxtran", stratSplit=stratSplit,
#'            stratFilter=stratFilter, stratScale=stratScale)
#'   vpcStats(project="RsmlxDemo1.mlxtran", xBinsSettings=xBinsSettings)
#'   vpcStats(project="RsmlxDemo1.mlxtran", time = "timeSinceLastDose")
#' }
#' 
#' @seealso \link{vpc} \link{plotVpc}
vpcStats <- function(project, time = "time", obsName = NULL, 
                     stratSplit = NULL, stratFilter = NULL, stratScale = NULL,
                     level = 90, higherPercentile = 90,
                     useCorrpred = FALSE, useCensored = TRUE, censoring = "simulated",
                     event.averageNumberEvents = FALSE,
                     xBinsSettings = NULL, yBinsSettings = NULL) {

  ## Check Arguments -----------------------------------------------------------
  params <- as.list(match.call(expand.dots = TRUE))[-1]

  if (missing(project)) stop("`project` is missing")
  # load project
  r <- .loadProject(project)
  project <- r$project
  
  time <- .check_time(time, "time")
  obsName <- .check_obs(obsName, "obsName")
  obsType <- .getObsType(obsName)
  dataType <- obsType$type
  dataSubType <- obsType$subType
  
  if ("censoring" %in% names(params) & ! "useCensored"  %in% names(params))
    useCensored <- TRUE
  .check_bool(useCorrpred, "useCorrpred")
  .check_bool(useCensored, "useCensored")
  if (useCensored)
    censoring <- .check_cens(censoring, "censoring")

  if (dataType == "event")
    .check_bool(event.averageNumberEvents, "event.averageNumberEvents")
  
  .check_perc(level, "level")
  .check_perc(higherPercentile, "higherPercentile")

  # bins settings
  xBinsSettings <- .check_bins(xBinsSettings, "xBinsSettings")
  yBinsSettings <- .check_bins(yBinsSettings, "yBinsSettings")

  # Read Data ------------------------------------------------------------------
  simName <- paste0("sim_", obsName)
  subjocc = .getSubjocc()
  
  # check Population parameters
  if (length(mlx.getEstimatedPopulationParameters()) == 0) {
    # Run population parameter estimation
    message("[INFO] Run population parameter estimation.")
    mlx.runPopulationParameterEstimation()
  }
  
  # get obs data
  obsData <- .getObservationData(obsName, dataType)
  simData <- .getSimulationData(obsName)
  obsData <- obsData[, !names(obsData) %in% c("split", "filter", "color")]
  simData <- simData[, !names(simData) %in% c("split", "filter", "color")]
  
  if (dataType == "continuous") {
    # censored data
    if (length(unique(obsData$censored)) > 1) {
      if (useCensored & censoring == "simulated") {
        obsName <- paste0(obsName, "_simBlq")
      }
    }
    else {
      useCensored <- FALSE
    }
  }
  
  if (dataType == "event") obsData <- .addTTECensoredColumn(obsData, obsName)
  
  # compute time since last dose if needed
  if (time == "timeSinceLastDose") {
    obsData <- .addDoseRelTime(obsData)
    simData <- .addDoseRelTime(simData)
    timeName <- "time_rel"
  } else {
    timeName <- "time"
  }

  # Stratify: split and filter -------------------------------------------------
  covStrat <- stratifyCovariates(split = stratSplit, filter = stratFilter, scaling = stratScale)
    
  if (!is.null(covStrat)) {
    obsData <- merge(obsData, subset(covStrat, select = c(subjocc, "split", "filter")), by = subjocc, sort = F)
    simData <- merge(simData, subset(covStrat, select = c(subjocc, "split", "filter")), by = subjocc, sort = F)
  } else {
    obsData$split <- rep("All", nrow(obsData))
    obsData$filter <- rep(FALSE, nrow(obsData))
    simData$split <- rep("All", nrow(simData))
    simData$filter <- rep(FALSE, nrow(simData))
  }
  filter <- FALSE
  obsData <- subset(obsData, filter == FALSE)
  simData <- subset(simData, filter == FALSE)
  
  # Bins -----------------------------------------------------------------------
  obsDataCensored <- obsData
  simDataCensored <- simData
  if (dataType == "continuous") {
    if (is.element("censored", names(obsData)) & !useCensored) {
      censored <- FALSE
      obsDataCensored <- subset(obsData, censored == 0)
      simDataCensored <- subset(simData, censored == 0)
    }
  }
  
  # x-axis bins (only if datatype is continuous or discrete)
  if (is.element(dataType, c("continuous", "discrete"))) {
    xBins <- computeVpcBins(obsDataCensored[[timeName]], split = obsDataCensored$split,
                            type = "continuous", binsSettings = xBinsSettings)
  } else {
    xBins <- NULL
  }
  obsDataCensored <- .addBinsIndex(obsDataCensored, timeName, xBins, "binIndex")
  simDataCensored <- .addBinsIndex(simDataCensored, timeName, xBins, "binIndex")
  
  # y-axis bins (only if datatype is discrete)
  if (dataType == "discrete") {
    type <- ifelse(dataSubType == "categorical", "categorical", "continuous")
    yBins <- computeVpcBins(obsDataCensored[[obsName]], split = obsDataCensored$split,
                            type = type, binsSettings = yBinsSettings,
                            is.binsName = TRUE, dataName = obsName)
    obsDataCensored <- .addBinsIndex(obsDataCensored, obsName, yBins, "category")
    simDataCensored <- .addBinsIndex(simDataCensored, simName, yBins, "category")
  } else {
    yBins <- NULL
  }

  # Corrected Prediction -------------------------------------------------------
  if (useCorrpred) {
    res <- .applyCorrectedPrediction(obsDataCensored, simDataCensored, obsName, simName)
    obsDataCensored <- res$obs
    simDataCensored <- res$sim
    obsName <- paste0(obsName, "_pc")
    simName <- paste0(simName, "_pc")
  }
  
  # Statistics -----------------------------------------------------------------
  # Empirical Data statistics per bin
  if (dataType == "continuous") {
    vpcData <- do.call(rbind, by(
      simDataCensored,
      simDataCensored$split,
      function(simSplit) {
        s <- unique(simSplit$split)
        obsSplit <- subset(obsDataCensored, split == s, select = c(obsName, timeName, "binIndex"))
        simSplit <- subset(simSplit, select = c("rep", simName, timeName, "binIndex"))
        res <- computeContinuousVPC(
          obsSplit, obsName, simSplit, simName,
          higherPercentile, level
        )
        res <- cbind(res, list(split = s))
      }
    ))
    vpcData <- merge(xBins, vpcData, by=c("bins", "split"))
    vpcData <- vpcData[order(vpcData$split, vpcData$bins_middles),]
    
  } else if (dataType == "discrete") {
    categories <- unique(obsDataCensored$category)
    vpcData <- do.call(rbind,by(
      subset(simDataCensored, select = c("rep", simName, "split", timeName, "binIndex", "category")),
      simData$split,
      function(simSplit) {
        s <- unique(simSplit$split)
        obsSplit <- subset(obsDataCensored, split == s, select = c(obsName, "split", timeName, "binIndex", "category"))
        res <- computeDiscreteVPC(
          obsSplit, obsName, simSplit, simName,
          higherPercentile, level, categories
        )
        res <- cbind(res, list(split = s))
      }
    ))
    vpcData <- merge(xBins, vpcData, by=c("bins", "split"))
    vpcData <- vpcData[order(vpcData$split, vpcData$category, vpcData$bins_middles),]
    vpcData$category <- sapply(vpcData$category, function(cat) yBins$binsName[yBins$bins == cat])
    
  } else if (dataType == "event") {
    vpcData <- do.call(rbind, by(
      simData,
      simData$split,
      function(simSplit) {
        s <- unique(simSplit$split)
        obsSplit <- subset(obsData, split == s)
        res <- computeEventVPC(
          obsSplit, obsName, simSplit, simName, timeName, dataSubType,
          xBinsSettings$nbDataPoints, level, event.averageNumberEvents
        )
        res <- cbind(list(split = s), res)
      }
    ))
  }

  result <- list(
    vpcPercentiles = vpcData,
    observations = obsData,
    obsName = obsName,
    timeName = timeName
  )
  return(result)
}

## Check -----------------------------------------------------------------------
.check_perc <- function(perc, argname) {
  if (!is.numeric(perc) | perc < 0 | perc > 100)
    stop("`", argname, "` must be in [0, 100]", call. = FALSE)
  return(perc)
}

.check_cens <- function(cens, argname) {
  if (is.null(cens)) cens <- "simulated"
  .check_in_vector(cens, argname, c("simulated", "blq"))
  return(cens)
}

.check_time <- function(time, argname) {
  if (is.null(time)) time <- "time"
  .check_in_vector(time, argname, c("time", "timeSinceLastDose"))
  # check if dose in dataset
  if (time == "timeSinceLastDose" & !is.element("amount", mlx.getData()$headerTypes))
    stop("Unexpected value encountered. `timeSinceLastDose` only valid when ",
         "dosing data is available in the dataset.", call. = FALSE)
  return(time)
}

.check_obs <- function(obs, argname) {
  obsnames <- mlx.getObservationInformation()$name
  if (is.null(obs)) {
    obs <- obsnames[1]
    if (length(obsnames) > 1)
      warning("`", argname, "` not specified. `", argname, "` set to ", obs, ".",
              call. = FALSE)
  }
  .check_in_vector(obs, argname, obsnames)
  return(obs)
}


## Covariates ------------------------------------------------------------------

#' Generate stratification by covariates
#'
#' @param split (\emph{vector(string)}) (\emph{optional}) vector of covariate names used to split dataset (no split by default).
#' @param filter (\emph{list}) (\emph{optional}) list of covariate names used to filter vpc plot
#' names are covariate names; values are category name in case of categorical covariate
#' and group id in case of continuous covariate.
#' (no filter by default).
#' @param scaling (\emph{list}) (\emph{optional}) list of continuous covariate scaling
#' Names are continuous covariate names; values are vector of break doubles.
#' If not scaling defined, by default break is defined at the median.
#' 
#' @return covariate dataframe with additional filter and split columns
stratifyCovariates <- function(split = c(), filter = list(), scaling = list()) {
  ## Get covariates information ------------------------------------------------
  covList <- mlx.getCovariateInformation()
  covData <- covList$covariate
  covtypes <- covList$type
  covNames <- c(covList$name, mlx.getData()$header[mlx.getData()$headerTypes == "occ"])
  if (length(.matchToData("occ")) > 0) {
    covtypes[.matchToData("occ")] <- "categorical"
  }
  contCovariates <- covList$name[covtypes == "continuous"]
  catCovariates <- covList$name[covtypes == "categorical"]

  ## Check Arguments -----------------------------------------------------------
  if (!is.null(split) & !is.vector(split)) stop("`split` must be a vector")
  if (!all(is.element(split, covNames))) {
    split <- split[split %in% covNames]
    warning(paste0("Covariates for split stratification filtered in", paste(covNames, collapse = ", ")))
  }
  if (!is.null(filter) & !is.list(filter)) stop("`filter` must be a list")
  if (!all(is.element(names(filter), covNames))) {
    filter <- filter[names(filter) %in% covNames]
    warning(paste0("Covariates for filter stratification filtered in", paste(covNames, collapse = ", ")))
  }
  for (i in seq_len(length(filter))) {
    name <- names(filter)[i]
    value <- filter[[name]]
    if (name %in% catCovariates & !value %in% unique(covData[[name]])) {
      stop(paste0("Unexepect value for filter `", name, "` . You must ",
                  "choose in (", paste(unique(covData[name]), collapse = ", "), ")"
      ))
    } else if (name %in% contCovariates & !value > 0 & !as.integer(value) == value) {
      stop(paste0(
        "Unexepect value for filter. Group id associated to ",
        "continous covariates `", name, "` must be a positive integer"
      ))
      isValid = FALSE
    }
  }
  if (!is.null(scaling) & !is.list(scaling)) stop("`scaling` must be a list")
  if (!all(is.element(names(scaling), contCovariates))) {
    scaling <- scaling[names(scaling) %in% contCovariates]
    warning(paste0("Covariates for scaling stratification filtered in", paste(contCovariates, collapse = ", ")))
  }
  
  if (!is.null(covList)) {
    # split column
    covData <- .addStratificationColumn(
      covData, covtypes,
      stratifiers = split,
      breaks = scaling,
      columnName = "split"
    )
    # filter column
    covData <- .addStratificationColumn(
      covData, covtypes,
      stratifiers = names(filter),
      breaks = scaling,
      filtering = TRUE,
      filterGroups = filter,
      columnName = "filter"
    )
    covData <- .addStratificationColumn(
      covData, covtypes,
      stratifiers = names(filter),
      breaks = scaling,
      columnName = "filterNames"
    )
    filterName <- unique(covData$filterNames[covData$filter == 0])
    
    # add filtering name to split column
    covData$split <- sapply(covData$split, function(s) paste0(s, filterName))
    covData$split <- sapply(covData$split, function(s) ifelse(s == "", "All", s))
  }
  return(covData)
}

# Read data --------------------------------------------------------------------
.getObservationData <- function(obsname, obsType) {
  if (obsType %in% c("discrete", "event")) {
    obsData <- mlx.getObservationInformation()[[obsname]]
  } else {
    obsFilename <- paste0(
      mlx.getProjectSettings()$directory, "/ChartsData/VisualPredictiveCheck/",
      obsname, "_observations.txt"
    )
    if (! file.exists(obsFilename)) {
      # Get VPC chart data files
      message("[INFO] Download Charts data.")
      s = mlx.getScenario()
      s$plotList = c(s$plotList, "vpc")
      mlx.setScenario(s)
      mlx.computeChartsData(exportVPCSimulations = TRUE)
    }
    obsData <- .renameColumns(.readDataset(obsFilename), "ID", "id")
  }
  return(obsData)
}

.getSimulationData <- function(obsname) {
  simFilename <- paste0(
    mlx.getProjectSettings()$directory, "/ChartsData/VisualPredictiveCheck/",
    obsname, "_simulations.txt"
  )
  if (! file.exists(simFilename)) {
    # Get VPC chart data files
    message("[INFO] Download Charts data.")
    s = mlx.getScenario()
    s$plotList = c(s$plotList, "vpc")
    mlx.setScenario(s)
    mlx.computeChartsData(exportVPCSimulations = TRUE)
  }
  simData <- .renameColumns(.readDataset(simFilename), "ID", "id")
  return(simData)
}

# Transform Times --------------------------------------------------------------
.addDoseRelTime <- function(data) {
  subjocc <- .getSubjocc()
  doseDf <- getDoseInformation()
  doseDfID <- .renameColumns(aggregate(doseDf$time, by = as.list(subset(doseDf, select = subjocc)), c), "x", "dose")
  data <- merge(data, doseDfID, by = subjocc, sort = F)
  data$time_rel <- apply(
    subset(data, select = c("dose", "time")),
    1,
    function(row) {
      if ("dose" %in% names(row)) {
        row <- as.list(row)
        return(min((row$time - row$dose)[row$time - row$dose >= 0]))
      } else {
        doses <- unname(row[names(row) != "time"])
        t <- row[["time"]]
        return(min((t - doses)[t - doses >= 0]))
      }
    }
  )
  data <- data[!names(data) %in% c("dose")]
  return(data)
}
