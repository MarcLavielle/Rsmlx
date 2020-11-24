#' Generate statistics for Visual Predictive Check
#' Creates data to further plot VPC from monolix project
#' 
#' @param project Monolix project
#' @param time emph{[Optional]} in {`time`, `timeSinceLastDose`} (default `time`).
#' `timeSinceLastDose` only possible when there is an `amount` column in the dataset
#' @param obsName emph{[Optional]} observation name. By default the first observation is considered.
#' 
#' \strong{parameters startification}
#' @param stratSplit emph{[Optional]} vector of covariate names used to split vpc plot.
#' (no split by default).
#' @param stratFilter emph{[Optional]} list of covariate names used to filter vpc plot
#' names are covariate names and values are category name in case of categorical covariate
#' and group id in case of continuous covariate.
#' (no filter by default).
#' @param stratScale emph{[Optional]} list of continuous covariate scaling
#' Names are continuous covariate names and values are vector of break doubles.
#' If not scaling defined, by default break is defined at the median.
#'
#' \strong{parameters for vpc computation}
#' @param level emph{[Optional]} (integer) level for prediction intervals computation (default 90).
#' @param higherPercentile emph{[Optional]}[continuous data] (integer)
#' Higher percentile for empirical and predicted percentiles computation (default 90). For continuous data only.
#' @param useCorrpred emph{[Optional]}[continuous data] (boolean) if TRUE, pcVPC are computed using Uppsala prediction correction (default FALSE). For continuous data only.
#' @param useCensored emph{[Optional]}[continuous data] (boolean) Choose to use BLQ data or to ignore it to compute the VPC (default TRUE). For continuous data only.
#' @param censoring emph{[Optional]}[continuous data] in {`simulated`, `loq`}, BLQ data can be simulated, or can be equal to the limit of quantification (LOQ) (default `simulated`). For continuous data only.
#' @param event.averageNumberEvents emph{[Optional]}[event data] (boolean) If TRUE compute average number events in addition to survival curve (default FALSE)
#' If you do not need it, leave it to FALSE, it will speed up computations
#' 
#' \strong{parameters for x-axis  bins}
#' @param xBinsSettings emph{[Optional]}[continuous and discrete data] a binSettingsClass object \link{binsSettings}.
#' a list of settings for time axis binning.
#' @param yBinsSettings emph{[Optional]}[countable discrete data] a binsSettingsClass object \link{binsSettings}.
#' a list of settings for y axis binning.
#' @export
#' @examples
#' \dontrun{
#'   stratSplit <- c("sex", "wt")
#'   stratFilter = list("age" = 1)
#'   stratScale = list(age = c(65), wt = c(70))
#'   xBinsSettings <- getBinsSettings(is.fixedBins = TRUE, fixedBins = c(0, 24, 48, 72, 96, 120))
#'   vpcStats(project="RsmlxDemo1.mlxtran", stratSplit=stratSplit, stratFilter=stratFilter, stratScale=stratScale)
#'   vpcStats(project="RsmlxDemo1.mlxtran", xBinsSettings=xBinsSettings)
#'   vpcStats(project="RsmlxDemo1.mlxtran", time = "timeSinceLastDose")
#' }
#' 
#' @seealso \link{vpc} \link{plotVPC}
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
  r <- prcheck(project)
  project <- r$project
  
  if (!is.element(time, c("time", "timeSinceLastDose"))) {
    time <- "time"
    warning("Invalid argument `time`. `time` set to \"time\"")
  }
  if (time == "timeSinceLastDose" & !is.element("amount", mlx.getData()$headerTypes))
    stop(paste0(
      "Unexpected value encountered. `timeSinceLastDose` not available when no 
      dosing data available in the dataset"
    ))

  obsnames <- mlx.getObservationInformation()$name
  if (is.null(obsName)) obsName <- obsnames[1]
  if (!is.element(obsName, obsnames)) {
    obsName <- obsnames[1]
    warning(paste0("Invalid argument `obsName`. `obsname` set to \"", obsName, "\""))
  }
  obsType <- .getObsType(obsName)
  dataType <- obsType$type
  dataSubType <- obsType$subType
  
  if ("censoring" %in% names(params) & ! "useCensored"  %in% names(params)) useCensored <- TRUE
  if (!is.logical(useCorrpred)) stop("`useCorrpred` must be logical.")
  useCorrpred <- useCorrpred[1]
  if (!is.logical(useCensored)) stop("`useCensored` must be logical.")
  useCensored <- useCensored[1]
  if (useCensored & !is.element(censoring, c("simulated", "blq"))) {
    censoring <- "simulated"
    warning(paste0("Invalid argument `censoring`. `censoring` set to \"", censoring, "\""))
  }
  
  if (dataType == "event" & !is.logical(event.averageNumberEvents))
    stop("`event.averageNumberEvents` must be logical.")
  
  if (!is.double(level) | level < 0 | level > 100)
    stop("`level` must be in [0, 100]")
  if (!is.double(higherPercentile) | higherPercentile < 0 | higherPercentile > 100)
    stop("`higherPercentile` must be in [0, 100]")
  
  
  # bins settings
  if (is.null(xBinsSettings)) xBinsSettings <- getBinsSettings()
  if(!class(xBinsSettings) == "binSettingsClass")
    stop("xBinsSettings must be a binSettingsClass object")
  
  if (is.null(yBinsSettings)) yBinsSettings <- getBinsSettings()
  if(!class(yBinsSettings) == "binSettingsClass")
    stop("yBinsSettings must be a binSettingsClass object")

  # Read Data ------------------------------------------------------------------
  simName <- paste0("sim_", obsName)
  subjocc = .getSubjocc()
  
  # check Population parameters
  if (length(mlx.getEstimatedPopulationParameters()) == 0) {
    # Run population parameter estimation
    message("Run population parameter estimation.")
    mlx.runPopulationParameterEstimation()
  }
  
  # get obs data
  obsData <- .getObservationData(obsName, dataType)
  simData <- .getSimulationData(obsName)
  obsData <- obsData[, !names(obsData) %in% c("split", "filter", "color")]
  simData <- simData[, !names(simData) %in% c("split", "filter", "color")]
  
  if (dataType == "continuous") {
    if (useCensored & length(unique(obsData$censored)) == 1)
      useCensored <- FALSE
    if (useCensored & !is.null(censoring) & length(unique(obsData$censored)) > 1){
      loq <- max(obsData[obsName][obsData$censored == 1,])
      simData$censored <- simData[simName] < loq
      if (censoring == "simulated")
        obsName <- paste0(obsName, "_simBlq")
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
  obsData <- subset(obsData, filter == FALSE)
  simData <- subset(simData, filter == FALSE)
  
  # Bins -----------------------------------------------------------------------
  obsDataCensored <- obsData
  simDataCensored <- simData
  if (dataType == "continuous") {
    if (is.element("censored", names(obsData)) & !useCensored) {
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
  } else {
    yBins <- NULL
  }
  obsDataCensored <- .addBinsIndex(obsDataCensored, obsName, yBins, "category")
  simDataCensored <- .addBinsIndex(simDataCensored, simName, yBins, "category")

  # Corrected Prediction -------------------------------------------------------
  if (useCorrpred) {
    res <- .applyCorrectedPrediction(obsDataCensored, simDataCensored)
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
    vpcData <- vpcData[order(vpcData$split, vpcData$bins_middle),]
    
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
    vpcData <- vpcData[order(vpcData$split, vpcData$category, vpcData$bins_middle),]
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

## Bins -------------------------------------------------------------------------

#' Bins settings
#'
#' @param is.fixedBins emph{[Optional]} (boolean) If TRUE, specify manually bin list (default FALSE).
#' @param fixedBins emph{[Optional]} (list[integer]) Define mannually a list of bins. To use when `useFixedBins` is set to TRUE.
#' @param criteria emph{[Optional]} Choose the bining criteria among `equalwidth`, `equalsize` or `leastsquare` (default leastsquare).
#' @param is.fixedNbBins emph{[Optional]} (boolean) If TRUE define a fixed number of bins, else define a range for automatic selection (default FALSE).
#' @param nbBins emph{[Optional]} (integer) Define a fixed number of bins (default 10).
#' @param binRange emph{[Optional]} (list[integer]) Define a range for the number of bins (default c(5, 20)).
#' @param nbBinData emph{[Optional]} (list[integer]) Define a range for the number of data points per bin (default c(10, 200)).
#' @param nbDataPoints emph{[Optional]}[event data] (integer) Number of data point in event time grid (default 100)
#' @return binsSettingClass object
#' @export
#' @examples
#' \dontrun{
#'   getBinsSettings(is.fixedBins = TRUE, fixedBins = c(1, 2, 3))
#' }
getBinsSettings <- function(is.fixedBins = FALSE, fixedBins = c(),
                            criteria = "leastsquare", is.fixedNbBins = FALSE,
                            nbBins = 10, binRange = c(5, 20),
                            nbBinData = c(10, 200), nbDataPoints = 100) {
  ## Check Arguments -----------------------------------------------------------
  params <- as.list(match.call(expand.dots = TRUE))[-1]
  if ("fixedBins" %in% names(params) & ! "is.fixedBins"  %in% names(params)) is.fixedBins <- TRUE
  if ("nbBins" %in% names(params) & ! "is.fixedNbBins"  %in% names(params)) is.fixedNbBins <- TRUE
  
  if (!is.logical(is.fixedBins)) stop("`is.fixedBins` must be logical")
  if (is.fixedBins) {
    if(!is.vector(fixedBins) | !all(is.double(fixedBins)) | !all(fixedBins >= 0))
      stop("`fixedBins` must be a vector of positive doubles")
  }
  if (!is.element(criteria, c("equalwidth", "equalsize", "leastsquare")))
    stop("`criteria` must be in {`equalwidth`, `equalsize`, `leastsquare`}")
  if (!is.logical(is.fixedNbBins)) stop("`is.fixedNbBins` must be logical")
  if (!is.double(nbBins) | !nbBins > 0) stop("`nbBins` must a positive integer")
  if (!as.integer(nbBins) == nbBins) stop("`nbBins` must a positive integer")
  if (length(binRange) != 2 | !all(is.double(binRange)) | !all(binRange >= 0))
    stop("`binRange` must a range of two positive integer")
  if (!all(as.integer(binRange) == binRange))
    stop("`binRange` must a range of two positive integer")
  if (length(nbBinData) != 2 | !all(is.double(nbBinData)) | !all(nbBinData >= 0))
    stop("`nbBinData` must a range of two positive integer")
  if (!all(as.integer(nbBinData) == nbBinData))
    stop("`nbBinData` must a range of two positive integer")
  if (!is.double(nbDataPoints) | !nbDataPoints > 0) stop("`nbDataPoints` must a positive integer")
  if (!as.integer(nbDataPoints) == nbDataPoints) stop("`nbDataPoints` must a positive integer")
  
  structure(list(
    is.fixedBins = is.fixedBins[1], fixedBins = fixedBins,
    criteria = criteria[1], is.fixedNbBins = is.fixedNbBins[1],
    nbBins = nbBins[1], binRange = binRange,
    nbBinData = nbBinData, nbDataPoints = nbDataPoints
  ), class = "binSettingsClass"
  )
}

#' Generate Bins
#'
#' @param data (vector) vector of data (times or obs) to compute bins.
#' @param split (vector) split category associated with data (no split by default).
#' @param type in {"continuous", "categorical"} type of data ("continuous" by default).
#' @param binsSettings emph{[Optional]} (binSettingsClass) Bins settings.
#' @param is.binsName emph{[Optional]} (boolean) if TRUE add a column binsName in P(x<dataName<y) (default FALSE).
#' @param dataName emph{[Optional]} (list[integer]) dataName when is.binsName is TRUE.
#' @return bins dataframe
#' @export
#' @examples
#' \dontrun{
#'   computeVpcBins(seq_len(100), binsSettings = getBinsSettings(nbBins = 3))
#' }
computeVpcBins <- function(data, split = NULL, type = "continuous",
                           binsSettings = NULL,
                           is.binsName = FALSE, dataName = "obs") {
  ## Check Arguments -----------------------------------------------------------
  params <- as.list(match.call(expand.dots = TRUE))[-1]
  if ("dataName" %in% names(params) & ! "is.binsName"  %in% names(params)) is.binsName <- TRUE
  if (missing(data)) stop("`data` is missing.")
  if (!is.vector(data)) stop("`data` must be a vector of doubles")
  if (!all(is.double(data) | is.integer(data))) stop("`data` must be a vector of doubles")
  if (is.null(split)) split <- rep("All", length(data))
  if (length(data) != length(split)) stop("`data` and `split` must have the same length.")
  if (!is.element(type, c("continuous", "categorical"))) {
    warning("type must be either `continuous` or `categorical`. Set type to continuous.")
    type <- "continuous"
  }
  if (is.null(binsSettings)) binsSettings <- getBinsSettings()
  if (!class(binsSettings) == "binSettingsClass")
    stop("`binsSettings` must be a binSettingsClass object.")
  if (!is.logical(is.binsName)) {
    warning("is.binsName must be logical. Set is.binsName to FALSE.")
    is.binsName <- FALSE
  }
  if (!is.character(dataName)) {
    warning("dataName must be a string. Set dataName to `obs`.")
    dataName <- "obs"
  }
  
  params <- list(
    useFixedBins = binsSettings$is.fixedBins, fixedBins = binsSettings$fixedBins,
    binsSettings$criteria, useFixedNbBins = binsSettings$is.fixedNbBins,
    nbBins = binsSettings$nbBins, binRange = binsSettings$binRange,
    nbBinData = binsSettings$nbBinData
  )

  bins <- NULL
  for (s in unique(split)) {
    if (type == "continuous") {
      binsSplit <- .getBins(data = data[split == s], settings = params)
      bins <- rbind(bins, cbind(binsSplit, list(split = s)))
    } else if (type == "categorical") {
      binsSplit <- .getBins(data = data[split == s], bycat = TRUE, categories = unique(data))
      bins <- rbind(bins, cbind(binsSplit, list(split = s)))
    }
  }
  if (is.binsName) {
    bins$binsName <- bins$bins
    bins[bins$bins_start == bins$bins_middle,]$binsName <- sapply(
      bins[bins$bins_start == bins$bins_middle,]$bins_start,
      function(x) paste0("P(", dataName, "=", x, ")")
    )
    bins[bins$bins_start < bins$bins_middle,]$binsName <- mapply(
      function(x, y) paste0("P(", ceiling(x), "<=", dataName, "<=", trunc(y), ")"),
      bins[bins$bins_start < bins$bins_middle,]$bins_start,
      bins[bins$bins_start < bins$bins_middle,]$bins_stop
    )
  }
  return(bins)
}

.getBins <- function(data, settings = NULL, bycat = FALSE, categories = c()) {
  if (bycat) {
    if (length(categories) == 0) categories <- unique(data)
    binsData <- data.frame(list(
      bins_start = categories, bins_stop = categories,
      bins_middle = categories, bins = seq(1, length(categories))
    ))
  } else if (!is.null(settings)) {
    if (settings$useFixedBins){
      # rescale bins
      b <- unique(sort(c(
        min(data), pmin(pmax(settings$fixedBins, min(data)), max(data)), max(data)
      )))
      # bins data
      factor <- cut(data, breaks = b, include.lowest=TRUE)
      binsData <- data.frame(list(
        bins_start = head(b, -1),
        bins_stop = tail(b, -1),
        bins_middle = as.vector(tapply(data, factor, mean)),
        bins = seq(1, length(b) - 1)
      ))
    } else {
      # get bins
      options <- list(
        criteria = settings$criteria, usefixednb = settings$useFixedNbBins,
        fixednb = settings$fixedNbBins,
        estimatednb = settings$binRange, nbbindata = settings$nbBinData
      )
      binsList <- mlx.computeBins(data = data, options = options)
      binsData <- data.frame(list(
        bins_start = head(binsList$values, -1),
        bins_stop = tail(binsList$values, -1),
        bins_middle = binsList$middles,
        bins = seq(1, length(binsList$values) - 1)
      ))
    }
  } else {
    binsData <- NULL
  }
  return(binsData)
}

.addBinsIndex <- function(data, columnName, binsValue = NULL, binName = "binIndex") {
  if (!columnName %in% names(data)) stop(paste0("Invalid columnName ", columnName, "."))
  if (!is.null(binsValue)) {
    if (!"split" %in% names(data)) data <- cbind(data, list(split =  "All"))
    if (!"split" %in% names(binsValue)) binsValue <- cbind(binsValue, list(split = "All"))
    data <- do.call(rbind, c(by(
      data,
      data$split,
      function(df) {
        bins <- subset(binsValue, split == unique(df$split))
        if (all(bins$bins_start == bins$bins_stop)) {
          df[binName] = bins$bins[match(df[[columnName]], bins$bins_start)]
        } else {
          df[binName] <- cut(
            df[[columnName]],
            breaks = c(bins$bins_start, tail(bins$bins_stop, 1)),
            include.lowest=TRUE, labels = F
          )
        }
        return(df)
      }
    ), make.row.names = FALSE))
  }
  # } else {
  #   data[binName] <- 1
  # }
  return(data)
}

## Covariates ------------------------------------------------------------------

#' Generate stratification by covariates
#'
#' @param split emph{[Optional]} vector of covariate names used to split dataset (no split by default).
#' @param filter emph{[Optional]} list of covariate names used to filter vpc plot
#' names are covariate names; values are category name in case of categorical covariate
#' and group id in case of continuous covariate.
#' (no filter by default).
#' @param scaling emph{[Optional]} list of continuous covariate scaling
#' Names are continuous covariate names; values are vector of break doubles.
#' If not scaling defined, by default break is defined at the median.
#' 
#' @return covStratClass object
#' @export
#' @examples
#' \dontrun{
#'   
#' }
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
      message("Download Charts data.")
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
    message("Download Charts data.")
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
