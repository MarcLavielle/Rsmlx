#' Visual Predictive Check
#'
#' Generate Visual predictive check plots
#' @param project Monolix project 
#' @param stratifier [optional] a list of settings for to stratify vpc plots:
#' \itemize{
#' \item \code{covariateSplit} vector of covariate names used to split vpc plot 
#' (no split by default).
#' \item \code{covariateFilter} list of covariate names used to filter vpc plot
#' names are covariate names and values are category name in case of categorical covariate
#' and group id in case of continuous covariate.
#' (no filter by default).
#' \item \code{covariateScaling} list of continuous covariate scaling
#' Names are continuous covariate names and values are vector of break doubles.
#' If not scaling defined, by default break is defined at the median.
#' }
#' @param settings [optional] a list of settings for to stratify vpc plots:
#' \itemize{
#' \item \code{grid} boolean, if TRUE, grid is displayed (default TRUE).
#' \item \code{legend} boolean, if TRUE, legend is displayed (default FALSE).
#' \item \code{displayObservedData} boolean, if TRUE, observed data is displayed (default FALSE).
#' \item \code{displayCensoredData} boolean, if TRUE, censored data is displayed (default FALSE).
#' \item \code{displayEmpiricalPercentiles} boolean, if TRUE, empirical percentiles are displayed (default TRUE).
#' \item \code{displayPredictedPercentiles} boolean, if TRUE, theoretical percentiles are displayed (default FALSE).
#' \item \code{level} integer, level for prediction intervals computation (default 90).
#' \item \code{higherPercentile} integer, higher percentile for empirical and predicted percentiles computation (default 90).
#' \item \code{outliersDots} boolean, if TRUE, red dots indicating empirical percentiles that are outside prediction intervals are displayed (default TRUE).
#' \item \code{outlierAreas} boolean, if TRUE, red areas indicating empirical percentiles that are outside prediction intervals are displayed (default TRUE).
#' \item \code{correctedPredictions} boolean, if TRUE, pcVPC are computed using Uppsala prediction correction (default FALSE).
#' \item \code{linearInterpolation} boolean, if FALSE, Set piece wise display for prediction intervals (default TRUE).
#' \item \code{useCensoredData} boolean, Choose to use BLQ data or to ignore it to compute the VPC (default TRUE).
#' \item \code{censoredData} in {`simulated`, `loq`}, BLQ data can be simulated, or can be equal to the limit of quantification (LOQ) (default `simulated`).
#' \item \code{binLimits} boolean, if TRUE, vertical lines added to the scatter plots to indicate the bins (default FALSE).
#' \item \code{binCriteria} Bining criteria in {`manual`, `equalwidth`, `equalsize`, `leastsquare`} (default `leastsquare`).
#' \item \code{useFixedNbBins} boolean, fixed number of bins. To define if `binCriteria` in {`equalwidth`, `equalsize`, `leastsquare`} (default FALSE).
#' \item \code{fixedNbBins} integer, fixed number of bins. To define if `useFixedNbBins` is TRUE (default 10).
#' \item \code{estimatedNbBins} integer vector of size 2, Bin range for automatic selection. To define if `useFixedNbBins` is FALSE (default `c(5, 30)`).
#' \item \code{nbBinData} integer vector of size 2, range for the number of data points per bin. To define if `binCriteria` in {`equalwidth`, `equalsize`, `leastsquare`} (default `c(10, 200)`).
#' \item \code{manualBins} vector, Manually define bins boundaries. To define if `binCriteria` is `manual`.
#' \item \code{xlogScale} boolean, if TRUE, plot in x log scale (default FALSE).
#' \item \code{ylogScale} boolean, if TRUE, plot in y log scale (default FALSE).
#' \item \code{xlabel} in {`time`, `timeSinceLastDose`} (default `time`).
#' }
#' @param preferences [optional] boolean to define if plot preferences
#' @importFrom ggplot2 ggplot element_rect element_line geom_ribbon geom_point
#'             geom_rect geom_line theme aes xlab ylab facet_wrap facet_grid
#'             alpha scale_color_manual scale_shape_manual scale_fill_manual
#'             guide_legend
#' @export
#' @examples
#' \dontrun{
#'   stratifier <- list(covariateSplit = c("sex", "wt"),
#'                      covariateFilter = list("age" = 1),
#'                      covariateScaling = list("age" = c(65), "wt" = c(70)))
#'   settings <- list(binCriteria = "manual", manualBins = c(0, 24, 48, 72, 96, 120))
#'   plotVPC(project="RsmlxDemo1.mlxtran", stratifier=stratifier)
#'   plotVPC(project="RsmlxDemo1.mlxtran", settings=settings)
#'   plotVPC(project="RsmlxDemo1.mlxtran", settings=list(xlabel = "timeSinceLastDose"))
#' } 
plotVPC <- function(project, stratifier = NULL, settings = NULL, preferences = NULL) {
  # TO DO: What happen when PKPD data ???

  # Check Arguments ------------------------------------------------------------
  # load project
  r <- prcheck(project)
  project <- r$project

  # Check and initialize the settings
  if(!is.null(settings)){
    if(!.checkVPCInput(inputName = "settings", inputValue = settings)){return(invisible(FALSE))}
  }
  if (is.null(settings$legend)) settings$legend <- FALSE 
  if (is.null(settings$grid)) settings$grid <- TRUE 
  if (is.null(settings$displayObservedData)) settings$displayObservedData <- FALSE 
  if (is.null(settings$displayCensoredData)) settings$displayCensoredData <- FALSE 
  if (is.null(settings$displayEmpiricalPercentiles)) settings$displayEmpiricalPercentiles <- TRUE 
  if (is.null(settings$displayPredictedPercentiles)) settings$displayPredictedPercentiles <- FALSE 
  if (is.null(settings$displayPredictionInterval)) settings$displayPredictionInterval <- TRUE 
  if (is.null(settings$level)) settings$level <- 90 
  if (is.null(settings$higherPercentile)) settings$higherPercentile <- 90 
  if (is.null(settings$outliersDots)) settings$outliersDots <- TRUE 
  if (is.null(settings$outlierAreas)) settings$outlierAreas <- TRUE 
  if (is.null(settings$correctedPredictions)) settings$correctedPredictions <- FALSE 
  if (is.null(settings$linearInterpolation)) settings$linearInterpolation <- TRUE
  if (is.null(settings$useCensoredData)) settings$useCensoredData <- TRUE 
  if (is.null(settings$censoredData)) settings$censoredData <- "simulated" 
  if (is.null(settings$binLimits)) settings$binLimits <- FALSE
  if (is.null(settings$binCriteria)) {
    settings$binCriteria <- ifelse(!is.null(settings$manualBins), "manual", "leastsquare")
  }
  if (!settings$binCriteria == "manual") {
    if (is.null(settings$useFixedNbBins)) settings$useFixedNbBins <- FALSE
    if (settings$useFixedNbBins & is.null(settings$fixedNbBins)) settings$fixedNbBins <- 10
    if (!settings$useFixedNbBins & is.null(settings$estimatedNbBins)) settings$estimatedNbBins <- c(5, 30)
    if (is.null(settings$nbBinData)) settings$nbBinData <- c(10, 200)
  }
  if (is.null(settings$xlogScale)) settings$xlogScale <- FALSE 
  if (is.null(settings$xlabel)) settings$xlabel <- "time" 
  if (is.null(settings$ylogScale)) settings$ylogScale <- FALSE
  # if (is.null(settings$ylabel)) settings$ylabel <- "" TO DO !!!!! Find an example
  
  # Check and initialize the stratifier
  if(!is.null(stratifier)){
    if(!.checkVPCInput(inputName = "stratifier", inputValue = stratifier)) return(invisible(FALSE))
  }
  if (is.null(stratifier$covariateSplit)) stratifier$covariateSplit <- c() 
  if (is.null(stratifier$covariateFilter)) stratifier$covariateFilter <- c()
  if (is.null(stratifier$covariateScaling)) stratifier$covariateScaling <- c()

  # Check and initialize the preferences
  # TO DO !!!

  # Read Data ------------------------------------------------------------------
  # TO DO: What if PKPD model< 2 VPC ? 
  # TO DO: what if not model run ? (prediction does not exist)
  obsName <- mlx.getObservationInformation()$name
  ylabel <- obsName
  simName <- paste0("sim_", obsName)
  subjocc = .getSubjocc()
  obsFilename <- paste0(
    mlx.getProjectSettings()$directory, "/ChartsData/VisualPredictiveCheck/",
    obsName, "_observations.txt"
  )
  vpcFilename <- paste0(
    mlx.getProjectSettings()$directory, "/ChartsData/VisualPredictiveCheck/",
    obsName, "_simulations.txt"
  )
  predFilename <- paste0(mlx.getProjectSettings()$directory, "/predictions.txt")
  if (! file.exists(predFilename)) {
    # Run population parameter estimation
    message(paste0("Run population parameter estimation."))
    mlx.runPopulationParameterEstimation()
  }
  if (! file.exists(vpcFilename) | ! file.exists(obsFilename)) {
    # Get VPC simulation file
    message(paste0("VPC file ", vpcFilename, " not found, Download Charts data."))
    mlx.computeChartsData(exportVPCSimulations = TRUE)
  }
  simDf <- .renameColumns(.readDataset(vpcFilename), "ID", "id")
  obsDf <- .renameColumns(.readDataset(obsFilename), "ID", "id")
  predDf <- subset(.readDataset(predFilename), select = c(subjocc, "time", "popPred"))

  if (settings$censoredData == "simulated" & length(unique(obsDf$censored)) > 1){
    obsName <- paste0(obsName, "_simBlq")
  }
  simDf <- subset(simDf, select = c("rep", subjocc, "time", simName, "censored"))
  obsDf <- subset(obsDf, select = c(subjocc, "time", obsName, "censored"))

  # compute time since last dose if needed
  if (settings$xlabel == "timeSinceLastDose") {
    doseDf <- getDoseInformation()
    doseDfID <- .renameColumns(aggregate(doseDf$time, by = as.list(subset(doseDf, select = .getSubjocc())), c), "x", "dose")
    obsDf <- merge(obsDf, doseDfID, by = subjocc, sort = F)
    obsDf$time_rel <- apply(
      obsDf,
      1,
      function(row) {
        dose_times <- as.double(row["dose"])
        t <- as.double(row["time"])
        return(min((t - dose_times)[t - dose_times >= 0]))
      }
    )
    obsDf <- obsDf[!names(obsDf) %in% c("dose")]
  }
  
  # Stratify: split and filter -------------------------------------------------
  covariates <- mlx.getCovariateInformation()
  covariatesDf <- covariates$covariate

  # split
  splitCol = rep("", nrow(covariatesDf))
  if (!is.null(stratifier$covariateSplit)) {
    stratifier$covariateSplit <- unique(stratifier$covariateSplit)
    for (cov in stratifier$covariateSplit) {
      if (covariates$type[covariates$name == cov] == "continuous") {
        extremes <- c(min(covariatesDf[[cov]]), max(covariatesDf[[cov]]))
        breaks = sort(unique(c(extremes, stratifier$covariateScaling[[cov]])))
        breaks <- breaks[breaks >= extremes[1] & breaks <= extremes[2]]
        factor = cut(covariatesDf[[cov]], breaks = breaks, include.lowest=TRUE)
        levels(factor) <- paste0("#", cov, ": ", levels(factor), " ")
        splitCol <- paste0(splitCol, factor)
      } else {
        splitCol <- paste0(splitCol, paste0("#", cov, ": ", covariatesDf[[cov]], " "))
      }
    }
  }

  # filter
  filterCol = rep(FALSE, nrow(covariatesDf))
  if (!is.null(stratifier$covariateFilter)) {
    stratifier$covariateSplit <- unique(stratifier$covariateSplit)
    for (i in seq(1, length(stratifier$covariateFilter))) {
      cov <- names(stratifier$covariateFilter)[i]
      group <- stratifier$covariateFilter[i]
      if (covariates$type[covariates$name == cov] == "continuous") {
        extremes <- c(min(covariatesDf[[cov]]), max(covariatesDf[[cov]]))
        breaks = sort(unique(c(extremes, stratifier$covariateScaling[[cov]])))
        breaks <- breaks[breaks >= extremes[1] & breaks <= extremes[2]]
        factor = cut(covariatesDf[[cov]], breaks = breaks, include.lowest=TRUE, labels = F)
        filterCol <- filterCol | !(factor == group)
        factor = cut(covariatesDf[[cov]], breaks = breaks, include.lowest=TRUE)
        levels(factor) <- paste0("#", cov, ": ", levels(factor), " ")
        splitCol <- paste0(splitCol, factor)
      } else {
        filterCol <- filterCol | !(covariatesDf[[cov]] == group)
        splitCol <- paste0(splitCol, paste0("#", cov, ": ", covariatesDf[[cov]], " "))
      }
    }
  }
  covariatesDf$split <- splitCol
  covariatesDf$filter <- filterCol
  obsDf <- merge(obsDf, subset(covariatesDf, select = c(subjocc, "split", "filter")), by = subjocc, sort = F)
  simDf <- merge(simDf, subset(covariatesDf, select = c(subjocc, "split", "filter")), by = subjocc, sort = F)
  predDf <- merge(predDf, subset(covariatesDf, select = c(subjocc, "split", "filter")), by = subjocc, sort = F)
  obsDf <- subset(obsDf, filter == FALSE)
  simDf <- subset(simDf, filter == FALSE)
  predDf <- subset(predDf, filter == FALSE)

  # Bins -----------------------------------------------------------------------
  if (! settings$useCensoredData) {
    obsDf <- subset(obsDf, censored == 0)
    simDf <- subset(simDf, censored == 0)
  }
  bins = NULL
  for (s in unique(obsDf$split)) {
    tbins <- subset(obsDf, split == s)$time
    if (!settings$binCriteria == "manual") {
      # get bins
      options <- list(
        criteria = settings$binCriteria, usefixednb = settings$useFixedNbBins,
        fixednb = settings$fixedNbBins, estimatednb = settings$estimatedNbBins,
        nbbindata = settings$nbBinData
      )
      binsList <- mlx.computeBins(data = tbins, options=options)
      # bins data
      bins <- rbind(bins, data.frame(list(
        bins_start = head(binsList$values, -1), bins_stop = binsList$values[-1],
        bins_middle = binsList$middles,
        bins = levels(cut(tbins, breaks = binsList$values, include.lowest=TRUE)),
        split = s
      )))
    } else {
      # rescale bins
      b <- settings$manualBins
      b <- unique(sort(c(min(tbins), b, max(tbins))))
      b <- b[b >= min(tbins) & b <= max(tbins)]
      # bins data
      factor <- cut(tbins, breaks = b, include.lowest=TRUE)
      bins <- rbind(bins, data.frame(list(
        bins_start = head(b, -1), bins_stop = b[-1],
        bins_middle = as.vector(tapply(tbins, factor, mean)), bins = levels(factor),
        split = rep(s, length(levels(factor)))
      )))
    }
  }

  # Corrected Prediction -------------------------------------------------------
  if (settings$correctedPredictions) {
    predDf <- do.call(rbind, by(
      predDf,
      predDf$split,
      function(df) {
        b <- bins[bins$split == unique(df$split),]
        factor <- cut(df$time, breaks = c(b$bins_start[1], b$bins_stop), include.lowest=TRUE)
        predictBin <- tapply(df$popPred, factor, function(x) median(x))
        predBin <- sapply(factor, function(f) predictBin[[f]])
        df$pcNorm <- ifelse(df$popPred == 0, 1, predBin / df$popPred)
        return(df)
      }
    ))
    obsDf <- merge(obsDf, predDf, by = c(subjocc, "time", "split", "filter"), sort = F)
    obsDf[paste0(obsName, "_pc")] <- obsDf[obsName] * obsDf$pcNorm
    obsName <- paste0(obsName, "_pc")
    ylabel <- paste0("Prediction corrected ", obsName)
    simDf <- merge(simDf, predDf, by = c(subjocc, "time", "split", "filter"), sort = F)
    simDf[paste0(simName, "_pc")] <- simDf[simName] * simDf$pcNorm
    simName <- paste0(simName, "_pc")
  }

  # Statistics -----------------------------------------------------------------
  # Empirical Data statistics per bin
  higherPercentile <- settings$higherPercentile
  level <- settings$level
  empiricalStats <- do.call(rbind, by(
    subset(obsDf, select = c(obsName, "split", "time")),
    obsDf$split,
    function(df) {
      b <- bins[bins$split == unique(df$split),]
      factor <- cut(df$time, breaks = c(b$bins_start[1], b$bins_stop), include.lowest=TRUE)
      s <- as.data.frame(do.call(
        rbind, tapply(df[[obsName]], factor, .stats, higherPercentile, header = "empirical")
      ))
      s$split <- unique(df$split)
      s <- cbind(bins = rownames(s), s)
      return(s)
    }
  ))

  simStats <- do.call(rbind, by (
    subset(simDf, select = c("rep", simName, "split", "time")),
    subset(simDf, select = c("rep", "split")),
    function(df) {
      b <- bins[bins$split == unique(df$split),]
      factor <- cut(df$time, breaks = c(b$bins_start[1], b$bins_stop), include.lowest=TRUE)
      st <- as.data.frame(do.call(rbind, tapply(df[[simName]], factor, .stats, higherPercentile, header = "theoretical")))
      st["rep"] <- unique(df$rep)
      st$bins <- rownames(st)
      st["split"] <- unique(df$split)
      return(st)
    }
  ))
  theoreticalStats <- do.call(rbind, by (
    simStats,
    subset(simStats, select = c("bins", "split")),
    function(df) {
      rlt <- as.data.frame(sapply(df[, !names(df) %in% c("bins", "split", "rep")], .stats, (100 - level) / 2))
      rlt$st <- rownames(rlt)
      rlt$bins <- unique(df$bins)
      rlt["split"] <- unique(df$split)
      rlt <- reshape(rlt, idvar = c("bins", "split"), timevar = "st", direction = "wide")
      return(rlt)
    }
  ))
  names(theoreticalStats) <- sapply(
    names(theoreticalStats), function(x) paste(strsplit(x, "\\.")[[1]], collapse = "_")
  )

  percStats <- merge(empiricalStats, theoreticalStats, by = c("bins", "split"))
  percStats <- merge(percStats, bins, by = c("bins", "split"))
  percStats <- percStats[order(percStats$split, percStats$bins_middle),]

  # Plot VPC -------------------------------------------------------------------
  p <- .plotVpc(percStats, obsDf, obsName, settings)
  return(invisible(p))
}

# Input Checks -----------------------------------------------------------------

.checkVPCInput = function(inputName, inputValue){
  isValid = TRUE
  inputName = inputName
  if (inputName == "settings") {
    if (is.list(inputValue) == FALSE) {
      message("ERROR: Unexpected type encountered. settings must be a list")
      isValid = FALSE
    } else {
      for (i in 1:length(inputValue)) {
        if (!.checkVPCSettings(settingName = names(inputValue)[i], settingValue = inputValue[[i]])) {
          isValid = FALSE
        }
      }
      if (!.checkVPCSettings(settingName = "settings", settingValue = inputValue)) {
        isValid = FALSE
      }
    }
  } else if (inputName == "stratifier") {
    if (is.list(inputValue) == FALSE) {
      message("ERROR: Unexpected type encountered. stratifier must be a list")
      isValid = FALSE
    } else {
      for (i in 1:length(inputValue)) {
        if (!.checkVPCStratifier(stratifierName = names(inputValue)[i], stratifierValue = inputValue[[i]])) {
          isValid = FALSE
        }
      }
      if (!.checkVPCStratifier(stratifierName = "stratifier", stratifier = inputValue)) {
        isValid = FALSE
      }
    }
  } else if (inputName == tolower("preferences")) {
    if (is.list(inputValue) == FALSE) {
      message("ERROR: Unexpected type encountered. preferences must be a list")
      isValid = FALSE
    } else {
      for (i in 1:length(inputValue)) {
        if (!.checkVPCPreferences(settingName = names(inputValue)[i], settingValue = inputValue[[i]])) {
          isValid = FALSE
        }
      }
    }
  }
  return(invisible(isValid))
}

.checkVPCSettings = function(settingName, settingValue){
  isValid = TRUE
  settingName = settingName
  if (settingName == "settings") {
    if (all(is.element(c("useCensoredData", "censoredData"), names(settingValue)))) {
      if (settingValue$useCensoredData == FALSE) {
        message(paste0(
          "ERROR: Unexpected setting encountered. when `useCensoredData` is set to FALSE,
          no `censoredData` method can be specified."
        ))
        isValid = FALSE
      } 
    } else if (is.element("usedFixedNbBins", names(settingValue))) {
      if (settingValue$usedFixedNbBins & !is.null(settingValue$estimatedNbBins)) {
        message(paste0(
          "ERROR: Unexpected setting encountered. when `usedFixedNbBins` is set to TRUE,
          no `estimatedNbBins` range can be specified, check `fixedNbBins` setting instead."
        ))
        isValid = FALSE
      } else if (!settingValue$usedFixedNbBins & !is.null(settingValue$fixedNbBins)) {
        message(paste0(
          "ERROR: Unexpected setting encountered. when usedFixedNbBins is set to TRUE,
          no `fixedNbBins` value can be specified, check `estimatedNbBins` setting instead."
        ))
        isValid = FALSE
      }
    } else if (is.element("binCriteria", names(settingValue))) {
      if (settingValue$binCriteria == "manual") {
        if (any(is.element(c("nbBinData", "useFixedNbBins", "fixedNbBins", "estimatedNbBins"),
                           names(settingValue)))) {
          message(paste0(
            "ERROR: Unexpected setting encountered. binCriteria set to manual,
            `nbBinData`, `useFixedNbBins`, `fixedNbBins`, estimatedNbBins` settings cannot be specified."
          ))
          isValid = FALSE
        } else if (!is.element("manualBins", names(settingValue))) {
          message(paste0(
            "ERROR: Unexpected setting encountered. binCriteria set to manual,
            `manualBins` must be specified."
          ))
          isValid = FALSE
        }
      } else {
        if (is.element("manualBins", names(settingValue))) {
          message(paste0(
            "ERROR: Unexpected setting encountered. binCriteria not manual,
            `manualBins` cannot be specified."
          ))
          isValid = FALSE
        }
      }
    }
  } else if (settingName %in% c(
    "legend", "grid", "displayObservedData", "displayCensoredData",
    "displayEmpiricalPercentiles", "displayPredictedPercentiles",
    "displayPredictionInterval", "outliersDots", "outlierAreas",
    "correctedPredictions", "linearInterpolation", "useCensoredData",
    "binLimits", "useFixedNbBins", "xlogScale", "ylogScale")) {
    if (!length(settingValue) == 1) {
      message(paste0(
        "ERROR: Unexpected type encountered. `", settingName, "` must be of size 1."
      ))
      isValid = FALSE
    } else if (!is.logical(settingValue)) {
      message(paste0(
        "ERROR: Unexpected type encountered. `", settingName, "` must be a boolean."
      ))
      isValid = FALSE
    }
  } else if (settingName %in% c("level", "higherPercentile")) {
    if (!length(settingValue) == 1) {
      message(paste0(
        "ERROR: Unexpected type encountered. `", settingName, "` must be of size 1."
      ))
      isValid = FALSE
    } else if (!is.double(settingValue) | settingValue < 0 | settingValue > 100)  {
      message(paste0(
        "ERROR: Unexpected type encountered. `", settingName, "` must be a double in [0, 100]."
      ))
      isValid = FALSE
    }
  } else if(settingName == "censoredData"){
    if (!length(settingValue) == 1) {
      message("ERROR: Unexpected type encountered. censoredData must be of size 1.")
      isValid = FALSE
    } else if (!is.element(tolower(settingValue), c("simulated", "loq"))) {
      message("ERROR: Unexpected type encountered. censoredData must be either `simulated` or `loq`.")
      isValid = FALSE
    }
  } else if(settingName == "binCriteria"){
    if (!length(settingValue) == 1) {
      message("ERROR: Unexpected type encountered. binCriteria must be of size 1.")
      isValid = FALSE
    } else if (!is.element(settingValue, c("manual", "equalwidth", "equalsize", "leastsquare"))) {
      message("ERROR: Unexpected type encountered. binCriteria must be `manual`, `equalwidth`, `equalsize` or `leastsquare`.")
      isValid = FALSE
    }
  } else if(settingName == "fixedNbBins"){
    if (length(settingValue)  != 1) {
      message(paste(
        "ERROR: Unexpected type encountered. fixedNbBins must be of size 1."
      ))
      isValid = FALSE
    } else if (! settingValue >= 0 | ! as.integer(settingValue) == settingValue) {
      message("ERROR: Unexpected type encountered. fixedNbBins must be positive integer.")
      isValid = FALSE
    }
  } else if(settingName %in% c("estimatedNbBins", "nbBinData")){
    if (length(settingValue)  != 2) {
      message(paste0(
        "ERROR: Unexpected type encountered. `", settingName, "` must be of size 2."
      ))
      isValid = FALSE
    } else if (!all(settingValue >= 0) | !all(as.integer(settingValue) == settingValue |
               settingValue[1] > settingValue[2])) {
      message(paste0(
        "ERROR: Unexpected type encountered. `", settingName, "` must be positive integers."))
      isValid = FALSE
    }
  } else if(settingName == "manualBins"){
    if (length(settingValue) < 1) {
      message("ERROR: Unexpected type encountered. manualBins must be a vector with at least 1 bin.")
      isValid = FALSE
    } else if (!all(settingValue >= 0) | !all(is.double(settingValue))) {
      message("ERROR: Unexpected type encountered. manualBins must be positive doubles")
      isValid = FALSE
    }
  } else if(settingName == "xlabel"){
    if (!length(settingValue) == 1) {
      message("ERROR: Unexpected type encountered. xlabel must be of size 1.")
      isValid = FALSE
    } else if (!is.element(settingValue, c("time", "timeSinceLastDose"))) {
      message("ERROR: Unexpected type encountered. xlabel must be either `time` or `timeSinceLastDose`")
      isValid = FALSE
    }
  } else {
    message(paste0("ERROR: Unknown setting ", settingName))
    isValid = FALSE
  }
  return(isValid)
}

.checkVPCStratifier = function(stratifierName, stratifierValue){
  isValid = TRUE
  allowedStratifiers <- c(
    mlx.getCovariateInformation()$name,
    mlx.getData()$header[mlx.getData()$headerTypes == "occ"]
  )
  if (stratifierName == "stratifier") {
    names(stratifierValue) <- names(stratifierValue)
    if (all(is.element(c("covariateSplit", "covariateFilter"), names(stratifierValue)))) {
      if (any(is.element(stratifierValue["covariateSplit"],
                         stratifierValue["covariateFilter"]))) {
        message("ERROR: Unexpected stratifier encountered. `covariateFilter` and `covariateFilter` list must be disjunct.")
        isValid = FALSE
      }
    } else if ("covariateFilter" %in% stratifierValue) {
      for (i in seq(1, length(stratifierValue))) {
        name <- names(stratifierValue)[i]
        value <- stratifierValue[i]
        if (
          (!is.null(covariateScaling[name]) & value >= length(covariateScaling[name])) |
          (is.null(covariateScaling[name]) & value > 2)
        ){
          nbGroup <- ifelse(is.null(covariateScaling[name]), 2, length(covariateScaling[name]))
          message(paste0(
            "ERROR: Unexpected stratifier encountered. Only ", nbGroup, "groups ",
            "are defined for `covariateFilter`", name, " You must choose among those groups."
          ))
        }
      }
      
    }
  } else if (stratifierName == "covariateSplit"){
    if (!all(is.element(stratifierValue, allowedStratifiers))) {
      message(paste0(
        "ERROR: Unknown `", stratifierName, "` names. You must choose in (",
        paste(allowedStratifiers, collapse = ", "), ")"
      ))
      isValid = FALSE
    }
  } else if (stratifierName == "covariateFilter") {
    if (!is.list(stratifierValue)) {
      message("ERROR: Unexpected type encountered. covariateFilter must be a list")
      isValid = FALSE
    } else {
      covList <- mlx.getCovariateInformation()
      covDf <- covList$covariate
      contCovariates <- covList$name[covList$type == "continuous"]
      catCovariates <- covList$name[covList$type == "categorical"]
      for (i in 1:length(stratifierValue)) {
        name <- names(stratifierValue)[i]
        value <- stratifierValue[[name]]
        if (!all(is.element(name, allowedStratifiers))) {
          message(paste0(
            "ERROR: Unknown `", stratifierName, "` names. You must choose in (",
            paste(allowedStratifiers, collapse = ", "), ")"
          ))
          isValid = FALSE
        } else if (name %in% catCovariates & !value %in% unique(covDf[[name]])) {
          message(paste0(
            "ERROR: Unexepect value for covariateFilter `", name, "` . You must ",
            "choose in (", paste(unique(covDf[name]), collapse = ", "), ")"
          ))
          isValid = FALSE
        } else if (name %in% contCovariates & !value > 0 & !as.integer(value) == value) {
          message(paste0(
            "ERROR: Unexpected type encountered in `covariateFilter`. Group id ",
            "associated to continous covariates `", name, "` must be a positive integer"
          ))
          isValid = FALSE
        }
      }
    }
  } else if (stratifierName == "covariateScaling") {
    if (is.list(stratifierValue) == FALSE) {
      message("ERROR: Unexpected type encountered. covariateScaling must be a list")
      isValid = FALSE
    } else {
      continousCovariates <- mlx.getCovariateInformation()$name[mlx.getCovariateInformation()$type == "continuous"]
      for (i in 1:length(stratifierValue)) {
        name <- names(stratifierValue)[i]
        value <- stratifierValue[[name]]
        covList <- mlx.getCovariateInformation()$covariate[name]
        if (!is.element(name, continousCovariates)) {
          message(paste0(
            "ERROR: Unknown `covariateScaling` name ", name,
            "You must choose continuous covariates."
          ))
          isValid = FALSE
        }
      }
    }
  } else {
    message(paste0("ERROR: Unknown stratifier ", stratifierName))
    isValid = FALSE
  }
  return(isValid)
}

# Stats on observations and simulations ----------------------------------------

.stats <- function(x, higherPerc, header = NULL) {
  percentiles = sort(c(1 - higherPerc /100, higherPerc /100))
  st <- c(
    median = median(x),
    lower = quantile(x, percentiles[1], type = 5)[[1]],
    upper = quantile(x, percentiles[2], type = 5)[[1]]
  )
  names(st) <- sapply(names(st), function(s) paste(c(header, s), collapse = "_"))
  return(st)
}

.interpolate <- function(data, xName, yNames) {
  x_range <- range(data[[xName]])
  x <- seq(x_range[1], x_range[2], length.out = 10000)
  dfOut <- NULL
  for (i in seq(1, length(yNames))) {
    df <- as.data.frame(approx(data[[xName]], data[[yNames[i]]], x))
    names(df) <- c(xName, yNames[i])
    if (is.null(dfOut)) dfOut <- df else dfOut <- merge(dfOut, df, by = c(xName))
  }
  return(dfOut)
}

# Plot -------------------------------------------------------------------------
.plotVpc <- function(percStats, obsDf, obsName, settings) {
  p <- ggplot(percStats) + xlab(settings$xlabel) + ylab(obsName) +
    theme(
      panel.background = element_rect(fill = "white", colour = "white",
                                      size = 2, linetype = "solid"),
      strip.background = element_rect(colour="white", fill="white", 
                                      size=1.5, linetype="solid")) +
    facet_wrap(.~split, ncol=3)

  colsName <- c("id", "label", "color", "linetype", "shape", "fill", "alpha")
  colsDf <- as.data.frame(matrix(ncol = 7, nrow = 0, dimnames = list(NULL, colsName)))

  if (settings$grid) {
    p <- p + theme(
      panel.grid.major = element_line(size = 0.1, linetype = 'solid',
                                      colour = "grey"),
      panel.grid.minor = element_line(size = 0.1, linetype = 'solid',
                                      colour = "grey")
    )
  }
  
  if (settings$displayObservedData) {
    color <- rgb(70, 130, 180, maxColorValue = 255)
    p <- p + geom_point(data = subset(obsDf, censored == 0),
                        aes(x = time, y = get(obsName), color = "c1", shape = "c1"))
    colsDf = rbind(colsDf, list(id = "c1", label = "Observed data", color = color,
                                linetype = NA, shape = 16, fill = NA, alpha = NA))
  }
  
  if (settings$displayCensoredData) {
    color <- rgb(212, 66, 66, maxColorValue = 255)
    p <- p + geom_point(data = subset(obsDf, censored != 0),
                        aes(x = time, y = get(obsName), color = "c2", shape = "c2"))

    colsDf = rbind(colsDf, list(id = "c2", label = "Censored data", color = color,
                                linetype = NA, shape = 16, fill = NA, alpha = NA))
  }
  
  if (settings$displayEmpiricalPercentiles) {
    color <- rgb(70, 130, 180, maxColorValue = 255)
    p <- p +
      geom_line(aes(x = bins_middle, y = empirical_lower, color = "c3")) +
      geom_line(aes(x = bins_middle, y = empirical_median, color = "c3")) +
      geom_line(aes(x = bins_middle, y = empirical_upper, color = "c3")) +
      geom_point(aes(x = bins_middle, y = empirical_lower), color = color) +
      geom_point(aes(x = bins_middle, y = empirical_median), color = color) +
      geom_point(aes(x = bins_middle, y = empirical_upper), color = color)

    colsDf = rbind(colsDf, list(id = "c3", label = "Empirical percentiles",
                                color = color, linetype = 1, shape = NA,
                                fill = NA, alpha = NA))
  }

  if (settings$displayPredictedPercentiles) {
    color <- rgb(0, 0, 0, maxColorValue = 255)
    p <- p +
      geom_line(aes(x = bins_middle, y = theoretical_lower_median, color = "c4"), linetype = "dashed") +
      geom_line(aes(x = bins_middle, y = theoretical_median_median, color = "c4"), linetype = "dashed") +
      geom_line(aes(x = bins_middle, y = theoretical_upper_median, color = "c4"), linetype = "dashed") +
      geom_point(aes(x = bins_middle, y = theoretical_lower_median), color = color) +
      geom_point(aes(x = bins_middle, y = theoretical_upper_median), color = color) +
      geom_point(aes(x = bins_middle, y = theoretical_upper_median), color = color)

    colsDf = rbind(colsDf, list(id = "c4", label = "Predicted median", color = color,
                                linetype = 2, shape = NA, fill = NA, alpha = NA))
  }

  if (settings$displayPredictionInterval) {
    color_extremes <- alpha(rgb(0, 128, 255, maxColorValue = 255), 0.3)
    color_median <- alpha(rgb(255, 0, 0, maxColorValue = 255), 0.2)
    colsDf = rbind(colsDf, list(id = "c5", label = "Prediction interval median",
                                color = color_median, linetype = NA, shape = NA,
                                fill = 1, alpha = 0.03))
    colsDf = rbind(colsDf, list(id = "c6", label = "Prediction interval percentiles",
                                color = color_extremes, linetype = NA, shape = NA,
                                fill = 1, alpha = 0.05))
    if (!settings$linearInterpolation) {
      p <- p +
        geom_rect(
          aes(xmin = bins_start, xmax = bins_stop,
              ymin = theoretical_lower_lower, ymax = theoretical_lower_upper,
              fill="c6")) +
        geom_rect(
          aes(xmin = bins_start, xmax = bins_stop,
              ymin = theoretical_upper_lower, ymax = theoretical_upper_upper,
              fill="c6")) +
        geom_rect(
          aes(xmin = bins_start, xmax = bins_stop,
              ymin = theoretical_median_lower, ymax = theoretical_median_upper,
              fill="c5"))
    } else {
      extremesDf <- do.call(rbind, by(
        percStats, percStats$split,
        function(df) {
          df <- df[order(df$bins_start),]
          rowf <- df[1,]
          rowl <- tail(df, n=1)
          rowf$bins_middle <- rowf$bins_start
          rowl$bins_middle <- rowl$bins_stop
          df <- rbind(rowf, df, rowl)
          return(df)
        }
      ))
      p <- p +
        geom_ribbon(
          data = extremesDf,
          aes(x = bins_middle, ymin = theoretical_lower_lower, ymax = theoretical_lower_upper, fill="c6")) +
        geom_ribbon(
          data = extremesDf,
          aes(x = bins_middle, ymin = theoretical_upper_lower, ymax = theoretical_upper_upper, fill="c6")) +
        geom_ribbon(
          data = extremesDf,
          aes(x = bins_middle, ymin = theoretical_median_lower, ymax = theoretical_median_upper, fill="c5"))

      if (settings$outlierAreas) {
        color <- rgb(255, 0, 0, maxColorValue = 255)
        percStatsInterp <- do.call(rbind, by(
          percStats,
          percStats$split,
          function(df) {
            yNames <- names(df)[!(names(df) == "split" | grepl("bins", names(df)))]
            r <- subset(df, select = c("bins_middle", "split", yNames))
            if (nrow(df) > 1)
              r <- .interpolate(df, "bins_middle", yNames)
              r["split"] <- unique(df$split)
            return(r)
          }
        ))
        p <- p +
          geom_ribbon(
            data = percStatsInterp,
            aes(x = bins_middle,
                ymin = pmin(theoretical_lower_upper, empirical_lower),
                ymax = empirical_lower, fill = "c7")
          ) +
          geom_ribbon(
            data = percStatsInterp,
            aes(x = bins_middle,
                ymin = pmin(empirical_lower, theoretical_lower_lower),
                ymax = theoretical_lower_lower, fill = "c7")) +
          geom_ribbon(
            data = percStatsInterp,
            aes(x = bins_middle,
                ymin = pmin(theoretical_median_upper, empirical_median),
                ymax = empirical_median, fill = "c7")) +
          geom_ribbon(
            data = percStatsInterp,
            aes(x = bins_middle,
                ymin = pmin(empirical_median, theoretical_median_lower),
                ymax = theoretical_median_lower, fill = "c7")) +
          geom_ribbon(
            data = percStatsInterp,
            aes(x = bins_middle,
                ymin = pmin(theoretical_upper_upper, empirical_upper),
                ymax = empirical_upper, fill = "c7")) +
          geom_ribbon(
            data = percStatsInterp,
            aes(x = bins_middle,
                ymin = pmin(empirical_upper, theoretical_upper_lower),
                ymax = theoretical_upper_lower, fill = "c6"))

        colsDf = rbind(colsDf, list(id = "c7", label = "Areas", color = color,
                                    linetype = NA, shape = NA, fill = 1, alpha = 1))
      }
    }
  }
  
  if (settings$outliersDots) {
    color <- rgb(255, 0, 0, maxColorValue = 255)
    p <- p +
      geom_point(
        data = subset(
          percStats,
          theoretical_lower_upper < empirical_lower | theoretical_lower_lower > empirical_lower
        ),
        aes(x = bins_middle, y = empirical_lower, color = "c8", shape = "c8"), size = 3) +
      geom_point(
        data = subset(
          percStats,
          theoretical_median_upper < empirical_median | theoretical_median_lower > empirical_median
        ),
        aes(x = bins_middle, y = empirical_median, color = "c8", shape = "c8"), size = 3) +
      geom_point(
        data = subset(
          percStats,
          theoretical_upper_upper < empirical_upper | theoretical_upper_lower > empirical_upper
        ),
        aes(x = bins_middle, y = empirical_upper, color = "c8", shape = "c8"), size = 3)

      colsDf = rbind(colsDf, list(id = "c8", label = "Dots", color = color,
                                  linetype = NA, shape = 1, fill = NA, alpha = NA))
  }
  
  if (settings$binLimits) {
    color <- rgb(255, 0, 0, maxColorValue = 255)
    datavlines = do.call(rbind, by(
      subset(percStats, select = c("bins_middle", "split")),
      percStats$split,
      function(df) data.frame(bins_middle = unique(df$bins_middle), split = unique(df$split))
    ))
    p <- p + geom_vline(data = datavlines, aes(xintercept = bins_middle), color = color, size = 0.2)
  }
  
  if (settings$ylogScale) p <- p + scale_y_log10()
  if (settings$xlogScale) p <- p + scale_x_log10()

  # legend
  cols <- setNames(colsDf[is.na(colsDf$fill),]$color, colsDf[is.na(colsDf$fill),]$id)
  shapes <- setNames(as.double(colsDf[is.na(colsDf$fill),]$shape), colsDf[is.na(colsDf$fill),]$id)
  linetypes <- as.double(colsDf[is.na(colsDf$fill),]$linetype)
  alphas <- as.double(colsDf[!is.na(colsDf$fill),]$alpha)
  p <- p + scale_color_manual(
    name = "legend", values = cols, labels = colsDf[is.na(colsDf$fill),]$label,
    guide=guide_legend(override.aes = list(shape = shapes, linetype = linetypes, fill = NA)))
  p <- p + scale_shape_manual(values = shapes, guide=F)
  fills <- setNames(colsDf[!is.na(colsDf$fill),]$color, colsDf[!is.na(colsDf$fill),]$id)
  p <- p + scale_fill_manual(
    name = "legend", values = fills, labels = colsDf[!is.na(colsDf$fill),]$label,
    guide = guide_legend(title = "legend", override.aes = list(alpha = alphas))
  )
  if (!settings$legend) {
    p <- p + theme(legend.position = "none")
  } else {
    p <- p + theme(
      legend.key = element_rect(fill = "white",colour = "white"),
      legend.title = element_blank(), legend.key.size = unit(0.5, "cm"),
      legend.position=c(.8,.75)
    )
  }
  return(p)
}
