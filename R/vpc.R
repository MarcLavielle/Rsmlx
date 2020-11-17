#' Visual Predictive Check
#'
#' Creates a VPC plot from monolix project
#' @param project Monolix project
#' @param time [optional] in {`time`, `timeSinceLastDose`} (default `time`).
#' `timeSinceLastDose` only possible when there is an `amount` column in the dataset
#' @param obsName [optional] observation name. By default the first observation is considered.
#' @param stratifier [optional] a list of settings to stratify vpc plots:
#' \itemize{
#' \item \code{covariateSplit} [optional] vector of covariate names used to split vpc plot 
#' (no split by default).
#' \item \code{covariateFilter} [optional] list of covariate names used to filter vpc plot
#' names are covariate names and values are category name in case of categorical covariate
#' and group id in case of continuous covariate.
#' (no filter by default).
#' \item \code{covariateScaling} [optional] list of continuous covariate scaling
#' Names are continuous covariate names and values are vector of break doubles.
#' If not scaling defined, by default break is defined at the median.
#' }
#' @param vpcSettings [optional] a list of settings for vpc computation:
#' \itemize{
#' \item \code{level} [optional] (integer) level for prediction intervals computation (default 90).
#' \item \code{higherPercentile} [optional][continuous data] (integer) Higher percentile for empirical and predicted percentiles computation (default 90). For continuous data only.
#' \item \code{correctedPredictions} [optional][continuous data] (boolean) if TRUE, pcVPC are computed using Uppsala prediction correction (default FALSE). For continuous data only.
#' \item \code{useCensoredData} [optional][continuous data] (boolean) Choose to use BLQ data or to ignore it to compute the VPC (default TRUE). For continuous data only.
#' \item \code{censoredData} [optional][continuous data] in {`simulated`, `loq`}, BLQ data can be simulated, or can be equal to the limit of quantification (LOQ) (default `simulated`). For continuous data only.
#' }
#' @param eventDataDisplay [optional][event data] a list of display settings for event data plots only.
#' \itemize{
#' \item \code{survivalCurve} [optional] (boolean) Add/remove plot for survival function (Kapan-Meier plot) (default TRUE).
#' \item \code{meanNumberEventsCurve} [optional] (boolean) Add/remove plot for mean number of events per individual (default FALSE).
#' \item \code{empiricalCurve} [optional] (boolean) Add/remove empirical curve (default TRUE).
#' \item \code{predictedMedian} [optional] (boolean) Add/remove predicted median (default FALSE).
#' \item \code{predictionInterval} [optional] (boolean) Add/remove prediction intervals given by the model (default TRUE).
#' }
#' @param discreteDataDisplay [optional][discrete data] a list of display settings for discrete data plots only.
#' \itemize{
#' \item \code{empiricalProbability} [optional] (boolean) Add/remove empirical probability (default TRUE).
#' \item \code{predictedMedian} [optional] (boolean) Add/remove theorical probability (default FALSE).
#' \item \code{predictionInterval} [optional] (boolean) Add/remove prediction intervals given by the model (default TRUE).
#' \item \code{outlierAreas} [optional] (boolean) Add/remove red areas indicating empirical percentiles that are outside prediction intervals (default FALSE).
#' }
#' @param continuousDataDisplay [optional][continuous data] a list of display settings for continuous data plots only.
#' \itemize{
#' \item \code{observedData} [optional] (boolean) Add/remove observed data (default FALSE).
#' \item \code{censoredData} [optional] (boolean) Add/remove BLQ data if present (default FALSE).
#' \item \code{empiricalPercentiles} [optional] (boolean) Add/remove empirical percentiles for the 10%, 50% and 90% quantiles (default TRUE).
#' \item \code{predictedPercentiles} [optional] (boolean) Add/remove theoretical percentiles for the 10% and 90% quantiles (default FALSE).
#' \item \code{predictionInterval} [optional] (boolean) Add/remove prediction intervals given by the model for the 10% and 90% quantiles (in blue) and the 50% quantile (in pink) (default TRUE).
#' \item \code{outliersDots} [optional] (boolean) Add/remove red dots indicating empirical percentiles that are outside prediction intervals (default TRUE).
#' \item \code{outlierAreas} [optional] (boolean) Add/remove red areas indicating empirical percentiles that are outside prediction intervals (default TRUE).
#' }
#' @param xBinsSettings [optional][continuous and discrete data] a list of settings for time axis binning.
#' \itemize{
#' \item \code{limits} [optional] (boolean) Add/remove vertical lines on the scatter plots to indicate the bins (default FALSE).
#' \item \code{useFixedBins} [optional] (boolean) If TRUE, specify manually bin list (default FALSE).
#' \item \code{fixedBins} [optional] (list[integer]) Define mannually a list of bins. To use when `useFixedBins` is set to TRUE.
#' \item \code{criteria} [optional] Choose the bining criteria among `equalwidth`, `equalsize` or `leastsquare` (default leastsquare).
#' \item \code{useFixedNbBins} [optional] (boolean) If TRUE define a fixed number of bins, else define a range for automatic selection (default FALSE).
#' \item \code{nbBins} [optional] (integer) Define a fixed number of bins (default 10).
#' \item \code{binRange} [optional] (list[integer]) Define a range for the number of bins (default c(5, 20)).
#' \item \code{nbBinData} [optional] (list[integer]) Define a range for the number of data points per bin (default c(10, 200)).
#' }
#' @param yBinsSettings [optional][countable discrete data] a list of settings for y axis binning.
#' \itemize{
#' \item \code{useFixedBins} [optional] (boolean) If TRUE, specify manually bin list (default FALSE).
#' \item \code{fixedBins} [optional] (list[integer]) Define mannually a list of bins. To use when `useFixedBins` is set to TRUE.
#' \item \code{criteria} [optional]  Choose the bining criteria among `equalwidth`, `equalsize` or `leastsquare` (default leastsquare).
#' \item \code{useFixedNbBins} [optional] (boolean) If TRUE define a fixed number of bins, else define a range for automatic selection (default FALSE).
#' \item \code{nbBins} [optional] (integer) Define a fixed number of bins (default 5).
#' \item \code{binRange} [optional] (list[integer]) Define a range for the number of bins (default c(3, 7)).
#' \item \code{nbBinData} [optional] (list[integer]) Define a range for the number of data points per bin (default c(10, 200)).
#' }
#' @param plotSettings [optional] a list of settings for plot
#' \itemize{
#' \item \code{legend} [optional] (boolean) Add/remove legend (default FALSE).
#' \item \code{grid} [optional] (boolean) Add/remove grid (default FALSE).
#' \item \code{xlogScale} [optional] (boolean) Add/remove log scale for x axis (default FALSE).
#' \item \code{ylogScale} [optional] (boolean) Add/remove log scale for x axis (default FALSE).
#' \item \code{linearInterpolation} [optional] (boolean) If TRUE set piece wise display for prediction intervals, else show bins as rectangular (default TRUE).
#' \item \code{xlab} [optional] (str) Time label (default "time" if time = "time", "time since last dose" if time = "timeSinceLastDose").
#' \item \code{ylab} [optional] (str) y label (default observation name).
#' }
#' @param vpcTheme [optional] theme to be used in VPC. Expects list of class vpc_theme created with function createVpcTheme()
#' @return a ggplot2 object
#' @importFrom ggplot2 ggplot element_rect element_line geom_ribbon geom_point
#'             geom_rect geom_line theme aes xlab ylab facet_wrap facet_grid
#'             alpha scale_color_manual scale_shape_manual scale_fill_manual
#'             guide_legend grid.arrange
#' @export
#' @examples
#' \dontrun{
#'   stratifier <- list(covariateSplit = c("sex", "wt"),
#'                      covariateFilter = list("age" = 1),
#'                      covariateScaling = list("age" = c(65), "wt" = c(70)))
#'   settings <- list(xBinCriteria = "manual", xManualBins = c(0, 24, 48, 72, 96, 120))
#'   plotVPC(project="RsmlxDemo1.mlxtran", stratifier=stratifier)
#'   plotVPC(project="RsmlxDemo1.mlxtran", settings=settings)
#'   plotVPC(project="RsmlxDemo1.mlxtran", settings=list(xlabel = "timeSinceLastDose"))
#' }
#' @seealso \link{createVpcTheme}

vpc <- function(project, time = NULL, obsName = NULL, 
                stratifier = NULL, vpcSettings = NULL,
                eventDataDisplay = NULL, discreteDataDisplay = NULL,
                continuousDataDisplay = NULL, xBinsSettings = NULL, yBinsSettings = NULL,
                plotSettings = NULL, vpcTheme = NULL) {

  # Check Arguments ------------------------------------------------------------
  # load project
  r <- prcheck(project)
  project <- r$project

  # Check and initialize time
  if(!is.null(time)){
    if(!.checkVPCInput(inputName = "time", inputValue = time))
      return(invisible(FALSE))
  } else {
    time <- "time"
  }
  if(!is.null(obsName)){
    if(!.checkVPCInput(inputName = "observationName", inputValue = obsName))
      return(invisible(FALSE))
  } else {
    obsName <- mlx.getObservationInformation()$name[1]
  }
  
  obsType <- .getObsType(obsName)
  dataType <- obsType$type
  dataSubType <- obsType$subType

  # Check and initialize vpc settings
  if(!is.null(vpcSettings)) {
    if(!.checkVPCInput(inputName = "vpcSettings", inputValue = vpcSettings))
      return(invisible(FALSE))
  } else {
    vpcSettings <- list()
  }
  defaultVpcSettings <- list(
    level = 90, higherPercentile = 90,
    correctedPredictions = FALSE,
    useCensoredData = ifelse(dataType == "continuous", TRUE, FALSE),
    censoredData = "simulated"
  )
  for (s in names(defaultVpcSettings)) {
    if (is.null(vpcSettings[[s]])) {
      vpcSettings[s] <- defaultVpcSettings[[s]]
    }
  }
  
  # Check and initialize curves display
  if(!is.null(eventDataDisplay)) {
    if(!.checkVPCInput(inputName = "eventDataDisplay", inputValue = eventDataDisplay, dataType = dataType)){return(invisible(FALSE))}
  } else {
    eventDataDisplay <- list()
  }
  defaultEventDataDisplay <- list(
    empiricalCurve = TRUE, predictedMedian = FALSE, predictionInterval = TRUE,
    survivalCurve = TRUE, meanNumberEventsCurve = FALSE
  )
  for (display in names(defaultEventDataDisplay)) {
    if (dataType != "event") {
      eventDataDisplay[display] <- FALSE
    } else if (is.null(eventDataDisplay[[display]])) {
      eventDataDisplay[display] <- defaultEventDataDisplay[[display]]
    }
  }
  
  if(!is.null(discreteDataDisplay)) {
    if(!.checkVPCInput(inputName = "discreteDataDisplay", inputValue = discreteDataDisplay, dataType = dataType)){return(invisible(FALSE))}
  } else {
    discreteDataDisplay <- list()
  }
  defaultDiscreteDataDisplay <- list(
    empiricalProbability = TRUE, theoricalProbability = FALSE,
    predictionInterval = TRUE, outlierAreas = FALSE
  )
  for (display in names(defaultDiscreteDataDisplay)) {
    if (dataType != "discrete") {
      discreteDataDisplay[display] <- FALSE
    } else if (is.null(discreteDataDisplay[[display]])) {
      discreteDataDisplay[display] <- defaultDiscreteDataDisplay[[display]]
    }
  }
  
  if(!is.null(continuousDataDisplay)) {
    if(!.checkVPCInput(inputName = "continuousDataDisplay", inputValue = continuousDataDisplay, dataType = dataType)){return(invisible(FALSE))}
  } else {
    continuousDataDisplay <- list()
  }
  defaultContinuousDataDisplay <- list(
    observedData = FALSE, censoredData = FALSE,
    empiricalPercentiles = TRUE, predictedPercentiles = FALSE,
    predictionInterval = TRUE, outliersDots = TRUE, outlierAreas = TRUE
  )
  for (display in names(defaultContinuousDataDisplay)) {
    if (dataType != "continuous") {
      continuousDataDisplay[display] <- FALSE
    } else if (is.null(continuousDataDisplay[[display]])) {
      continuousDataDisplay[display] <- defaultContinuousDataDisplay[[display]]
    }
  }

  # bins settings
  if(!is.null(xBinsSettings)) {
    if(!.checkVPCInput(inputName = "xBinsSettings", inputValue = xBinsSettings))
      return(invisible(FALSE))
  }
  xBinsSettings$bins <- ifelse(dataType != "event", TRUE, FALSE)
  if (is.null(xBinsSettings$nbDataPoints)) xBinsSettings$nbDataPoints <- 100
  if (is.null(xBinsSettings$limits)) xBinsSettings$limits <- FALSE
  if (is.null(xBinsSettings$useFixedBins))
    xBinsSettings$useFixedBins <- ifelse(!is.null(xBinsSettings$fixedBins), TRUE, FALSE)
  if (!xBinsSettings$useFixedBins) {
    if (is.null(xBinsSettings$criteria)) xBinsSettings$criteria <- "leastsquare"
    if (is.null(xBinsSettings$useFixedNbBins)) xBinsSettings$useFixedNbBins <- FALSE
    if (xBinsSettings$useFixedNbBins & is.null(xBinsSettings$nbBins)) xBinsSettings$nbBins <- 10
    if (!xBinsSettings$useFixedNbBins & is.null(xBinsSettings$binRange)) {
      if (dataType == "continuous") {
        xBinsSettings$binRange <- c(5, 20)
      } else {
        xBinsSettings$binRange <- c(5, 30)
      }
    }
    if (is.null(xBinsSettings$nbBinData)) xBinsSettings$nbBinData <- c(10, 200)
  }
  
  if(!is.null(yBinsSettings)) {
    if(!.checkVPCInput(inputName = "yBinsSettings", inputValue = yBinsSettings, dataSubType))
      return(invisible(FALSE))
  }
  yBinsSettings$bins = ifelse(dataType == "discrete" & dataSubType == "count", TRUE, FALSE)
  if (is.null(yBinsSettings$useFixedBins))
    yBinsSettings$useFixedBins <- ifelse(!is.null(yBinsSettings$fixedBins), TRUE, FALSE)
  if (!yBinsSettings$useFixedBins) {
    if (is.null(yBinsSettings$criteria)) yBinsSettings$criteria <- "leastsquare"
    if (is.null(yBinsSettings$useFixedNbBins)) yBinsSettings$useFixedNbBins <- FALSE
    if (yBinsSettings$useFixedNbBins & is.null(yBinsSettings$nbBins)) yBinsSettings$nbBins <- 5
    if (!yBinsSettings$useFixedNbBins & is.null(yBinsSettings$binRange)) yBinsSettings$binRange <- c(3, 7)
    if (is.null(yBinsSettings$nbBinData)) yBinsSettings$nbBinData <- c(10, 200)
  }
  
  # Check and initialize plot settings
  if(!is.null(plotSettings)){
    if(!.checkVPCInput(inputName = "plotSettings", inputValue = plotSettings))
      return(invisible(FALSE))
  }
  plotSettings <- list(
    legend = ifelse(is.null(plotSettings$legend), FALSE, plotSettings$legend),
    grid = ifelse(is.null(plotSettings$grid), TRUE, plotSettings$grid),
    xlogScale = ifelse(is.null(plotSettings$xlogScale), FALSE, plotSettings$xlogScale),
    ylogScale = ifelse(is.null(plotSettings$ylogScale), FALSE, plotSettings$ylogScale),
    linearInterpolation = ifelse(is.null(plotSettings$linearInterpolation), TRUE, plotSettings$linearInterpolation),
    xlab = ifelse(
      !is.null(plotSettings$xlab),
      plotSettings$xlab,
      "time"
    ),
    ylab = ifelse(
      !is.null(plotSettings$ylab),
      plotSettings$ylab,
      ifelse(!vpcSettings$correctedPredictions, obsName, paste0("Prediction corrected", obsName))
    )
  )

  # Check and initialize vpc theme
  if(!is.null(vpcTheme)){
    if(!.checkVPCInput(inputName = "vpcTheme", inputValue = vpcTheme))
      return(invisible(FALSE))
  }
  if (is.null(vpcTheme)) {
    if (dataType == "continuous") {
      vpcTheme <- continuous_theme()
    } else if (dataType == "discrete") {
      vpcTheme <- discrete_theme()
    } else {
      vpcTheme <- event_theme()
    }
  }
  
  displaySettings <- list(
    binLimits = xBinsSettings$limits,
    grid = plotSettings$grid,
    legend = plotSettings$legend,
    xlogScale = plotSettings$xlogScale,
    ylogScale = plotSettings$ylogScale,
    linearInterpolation = plotSettings$linearInterpolation,
    xlab = plotSettings$xlab,
    ylab = plotSettings$ylab,
    observedData = continuousDataDisplay$observedData,
    censoredData = continuousDataDisplay$censoredData,
    empiricalStats = any(
      continuousDataDisplay$empiricalPercentiles,
      discreteDataDisplay$empiricalProbability,
      eventDataDisplay$empiricalCurve
    ),
    theoricalMedian = any(
      continuousDataDisplay$predictedPercentiles,
      discreteDataDisplay$theoricalProbability,
      eventDataDisplay$predictedMedian
    ),
    predictionInterval = any(
      continuousDataDisplay$predictionInterval,
      discreteDataDisplay$predictionInterval,
      eventDataDisplay$predictionInterval
    ),
    outliersDots = continuousDataDisplay$outliersDots,
    outlierAreas = any(
      continuousDataDisplay$outlierAreas, discreteDataDisplay$outlierAreas
    ),
    survivalCurve = eventDataDisplay$survivalCurve,
    meanNumberEventsCurve = eventDataDisplay$meanNumberEventsCurve
  )
  
  # Check and initialize stratifier
  if(!is.null(stratifier)){
    if(!.checkVPCInput(inputName = "stratifier", inputValue = stratifier)) return(invisible(FALSE))
  }
  if (is.null(stratifier$covariateSplit)) stratifier$covariateSplit <- c() 
  if (is.null(stratifier$covariateFilter)) stratifier$covariateFilter <- c()
  if (is.null(stratifier$covariateScaling)) stratifier$covariateScaling <- c()
  
  # Read Data ------------------------------------------------------------------
  simName <- paste0("sim_", obsName)
  subjocc = .getSubjocc()
  
  # check Population parameters
  if (is.null(mlx.getEstimatedPopulationParameters())) {
    # Run population parameter estimation
    message(paste0("Run population parameter estimation."))
    mlx.runPopulationParameterEstimation()
  }
  
  # get obs data
  obsData <- .getObservationData(obsName, dataType)
  simData <- .getSimulationData(obsName)
  obsData <- obsData[, !names(obsData) %in% c("split", "filter", "color")]
  simData <- simData[, !names(simData) %in% c("split", "filter", "color")]
  
  if (dataType == "continuous") {
    if (vpcSettings$useCensoredData & !is.null(vpcSettings$censoredData)){
      if (vpcSettings$censoredData == "simulated" & length(unique(obsData$censored)) > 1)
        obsName <- paste0(obsName, "_simBlq")
    }
  }
  
  if (dataType == "event") obsData <- .addTTECensoredColumn(obsData)

  # compute time since last dose if needed
  if (time == "timeSinceLastDose") {
    obsData <- .addDoseRelTime(obsData)
    simData <- .addDoseRelTime(simData)
    if (is.null(displaySettings$xlab)) displaySettings$xlab <- "time since last dose"
    time <- "time_rel"
  } else {
    if (is.null(displaySettings$xlab)) displaySettings$xlab <- "time"
    time <- "time"
  }
  
  # Stratify: split and filter -------------------------------------------------
  if (!is.null(mlx.getCovariateInformation())) {
    covariates <- mlx.getCovariateInformation()
    covData <- covariates$covariate
    covTypes <- covariates$type
    if (length(.matchToData("occ")) > 0) {
      covTypes[.matchToData("occ")] <- "categorical"
    }
    # split column
    covData <- .addStratificationColumn(
      covData, covtypes,
      stratifiers = stratifier$covariateSplit,
      breaks = stratifier$covariateScaling,
      columnName = "split"
    )
    # filter column
    covData <- .addStratificationColumn(
      covData, covtypes,
      stratifiers = names(stratifier$covariateFilter),
      breaks = stratifier$covariateScaling,
      filtering = TRUE,
      filterGroups = stratifier$covariateFilter,
      columnName = "filter"
    )
    covData <- .addStratificationColumn(
      covData, covtypes,
      stratifiers = names(stratifier$covariateFilter),
      breaks = stratifier$covariateScaling,
      columnName = "filterNames"
    )
    filterName <- unique(covData$filterNames[covData$filter == 0])
    
    # add filtering name to split column
    covData$split <- sapply(covData$split, function(s) paste0(s, filterName))
    covData$split <- sapply(covData$split, function(s) ifelse(s == "", "All", s))
    obsData <- merge(obsData, subset(covData, select = c(subjocc, "split", "filter")), by = subjocc, sort = F)
    simData <- merge(simData, subset(covData, select = c(subjocc, "split", "filter")), by = subjocc, sort = F)
  } else {
    obsData$split <- rep("All", nrow(obsData))
    obsData$filter <- rep(FALSE, nrow(obsData))
    simData$split <- rep("All", nrow(simData))
    simData$filter <- rep(FALSE, nrow(simData))
  }
  obsData <- subset(obsData, filter == FALSE)
  simData <- subset(simData, filter == FALSE)
  
  # Bins -----------------------------------------------------------------------
  if (dataType == "continuous") {
    if (is.element("censored", names(obsData)) & !vpcSettings$useCensoredData) {
      obsData <- subset(obsData, censored == 0)
      simData <- subset(simData, censored == 0)
    }
  }
  xBinsValue <- NULL
  yBinsValue <- NULL
  for (s in unique(obsData$split)) {
    # x bins
    if (xBinsSettings$bins) {
      binsSplit <- .getBins(
        dataframe = subset(obsData, split == s), columnName = timeName,
        settings = xBinsSettings
      )
      xBinsValue <- rbind(xBinsValue, cbind(binsSplit, list(split = s)))
    }
    # y bins
    if (dataType == "discrete") {
      if (dataSubType == "categorical") {
        binsSplit <- .getBins(
          dataframe = subset(obsData, split == s), columnName = obsName,
          bycat = TRUE, categories = unique(obsData[[obsName]])
        )
        yBinsValue <- rbind(yBinsValue, cbind(binsSplit, list(split = s)))
      } else if (dataSubType == "count") {
        binsSplit <- .getBins(
          dataframe = subset(obsData, split == s), columnName = obsName,
          settings = yBinsSettings
        )
        yBinsValue <- rbind(yBinsValue, cbind(binsSplit, list(split = s)))
      }
      yBinsValue$binsName <- yBinsValue$bins
      yBinsValue[yBinsValue$bins_start == yBinsValue$bins_middle,]$binsName <- sapply(
        yBinsValue[yBinsValue$bins_start == yBinsValue$bins_middle,]$bins_start,
        function(x) paste0("P(", obsName, "=", x, ")")
      )
      yBinsValue[yBinsValue$bins_start < yBinsValue$bins_middle,]$binsName <- mapply(
        function(x, y) paste0("P(", ceiling(x), "<=", obsName, "<=", trunc(y), ")"),
        yBinsValue[yBinsValue$bins_start < yBinsValue$bins_middle,]$bins_start,
        yBinsValue[yBinsValue$bins_start < yBinsValue$bins_middle,]$bins_stop
      ) 
    }
  }
  obsData <- .addBinsIndex(obsData, xBinsValue, timeName, "binIndex")
  simData <- .addBinsIndex(simData, xBinsValue, timeName, "binIndex")
  obsData <- .addBinsIndex(obsData, yBinsValue, obsName, "category")
  simData <- .addBinsIndex(simData, yBinsValue, simName, "category")
  
  # Corrected Prediction -------------------------------------------------------
  if (vpcSettings$correctedPredictions) {
    res <- .applyCorrectedPrediction(obsData, simData)
    obsData <- res$obs
    simData <- res$sim
    obsName <- paste0(obsName, "_pc")
    simName <- paste0(simName, "_pc")
  }
  
  # Statistics -----------------------------------------------------------------
  # Empirical Data statistics per bin
  if (dataType == "continuous") {
    vpcData <- do.call(rbind, by(
      simData,
      simData$split,
      function(simSplit) {
        s <- unique(simSplit$split)
        obsSplit <- subset(obsData, split == s, select = c(obsName, timeName, "binIndex"))
        simSplit <- subset(simSplit, select = c("rep", simName, timeName, "binIndex"))
        res <- .computeContinuousVPC(
          obsSplit, obsName, simSplit, simName,
          vpcSettings$higherPercentile, vpcSettings$level
        )
        res <- cbind(list(split = s), res)
      }
    ))
    vpcData <- merge(vpcData, xBinsValue, by=c("bins", "split"))
    vpcData <- vpcData[order(vpcData$split, vpcData$bins_middle),]

  } else if (dataType == "discrete") {
    categories <- unique(obsData$category)
    vpcData <- do.call(rbind, by(
      subset(simData, select = c("rep", simName, "split", timeName, "binIndex", "category")),
      simData$split,
      function(simSplit) {
        s <- unique(simSplit$split)
        obsSplit <- subset(obsData, split == s, select = c(obsName, "split", timeName, "binIndex", "category"))
        res <- .computeDiscreteVPC(
          obsData, obsName, simData, simName,
          vpcSettings$higherPercentile, vpcSettings$level, categories
        )
        res <- cbind(list(split = s), res)
      }
    ))
    vpcData <- merge(vpcData, xBinsValue, by=c("bins", "split"))
    vpcData <- vpcData[order(vpcData$split, vpcData$category, vpcData$bins_middle),]
    vpcData$category <- sapply(vpcData$category, function(cat) yBinsValue$binsName[yBinsValue$bins == cat])

  } else if (dataType == "event") {
    vpcData <- do.call(rbind, by(
      simData,
      simData$split,
      function(simSplit) {
        s <- unique(simSplit$split)
        obsSplit <- subset(obsData, split == s)
        res <- .computeEventVPC(
          obsData, obsName, simData, simName,
          xBinsSettings$nbDataPoints, vpcSettings$level
        )
        res <- cbind(list(split = s), res)
      }
    ))
  }

  # Plot VPC -------------------------------------------------------------------
  if (dataType %in% c("discrete", "continuous")) {
    p <- plotVpc(vpcData, obsData, obsName, time, displaySettings, vpcTheme)
  } else {
    if (displaySettings$survivalCurve) {
      displaySettings$ylab <- "Survival Function"
      data <- subset(vpcData, select = c("split", "time", names(vpcData)[grepl("survivalFunction", names(vpcData))]))
      p1 <- plotVpc(data, obsData, obsName, time, displaySettings, vpcTheme)
      if (displaySettings$meanNumberEventsCurve) {
        displaySettings$ylab <- "Mean number of events per subject"
        data <- subset(vpcData, select = c("split", "time", names(vpcData)[grepl("averageEventNumber", names(vpcData))]))
        p2 <- plotVpc(data, obsData, obsName, time, displaySettings, vpcTheme)
        p <- grid.arrange(p1, p2, vp, sc, ncol=2)
      } else {
        p <- p1
      }
    } else if (displaySettings$meanNumberEventsCurve) {
      displaySettings$ylab <- "Mean number of events per subject"
      data <- subset(vpcData, select = c("split", "time", names(vpcData)[grepl("averageEventNumber", names(vpcData))]))
      p <- plotVpc(data, obsData, "meanNumberEventsCurve", time, displaySettings, vpcTheme)
    }
  }
  return(invisible(p))
}

# Input Checks -----------------------------------------------------------------
.checkVPCInput = function(inputName, inputValue, dataType = NULL){
  isValid = TRUE
  inputName = inputName
  
  if (inputName == "obsName") {
    obsnames <- mlx.getObservationInformation()$name
    if (!length(inputValue) == 1) {
      message("ERROR: Unexpected type encountered. obsName must be of size 1.")
      isValid = FALSE
    } else if (!is.element(inputValue, obsnames)) {
      message(paste0(
        "ERROR: Unexpected type encountered. obsName must be either in ",
        paste(obsName, collapse = ", ")
      ))
      isValid = FALSE
    }
  } else if (inputName == "time") {
    if (!length(inputValue) == 1) {
      message("ERROR: Unexpected type encountered. time must be of size 1.")
      isValid = FALSE
    } else if (!is.element(inputValue, c("time", "timeSinceLastDose"))) {
      message("ERROR: Unexpected type encountered. time must be either `time` or `timeSinceLastDose`")
      isValid = FALSE
    } else if (inputValue == "timeSinceLastDose") {
      if (! is.element("amount", mlx.getData()$headerTypes)) {
        message(paste0(
          "ERROR: Unexpected value encountered. `timeSinceLastDose` not available ",
          "for `time` setting when no amount column in the dataset"
        ))
        isValid = FALSE
      }
    }
  } else if (inputName == "vpcTheme") {
    if (class(vpc_theme) != "vpc_theme") {
      message("ERROR: Unexpected type encountered. vpcTheme must be a `vpc_theme` object")
    }
  } else if (inputName == "plotSettings") {
    if (is.list(inputValue) == FALSE) {
      message("ERROR: Unexpected type encountered. plotSettings must be a list")
      isValid = FALSE
    } else {
      for (i in 1:length(inputValue)) {
        if (!.checkVpcPlotSettings(settingName = names(inputValue)[i], settingValue = inputValue[[i]])) {
          isValid = FALSE
        }
      }
    }
  } else if (inputName == "vpcSettings") {
    if (is.list(inputValue) == FALSE) {
      message("ERROR: Unexpected type encountered. vpcSettings must be a list")
      isValid = FALSE
    } else {
      for (i in 1:length(inputValue)) {
        if (!.checkVpcSettings(settingName = names(inputValue)[i], settingValue = inputValue[[i]])) {
          isValid = FALSE
        }
      }
    }
  } else if (inputName == "eventDataDisplay") {
    if (dataType != "event") {
      message(paste0(
        "WARNING: Unexpected setting encountered. eventDataDisplay ignored when ",
        dataType, " data."
      ))
    } else if (is.list(inputValue) == FALSE) {
      message("ERROR: Unexpected type encountered. eventDataDisplay must be a list")
      isValid = FALSE
    } else {
      for (i in 1:length(inputValue)) {
        if (!.checkVpcDisplay(settingName = names(inputValue)[i], settingValue = inputValue[[i]])) {
          isValid = FALSE
        }
      }
    }
  } else if (inputName == "discreteDataDisplay") {
    if (dataType != "discrete") {
      message(paste0(
        "WARNING: Unexpected setting encountered. discreteDataDisplay ignored when ",
        dataType, " data."
      ))
    } else if (is.list(inputValue) == FALSE) {
      message("ERROR: Unexpected type encountered. discreteDataDisplay must be a list")
      isValid = FALSE
    } else {
      for (i in 1:length(inputValue)) {
        if (!.checkVpcDisplay(settingName = names(inputValue)[i], settingValue = inputValue[[i]])) {
          isValid = FALSE
        }
      }
    }
  } else if (inputName == "continuousDataDisplay") {
    if (dataType != "continuous") {
      message(paste0(
        "WARNING: Unexpected setting encountered. continuousDataDisplay ignored when ",
        dataType, " data."
      ))
    } else if (is.list(inputValue) == FALSE) {
      message("ERROR: Unexpected type encountered. continuousDataDisplay must be a list")
      isValid = FALSE
    } else {
      for (i in 1:length(inputValue)) {
        if (!.checkVpcDisplay(settingName = names(inputValue)[i], settingValue = inputValue[[i]])) {
          isValid = FALSE
        }
      }
    }
  } else if (inputName == "xBinsSettings") {
    if (is.list(inputValue) == FALSE) {
      message("ERROR: Unexpected type encountered. xBinsSettings must be a list")
      isValid = FALSE
    } else {
      for (i in 1:length(inputValue)) {
        if (!.checkVpcBins(settingName = names(inputValue)[i], settingValue = inputValue[[i]])) {
          isValid = FALSE
        }
      }
    }
  } else if (inputName == "yBinsSettings") {
    if (dataType != "count") {
      message(paste0(
        "WARNING: Unexpected setting encountered. yBinsSettings ignored when ",
        dataType, " data."
      ))
    } else if (is.list(inputValue) == FALSE) {
      message("ERROR: Unexpected type encountered. yBinsSettings must be a list")
      isValid = FALSE
    } else {
      for (i in 1:length(inputValue)) {
        if (!.checkVpcBins(settingName = names(inputValue)[i], settingValue = inputValue[[i]])) {
          isValid = FALSE
        }
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

.checkVpcPlotSettings <- function(settingName, settingValue){
  isValid = TRUE
  settingName = settingName
  if (settingName %in% c("legend", "grid", "xlogScale", "ylogScale", "linearInterpolation")) {
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
  } else if (settingName %in% c("xlab", "ylab")) {
    if (!length(settingValue) == 1) {
      message(paste0(
        "ERROR: Unexpected type encountered. `", settingName, "` must be of size 1."
      ))
      isValid = FALSE
    } else if (!is.character(settingValue)) {
      message(paste0(
        "ERROR: Unexpected type encountered. `", settingName, "` must be a string"
      ))
      isValid = FALSE
    }
  } else {
    message(paste0("ERROR: Unexpected value encountered. Unknown setting ", settingName))
    isValid = FALSE
  }
  return(isValid)
}

.checkVpcSettings <- function(settingName, settingValue){
  isValid = TRUE
  settingName = settingName
  if (settingName %in% c("correctedPredictions", "useCensoredData")) {
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
  } else {
    message(paste0("ERROR: Unexpected value encountered. Unknown setting ", settingName))
    isValid = FALSE
  }
  return(isValid)
}

.checkVpcDisplay <- function(settingName, settingValue){
  isValid = TRUE
  settingName = settingName
  if (settingName %in% c(
    "observedData", "censoredData", "empiricalPercentiles", "predictedPercentiles",
    "empiricalProbability", "theoricalProbability", "empiricalCurve",
    "predictedMedian", "predictionInterval",
    "survivalCurve", "meanNumberEventsCurve", "outliersDots", "outlierAreas")) {
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
  } else {
    message(paste0("ERROR: Unexpected value encountered. Unknown setting ", settingName))
    isValid = FALSE
  }
  return(isValid)
}

.checkVpcBins <- function(settingName, settingValue){
  isValid = TRUE
  settingName = settingName
  if (settingName %in% c("bins", "limits", "useFixedBins", "useFixedNbBins")) {
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
  } else if (settingName == "criteria") {
    if (!length(settingValue) == 1) {
      message("ERROR: Unexpected type encountered. criteria must be of size 1.")
      isValid = FALSE
    } else if (!is.element(settingValue, c("equalwidth", "equalsize", "leastsquare"))) {
      message(paste0(
        "ERROR: Unexpected type encountered. criteria must be `equalwidth`, ",
        "`equalsize` or `leastsquare`."
      ))
      isValid = FALSE
    }
  } else if (settingName %in% c("nbBins", "nbDataPoints")){
    if (length(settingValue)  != 1) {
      message(paste0(
        "ERROR: Unexpected type encountered. `", settingName, "` must be of size 1."
      ))
      isValid = FALSE
    } else if (! settingValue >= 0 | ! as.integer(settingValue) == settingValue) {
      message(paste0(
        "ERROR: Unexpected type encountered. `", settingName, "` must be positive integer."
      ))
      isValid = FALSE
    }
  } else if(settingName %in% c("binRange", "nbBinData")){
    if (length(settingValue)  != 2) {
      message(paste0(
        "ERROR: Unexpected type encountered. `", settingName, "` must be of size 2."
      ))
      isValid = FALSE
    } else if (!all(settingValue >= 0) | !all(as.integer(settingValue) == settingValue |
                                              settingValue[1] > settingValue[2])) {
      message(paste0(
        "ERROR: Unexpected type encountered. `", settingName, "` must be positive integers."
      ))
      isValid = FALSE
    }
  } else if(settingName == "fixedBins") {
    if (length(settingValue) < 1) {
      message(
        "ERROR: Unexpected type encountered. fixedBins must be a vector with at least 1 bin."
      )
      isValid = FALSE
    } else if (!all(settingValue >= 0) | !all(is.double(settingValue))) {
      message(paste0(
        "ERROR: Unexpected type encountered. fixedBins must be positive doubles"
      ))
      isValid = FALSE
    }
  } else {
    message(paste0("ERROR: Unexpected value encountered. Unknown setting ", settingName))
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

# Read data --------------------------------------------------------------------
.getObservationData <- function(name, obsType) {
  if (obsType %in% c("discrete", "event")) {
    obsData <- mlx.getObservationInformation()[[name]]
  } else {
    obsFilename <- paste0(
      mlx.getProjectSettings()$directory, "/ChartsData/VisualPredictiveCheck/",
      name, "_observations.txt"
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

.getSimulationData <- function(name) {
  simFilename <- paste0(
    mlx.getProjectSettings()$directory, "/ChartsData/VisualPredictiveCheck/",
    name, "_simulations.txt"
  )
  if (! file.exists(simFilename)) {
    # Get VPC chart data files
    message("Download Charts data.")
    mlx.computeChartsData(exportVPCSimulations = TRUE)
  }
  simData <- .renameColumns(.readDataset(simFilename), "ID", "id")
  return(simData)
}

# Stats on observations and simulations ----------------------------------------
.getBins <- function(dataframe, columnName, bycat = FALSE, settings = NULL, categories = c()) {
  data <- dataframe[[columnName]]
  if (bycat & length(categories) > 0) {
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

.interpolate <- function(data, xName, yNames, xInterp = NULL) {
  if (is.null(xInterp)) {
    x_range <- range(data[[xName]])
    x <- seq(x_range[1], x_range[2], length.out = 10000)
  } else {
    x <- xInterp
  }
  dfOut <- NULL
  for (i in seq(1, length(yNames))) {
    df <- as.data.frame(approx(data[[xName]], data[[yNames[i]]], x))
    names(df) <- c(xName, yNames[i])
    if (is.null(dfOut)) dfOut <- df else dfOut <- merge(dfOut, df, by = c(xName))
  }
  return(dfOut)
}

.addDoseRelTime <- function(data) {
  subjocc <- .getSubjocc()
  doseDf <- getDoseInformation()
  doseDfID <- .renameColumns(aggregate(doseDf$time, by = as.list(subset(doseDf, select = .getSubjocc())), c), "x", "dose")
  data <- merge(data, doseDfID, by = subjocc, sort = F)
  data$time_rel <- apply(
    subset(data, select = c("dose", "time")),
    1,
    function(row) min((row$time - row$dose)[row$time - row$dose >= 0])
  )
  data <- data[!names(data) %in% c("dose")]
  return(data)
}

.addBinsIndex <- function(data, binsValue, columnName, binName = "binIndex") {
  if (!is.null(binsValue)) {
    data <- do.call(rbind, by(
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
    ))
  }
  # } else {
  #   data[binName] <- 1
  # }
  return(data)
}
