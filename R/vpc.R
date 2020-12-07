#' Visual Predictive Check
#'
#' Creates a VPC plot from monolix project
#' @param project (\emph{string}) Monolix project
#' @param time (\emph{string}) (\emph{optional}) in {`time`, `timeSinceLastDose`} (default `time`).
#' `timeSinceLastDose` only possible when there is an `amount` column in the dataset
#' @param obsName (\emph{string}) (\emph{optional}) observation name.
#' By default when several observations in the dataset, the first observation is considered.
#' @param plot (\emph{bool}) (\emph{optional}) if TRUE return plot (ggplot object), else return vpc stat dataframe (default TRUE)
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
#' @param higherPercentile (\emph{int}) (\emph{optional}) (\emph{continuous data}) Higher percentile for empirical and predicted percentiles computation.
#' (default 90) For continuous data only.
#' @param useCorrpred (\emph{bool}) (\emph{optional}) (\emph{continuous data}) if TRUE, pcVPC are computed using Uppsala prediction correction (default FALSE). For continuous data only.
#' @param useCensored (\emph{bool}) (\emph{optional}) (\emph{continuous data}) Choose to use BLQ data or to ignore it to compute the VPC (default TRUE). For continuous data only.
#' @param censoring (\emph{string}) (\emph{optional}) (\emph{continuous data}) in {`simulated`, `loq`}, BLQ data can be simulated, or can be equal to the limit of quantification (LOQ) (default `simulated`). For continuous data only.
#' 
#' \strong{display parameters for plot}
#' @param eventDisplay (\emph{vector(string)}) (\emph{optional}) (\emph{event data}) List of event data curves to display in VPC plots.
#' (by default eventDisplay = c("survivalCurve", "empiricalCurve", "predictionInterval"))
#' You must choose among the following curves:
#' \itemize{
#' \item \code{survivalCurve} Add plot for survival function (Kapan-Meier plot).
#' \item \code{meanNumberEventsCurve} Add plot for mean number of events per individual.
#' \item \code{empiricalCurve} Add empirical curve.
#' \item \code{predictedMedian} Add predicted median.
#' \item \code{predictionInterval} Add prediction intervals given by the model.
#' }
#' @param discreteDisplay (\emph{vector(string)}) (\emph{optional}) (\emph{discrete data}) List of discrete data curves to display in VPC plots.
#' (by default discreteDisplay = c("empiricalProbability", "predictionInterval"))
#' You must choose among the following curves:
#' \itemize{
#' \item \code{empiricalProbability} Add empirical probability.
#' \item \code{predictedMedian} Add theorical probability.
#' \item \code{predictionInterval} Add prediction intervals given by the model.
#' \item \code{outlierAreas} Add red areas indicating empirical percentiles that are outside prediction intervals.
#' }
#' @param continuousDisplay (\emph{vector(string)}) (\emph{optional}) (\emph{continuous data}) List of continuous data curves to display in VPC plots.
#' (by default continuousDisplay = c("empiricalPercentiles", "predictionInterval", "outliersDots", "outlierAreas"))
#' You must choose among the following curves:
#' \itemize{
#' \item \code{observedData} Add observed data.
#' \item \code{censoredData} Add BLQ data if present.
#' \item \code{empiricalPercentiles} Add empirical percentiles for the 10%, 50% and 90% quantiles.
#' \item \code{predictedPercentiles} Add theoretical percentiles for the 10% and 90% quantiles.
#' \item \code{predictionInterval} Add prediction intervals given by the model for the 10% and 90% quantiles (in blue) and the 50% quantile (in pink).
#' \item \code{outliersDots} Add red dots indicating empirical percentiles that are outside prediction intervals.
#' \item \code{outlierAreas} Add red areas indicating empirical percentiles that are outside prediction intervals.
#' }
#' 
#' \strong{parameters for x-axis bins - only for continuous and discrete data}
#' @param xBins.useFixedBins (\emph{bool}) (\emph{optional}) If TRUE, specify manually bin list (default FALSE).
#' @param xBins.fixedBins (\emph{vector(double)}) (\emph{optional}) Define mannually a list of bins. To use when `useFixedBins` is set to TRUE.
#' @param xBins.criteria (\emph{vector(string)}) (\emph{optional}) Choose the bining criteria among `equalwidth`, `equalsize` or `leastsquare` (default leastsquare).
#' @param xBins.useFixedNbBins (\emph{vector(bool)}) (\emph{optional}) If TRUE define a fixed number of bins, else define a range for automatic selection (default FALSE).
#' @param xBins.nbBins (\emph{int}) (\emph{optional}) Define a fixed number of bins (default 10).
#' @param xBins.binRange (\emph{vector(int, int)}) (\emph{optional}) Define a range for the number of bins (default c(5, 20)).
#' @param xBins.nbBinData (\emph{vector(int, int)}) (\emph{optional}) Define a range for the number of data points per bin (default c(10, 200)).
#' @param xBins.nbDataPoints (\emph{int}) (\emph{optional}) (\emph{event data}) Number of data point in event time grid (default 100)
#' 
#' \strong{parameters for y-axis bins - only for discrete countable data}
#' @param yBins.useFixedBins (\emph{bool}) (\emph{optional}) If TRUE, specify manually bin list (default FALSE).
#' @param yBins.fixedBins (\emph{vector(double)}) (\emph{optional}) Define mannually a list of bins. To use when `useFixedBins` is set to TRUE.
#' @param yBins.criteria (\emph{string}) (\emph{optional}) Choose the bining criteria among `equalwidth`, `equalsize` or `leastsquare` (default leastsquare).
#' @param yBins.useFixedNbBins (\emph{vector(bool)}) (\emph{optional}) If TRUE define a fixed number of bins, else define a range for automatic selection (default FALSE).
#' @param yBins.nbBins (\emph{int}) (\emph{optional}) Define a fixed number of bins (default 5).
#' @param yBins.binRange (\emph{vector(int, int)}) (\emph{optional}) Define a range for the number of bins (default c(3, 7)).
#' @param yBins.nbBinData (\emph{vector(int, int)}) (\emph{optional}) Define a range for the number of data points per bin (default c(10, 200)).
#' 
#' \strong{General settings for plot}
#' @param legend (\emph{bool}) (\emph{optional}) Add/remove legend (default FALSE).
#' @param grid (\emph{bool}) (\emph{optional}) Add/remove grid (default FALSE).
#' @param xlogScale (\emph{bool}) (\emph{optional}) Add/remove log scale for x axis (default FALSE).
#' @param ylogScale (\emph{bool}) (\emph{optional}) Add/remove log scale for x axis (default FALSE).
#' @param linearInterpolation (\emph{bool}) (\emph{optional}) If TRUE set piece wise display for prediction intervals, else show bins as rectangular (default TRUE).
#' @param xlab (\emph{string}) (\emph{optional}) Time label (default "time" if time = "time", "time since last dose" if time = "timeSinceLastDose").
#' @param ylab (\emph{string}) (\emph{optional}) y label (default observation name).
#' @param binLimits (\emph{bool}) (\emph{optional}) Add/remove vertical lines on the scatter plots to indicate the bins (default FALSE).
#' 
#' @param vpcTheme (\emph{vpc_theme}) (\emph{optional}) theme to be used in VPC. Expects list of class vpc_theme created with function createVpcTheme()
#' @return If plot set to TRUE, return a ggplot2 object, else return a list with the following info
#' \itemize{
#'   \item \code{vpcPercentiles}: a dataframe with vpc percentiles (bins, empirical and theorical vpc) 
#'   \item \code{observations}: a dataframe with observations
#'   \item \code{obsName}: Name of observation
#'   \item \code{timeName}: Name of time
#' }
#' @importFrom gridExtra grid.arrange
#' @export
#' @examples
#' \dontrun{
#'   p <- vpc(
#'     project = "RsmlxDemo1.mlxtran",
#'     stratSplit = c("sex", "wt"),
#'     stratFilter = list("age" = 1),
#'     stratScale = list("age" = c(65), "wt" = c(70))
#'   )
#'   p <- vpc(project = "RsmlxDemo1.mlxtran",
#'       xBins.useFixedBins = TRUE,
#'       xBins.fixedBins = c(0, 24, 48, 72, 96, 120)
#'   )
#'   p <- vpc(project = "RsmlxDemo1.mlxtran", time = "timeSinceLastDose")
#'   vpcPerc <- vpc(project = "RsmlxDemo1.mlxtran", plot = FALSE)
#'   
#'   continuousDisplay = c("censoredData", "empiricalPercentiles",
#'                         "predictionInterval", "outlierAreas")
#'   p <- vpc(project = "RsmlxDemo1.mlxtran", plot = TRUE,
#'       continuousDisplay = continuousDisplay
#'   )
#' }
#' @seealso \link{vpcStats} \link{createVpcTheme} \link{plotVpc}
vpc <- function(
  project, time = "time", obsName = NULL, plot = TRUE,
  stratSplit = c(), stratFilter = list(), stratScale = list(),
  level = 90, higherPercentile = 90, useCorrpred = FALSE,
  useCensored = TRUE, censoring = "simulated",
  eventDisplay = c("survivalCurve", "empiricalCurve", "predictionInterval"),
  discreteDisplay = c("empiricalProbability", "predictionInterval"),
  continuousDisplay = c("empiricalPercentiles", "predictionInterval", "outliersDots", "outlierAreas"),
  xBins.useFixedBins = FALSE, xBins.fixedBins = c(), xBins.criteria = "leastsquare",
  xBins.useFixedNbBins = FALSE, xBins.nbBins = 10, xBins.binRange = c(5, 20),
  xBins.nbBinData = c(10, 200), xBins.nbDataPoints = 100,
  yBins.useFixedBins = FALSE, yBins.fixedBins = c(), yBins.criteria = "leastsquare",
  yBins.useFixedNbBins = FALSE, yBins.nbBins = 5, yBins.binRange = c(3, 7),
  yBins.nbBinData = c(10, 200),
  legend = FALSE, grid = FALSE, xlogScale = FALSE, ylogScale = FALSE,
  linearInterpolation = TRUE, xlab = "time", ylab = "y", binLimits = FALSE,
  vpcTheme = NULL) {

  # Check Arguments ------------------------------------------------------------
  params <- as.list(match.call(expand.dots = TRUE))[-1]
  
  # load project
  r <- prcheck(project)
  project <- r$project
  
  plot <- check_bool(plot, "plot")
  
  obsName <- check_obs(obsName, "obsName")
  obsType <- .getObsType(obsName)
  dataType <- obsType$type
  dataSubType <- obsType$subType

  # Check and initialize curves display
  if (dataType == "event") {
    eventCurves <- c("survivalCurve", "meanNumberEventsCurve", "empiricalCurve",
                     "predictedMedian", "predictionInterval")
    eventDisplay <- check_display(eventDisplay, "eventDisplay", eventCurves)
  } else {
    check_ignoring_arg_on_condition("eventDisplay", params, paste0(dataType, " data"))
    eventDisplay <- c()
  }
  
  if (dataType == "discrete") {
    discreteCurves <- c("empiricalProbability", "predictedMedian",
                        "predictionInterval", "outlierAreas")
    discreteDisplay <- check_display(discreteDisplay, "discreteDisplay", discreteCurves)
  } else {
    check_ignoring_arg_on_condition("discreteDisplay", params, paste0(dataType, " data"))
    discreteDisplay <- c()
  }

  if (dataType == "continuous") {
    continuousCurves <- c(
      "observedData", "censoredData", "empiricalPercentiles", "predictedPercentiles",
      "predictionInterval", "outliersDots", "outlierAreas"
    )
    continuousDisplay <- check_display(continuousDisplay, "continuousDisplay", continuousCurves)
  } else {
    check_ignoring_arg_on_condition("continuousDisplay", params, paste0(dataType, " data"))
    continuousDisplay <- c()
  }

  # Check and initialize plot settings
  legend <- check_bool(legend, "legend")
  grid <- check_bool(grid, "grid")
  xlogScale <- check_bool(xlogScale, "xlogScale")
  ylogScale <- check_bool(ylogScale, "ylogScale")
  linearInterpolation <- check_bool(linearInterpolation, "linearInterpolation")
  binLimits <- check_bool(binLimits, "binLimits")
  xlab <- check_char(xlab, "xlab")
  if (!"xlab" %in% names(params))
    xlab <- ifelse(time == "timeSinceLastDose", "time since last dose", "time")
  ylab <- check_char(ylab, "ylab")
  if (!"ylab" %in% names(params))
    ylab <- obsName

  # Check and initialize vpc theme
  if(is.null(vpcTheme)) {
    defaultTheme <- list(
      continuous = continuous_theme(),
      discrete = discrete_theme(),
      event = event_theme()
    )
    vpcTheme <- defaultTheme[[dataType]]
  }
  check_object_class(vpcTheme, "vpcTheme", "vpc_theme")

  if (!"xBins.nbBinData" %in% params & dataType == "discrete")
    xBins.nbBinData <- c(5, 30)

  displaySettings <- list(
    observedData = "observedData" %in% continuousDisplay,
    censoredData = "censoredData" %in% continuousDisplay,
    empiricalData = any(
      "empiricalPercentiles" %in% continuousDisplay,
      "empiricalProbability" %in% discreteDisplay,
      "empiricalCurve" %in% eventDisplay
    ),
    theoreticalData = any(
      "predictedPercentiles" %in% continuousDisplay,
      "theoreticalProbability" %in% discreteDisplay,
      "predictedMedian" %in% eventDisplay
    ),
    predictionInterval = any(
      "predictionInterval" %in% continuousDisplay,
      "predictionInterval" %in% discreteDisplay,
      "predictionInterval" %in% eventDisplay
    ),
    outliersDots = "outliersDots" %in% continuousDisplay,
    outlierAreas = any(
      "outlierAreas" %in% continuousDisplay, "outlierAreas" %in% discreteDisplay
    )
  )
  curvesDisplay <- names(displaySettings[displaySettings == TRUE])

  ## Generate Chart data -------------------------------------------------------
  
  xBinsSettings <- getBinsSettings(
    is.fixedBins = xBins.useFixedBins, fixedBins = xBins.fixedBins,
    criteria = xBins.criteria, is.fixedNbBins = xBins.useFixedNbBins,
    nbBins = xBins.nbBins, binRange = xBins.binRange, nbBinData = xBins.nbBinData,
    nbDataPoints = xBins.nbDataPoints
  )
  yBinsSettings <- getBinsSettings(
    is.fixedBins = yBins.useFixedBins, fixedBins = yBins.fixedBins,
    criteria = yBins.criteria, is.fixedNbBins = yBins.useFixedNbBins,
    nbBins = yBins.nbBins, binRange = yBins.binRange, nbBinData = yBins.nbBinData
  )
  
  vpcStats <- vpcStats(
    project, time = time, obsName = obsName,
    stratSplit = stratSplit, stratFilter = stratFilter, stratScale = stratScale,
    level = level, higherPercentile = higherPercentile,
    useCorrpred = useCorrpred, useCensored = useCensored, censoring = censoring,
    event.averageNumberEvents = "averageEventNumber" %in% eventDisplay,
    xBinsSettings = xBinsSettings, yBinsSettings = yBinsSettings
  )

  ## Plot VPC -------------------------------------------------------------------
  if (!plot) return(vpcStats)

  if (dataType %in% c("discrete", "continuous")) {
    p <- plotVpc(
      vpcStats$vpcPercentiles, vpcStats$observations, vpcStats$obsName,
      vpcStats$timeName, curvesDisplay = curvesDisplay, legend = legend,
      grid = grid, xlogScale = xlogScale, ylogScale = ylogScale,
      binLimits = binLimits, linearInterpolation = linearInterpolation,
      xlab = xlab, ylab = ylab, theme = vpcTheme
    )
  } else {
    survivalCurve <- "survivalCurve" %in% eventDisplay
    averageEventNumber <- "meanNumberEventsCurve" %in% eventDisplay
    p <- plotEvent(vpcStats, survivalCurve = survivalCurve, averageEventNumber = averageEventNumber,
                   curvesDisplay = curvesDisplay, legend = legend,
                   grid = grid, xlogScale = xlogScale, ylogScale = ylogScale,
                   binLimits = binLimits, linearInterpolation = linearInterpolation,
                   xlab = xlab, theme = vpcTheme)
  }
  return(p)
}

plotEvent <- function(vpcStats, survivalCurve, averageEventNumber, ...) {
  vpcData <- vpcStats$vpcPercentiles
  dataSF <- subset(
    vpcData,
    select = c("split", vpcStats$timeName,
               names(vpcData)[grepl("survivalFunction", names(vpcData))])
  )
  dataAEN <- subset(
    vpcData,
    select = c("split", vpcStats$timeName,
               names(vpcData)[grepl("averageEventNumber", names(vpcData))])
  )
  print(survivalCurve)
  print(averageEventNumber)
  if (survivalCurve) {
    p1 <- plotVpc(
      dataSF, vpcStats$observations, vpcStats$obsName, vpcStats$timeName,
      ylab = "Survival Function", ...
    )
    if (averageEventNumber) {
      p2 <- plotVpc(
        dataAEN, vpcStats$observations, vpcStats$obsName, vpcStats$timeName,
        ylab = "Mean number of events per subject", ...
      )
      p <- grid.arrange(p1, p2, ncol=2)
    } else {
      p <- p1
    }
  } else if (averageEventNumber) {
    p <- plotVpc(
      dataAEN, vpcStats$observations, vpcStats$obsName, vpcStats$timeName,
      ylab = "Mean number of events per subject", ...
    )
  }
  return(p)
}

check_ignoring_arg_on_condition <- function(argName, params, condition) {
  if ("argName" %in% params)
    warning("When ", condition, ", `", argName, "` is ignored.", call. = FALSE)
  return(c())
}

check_arg_on_condition <- function(arg, argName, condition) {
  if (!is.null(arg)) {
    warning("When ", condition, " `", argName, "` is ignored.", call. = FALSE)
    arg <- NULL
  }
  return(arg)
}
