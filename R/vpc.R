#' Visual Predictive Check
#'
#' Creates a VPC plot from monolix project
#' @param project Monolix project
#' @param time emph{[Optional]} in {`time`, `timeSinceLastDose`} (default `time`).
#' `timeSinceLastDose` only possible when there is an `amount` column in the dataset
#' @param obsName emph{[Optional]} observation name. By default the first observation is considered.
#' @param plot emph{[Optional]} if TRUE return plot (ggplot object), else return vpc stat dataframe (default TRUE)
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
#' 
#' \strong{display parameters for plot}
#' @param eventDisplay emph{[Optional]}[event data] (vector[string]) List of event data curves to display in VPC plots.
#' (by default eventDisplay = c("survivalCurve", "empiricalCurve", "predictionInterval"))
#' You must choose among the following curves:
#' \itemize{
#' \item \code{survivalCurve} Add plot for survival function (Kapan-Meier plot) (default TRUE).
#' \item \code{meanNumberEventsCurve} Add plot for mean number of events per individual (default FALSE).
#' \item \code{empiricalCurve} Add empirical curve (default TRUE).
#' \item \code{predictedMedian} Add predicted median (default FALSE).
#' \item \code{predictionInterval} Add prediction intervals given by the model (default TRUE).
#' }
#' @param discreteDisplay emph{[Optional]}[discrete data] (vector[string]) List of discrete data curves to display in VPC plots.
#' (by default discreteDisplay = c("empiricalProbability", "predictionInterval"))
#' You must choose among the following curves:
#' \itemize{
#' \item \code{empiricalProbability} Add empirical probability.
#' \item \code{predictedMedian} Add theorical probability.
#' \item \code{predictionInterval} Add prediction intervals given by the model.
#' \item \code{outlierAreas} Add red areas indicating empirical percentiles that are outside prediction intervals.
#' }
#' @param continuousDisplay emph{[Optional]}[continuous data] (vector[string]) List of continuous data curves to display in VPC plots.
#' (by default continuousDisplay = c("empiricalPercentiles", "predictionInterval", "outliersDots", "outlierAreas"))
#' You must choose among the following curves:
#' \itemize{
#' \item \code{observedData} Add observed data (default FALSE).
#' \item \code{censoredData} Add BLQ data if present (default FALSE).
#' \item \code{empiricalPercentiles} Add empirical percentiles for the 10%, 50% and 90% quantiles.
#' \item \code{predictedPercentiles} Add theoretical percentiles for the 10% and 90% quantiles.
#' \item \code{predictionInterval} Add prediction intervals given by the model for the 10% and 90% quantiles (in blue) and the 50% quantile (in pink).
#' \item \code{outliersDots} Add red dots indicating empirical percentiles that are outside prediction intervals.
#' \item \code{outlierAreas} Add red areas indicating empirical percentiles that are outside prediction intervals.
#' }
#' 
#' \strong{parameters for x-axis bins - only for continuous and discrete data}
#' @param xBins.useFixedBins emph{[Optional]} (boolean) If TRUE, specify manually bin list (default FALSE).
#' @param xBins.fixedBins emph{[Optional]} (list[integer]) Define mannually a list of bins. To use when `useFixedBins` is set to TRUE.
#' @param xBins.criteria emph{[Optional]} Choose the bining criteria among `equalwidth`, `equalsize` or `leastsquare` (default leastsquare).
#' @param xBins.useFixedNbBins emph{[Optional]} (boolean) If TRUE define a fixed number of bins, else define a range for automatic selection (default FALSE).
#' @param xBins.nbBins emph{[Optional]} (integer) Define a fixed number of bins (default 10).
#' @param xBins.binRange emph{[Optional]} (list[integer]) Define a range for the number of bins (default c(5, 20)).
#' @param xBins.nbBinData emph{[Optional]} (list[integer]) Define a range for the number of data points per bin (default c(10, 200)).
#' @param xBins.nbDataPoints emph{[Optional]}[event data] (integer) Number of data point in event time grid (default 100)
#' 
#' \strong{parameters for y-axis bins - only for discrete countable data}
#' @param yBins.useFixedBins emph{[Optional]} (boolean) If TRUE, specify manually bin list (default FALSE).
#' @param yBins.fixedBins emph{[Optional]} (list[integer]) Define mannually a list of bins. To use when `useFixedBins` is set to TRUE.
#' @param yBins.criteria emph{[Optional]} Choose the bining criteria among `equalwidth`, `equalsize` or `leastsquare` (default leastsquare).
#' @param yBins.useFixedNbBins emph{[Optional]} (boolean) If TRUE define a fixed number of bins, else define a range for automatic selection (default FALSE).
#' @param yBins.nbBins emph{[Optional]} (integer) Define a fixed number of bins (default 5).
#' @param yBins.binRange emph{[Optional]} (list[integer]) Define a range for the number of bins (default c(3, 7)).
#' @param yBins.nbBinData emph{[Optional]} (list[integer]) Define a range for the number of data points per bin (default c(10, 200)).
#' 
#' \strong{General settings for plot}
#' @param legend emph{[Optional]} (boolean) Add/remove legend (default FALSE).
#' @param grid emph{[Optional]} (boolean) Add/remove grid (default FALSE).
#' @param xlogScale emph{[Optional]} (boolean) Add/remove log scale for x axis (default FALSE).
#' @param ylogScale emph{[Optional]} (boolean) Add/remove log scale for x axis (default FALSE).
#' @param linearInterpolation emph{[Optional]} (boolean) If TRUE set piece wise display for prediction intervals, else show bins as rectangular (default TRUE).
#' @param xlab emph{[Optional]} (str) Time label (default "time" if time = "time", "time since last dose" if time = "timeSinceLastDose").
#' @param ylab emph{[Optional]} (str) y label (default observation name).
#' @param binLimits emph{[Optional]} (boolean) Add/remove vertical lines on the scatter plots to indicate the bins (default FALSE).
#' 
#' @param vpcTheme emph{[Optional]} theme to be used in VPC. Expects list of class vpc_theme created with function createVpcTheme()
#' @return If plot set to TRUE, return a ggplot2 object, else return chart Data

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
#'   p <- vpc(project = "RsmlxDemo1.mlxtran", plot = TRUE,
#'       continuousDisplay = c("censoredData", "empiricalPercentiles", "predictionInterval", "outlierAreas")
#'   )
#' }
#' @seealso \link{vpcStats} \link{createVpcTheme} \link{plotVPC}
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
  
  if (missing(project)) stop("`project` is missing")
  # load project
  r <- prcheck(project)
  project <- r$project
  
  if (is.null(plot) | !is.logical(plot)) plot = TRUE

  obsnames <- mlx.getObservationInformation()$name
  if (is.null(obsName)) obsName <- obsnames[1]
  if (!is.element(obsName, obsnames)) {
    obsName <- obsnames[1]
    warning(paste0("Invalid argument `obsName`. `obsname` set to \"", obsName, "\""))
  }
  obsType <- .getObsType(obsName)
  dataType <- obsType$type
  dataSubType <- obsType$subType

  # Check and initialize curves display
  if (dataType != "event") {
    if ("eventDisplay" %in% params)
      warning(paste0("`eventDisplay` ignored in in case of ", dataType, " data."))
    eventDisplay <- c()
  } else {
    eventCurves <- c("survivalCurve", "meanNumberEventsCurve", "empiricalCurve", "predictedMedian", "predictionInterval")
    if (!is.vector(eventDisplay)) stop("`eventDisplay` must be a vector.")
    if (!all(is.element(eventDisplay, eventCurves)))
      warning(paste0(
        "`eventDisplay` values must be in ",
        paste(eventCurves, collapse = ", "), " list. Ignore invalid curves."
      ))
    eventDisplay <- eventDisplay[eventDisplay %in% eventCurves]
  }
  if (dataType != "discrete") {
    if ("discreteDisplay" %in% params)
      warning(paste0("`discreteDisplay` ignored in in case of ", dataType, " data."))
    discreteDisplay <- c()
  } else {
    discreteCurves <- c("empiricalProbability", "predictedMedian", "predictionInterval", "outlierAreas")
    if (!is.vector(discreteDisplay)) stop("`discreteDisplay` must be a vector.")
    if (!all(is.element(discreteDisplay, discreteCurves)))
        warning(paste0(
          "`discreteDisplay` values must be in ",
          paste(discreteCurves, collapse = ", "), " list. Ignore invalid curves."
        ))
    discreteDisplay <- discreteDisplay[discreteDisplay %in% discreteCurves]
  }
  if (dataType != "continuous") {
    if ("continuousDisplay" %in% params)
      warning(paste0("`continuousDisplay` ignored in in case of ", dataType, " data."))
    continuousDisplay <- c()
  } else {
    continuousCurves <- c(
      "observedData", "censoredData", "empiricalPercentiles", "predictedPercentiles",
      "predictionInterval", "outliersDots", "outlierAreas"
    )
    if (!is.vector(continuousDisplay)) stop("`continuousDisplay` must be a vector.")
    if (!all(is.element(continuousDisplay, continuousCurves)))
      warning(paste0(
        "`continuousDisplay` values must be in ",
        paste(continuousCurves, collapse = ", "), " list. Ignore invalid curves."
      ))
    continuousDisplay <- continuousDisplay[continuousDisplay %in% continuousCurves]
  }

  # Check and initialize plot settings
  if (!is.logical(legend)) {
    warning("`legend` must be a boolean. Set `legend` to FALSE.")
    legend <- FALSE
  }
  if (!is.logical(grid)) {
    warning("`grid` must be a boolean. Set `grid` to FALSE.")
    grid <- FALSE
  }
  if (!is.logical(xlogScale)) {
    warning("`xlogScale` must be a boolean. Set `xlogScale` to FALSE.")
    xlogScale <- FALSE
  }
  if (!is.logical(ylogScale)) {
    warning("`ylogScale` must be a boolean. Set `ylogScale` to FALSE.")
    ylogScale <- FALSE
  }
  if (!is.logical(linearInterpolation)) {
    warning("`linearInterpolation` must be a boolean. Set `linearInterpolation` to TRUE")
    linearInterpolation <- TRUE
  }
  if (!is.logical(binLimits)) {
    warning("`binLimits` must be a boolean. Set `binLimits` to FALSE")
    binLimits <- FALSE
  }
  if (!is.character(xlab)) {
    warning("`xlab` must be a string Set `xlab` to time.")
  }
  if ("xlab" %in% names(params)) {
    if (time == "timeSinceLastDose") {
      xlab <- "time since last dose"
    } else {
      xlab <- "time"
    }
  }
  if (!is.character(ylab)) {
    warning(paste0("`ylab` must be a string Set `ylab` to ", obsName, "."))
    ylab <- obsName
  }

  # Check and initialize vpc theme
  if(is.null(vpcTheme)) {
    if (dataType == "continuous") {
      vpcTheme <- continuous_theme()
    } else if (dataType == "discrete") {
      vpcTheme <- discrete_theme()
    } else {
      vpcTheme <- event_theme()
    }
  }
  if(class(vpcTheme) != "vpc_theme") stop("`vpcTheme` must be a `vpc_theme` object")
  
  if (!"xBins.nbBinData" %in% params & dataType == "discrete")
    xBins.nbBinData <- c(5, 30)

  displaySettings <- list(
    binLimits = binLimits,
    grid = grid,
    legend = legend,
    xlogScale = xlogScale,
    ylogScale = ylogScale,
    linearInterpolation = linearInterpolation,
    xlab = xlab,
    ylab = ylab,
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
    ),
    survivalCurve = "survivalCurve" %in% eventDisplay,
    meanNumberEventsCurve = "meanNumberEventsCurve" %in% eventDisplay
  )

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
    event.averageNumberEvents = displaySettings$meanNumberEventsCurve,
    xBinsSettings = xBinsSettings, yBinsSettings = yBinsSettings
  )

  # Plot VPC -------------------------------------------------------------------
  if (!plot) return(vpcStats)

  if (dataType %in% c("discrete", "continuous")) {
    p <- plotVpc(
      vpcStats$vpcPercentiles, vpcStats$observations, vpcStats$obsName,
      vpcStats$timeName, displaySettings, vpcTheme
    )
  } else {
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
    if (displaySettings$survivalCurve) {
      displaySettings$ylab <- "Survival Function"
      p1 <- plotVpc(
        dataSF, vpcStats$observations, vpcStats$obsName, vpcStats$timeName,
        displaySettings, vpcTheme
      )
      if (displaySettings$meanNumberEventsCurve) {
        displaySettings$ylab <- "Mean number of events per subject"
         p2 <- plotVpc(
          dataAEN, vpcStats$observations, vpcStats$obsName, vpcStats$timeName,
          displaySettings, vpcTheme
        )
        p <- grid.arrange(p1, p2, vp, sc, ncol=2)
      } else {
        p <- p1
      }
    } else if (displaySettings$meanNumberEventsCurve) {
      displaySettings$ylab <- "Mean number of events per subject"
      p <- plotVpc(
        dataAEN, vpcStats$observations, vpcStats$obsName, vpcStats$timeName,
        displaySettings, vpcTheme
      )
    }
  }
  return(invisible(p))
}
