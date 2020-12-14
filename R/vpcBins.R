#' Bins settings
#'
#' @param is.fixedBins (\emph{bool}) (\emph{optional}) If TRUE, specify manually bin vector (default FALSE).
#' @param fixedBins (\emph{vector(double)}) (\emph{optional}) Define mannually a vector of bins. To use when `useFixedBins` is set to TRUE.
#' @param criteria (\emph{string}) (\emph{optional}) Choose the bining criteria among `equalwidth`, `equalsize` or `leastsquare` (default leastsquare).
#' @param is.fixedNbBins (\emph{bool}) (\emph{optional}) If TRUE define a fixed number of bins, else define a range for automatic selection (default FALSE).
#' @param nbBins (\emph{int}) (\emph{optional}) Define a fixed number of bins (default 10).
#' @param binRange (\emph{vector(int, int)}) (\emph{optional}) Define a range for the number of bins (default c(5, 20)).
#' @param nbBinData (\emph{vector(int, int)}) (\emph{optional}) Define a range for the number of data points per bin (default c(10, 200)).
#' @param nbDataPoints (\emph{int}) (\emph{optional}) (\emph{event data}) Number of data point in event time grid (default 100)
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
  
  .check_bool(is.fixedBins, "is.fixedBins")
  if (is.fixedBins) .check_pos_double(fixedBins, "fixedBins")
  .check_in_vector(criteria, "criteria", c("equalwidth", "equalsize", "leastsquare"))
  .check_bool(is.fixedNbBins, "is.fixedNbBins")
  .check_integer(nbBins, "nbBins")
  .check_range(binRange, "binRange")
  .check_integer(binRange, "binRange")
  .check_range(nbBinData, "nbBinData")
  .check_integer(nbBinData, "nbBinData")
  .check_integer(nbDataPoints, "nbDataPoints")

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
#' @param data (\emph{vector(double)}) vector of data (times or obs) to compute bins.
#' @param split (\emph{vector}) split category associated with data (no split by default).
#' @param type in {"continuous", "categorical"} type of data ("continuous" by default).
#' @param binsSettings (\emph{binSettingsClass}) (\emph{optional}) Bins settings.
#' @param is.binsName (\emph{bool}) (\emph{optional}) if TRUE add a column binsName in P(x<dataName<y) (default FALSE).
#' @param dataName (\emph{string}) (\emph{optional}) dataName when is.binsName is TRUE.
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
  .check_double(data, "data")
  if (is.null(split)) split <- rep("All", length(data))
  .check_length(data, "data", split, "split")
  .check_in_vector(type, "type", c("continuous", "categorical"), type = "warning")
  if (is.null(type)) type <- "continuous"
  binsSettings <- .check_bins(binsSettings, "binsSettings")
  .check_bool(is.binsName, "is.binsName")
  .check_char(dataName, "dataName")

  params <- list(
    useFixedBins = binsSettings$is.fixedBins, fixedBins = binsSettings$fixedBins,
    criteria = binsSettings$criteria, useFixedNbBins = binsSettings$is.fixedNbBins,
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
    bins[bins$bins_start == bins$bins_middles,]$binsName <- sapply(
      bins[bins$bins_start == bins$bins_middles,]$bins_start,
      function(x) paste0("P(", dataName, "=", x, ")")
    )
    bins[bins$bins_start < bins$bins_middles,]$binsName <- mapply(
      function(x, y) paste0("P(", ceiling(x), "<=", dataName, "<=", trunc(y), ")"),
      bins[bins$bins_start < bins$bins_middles,]$bins_start,
      bins[bins$bins_start < bins$bins_middles,]$bins_stop
    )
  }
  return(bins)
}

.getBins <- function(data, settings = NULL, bycat = FALSE, categories = c()) {
  if (bycat) {
    if (length(categories) == 0) categories <- unique(data)
    binsData <- data.frame(list(
      bins_start = categories, bins_stop = categories,
      bins_middles = categories, bins = seq(1, length(categories))
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
        bins_start = utils::head(b, -1),
        bins_stop = utils::tail(b, -1),
        bins_middles = as.vector(tapply(data, factor, mean)),
        bins = seq(1, length(b) - 1)
      ))
    } else {
      # get bins
      options <- list(
        criteria = settings$criteria, usefixednb = settings$useFixedNbBins,
        fixednb = settings$nbBins,
        estimatednb = settings$binRange, nbbindata = settings$nbBinData
      )
      binsList <- mlx.computeBins(data = data, options = options)
      binsData <- data.frame(list(
        bins_start = utils::head(binsList$values, -1),
        bins_stop = utils::tail(binsList$values, -1),
        bins_middles = binsList$middles,
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
            breaks = c(bins$bins_start, utils::tail(bins$bins_stop, 1)),
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

.check_bins <- function(binsSettings, argname) {
  if (is.null(binsSettings)) {
    binsSettings <- getBinsSettings()
  }
  .check_object_class(binsSettings, "binsSettings", "binSettingsClass")
  return(binsSettings)
}

.check_length <- function(arg1, arg1name, arg2, arg2name) {
  if (length(arg1) != length(arg2))
    stop("`", arg1name, "` and `", arg2name, "` must have the same length.", .call = FALSE)
}
