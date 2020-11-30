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
