.computeContinuousVPC <- function(obsDf, obsName, simDf, simName, higherPercentile=90, level=90) {
  empiricalData <- as.data.frame(do.call(
        rbind, tapply(obsDf[[obsName]], obsDf$binIndex, .stats, higherPercentile, header = "empirical")
      ))
  empiricalData <- cbind(list(bins = rownames(empiricalData)), empiricalData)
  simData <- do.call(rbind, by (
    simDf,
    subset(simDf, select = "rep"),
    function(df) {
      res <- as.data.frame(do.call(rbind, tapply(df[[simName]], df$binIndex, .stats, higherPercentile, header = "theoretical")))
      res <- cbind(list(bins = rownames(res), rep = unique(df$rep)), res)
      return(res)
    }
  ))
  theoreticalData <- do.call(rbind, by (
    simData,
    subset(simData, select = c("bins")),
    function(df) {
      res <- as.data.frame(sapply(
        df[!names(df) %in% c("bins", "rep")],
        .stats, (100 - level) / 2
      ))
      res <- cbind(list(stat = rownames(res), bins = unique(df$bins)), res)
      res <- reshape(res, idvar = c("bins"), timevar = "stat", direction = "wide")
      return(res)
    }
  ))
  names(theoreticalData) <- sapply(
    names(theoreticalData), function(x) paste(strsplit(x, "\\.")[[1]], collapse = "_")
  )
  res <- merge(empiricalData, theoreticalData, by = c("bins"))
  return(res)
}

.computeDiscreteVPC <- function(obsDf, obsName, simDf, simName, higherPercentile=90, level=90, categories = NULL) {
  if (is.null(categories)) categories <- unique(obsDf$category)
  empiricalData <- .discreteDataProba(obsDf, obsName, categories = categories, name = "empirical")
  simData <- do.call(rbind, by (
    simDf,
    subset(simDf, select = "rep"),
    function(df) {
      res <- .discreteDataProba(df, simName, categories = categories)
      res <- cbind(list(rep = unique(df$rep)) , res)
      return(res)
    }
  ))
  theoreticalData <- do.call(rbind, by (
    simData,
    subset(simData, select = c("bins", "category")),
    function(df) {
      res <- as.data.frame(sapply(
        df[!names(df) %in% c("bins", "rep", "category")],
        .stats, (100 - level) / 2
      ))
      res <- cbind(list(stat = rownames(res), bins = unique(df$bins)), res)
      res <- reshape(res, idvar = c("bins"), timevar = "stat", direction = "wide")
      res$category <- unique(df$category)
      return(res)
    }
  ))
  names(theoreticalData) <- sapply(
    names(theoreticalData), function(x) paste(strsplit(x, "\\.")[[1]], collapse = "_")
  )
  res <- merge(empiricalData, theoreticalData, by = c("bins", "category"))
}

.computeEventVPC <- function(obsDf, obsName, simDf, simName, nbDataPoints = 100, level = 90) {
  subjocc <- .getSubjocc()
  tRange <- range(obsDf$time)
  timeGrid <- sort(unique(c(
    obsDf[obsDf[obsName] == 1, "time"],  # obs events
    seq(tRange[1], tRange[2], (tRange[2] - tRange[1]) / (nbDataPoints - 1))  # grid
  )))
  
  empiricalData <- .computeSurvivalCurve(obsDf, obsName, timeGrid) 
  
  simData <- do.call(rbind, by (
    subset(simDf, select = c(subjocc, "rep", simName, "time", "censored")),
    subset(simDf, select = c("rep")),
    function(df) {
      res <- .computeSurvivalCurve(df, simName, timeGrid)
      res <- cbind(rep = unique(df$rep), res)
      return(res)
    }
  ))
  
  theoreticalData <- do.call(rbind, by (
    simData,
    subset(simData, select = c("time")),
    function(df) {
      res <- as.data.frame(sapply(
        df[!names(df) %in% c("time", "rep")],
        .stats, (100 - level) / 2
      ))
      res <- cbind(list(stat = rownames(res), time = unique(df$time)), res)
      res <- reshape(res, idvar = c("time"), timevar = "stat", direction = "wide")
      return(res)
    }
  ))
  names(theoreticalData) <- sapply(
    names(theoreticalData), function(x) paste(strsplit(x, "\\.")[[1]], collapse = "_")
  )
  
  res <- merge(empiricalData, theoreticalData, by = c("time"))
  return(res)
}

.computeSurvivalCurve <- function(data, eventName, timeGrid) {
  subjocc <- .getSubjocc()
  data <- data[order(subset(data, select = time)),]
  data$kEvent <- 0
  subjoccs <- unique(subset(data, select = subjocc))
  for (i in seq_len(nrow(subjoccs))) {
    so <- subjoccs[i, ]
    df <- subset(
      data,
      apply(subset(data, select = subjocc), 1, function(x) all(x==so))
    )
    df$kEvent[df[eventName] == 1] <- seq_len(nrow(df[df[eventName] == 1,]))
    df$kEvent[df$censored] <- - (max(df$kEvent) + 1)
    data[apply(subset(data, select = subjocc), 1, function(x) all(x==so)),] <- df
  }

  times <- sort(unique(c(
    c(0, data[data[eventName] == 1, "time"]),
    timeGrid
  )))
  nbMaxEvents <- max(data$kEvent)
  S <- sapply(
    seq_len(nbMaxEvents),
    function(k) {
      cumprod(sapply(
        times,
        function(t) {
          nbSubject <- nrow(unique(subset(
            data,
            time >= t & (kEvent == k | (censored == 1 & kEvent >= -k)),
            select = subjocc
          )))
          res <- 1 - ifelse(
            nbSubject > 0,
            nrow(subset(data, time == t & kEvent == k)) / nbSubject,
            0
          )
          return(res)
        }
      ))
    }
  )
  survivalFunction <- S[, 1]
  averageEventNumber <- rowSums(1 - S)
  res <- as.data.frame(list(
    time = timeGrid,
    survivalFunction = survivalFunction[times %in% timeGrid],
    averageEventNumber = averageEventNumber[times %in% timeGrid]
  ))
  return(res)
  
}

.addTTECensoredColumn <- function(data) {
  subjocc <- .getSubjocc()
  data <- data[order(subset(data, select = time)),]
  data$censored <- 0
  subjoccs <- unique(subset(data, select = subjocc))
  for (i in seq_len(nrow(subjoccs))) {
    so <- subjoccs[i, ]
    df <- subset(
      data,
      apply(subset(data, select = subjocc), 1, function(x) all(x==so))
    )
    df$censored[df[eventName] == 0 & df$time == max(df$time)] <- 1
    data[apply(subset(data, select = subjocc), 1, function(x) all(x==so)),] <- df
  }
  return(data)
  
}

.discreteDataProba <- function(data, discreteName, categories = NULL, name = NULL) {
  res <- as.data.frame(tapply(data[[discreteName]], list(data$binIndex, data$category), length))
  if (!is.null(categories))
    for (c in setdiff(categories, names(res))) res[c] <- 0
  res[is.na(res)] <- 0
  res <- res / rowSums(res)
  name <- ifelse(is.null(name), "propCategory", paste0("propCategory_", name))
  res <- reshape(
    res, idvar = "bins", ids = row.names(res), times = names(res),
    timevar = "category", varying = names(res),
    v.names = name, direction = "long"
  )
  return(res)
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

.applyCorrectedPrediction <- function(obsData, simData) {
  subjocc <- .getSubjocc()
  predFilename <- paste0(mlx.getProjectSettings()$directory, "/predictions.txt")
  predData <- subset(.readDataset(predFilename), select = c(subjocc, "time", "popPred"))
  # merge popPred
  obsData <- merge(obsData, predData, by = c(subjocc, "time"), sort = F)
  obsData <- do.call(rbind, by(
    obsData,
    data$split,
    function(df) {
      predictBin <- tapply(df$popPred, df$binIndex, median)
      obsData$predBin <- sapply(df$binIndex, function(index) predictBin[[index]])
      obsData$pcNorm <- ifelse(df$popPred == 0, 1, df$predBin / df$popPred)
      obsData[paste0(obsName, "_pc")] <- df[[obsName]] * df$pcNorm
      return(df)
    }
  ))
  
  simData <- merge(
    simData,
    subset(obsData, select = c(subjocc, "time", "pcNorm")),
    by = c(subjocc, "time"), sort = F
  )
  simData[paste0(simName, "_pc")] <- simData[simName] * simData$pcNorm
  return(sim = simData, obs = obsData)
}
