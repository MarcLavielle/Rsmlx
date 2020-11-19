computeContinuousVPC <- function(obsDf, obsName, simDf, simName, higherPercentile=90, level=90) {
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

computeDiscreteVPC <- function(obsDf, obsName, simDf, simName, higherPercentile=90, level=90, categories = NULL) {
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

computeEventVPC <- function(obsDf, obsName, simDf, simName, timeName, eventType, nbDataPoints = 100, level = 90, averageEventNumber = FALSE) {
  exact <- (eventType == "exactEvent")
  
  tRange <- range(obsDf$time)
  timeGrid <- seq(tRange[1], tRange[2], (tRange[2] - tRange[1]) / (nbDataPoints - 1))
  if (exact) timeGrid <- c(timeGrid, obsDf[obsDf[obsName] == 1, timeName])
  timeGrid <- sort(unique(timeGrid))
  obsDf <- .addEventIdColumn(obsDf, obsName)
  simDf <- .addEventIdColumn(simDf, simName)
  empiricalData <- computeSurvivalCurves(obsDf, timeName, obsName, timeGrid, exact, averageEventNumber) 

  simData <- do.call(rbind, by (
    simDf,
    subset(simDf, select = c("rep")),
    function(df) {
      res <- computeSurvivalCurves(df, timeName, simName, timeGrid, exact, averageEventNumber)
      res <- cbind(rep = unique(df$rep), res)
      return(res)
    }
  ))
  
  theoreticalData <- do.call(rbind, by (
    simData,
    subset(simData, select = c(timeName)),
    function(df) {
      res <- as.data.frame(sapply(
        df[!names(df) %in% c(timeName, "rep")],
        .stats, (100 - level) / 2
      ))
      res <- cbind(list(stat = rownames(res), time = unique(df$time)), res)
      res <- reshape(res, idvar = c(timeName), timevar = "stat", direction = "wide")
      return(res)
    }
  ))
  names(theoreticalData) <- sapply(
    names(theoreticalData), function(x) paste(strsplit(x, "\\.")[[1]], collapse = "_")
  )
  
  res <- merge(empiricalData, theoreticalData, by = c("time"))
  return(res)
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
