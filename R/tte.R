computeSurvivalCurves <- function(data, timeName, eventName, timeGrid, exact, averageEventNumber = FALSE) {
  nbMaxEvents <- max(data$nEvent)
  # first event
  k <- 1
  data1 = data[
    (data[eventName] == 1 & data$nEvent == k) |
      (data[eventName] == 0 & data$nEvent < k ),
  ]
  S <- .computeSurvivalCurve(data1, timeName, eventName, timeGrid, exact)
  S <- .renameColumns(S, "survivalFunction", k)
  
  if (averageEventNumber & nbMaxEvents > 1) {
    for (k in seq(2, nbMaxEvents)) {
      datak = data[
        (data[eventName] == 1 & data$nEvent == k) |
        (data[eventName] == 0 & data$nEvent < k ),
      ]
      S <- merge(
        S,
        .renameColumns(
          .computeSurvivalCurve(datak, timeName, timeGrid, exact),
          "survivalFunction",
          k
        ), by = timeName, all = TRUE)
    }
  }

  eventNumber <- rowSums(1 - S[! names(S) == timeName])
  res <- as.data.frame(list(
    time = sort(timeGrid),
    survivalFunction = S["1"][S[[timeName]] %in% timeGrid,],
    averageEventNumber = eventNumber[S[[timeName]] %in% timeGrid]
  ))
  return(res)
}

.computeSurvivalCurve <- function(data, timeName, eventName, timeGrid, exact) {
  if (exact) {
    res <- .computeKaplanMeier(data, timeName, eventName, timeGrid)
  } else {
    res <- .computeTurnbull(data, timeName, eventName, timeGrid)
  }
  return(res)
}

.computeKaplanMeier <- function(data, timeName, eventName, timeGrid) {
  subjocc <- .getSubjocc()
  times <- sort(unique(c(timeGrid, 0, data[data[eventName] == 1, timeName])))
  survFunc <- cumprod(sapply(
    times,
    function(t) {
      nbEvents <- nrow(data[data[timeName] == t & data[eventName] == 1,])
      nbSubjects <- nrow(unique(subset(data[data$time >= t,], select = subjocc)))
      res <- 1 - ifelse(nbSubjects > 0, nbEvents / nbSubjects, 0)
      return(res)
    }
  ))
  res <- list(sort(timeGrid), survFunc[times %in% timeGrid])
  names(res) <- c(timeName, "survivalFunction")
  return(as.data.frame(res))
}

.computeTurnbull <- function(data, timeName, eventName, timeGrid) {
  bounds <- .getBounds(data, timeName, eventName)
  tau <- sort(unique(c(0, timeGrid, bounds$lb, bounds$ub, Inf)))

  subjocc <- .getSubjocc()
  subjoccs <- unique(subset(data, select = subjocc))
  
  weights <- matrix(0, nrow=nrow(subjoccs), ncol=length(tau) - 1)
  intervalProbabilities <- matrix(1 / length(tau), nrow = length(tau) - 1)
  Q <- matrix(1, nrow = length(tau) - 1)
  for (i in seq_len(nrow(subjoccs))) {
    so <- subjoccs[i, ]
    df <- subset(
      bounds,
      apply(subset(bounds, select = c(subjocc)), 1, function(x) all(x==c(so)))
    )
    weights[i, head(tau,- 1) >= df$lb & tail(tau, -1) <= df$ub] <- 1
  }
  difference <- abs(Q - intervalProbabilities)
  
  iterMax <- 200
  eps <- 1e-4;
  maxDifference <- max(difference)

  iteration <- 1

  while (maxDifference > eps & iteration < iterMax) {
    iteration <- iteration + 1
    Q <- intervalProbabilities
    C <- weights %*% intervalProbabilities
    inverseC <- 1 / C
    tempP <- t(weights) %*% inverseC
    intervalProbabilities <- 1. / nrow(subjoccs) * intervalProbabilities * tempP

    difference <- abs(Q - intervalProbabilities)
    maxDifference <- sum(difference ** 2) / sum(intervalProbabilities ** 2)
    
  }

  # filling survival curve
  survivalCurve <- 1 - c(0, cumsum(intervalProbabilities))
  res <- list(sort(timeGrid), survivalCurve[tau %in% timeGrid])
  names(res) <- c(timeName, "survivalFunction")
  return(as.data.frame(res))
}

.getBounds <- function(data, timeName, eventName) {
  subjocc <- .getSubjocc()
  bounds <- do.call(rbind, by(
    data,
    subset(data, select = subjocc),
    function(df) {
      so <- unname(unique(subset(df, select = subjocc)))
      df <- df[order(df$time),]
      eventidx <- which(df[[eventName]] == 1)
      lb <- ifelse (length(eventidx) > 0, df[eventidx - 1, timeName], df[df$censored == 1, timeName])
      ub <- ifelse(length(eventidx) > 0, df[eventidx, timeName], Inf)
      b <- list(so, lb, ub)
      names(b) <- c(subjocc, "lb", "ub")
      return(as.data.frame(b))
    }
  ))
  return(bounds)
}

.addEventIdColumn <- function(data, eventName) {
  subjocc <- intersect(c(.getSubjocc(), "rep"), names(data))
  data$nEvent <- 0
  data <- data[order(subset(data, select = time)),]
  data <- do.call(rbind, by(
    data,
    subset(data, select = subjocc),
    function(df) {
      df$nEvent <- cumsum(df[[eventName]])
      return(df)
    }
  ))
  return(data)
}

.addTTECensoredColumn <- function(data, eventName) {
  subjocc <- .getSubjocc()
  data <- data[order(subset(data, select = time)),]
  data$censored <- 0
  data <- do.call(rbind, by(
    data,
    subset(data, select = subjocc),
    function(df) {
      df$censored[df[eventName] == 0 & df$time == max(df$time)] <- 1
      return(df)
    }
  ))
  return(data)
}
