#' Read pkanalix input dataset
#' @return Input dataset in dataframe format. Headers are renamed to match mlx names
#' @examples
#' \dontrun{
#' readmlxDataset()
#' }
readmlxDataset <-function() {
  filename <- mlx.getData()$dataFile
  df <- .readDataset(filename)
  names(df) <- mlx.getData()$header

  renameHeaders <- .filterHeaders(c("id", "time", "observation", "amount", "evid",
                                    "obsid", "cens", "limit", "regressor", "admid",
                                    "rate", "tinf", "ss", "addl", "mdv", "ii"))
  df <- .renameColumns(df, .matchToData(renameHeaders), renameHeaders)
  return(df)
}

#' Get dose rows from input dataset
#' @return Dose rows of input dataset (dataframe)
#' @examples
#' \dontrun{
#' getDoseInformation()
#' }
getDoseInformation <- function() {
  df <- readmlxDataset()
  if (!is.element("evid", names(df))) {
    subDf <- subset(df, df$amount != ".")
  } else {
    subDf <- subset(df, (df$amount != ".") & (df$evid %in% c(1, 3, 4)))
  }
  return(subDf)
}

.getDelim <- function(filename, seps = c("\t"," ",";",","), nrows = 10) {
  if (!file.exists(filename)){
    message(paste0("ERROR: File ", filename, " does not exist !"))
    return(invisible(NULL))
  }
  nbColumns <- 1
  i <- 1
  while ((nbColumns == 1) & (i <= length(seps))) {
    sep <- seps[i]
    df <- read.csv(filename, sep = sep, nrows = nrows)
    nbColumns <- length(df)
    i <- i + 1
  }
  if (nbColumns == 1) {
    message(paste0("ERROR: Cannot read file ", filename, "."))
    return(invisible(NULL))
  }
  return(invisible(sep))
}

.readDataset <- function(filename) {
  if (!file.exists(filename)){
    message(paste0("ERROR: File ", filename, " does not exist !"))
    return(invisible(NULL))
  }
  sep <- .getDelim(filename)
  if (is.null(sep)){
    message(paste0("ERROR: Cannot read file ", filename, "."))
    return(invisible(NULL))
  }
  df <- read.csv(filename, sep = sep)
  return(df)
}

.filterHeaders <- function(headersList) {
  outputHeadersList <- c()
  for (h in headersList) {
    if (h %in% mlx.getData()$headerTypes) {
      outputHeadersList <- c(outputHeadersList, h)
    } else if (is.element(h, mlx.getData()$header)) {
      outputHeadersList <- c(outputHeadersList, h)
    }
  }
  return(outputHeadersList)
}

.renameColumns <- function(df, oldName, newName){
  if (length(oldName) != length(newName)) {
    message("ERROR: vector of old names and new names must match in size")
    return(invisible(df))
  }
  for (i in seq(1, length(oldName))) {
    old <- oldName[i]
    new <- newName[i]
    if (old %in% names(df)) {
      names(df)[names(df) == old] <- new
    }
  }
  return(invisible(df))
}

.matchToData <- function(headermlx){
  headersmlx <- mlx.getData()$headerTypes
  headersDataset <- mlx.getData()$header
  
  headerDataset <- c()
  for (h in headermlx) {
    if (is.element(h, headersmlx)) {
      idxHeader <- which(headersmlx %in% h)
      headerDataset <- c(headerDataset, headersDataset[idxHeader])
    }
  }
  return(headerDataset)
}

.isHeader <- function(headers) {
  return(invisible(headers %in% mlx.getData()$headerTypes))
}

.getObsType <- function(obsName) {
  obsInformation <- mlx.getObservationInformation()
  obsType <- obsInformation$type[[obsInformation$mapping[[obsName]]]]
  
  subType = obsType
  modelFile <- mlx.getStructuralModel()
  if (grepl( "^lib:", modelFile)) {
    lixoftConnectorsState <- mlx.getLixoftConnectorsState(quietly = TRUE)
    libPath <- paste(lixoftConnectorsState$path, "factory", "library", sep = "/")
    modelFile <- paste(
      libPath,
      list.files(
        path = libPath,
        pattern = paste0("*", gsub("^lib:", "", modelFile)),
        recursive = TRUE
      ),
      sep = "/"
    )
  }
  lines  <- paste(readLines(modelFile, warn = FALSE), collapse = " ")
  if (obsType == "discrete") {
    m <- regmatches(lines, regexpr(paste0(obsName, "\\s*=\\s*\\{([a-zA-Z0-9]|[\\(\\)\\+\\-\\*\\./\"_.=,;]|\\s)*}"), lines, perl=TRUE))
    m2 <- regmatches(lines, regexpr("type\\s*=\\s*[a-zA-Z]*", lines, perl=TRUE))
    subType <- gsub("type\\s*=\\s*", "", m2)
  } else if (obsType == "event") {
    subType <- "exactEvent"
    m <- regmatches(lines, regexpr(paste0(obsName, "\\s*=\\s*\\{([a-zA-Z0-9]|[\\(\\)\\+\\-\\*\\./\"_.=,;]|\\s)*}"), lines, perl=TRUE))
    m2 <- regmatches(m, regexpr("eventType\\s*=\\s*[a-zA-Z]*", m, perl=TRUE))
    if (length(m2) > 0) {
      subType <- gsub("eventType\\s*=\\s*", "", m2)
    }
  }
  return(list(type = obsType, subType = subType))
}

.getSubjocc <- function(){
  if (.isHeader("occ")) {
    subjocc <- c("id", .matchToData("occ"))
  } else {
    subjocc <- "id"
  }
  return(subjocc)
}

.stratify <- function(stratifierData, type, covName, breaks = NULL, filtering = FALSE, filterGroup = NULL) {
  if (type == "continuous") {
    rangeCov <- range(stratifierData)
    if (is.null(breaks)) breaks <- median(stratifierData)
    breaks <- sort(unique(c(rangeCov, pmin(pmax(rangeCov[1], breaks), rangeCov[2]))))
  
    if (!filtering) {
      factor <- cut(stratifierData, breaks = breaks, include.lowest=TRUE)
      levels(factor) <- paste0("#", covName, ": ", levels(factor), " ")
      res <- factor
    } else {
      factor <- cut(stratifierData, breaks = breaks, include.lowest=TRUE, labels = F)
      res <- !(factor == filterGroup)
    }
  } else {
    if (!filtering) {
      res <- paste0("#", covName, ": ", stratifierData, " ")
    } else {
      res <- !(stratifierData == filterGroup)
    }
  }
  return(res)
}

.addStratificationColumn <- function(covData, covtypes, stratifiers, breaks,
                                     filtering = FALSE, filterGroups = NULL,
                                     columnName = "split") {
  if (!filtering) {
    stratCol <- rep("", nrow(covData))
  } else {
    stratCol <- rep(FALSE, nrow(covData))
  }
  for (cov in stratifiers) {
    strat <- .stratify(
      covData[[cov]], covtypes[[cov]], breaks[[cov]], covName = cov,
      filtering = filtering, filterGroup = filterGroups[[cov]]
    )
    if (!filtering) {
      stratCol <- paste0(stratCol, strat)
    } else {
      stratCol <- stratCol | strat
    }
  }
  covData[columnName] <- stratCol
  return(covData)
}

fill.na <- function(df, columns) {
  for(col in columns) {
    naidx <- which(is.na(df[col]))
    valuesidx <- which(!is.na(df[col]))
    closestidx <- sapply(naidx, function(x) tail(valuesidx[valuesidx < x], n = 1))
    df[col][naidx,] <- df[col][closestidx,]
  }
  return(df)
}
