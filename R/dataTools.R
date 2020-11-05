#' Read pkanalix input dataset
#' @return Input dataset in dataframe format. Headers are renamed to match mlx names
#' @examples
#' \dontrun{
#' readmlxDataset()
#' }
readmlxDataset <-function(obsid = NULL) {
  filename <- mlx.getData()$dataFile
  df <- .readDataset(filename)
  names(df) <- mlx.getData()$header
  if (.isHeader("evid")) {
    message("WARNING: EVID dataset ! ")
  }
  renameHeaders <- .filterHeaders(c("id", "time", "observation", "amount", "evid",
                                    "obsid", "cens", "limit", "regressor", "admid",
                                    "rate", "tinf", "ss", "addl", "mdv", "ii"))
  df <- .renameColumns(df, .matchToData(renameHeaders), renameHeaders)
  return(df)
}

#' Get censoring information from input dataset
#' @return Censoring columns of input dataset (dataframe).
#' @examples
#' \dontrun{
#' getCensoringInformation()
#' }
getCensoringInformation <- function() {
  subjocc <- .getSubjocc()
  censoring <- getPKObservations()
  censoring <- subset(censoring, select = .filterHeaders(c(subjocc, "time", "cens", "limit")))
  return(censoring)
}

#' Get dose rows from input dataset
#' @return Dose rows of input dataset (dataframe)
#' @examples
#' \dontrun{
#' getDoseInformation()
#' }
getDoseInformation <- function() {
  df <- readmlxDataset()
  if (!is.element("evid", df)) {
    subDf <- subset(df, df$amount != ".")
  } else {
    subDf <- subset(df, (df$amount != ".") & (df$evid %in% c(1, 3, 4)))
  }
  return(subDf)
}

#' Get observations rows from input dataset
#' @return PK Observations rows of input dataset (dataframe)
#' @examples
#' \dontrun{
#' getPKObservations()
#' }
getPKObservations <- function() {
  df <- readmlxDataset()
  if (!is.element("evid", df)) {
    subDf <- subset(df, df$observation != ".")
  } else {
    subDf <- subset(df, (df$observation != ".") & (df$evid %in% c(0)))
  }
  if (is.element("mdv", df)){
    subDf <- subset(df, mdv == 1)
  }
  if (is.element("obsid", df)) {
    subDf <- subset(subDf, df$obsid == mlx.getGlobalObsIdToUse())
  }
  return(subDf)
}

#' Get dataset to get only observations after the last dose
#' @param df Input dataframe.
#' @return observations after last dose (dataframe)
#' @examples
#' \dontrun{
#' cutData()
#' }
cutData <- function(df) {
  starttime =  unique(as.numeric(df$startTime))
  df <- subset(df, df$time >= starttime)
  df["time"] <- sapply(df$time, FUN=function(time) {time - starttime})
  df <-  df[,!(names(df) %in% c("startTime"))]
  return(df)
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

.getStartTime <- function(df) {
  subjocc <- .getSubjocc()
  dfStart <- subset(getDoseInformation(), select = c(subjocc, "time"))
  dfStart <- aggregate(
    dfStart$time,
    by = as.list(subset(dfStart, select = subjocc)),
    FUN=max
  )
  dfStart <- .renameColumns(dfStart, "x", "startTime")
  return(dfStart)
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

.getSubjocc <- function(){
  if (.isHeader("occ")) {
    subjocc <- c("id", .matchToData("occ"))
  } else {
    subjocc <- "id"
  }
  return(subjocc)
}

.isHeader <- function(headers) {
  return(invisible(headers %in% mlx.getData()$headerTypes))
}
