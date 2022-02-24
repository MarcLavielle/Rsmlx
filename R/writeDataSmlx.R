#' Write Simulx Dataset
#'
#' Format outputs of simulx simulations and write datasets in monolix
#' and pkanalix project format.
#' 
#' WARNING: `writeData` function is not implemented for simulx project with regressors in MonolixSuite version 2020R1
#' 
#' @param filename (\emph{string}) (\emph{optional}) file path to dataset.
#' (default "simulated_dataset.csv")
#' In case of multiple replicates, the function creates one dataset per replicate with name $filename_repi
#' If filename contains an extension, it must be "csv" or "txt". If it does not, extension is defined by `ext` argument.
#' @param sep (\emph{string}) (\emph{optional}) Separator used to write dataset file. (default ",")
#' It must be one of {"\\t", " ", ";", ","}
#' @param ext (\emph{bool}) (\emph{optional}) Extension used to write dataset file. (default "csv")
#' It must be one of {"csv", "txt"}
#' To defined only if filename with no extension
#' @param nbdigits (\emph{integer}) (\emph{optional}) number of decimal digits in output file.
#' (default = 5)
#' @param mapObservation (\emph{name vector}) (\emph{optional}) mapping of observation name
#'
#' @return a dataframe if one single simulation, a list of dataframe if multiple replicates.
writeDataSmlx <- function(filename = "simulated_dataset.csv", sep = ",",
                          ext = "csv", nbdigits = 5, mapObservation = NULL) {
  if (is.null(filename)) filename = "simulated_dataset.csv"
  if (! is.null(.getFileExt(filename))) ext <- .getFileExt(filename)

  .checkDelimiter(sep, "sep")
  .checkExtension(ext, "ext")
  .check_pos_integer(nbdigits, "nbdigits")
  
  # check input smlx project and load it
  if (!.isProjectLoaded()) {
      stop("No simulx project loaded", call. = FALSE)
  }

  # run simulation if needed
  if (is.null(smlx.getSimulationResults())) smlx.runSimulation()

  filename <- paste0(tools::file_path_sans_ext(filename), ".", ext)

  # get Simulx data
  groups <- smlx.getGroups()
  treatment <- .getTreatment()
  regressor <- .getRegressor()
  covariates <- smlx.getCovariateElements()
  simulation <- smlx.getSimulationResults()
  
  # Write data for a simulx project with Regressors not implemented
  version <- mlx.getLixoftConnectorsState()$version
  if (!is.null(regressor) & (version == "2020R1")) {
    stop("`writeData` function is not implemented for simulx project with regressors in MonolixSuite version 2020R1", call. = F)
  }
  
  nbRep <- smlx.getNbReplicates()

  # create a dataset for each replicate
  res <- list()
  for (r in seq_len(nbRep)) {
    simrep <- .filterReplicate(simulation, r)
    if (nbRep == 1) {
      treatrep <- treatment
      regrep <- regressor
    } else {
      treatrep <- NULL
      regrep <- NULL
      if (!is.null(treatment)) treatrep <- subset(treatment, rep == r)
      if (!is.null(regressor)) regrep <- subset(regressor, rep == r)
    }
    rdata <- .compileData(simrep, treatrep, regrep, covariates, groups, mapObservation)
    rdata <- rdata[names(rdata) != "rep"]

    # rename observation name
    # rdata <- .renameColumns(rdata, names(mapObservation), unname(mapObservation))
    res[[paste0("rep", r)]] <- rdata
  }
  
  # write datasets
  if (nbRep == 1) {
    res <- res$rep1
    # cut digits
    res <- .roundDataframe(res, nbdigits)
    utils::write.table(
      x = res, file = filename,
      row.names = FALSE, sep = sep,
      quote = FALSE
    )
  } else {
    for (r in seq_along(res)) {
      rfile <- paste0(
        tools::file_path_sans_ext(filename),
        "_", r, ".", .getFileExt(filename)
      )
      # cut digits
      res[[paste0("rep", r)]] <- .roundDataframe(res[[paste0("rep", r)]], nbdigits)

      # write file
      utils::write.table(
        x = res[[paste0("rep", r)]], file = rfile,
        row.names = FALSE, sep = sep,
        quote = FALSE
      )
    }
  }
  return(res)
}

.compileData <- function(simulation, treatment, regressor, covariates, groups, mapObservation) {
  # Create simulation dataset
  simData <- .convertSimulations(simulation, mapObservation)
  if (length(groups) == 1) simData$group <- groups[[1]]$name

  # Create covariates dataset
  covData <- .convertCovariates(simulation, covariates, groups)

  # create dose dataset
  treatmentData <- .convertTreatment(treatment)
  
  # create regressor dataset
  sos <- unique(subset(simData, select = .getSubjocc()))
  regressor <- .convertRegressor(regressor, sos)

  # Create compiled dataset
  if (!is.null(treatmentData)) {
    obsNames <- .getObservationNames()
    obsheader <- ifelse (length(obsNames) == 1, obsNames, "y")
    treatmentData[obsheader] <- "."
    if (is.element("evid", names(treatmentData)))
      simData$evid <- 0
    simData <- merge(simData, treatmentData, all = TRUE)
  }
  if (!is.null(regressor)) {
    simData <- merge(simData, regressor, all = TRUE)
  }
  if (!is.null(covData)) {
    simData <- merge(simData, covData)
  }
  # fill NAN by "."
  for (h in names(simData)) {
    simData[h][is.na(simData[h]),] <- "."
  }
  
  # sort dataframe
  tosort <- intersect(names(simData), c(.getSubjocc(), "time"))
  simData <- simData[do.call(order, simData[,tosort]),]
  
  if (length(unique(simData$group)) == 1) {
    group <- NULL
    simData <- subset(simData, select = - group)
  }
  return(simData)
}

.getSubjocc <- function() {
  subjocc <- c("id")
  occ <- smlx.getOccasionElements()$names
  if (length(occ) > 0) {
    subjocc <- c("id", occ)
  }
  return(subjocc)
}

.filterReplicate <- function(simulationList, r) {
  nbRep <- smlx.getNbReplicates()
  if (nbRep == 1) return(simulationList)
  res <- simulationList
  for (n in names(simulationList)) {
    for (m in names(simulationList[[n]])) {
      res[[n]][[m]] <- subset(simulationList[[n]][[m]], rep == r)
      res[[n]][[m]] <- subset(res[[n]][[m]], select = -rep)
    }
  }
  return(res)
}

.getObservationNames <- function() {
  obsnames <- names(smlx.getSimulationResults()$res)
  return(obsnames)
}

# Create dataset with simulations associated to each group
.convertSimulations <- function(simulation, mapObservation) { 
  obsnames <- .getObservationNames()
  simData <- simulation$res[[obsnames[1]]]
  if (length(obsnames) > 1) {
    simData$ytype <- mapObservation[[obsnames[1]]]
    simData <- .renameColumns(simData, obsnames[1], "y")
    for (i in seq(2, length(obsnames))) {
      dataobs <- simulation$res[[obsnames[i]]]
      dataobs$ytype <- mapObservation[[obsnames[i]]]
      dataobs <- .renameColumns(dataobs, obsnames[i], "y")
      simData <- rbind(simData, dataobs)
    }
  }
  return(simData)
}

# Create dataset with covariates associated to each group
.convertCovariates <- function(simulation, covariates, groups) { 
  covData <- NULL
  subjocc <- .getSubjocc()
  if (!length(covariates)) return(covData)
  for (g in seq_along(groups)) {
    gname <- groups[[g]]$name
    gcov <- groups[[g]]$covariate
    if (is.null(gcov)) next
    cov <- covariates[[gcov]]$data
    if (is.element("name", names(cov))) {
      gcovnames <- unique(cov$name)
    } else {
      gcovnames <- setdiff(names(cov), c("ID", subjocc))
    }
    gsim <- simulation$IndividualParameters[[gname]]
    gsim <- subset(gsim, select = c(subjocc, gcovnames))
    gsim$group <- gname
    covData <- rbind(covData, gsim)
  }
  return(covData)
}

# Create dataset with treatments associated to each group
.convertTreatment <- function(treatment) {
  if (is.null(treatment)) return(treatment)
  if (!all(treatment$washout == FALSE))
    treatment$evid <- sapply(treatment$washout, function(x) ifelse(x, 3, 1))
  # remove unused columns
  washout <- tinf <- admtype <- NULL
  treatment <- subset(treatment, select = - washout)
  if (all(is.na(treatment$tinf)))
    treatment <- subset(treatment, select = - tinf)
  if (all(treatment$admtype == 1))
    treatment <- subset(treatment, select = - admtype)
  return(treatment)
}

# Temporary treatment of regressor - when only 1 id in regressor file, spread value of the id to others
.convertRegressor <- function(regressor, sos) {
  if (is.null(regressor)) return(regressor)
  if (all(regressor$id == 1)) {
    ntimes <- nrow(regressor)
    regressor <- do.call("rbind", replicate(nrow(sos), regressor, simplify = FALSE))
    regressor[names(sos)] <- as.data.frame(lapply(sos, rep, each = ntimes))
  }
  return(regressor)
}

.getTreatment <- function() {
  filename <- file.path(smlx.getProjectSettings()$directory, "Simulation", "doses.txt")
  if (file.exists(filename)) {
    treatment <- utils::read.csv(filename)
    if (nrow(treatment) == 0) treatment <- NULL
  } else {
    treatment <- NULL
  }
  return(treatment)
}

.getRegressor <- function() {
  filename <- file.path(smlx.getProjectSettings()$directory, "Simulation", "regressors.txt")
  if (file.exists(filename)) {
    regressor <- utils::read.csv(filename)
    if (nrow(regressor) == 0) regressor <- NULL
  } else {
    regressor <- NULL
  }
  return(regressor)
}

################################################################################
# Checks
################################################################################
# Checks if a delimiter is valid (valid delimiters "\t", " ", ";", ",")
.checkDelimiter <- function(delimiter, argname) {
  if (is.null(delimiter)) return(invisible(TRUE))
  if (! is.element(delimiter, c("\t", " ", ";", ",")))
    stop("`", argname, "` must be one of: {\"",
         paste(c("\\t", " ", ";", ","), collapse="\", \""), "\"}.",
         call. = FALSE)
  return(invisible(TRUE))
}

# Check if an extension is valid (Valid extensions are "csv", "txt")
.checkExtension <- function(ext, argname) {
  if (is.null(ext)) return(invisible(TRUE))
  if (! is.element(ext, c("csv", "txt")))
    stop("`", argname, "` must be one of: {\"",
         paste(c("csv", "txt"), collapse="\", \""), "\"}.",
         call. = FALSE)
  return(invisible(TRUE))
}

# Checks if an input is a strictly positive integer and returns an error -------
.check_pos_integer <- function(int, argname) {
  if(!(is.double(int)||is.integer(int)))
    stop("Invalid ", argname, ". It must be a positive integer.", call. = F)
  if ((int < 0) || (!as.integer(int) == int))
    stop("Invalid ", argname, ". It must be a positive integer.", call. = F)
  return(int)
}

# Check is a project is loaded -------------------------------------------------
.isProjectLoaded <- function() {
  # disable monolix errors
  op <- getOption("lixoft_notificationOptions")
  op$errors <- 1
  options(lixoft_notificationOptions = op)
  if (is.null(smlx.getProjectSettings())) {
    isloaded <- FALSE
  } else {
    isloaded <- TRUE
  }
  # reactivate monolix errors
  op <- getOption("lixoft_notificationOptions")
  op$errors <- 0
  options(lixoft_notificationOptions = op)
  return(isloaded)
}

