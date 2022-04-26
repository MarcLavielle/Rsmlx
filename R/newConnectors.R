#' Get estimated individual and population parameters
#' 
#' Get the individual individual parameters, the population parameters with the population covariates 
#' and the population parameters with the individual covariates.
#' @return a list of data frames.
#' @examples
#' \dontrun{
#' # Assume that the Monolix project "warfarinPKPD.mlxtran" has been loaded
#' r = getEstimatedIndividualParameters2() 
#' 
#' # r is a list with elements "saem", "conditionalMean", "conditionalSD",   "conditionalMode",
#' # "popPopCov" and "popIndCov"
#' 
#' # See http://rsmlx.webpopix.org/userguide/newconnectors/ for more detailed examples
#' # Download the demo examples here: http://rsmlx.webpopix.org/installation

#' }
#' @export
getEstimatedIndividualParameters2 <- function() {
  
  if (!initRsmlx()$status)
    return()
  
  ind.param <- mlx.getEstimatedIndividualParameters()
  N <- nrow(ind.param$saem)
  pop.param <- mlx.getEstimatedPopulationParameters()
  pop.param <- pop.param[grep("_pop",names(pop.param))]
  df <- as.data.frame(matrix(pop.param,nrow=N,ncol=length(pop.param),byrow=TRUE))
  names(df) <- gsub("_pop","",names(pop.param))
  ind.param$popPopCov <- data.frame(id=ind.param$saem["id"],df)
  
  rand.eff <- mlx.getEstimatedRandomEffects()
  
  if (!is.null(ind.param$conditionalMean)) {
    ip <- ind.param$conditionalMean
    re <- rand.eff$conditionalMean
  } else if (!is.null(ind.param$conditionalMode)) {
    ip <- ind.param$conditionalMode
    re <- rand.eff$conditionalMode
  } else
    stop("The conditional mean or the conditional model should have been computed", call.=FALSE)
  
  ind.dist <- mlx.getIndividualParameterModel()$distribution
  var.param <- names(ind.dist)
  ind.param$popIndCov <- ind.param$popPopCov
  for (nj in var.param) {
    dj <- tolower(ind.dist[nj])
    yj <- ip[[nj]]
    rj <- re[[paste0("eta_",nj)]]
    if (dj == "normal") {
      yjc <- yj-rj
    } else if (dj == "lognormal") {
      yjc <- exp(log(yj)-rj)
    } else if (dj == "logitnormal") {
      yjc <- 1/(1+exp(-log(yj/(1-yj))+rj))
    } else if (dj == "probitnormal") {
      yjc <- pnorm(qnorm(yj) - rj)
    } 
    ind.param$popIndCov[nj] <- yjc
  }
  return(ind.param)
}


#--------------------------------------------------------------------

#' Get estimated predictions
#' 
#' Get the individual predictions obtained with the estimated individual parameters :
#' @return a list of data frames (one data frame per output).
#' @examples
#' \dontrun{
#' # Assume that the Monolix project "warfarinPKPD.mlxtran" has been loaded
#' r = getEstimatedPredictions() # r is a list with elements "y1" and "y2"
#' 
#' # See http://rsmlx.webpopix.org/userguide/newconnectors/ for more detailed examples
#' # Download the demo examples here: http://rsmlx.webpopix.org/installation

#' }
#' @export
getEstimatedPredictions <- function() {
  
  if (!initRsmlx()$status)
    return()
  
  ip <- getEstimatedIndividualParameters2()
  
  obs.info <- mlx.getObservationInformation()
  
  df <- list()
  obsNames <- names(obs.info$mapping)
  contNames <- obsNames[obs.info$type[obsNames]=="continuous"]
  for (j in seq_along(contNames)) {
    name <- contNames[j]
    df[j] <- obs.info[name]
    df[[j]][name] <- NULL
  }
  
  f.pop1 <- mlx.computePredictions(ip$popPopCov)
  for (j in seq_along(contNames)) {df[[j]]$popPopCov <- f.pop1[[j]]}
  f.pop2 <- mlx.computePredictions(ip$popIndCov)
  for (j in seq_along(contNames)) {df[[j]]$popIndCov <- f.pop2[[j]]}
  
  if (!is.null(ip$conditionalMean)) {
    f.mean <- mlx.computePredictions(ip$conditionalMean)
    for (j in seq_along(contNames)) {df[[j]]$conditionalMean <- f.mean[[j]]}
  }
  if (!is.null(ip$conditionalMode)) {
    f.mode <- mlx.computePredictions(ip$conditionalMode)
    for (j in seq_along(contNames)) {df[[j]]$conditionalMode <- f.mode[[j]]}
  }
  names(df) <- names(f.pop1)
  return(df)
}

#--------------------------------------------------------------------

#' Get estimated residuals
#' 
#' Get the residuals computed from the individual predictions obtained 
#' with the estimated individual parameters:
#' @return a list of data frames (one data frame per output).
#' @examples
#' \dontrun{
#' # Assume that the Monolix project "warfarinPKPD.mlxtran" has been loaded
#' r = getEstimatedResiduals()  # r is a list with elements "y1" and "y2" 
#' 
#' # See http://rsmlx.webpopix.org/userguide/newconnectors/ for more detailed examples
#' # Download the demo examples here: http://rsmlx.webpopix.org/installation

#' }
#' @export
getEstimatedResiduals <- function() {
  
  if (!initRsmlx()$status)
    return()
  
  df <- getEstimatedPredictions()
  obs.info <- mlx.getObservationInformation()
  nip <- c("popPopCov", "popIndCov", "conditionalMean", "conditionalMode")
  
  obsNames <- names(obs.info$mapping)
  contNames <- obsNames[obs.info$type[obsNames]=="continuous"]
  nobs <- length(contNames)
  error.model <- mlx.getContinuousObservationModel()$errorModel
  error.dist <- mlx.getContinuousObservationModel()$distribution
  ep <- error.parameter()
  pop.param <- mlx.getEstimatedPopulationParameters()
  param.error <- list()
  for (j in seq_along(contNames)) {
    name <- contNames[j]
    dfj <- df[[j]]
    ij <- which(names(dfj) %in% nip)
    yoj <- replicate(length(ij), obs.info[[name]][[name]])
    erj <- tolower(error.dist[j])
    if (erj=="normal") {
      ypj <- dfj[,ij]
    } else if (erj=="lognormal") {
      ypj <- log(dfj[,ij])
      yoj <- log(yoj)
    } else if (erj=="logitnormal") {
      limiti <- mlx.getContinuousObservationModel()$limits[[names(error.dist)[j]]]
      ypj <- log((dfj[,ij]-limiti[1])/(limiti[2]-dfj[,ij]))
      yoj <- log((yoj-limiti[1])/(limiti[2]-yoj))
    }
    epj <- pop.param[ep[[j]]]
    a <- epj[grep("a", names(epj))]
    if (length(a)==0)  
      a <- 0
    b <- epj[grep("b", names(epj))]
    if (length(b)==0)  
      b <- 0
    c <- epj[grep("c", names(epj))]
    if (length(c)==0)  
      c <- 1
    pei <- c(a, b, c)
    if (error.model[j]=="combined2")
      dfj <- (yoj - ypj)/sqrt(pei[1]^2 + (pei[2]*ypj^pei[3])^2)
    else
      dfj <- (yoj - ypj)/(pei[1] + pei[2]*ypj^pei[3])
    df[[j]][names(dfj)] <- dfj
  }
  names(df) <- contNames
  return(df)
}

#--------------------------------------------------------------------
#' Get simulated predictions
#' 
#' Get the individual predictions obtained with the simulated individual parameters :
#' @return a list of data frames (one data frame per output).
#' @examples
#' \dontrun{
#' # Assume that the Monolix project "warfarinPKPD.mlxtran" has been loaded
#' r = getSimulatedPredictions()  # r is a list with elements "Cc" and "E" 
#' 
#' # See http://rsmlx.webpopix.org/userguide/newconnectors/ for more detailed examples
#' # Download the demo examples here: http://rsmlx.webpopix.org/installation

#' }
#' @export
getSimulatedPredictions <- function() {
  
  if (!initRsmlx()$status)
    return()
  
  sip <- mlx.getSimulatedIndividualParameters()
  if (is.null(sip$rep)) 
    sip$rep <- 1
  nrep <- max(sip$rep)
  
  obs.info <- mlx.getObservationInformation()
  df <- list()
  obsNames <- names(obs.info$mapping)
  contNames <- obsNames[obs.info$type[obsNames]=="continuous"]
  pred <- mlx.getContinuousObservationModel()$prediction
  for (j in seq_along(contNames)) {
    name <- contNames[j]
    df[j] <- obs.info[name]
    df[[j]][name] <- NULL
    df[[j]][pred[name]] <- 0
    df[[j]] <- cbind(rep=1, df[[j]])
  }
  
  col.el <- which(!(names(sip) %in% c("rep","id")))
  res <- list()
  for (irep in seq_len(nrep)) {
    parami <- subset(sip, rep==irep)[,col.el]
    fi <- mlx.computePredictions(parami)
    for (j in seq_along(contNames)) {
      name <- contNames[j]
      df[[j]][pred[name]] <- fi[[pred[name]]]
      df[[j]]["rep"] <- irep
      if (irep==1)
        res[[j]] <- df[[j]]
      else
        res[[j]] <- rbind(res[[j]], df[[j]])
    }
    
  }
  names(res) <- pred[contNames]
  return(res)
}

#--------------------------------------------------------------------
#' Get simulated residuals
#' 
#' Get the residuals computed from the individual predictions obtained 
#' with the simulated individual parameters:
#' @return a list of data frames (one data frame per output).
#' @examples
#' \dontrun{
#' # Assume that the Monolix project "warfarinPKPD.mlxtran" has been loaded
#' r = getSimulatedResiduals()  # r is a list with elements "y1" and "y2" 
#' 
#' # See http://rsmlx.webpopix.org/userguide/newconnectors/ for more detailed examples
#' # Download the demo examples here: http://rsmlx.webpopix.org/installation

#' }
#' @export
getSimulatedResiduals <- function() {
  
  if (!initRsmlx()$status)
    return()
  
  df <- getSimulatedPredictions()
  obs.info <- mlx.getObservationInformation()
  obsNames <- names(obs.info$mapping)
  contNames <- obsNames[obs.info$type[obsNames]=="continuous"]
  
  error.model <- mlx.getContinuousObservationModel()$errorModel
  error.dist <- mlx.getContinuousObservationModel()$distribution
  ep <- error.parameter()
  pop.param <- mlx.getEstimatedPopulationParameters()
  param.error <- list()
  nrep <- max(df[[1]]["rep"])
  for (j in seq_along(contNames)) {
    name <- contNames[j]
    dfj <- df[[j]]
    ij <- which(names(dfj) == names(df)[j])
    yoj <- rep(obs.info[[name]][[name]],nrep)
    erj <- tolower(error.dist[j])
    if (erj=="normal") {
      ypj <- dfj[,ij]
    } else if (erj=="lognormal") {
      ypj <- log(dfj[,ij])
      yoj <- log(yoj)
    } else if (erj=="logitnormal") {
      limiti <- mlx.getContinuousObservationModel()$limits[[names(error.dist)[j]]]
      ypj <- log((dfj[,ij]-limiti[1])/(limiti[2]-dfj[,ij]))
      yoj <- log((yoj-limiti[1])/(limiti[2]-yoj))
    }
    epj <- pop.param[ep[[j]]]
    a <- epj[grep("a", names(epj))]
    if (length(a)==0)  
      a <- 0
    b <- epj[grep("b", names(epj))]
    if (length(b)==0)  
      b <- 0
    c <- epj[grep("c", names(epj))]
    if (length(c)==0)  
      c <- 1
    pei <- c(a, b, c)
    if (error.model[j]=="combined2")
      dfj <- (yoj - ypj)/sqrt(pei[1]^2 + (pei[2]*ypj^pei[3])^2)
    else
      dfj <- (yoj - ypj)/(pei[1] + pei[2]*ypj^pei[3])
    df[[j]][,ij] <- dfj
    names(df[[j]])[ij] <- "residual"
  }
  names(df) <- contNames
  return(df)
}

#--------------------------------------------------------------------
#' Get estimated covariance and correlation matrices
#' 
#' @return a list of two matrices.
#' @examples
#' \dontrun{
#' # Assume that the Monolix project "warfarinPKPD.mlxtran" has been loaded
#' r = getEstimatedCovarianceMatrix()  # r is a list with elements "cor.matrix" and "cov.matrix"
#' 
#' # See http://rsmlx.webpopix.org/userguide/newconnectors/ for more detailed examples
#' # Download the demo examples here: http://rsmlx.webpopix.org/installation

#' }
#' @export
getEstimatedCovarianceMatrix <- function() {
  
  if (!initRsmlx()$status)
    return()
  
  param <- mlx.getEstimatedPopulationParameters()
  pname <- names(param)
  i.omega <- grep("^omega_",pname)
  if (length(i.omega)>0) {
    oest <- 1
    oname <- gsub("^omega_","",pname[i.omega])
    omega <- param[i.omega]
  } else {
    i.omega <- grep("^omega2_",pname)
    oname <- gsub("^omega2_","",pname[i.omega])
    omega <- sqrt(param[i.omega])
    oest <- 2
  }
  i.corr <- grep("^corr_",pname)
  d <- length(i.omega)
  c <- param[i.corr]
  R <- diag(rep(1,d))
  rownames(R) <- colnames(R) <- oname
  if (length(c)>0) {
    for (j in 1:length(c)) {
      cj <- names(c)[j]
      sj <- strsplit(cj,"_")[[1]]
      R[sj[3],sj[2]] <- R[sj[2],sj[3]] <- c[j]
    }
  }
  C <- diag(omega)%*%R%*%diag(omega)
  return(list(cor.matrix=R, cov.matrix=C))
}

#--------------------------------------------------------------------
error.parameter <- function(project=NULL) {
  if (is.null(project)) {
    dp <- mlx.getProjectSettings()$directory
    if (!is.null(dp))
      project <- paste0(dp,".mlxtran")
    #    project <- paste0(basename(dp),".mlxtran")
  }
  if (!file.exists(project)) 
    stop("Enter a valid project", call.=FALSE)
  
  con        = file(project, open = "r")
  lines      = readLines(con, warn=FALSE)
  close(con)
  lines <- lines[grep("errorModel", lines)]
  r <- list()
  if (length(lines)>0) {
    lines <- gsub(" ","",lines)
    for (k in (1: length(lines))) {
      lk <- lines[k]
      i1 <- regexpr("\\(",lk)
      i2 <- regexpr("\\)",lk)
      r[[k]] <- strsplit(substr(lk,(i1+1),(i2-1)),",")[[1]]
    }
  }
  return(r)
}

#-------------------------------------------------
mlx.getFit <- function() {
  project <- mlx.getProjectSettings()$project
  con     <- file(project, open = "r")
  lines   <- readLines(con, warn=FALSE)
  close(con)
  lines <- gsub("\\;.*","",lines)
  l.fit <- grep("<FIT>", lines)
  lines <- lines[l.fit:length(lines)]
  l.data <- lines[grep("data", lines)[1]]
  l.model <- lines[grep("model", lines)[1]]
  
  ll <- gsub(" ","",l.data)
  ll <- gsub("data=","",ll)
  if (grepl("\\{",ll)) {
    i1 <- regexpr("\\{",ll)
    i2 <- regexpr("\\}",ll)
    ll <- substr(ll,i1+1,i2-1)
    data.names  <- strsplit(ll,",")[[1]]
  } else {
    data.names <- ll
  }
  ll <- gsub(" ","",l.model)
  ll <- gsub("model=","",ll)
  if (grepl("\\{",ll)) {
    i1 <- regexpr("\\{",ll)
    i2 <- regexpr("\\}",ll)
    ll <- substr(ll,i1+1,i2-1)
    model.names  <- strsplit(ll,",")[[1]]
  } else {
    model.names <- ll
  }
  
  return(list(data=data.names, model=model.names))
}
  
