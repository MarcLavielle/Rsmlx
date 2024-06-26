prcheck <- function(project, f=NULL, settings=NULL, model=NULL, paramToUse=NULL,
                    parameters=NULL, level=NULL, tests=NULL, nboot=NULL, method=NULL) {
  #prcheck <- function(project) {
  RsmlxDemo1.project <- RsmlxDemo2.project <- warfarin.data  <- resMonolix <- NULL
  if (identical(substr(project,1,9),"RsmlxDemo")) {
    RsmlxDemo1.project <- RsmlxDemo2.project <- warfarin.data  <- resMonolix <- NULL
    rm(RsmlxDemo1.project, RsmlxDemo2.project, warfarin.data, resMonolix)
    eval(parse(text="data(RsmlxDemo)"))
    tmp.dir <- tempdir()
    write(RsmlxDemo1.project, file=file.path(tmp.dir,"RsmlxDemo1.mlxtran"))
    write(RsmlxDemo2.project, file=file.path(tmp.dir,"RsmlxDemo2.mlxtran"))
    write.csv(warfarin.data, file=file.path(tmp.dir,"warfarin_data.csv"), quote=FALSE, row.names = FALSE)
    project <- file.path(tmp.dir,project)
    demo <- TRUE
    if (!is.null(f)) {
      if (f=="boot") {
        if (is.null(settings))
          res <- resMonolix$r1.boot
        else if (!is.null(settings$N) & is.null(settings$covStrat))
          res <- resMonolix$r2.boot
        else
          res <- resMonolix$r3.boot
      } else if (f=="build") {
        if (identical(model,"all") & identical(paramToUse,"all")) 
          res <- resMonolix$r1.build
        else if (identical(model,"all")) 
          res <- resMonolix$r2.build
        else 
          res <- resMonolix$r3.build
      } else if (f=="conf") {
        if (method == "fim" & level==0.90)
          res <- resMonolix$r1.conf
        else if (method == "fim" & level==0.95)
          res <- resMonolix$r2.conf
        else if (method == "proflike")
          res <- resMonolix$r3.conf
        else
          res <- resMonolix$r4.conf
      } else if (f=="cov") {
        if (identical(method,"COSSAC") & identical(paramToUse,"all")) 
          res <- resMonolix$r1.cov
        else if (identical(method,"SCM")) 
          res <- resMonolix$r2.cov
        else 
          res <- resMonolix$r3.cov
      } else if (f=="test") {
        if (length(tests)==4) 
          res <- resMonolix$r1.test
        else 
          res <- resMonolix$r2.test
      } else if (f=="set")
        res="foo"
    }
    
  } else {
    
    if (!grepl("\\.",project))
      project <- paste0(project,".mlxtran")
    
    if(!file.exists(project))
      stop(paste0("Project '", project, "' does not exist"), call.=FALSE)
    
    lp <- mlx.loadProject(project) 
    if (!lp) 
      stop(paste0("Could not load project '", project, "'"), call.=FALSE)
    
    demo <- FALSE
    res <- NULL
  }
  
  return(list(project=project, demo=demo, res=res))
  #  return(project)
}







#' Initialize Rsmlx library
#' 
#' Initialize Rsmlx library and lixoftConnectors. Prints information about the versions of Monolix and lixoftConnectors used.
#' @param path Monolix path 
#' @return A list:
#' \itemize{
#'   \item \code{software}: the software that is used (should be monolix with Rsmlx)
#'   \item \code{path}: the path to MonolixSuite
#'   \item \code{version}: the version of MonolixSuite that is used
#'   \item \code{status}: boolean equaling TRUE if the initialization has been successful.
#' }
#' @examples
#' \dontrun{
#' initRsmlx() 
#' initRsmlx(path="C:/ProgramData/Lixoft/MonolixSuite2024R1")  # specifiy a specific path
#' }
#' @export
initRsmlx <- function(path=NULL){
  if (system.file(package = "lixoftConnectors") == "")
    stop("You need to install the lixoftConnectors package in order to use Rsmlx", call. = FALSE)
  
  ver_Connectors <- packageVersion("lixoftConnectors")
  ver_Rsmlx <- packageVersion("Rsmlx")
  
  if (ver_Rsmlx$major != ver_Connectors$major) {
    stop(paste0("The major version number for the lixoftConnectors package and the Rsmlx package must be the same:\nlixoftConnectors package version -> ",
                ver_Connectors, "\nRsmlx package version -> ", ver_Rsmlx), call. = FALSE)
  }
  
  lixoftConnectorsState <- mlx.getLixoftConnectorsState(quietly = TRUE)
  
  if (!is.null(lixoftConnectorsState)){
    
    if (lixoftConnectorsState$software == "monolix"  && is.null(path)) {
      status=TRUE
    } else {
      status = mlx.initializeLixoftConnectors(path=path)
    }
    
  } else {
    status = mlx.initializeLixoftConnectors(path=path)
  }
  
  lixoftConnectorsState <- mlx.getLixoftConnectorsState(quietly = TRUE)
  lixoftConnectorsState$status <- status
  ver_lixoft <- lixoftConnectorsState$version
  ver_lixoft_major <- as.numeric(sub("([0-9]+).*$", "\\1", ver_lixoft))
  
  if (ver_Rsmlx$major != ver_lixoft_major) {
    stop(paste0("The major version number for MonolixSuite and the Rsmlx package must be the same:\nMonolixSuite version -> ",
                ver_lixoft, "\nRsmlx package version -> ", ver_Rsmlx), call. = FALSE)
  }
  
  if (!status)
    stop("Failed to initialize lixoftConnectors package.", call. = FALSE)
  
  if (is.null(path))
    return(lixoftConnectorsState)
  else
    return(invisible(lixoftConnectorsState))
}


#-------------------------------------------------
mlx.generatePKmodel <-  function(parameter, model="pk_model.txt", output=NULL) {
  str1 <- paste(parameter,collapse=",")
  str2 <- str1
  model_txt="
[LONGITUDINAL]
input = {param.input}
EQUATION:"
  
  if ("Q" %in% parameter) {
    model_txt <- paste0(model_txt,"\nV = V1\nk12 = Q/V1\nk21 = Q/V2")
    str2 <- gsub("V1","V", str2)
    str2 <- gsub("Q","k12", str2)
    str2 <- gsub("V2","k21", str2)
  }
  if ("Q2" %in% parameter) {
    model_txt <- paste0(model_txt,"\nV = V1\nk12 = Q2/V1\nk21 = Q2/V2\nk13 = Q3/V1\nk31 = Q3/V3")
    str2 <- gsub("V1","V", str2)
    str2 <- gsub("Q2","k12", str2)
    str2 <- gsub("V2","k21", str2)
    str2 <- gsub("Q3","k13", str2)
    str2 <- gsub("V3","k31", str2)
  }
  if ("Cl" %in% parameter) {
    model_txt <- paste0(model_txt,"\nk = Cl/V")
    str2 <- gsub("Cl","k", str2)
  }
  model_txt <- paste0(model_txt,"\nCc = pkmodel(param.pkmodel)
OUTPUT:
output = Cc
")
  model_txt = gsub("param.input",str1,model_txt)
  model_txt = gsub("param.pkmodel",str2,model_txt)
  if(!is.null(output))
    model_txt = gsub("Cc",output,model_txt)
  write(model_txt, model)
  #  return(invisible())
}


#-------------------------------------------------
setPKproject <- function(parameter, project="pk_project.mlxtran", model="pk_model.txt") {
  if (is.character(parameter)) 
    param.list <- parameter
  else
    param.list <- names(parameter)
  mlx.generatePKmodel(parameter=param.list, model=model)
  mlx.setStructuralModel(modelFile=model)
  if (is.numeric(parameter)) {
    pop.ini <- mlx.getPopulationParameterInformation()
    j.pop <- which(pop.ini$name %in% paste0(param.list,"_pop"))
    pop.ini$initialValue[j.pop] <- parameter
    mlx.setPopulationParameterInformation(pop.ini)
  }
}


#-------------------------------------------------
pk.estim <- function(r, admin) {  
  time <- NULL
  g=mlx.getObservationInformation()
  gn <- g$name[[1]]
  gy <- g[[gn]]
  
  treat <- r$treatment
  if (!is.null(treat$rate))
    treat$tinf <- treat$amount/treat$rate
  else
    treat$tinf <- 0
  abs <- elim <- max1 <- trid <- amount2 <- NULL
  for (id in r$id) {
    ji <- which(treat$id==id)
    tri <- treat[ji,]
    tri <- tri[order(tri$time),]
    
    ji <- which(gy$id==id)
    yi <- gy[ji,]
    yi <- yi[yi[['time']]>= min(tri[['time']]),]
    
    jty1 <- which(tri[['time']]<=min(yi[['time']]))
    tri <- tri[max(jty1):nrow(tri), ]
    jty2 <- which(tri[['time']]>max(yi[['time']]))
    if (length(jty2)>1)
      tri <- tri[1:(min(jty2)-1), ]
    
    ndi <- nrow(tri)
    tri.inf <- tri[ndi,]
    tri.inf['time'] <- Inf
    tri <- rbind(tri, tri.inf)
    
    tri1 <- tri[1:2,]
    tri2 <- tri[ndi:(ndi+1),]
    #trid <- rbind(trid, tri1[1,])
    
    yi1 <- subset(yi, time>=tri1$time[1] & time<tri1$time[2] )
    if (nrow(yi1)>0) {
      yi1[gn] <- yi1[gn]/tri1[['amount']][1]
      yi1['time'] <- yi1['time'] - tri1$time[1]
      if (admin=="oral") {
        j.max1 <- which.max(yi1[[gn]])
        if (j.max1 == 1) {
          yi1 <- rbind(yi1[1,], yi1)
          yi1[1,'time'] <- yi1[1,'time']/4
          yi1[1,gn] <- yi1[1,gn]/2
          j.max1 <- 2
        }
        if (length(j.max1)>0 && j.max1>1) {
          abs <- rbind(abs, yi1[1:(j.max1-1),])
          # if (length(j.max1)>0) {
          #   abs <- rbind(abs, yi1[1:(j.max1),])
          if (j.max1<nrow(yi1)) 
            max1 <- rbind(max1,yi1[j.max1,])
        }
      } else {
        tinfi <- tri1[1,"tinf"]
        #yi1 <- yi1[yi1$time <= tinfi,]
        yi1$tinf <- tinfi
        abs <- rbind(abs, yi1)
      }
    }
    
    yi2 <- subset(yi, time>=tri2$time[1] & time<tri2$time[2] )
    if (nrow(yi2)>0) {
      amti2 <- tri2[['amount']][1]
      yi2[gn] <- yi2[gn]/amti2
      amount2 <- c(amount2, amti2)
      yi2['time'] <- yi2['time'] - tri2$time[1]
      yi2['amount'] <- amti2
      j.max2 <- which.max(yi2[[gn]])
      ni <- nrow(yi2)
      if (length(j.max2)>0 && j.max2<ni) 
        elim <- rbind(elim, yi2[(j.max2):ni,])
    }
  }
  if (!is.null(abs)) names(abs)[which(names(abs)==gn)] <- "y"
  if (!is.null(elim)) names(elim)[which(names(elim)==gn)] <- "y"
  if (!is.null(max1)) names(max1)[which(names(max1)==gn)] <- "y"
  return(list(abs=abs, elim=elim, max1=max1, amount2=amount2))
}


#-------------------------------------------------
compute.ini <- function(r, parameter) {
  
  y <- NULL
  if (("ka" %in% parameter) | ("Tk0" %in% parameter))
    admin <- "oral"
  else
    admin <- "iv"
  
  th <- pk.estim(r, admin)
  abs <- th$abs
  elim <- th$elim
  
  k_ini <- -lm(log(y) ~ time, data=subset(elim, y>0))$coefficients[[2]]
  
  if (admin=="oral") {
    ymax <- th$max1$y
    tmax <- th$max1$time
    ka_ini <- lm(log(y) ~  time, data=subset(abs, y>0))$coefficients[[2]]
    if (is.na(ka_ini) || ka_ini < 0)
      ka_ini <- 1
    #    ka_ini <- lm(y ~ -1 + time, data=subset(abs, y>0))$coefficients[[1]]
    Tk0_ini <- mean(tmax)
    if (ka_ini>0)
      V_ini <- 1/mean(ymax)*ka_ini/abs(ka_ini-k_ini)
    else
      V_ini <- 1/(Tk0_ini*k_ini*mean(ymax))*(1-exp(-k_ini*Tk0_ini))
    Tlag_ini <- Tk0_ini/5
    Mtt_ini <- Tlag_ini
    Ktr_ini <- ka_ini*5
    list.ini <- c(ka=ka_ini, V=V_ini, k=k_ini, Tk0=Tk0_ini, Tlag=Tlag_ini, Mtt=Mtt_ini, Ktr=Ktr_ini)
  } else {
    rV <- (1 - exp(-k_ini*abs$tinf)) /(abs$tinf*k_ini)
    rV[abs$tinf==0] <- 1
    dt <- pmax(abs$time - abs$tinf, 0)
    rV <- rV*exp(-k_ini*dt)
    #V_ini <- exp(mean(log(rV/abs$y)))
    V_ini <- mean(rV^2)/mean(rV*abs$y)
    list.ini <- c(V=V_ini, k=k_ini)
  }
  Cl_ini <- k_ini*V_ini
  Cmax <- aggregate(elim$y*elim$amount, by=list(elim$id), FUN=max)
  Km_ini <- mean(Cmax[,2])
  Vm_ini <- Cl_ini*(2*Km_ini)
  k12_ini <- k_ini/2
  k21_ini <- k_ini/2
  k13_ini <- k_ini/4
  k31_ini <- k_ini/4
  # k13_ini <- k_ini/2
  # k31_ini <- k_ini/2
  
  list.ini <- c(list.ini, Cl=Cl_ini, Km=Km_ini, Vm=Vm_ini)
  list.ini <- c(list.ini, k12=k12_ini, k21=k21_ini, k13=k13_ini, k31=k31_ini)
  list.ini <- c(list.ini, V1=V_ini, Q=k12_ini*V_ini, Q2=k12_ini*V_ini, V2=k12_ini/k21_ini*V_ini, 
                Q3=k13_ini*V_ini, V3=k13_ini/k31_ini*V_ini)
  return(list.ini[parameter])
}

#-------------------------------------------------
err <-  function(parameter, y, p.ind, N, a) {
  p.ind[,] <- matrix(exp(parameter),nrow=N,ncol=length(parameter),byrow=TRUE)
  f <- as.numeric(mlx.computePredictions(p.ind)[[1]])
  if (any(is.nan(f)) | any(is.infinite(f)))
    e <- Inf
  else
    e <- mean((log(f+a)-log(y+a))^2)
  #    e <- mean((f^a-y^a)^2)
  return(e)
}

#-------------------------------------------------
pop.opt <- function(p0) {
  #setPKproject(parameter=p0)
  g=mlx.getObservationInformation()
  gn <- g$name[[1]]
  gy <- g[[gn]]
  N <- length(unique(gy[['id']]))
  y <- gy[[gn]]
  if (N>1)
    p.ind <- as.data.frame(t(p0)[rep(1,N),])
  else
    p.ind <- p0
  a <- max(-min(y) + 0.5, 0.5)
  #  if ("k12" %in% names(p0))  browser()
  r <- optim(log(p0), err, y=y, p.ind=p.ind, N=N, a=a)
  return(exp(r$par))
}


#-------------------------------------------------
compute.bic <- function(parameter, data, new.dir=NULL, level=NULL, par.ini=NULL, linearization=F) {
  cat("\n")
  r <- pkpopini(parameter=parameter, data=data, new.dir=new.dir, par.ini=par.ini) 
  print(r$project)
  mlx.loadProject(projectFile = r$project)
  g=mlx.getObservationInformation()
  gn <- g$name[[1]]
  gy <- g[[gn]]
  N <- length(unique(gy[['id']]))
  n <- nrow(gy)
  scenario <- mlx.getScenario()
  scenario$tasks[seq_along(scenario$tasks)] <- FALSE
  scenario$tasks[c("populationParameterEstimation", "logLikelihoodEstimation")] <- TRUE
  if (linearization)
    scenario$tasks[c("conditionalModeEstimation")] <- TRUE
  else
    scenario$tasks[c("conditionalDistributionSampling")] <- TRUE
  mlx.setScenario(scenario)
  
  if (!is.null(level))
    setSettings(level=level)
  mlx.saveProject(r$project)
  w.dir <- getwd()
  setwd(new.dir)
  
  launched.tasks <- mlx.getLaunchedTasks()
  if (!launched.tasks[["populationParameterEstimation"]]) {
    cat("Estimation of the population parameters...\n")
    mlx.runPopulationParameterEstimation()
  }
  if (linearization) {
    if (!launched.tasks[["conditionalModeEstimation"]])
      cat("Estimation of the conditional modes... \n")
    mlx.runConditionalModeEstimation()
    if (!("linearization" %in% launched.tasks[["logLikelihoodEstimation"]])) { 
      cat("Estimation of the log-likelihood (linearization)... \n")
      mlx.runLogLikelihoodEstimation(linearization=T)
      g <- mlx.getEstimatedLogLikelihood()[["linearization"]]
    }
  } else {
    if (!launched.tasks[["conditionalDistributionSampling"]])
      cat("Estimation of the conditional distributions... \n")
    mlx.runConditionalDistributionSampling()
    if (!("importanceSampling" %in% launched.tasks[["logLikelihoodEstimation"]])) { 
      cat("Estimation of the log-likelihood (importance sampling)... \n")
      mlx.runLogLikelihoodEstimation(linearization=F)
      g <- mlx.getEstimatedLogLikelihood()[["importanceSampling"]]
    }
  }
  
  setwd(w.dir)
  r$ofv <- g[['OFV']]
  r$bicc <- g[['BICc']]
  r$bic <- g[['BIC']]
  r$aic <- g[['AIC']]
  r$pop.est <- mlx.getEstimatedPopulationParameters()
  print(g)
  return(r)
}


#-------------------------------------------------
read.res <- function(file) {
  d <- read.csv(file, sep=",")
  if (ncol(d)==1)
    d <- read.csv(file, sep=";")
  if (ncol(d)==1)
    d <- read.csv(file, sep=" ")
  if (ncol(d)==1)
    d <- read.csv(file, sep="\t")
  return(d)
}

# Get the extension of a file --------------------------------------------------
.getFileExt <- function(filename) {
  ext <- utils::tail(strsplit(filename, "\\.")[[1]], n=1)
  if (ext == filename) ext <- NULL
  return(ext)
}

# Round dataframe --------------------------------------------------------------
.roundDataframe <- function(df, nbDigits) {
  for (n in names(df)) {
    df[n] <- sapply(
      df[[n]],
      function(x, nbDigits) {
        if (suppressWarnings(!is.na(as.numeric(x)))) x <- round(as.numeric(x), nbDigits)
        return(x)
      },
      nbDigits
    )
  }
  return(df)
}

# Get the separator of a file --------------------------------------------------
.getDelimiter <- function(fileName, sep = NULL){
  L <- suppressMessages(suppressWarnings(readLines(fileName, n = 1)))
  L <- gsub("#","", L)
  sepToCheck = c(' ', '\t', ',', ';')
  if (!is.null(sep)) sepToCheck <- unique(c(sep, sepToCheck))
  nSepToCheck <- length(sepToCheck)
  numfields <- vector(length = nSepToCheck)
  for(index in 1:nSepToCheck){
    numfields[index] <- utils::count.fields(textConnection(L), sep = sepToCheck[index])
  }
  if (all(numfields == 0)) {
    sep <- NULL
  } else {
    sep <- sepToCheck[min(which.max(numfields))]
  }
  return(sep)
}

# Rename dataframe column ------------------------------------------------------
.renameColumns <- function(df, oldName, newName){
  if (length(oldName) != length(newName)) {
    message("[ERROR] vector of old names and new names must match in size")
    return(df)
  }
  for (i in seq_along(oldName)) {
    old <- oldName[i]
    new <- newName[i]
    if (old %in% names(df)) {
      names(df)[names(df) == old] <- new
    }
  }
  return(df)
}

################################################################################
# Set lixoft connectors options
################################################################################
set_options <- function(errors = NULL, warnings = NULL, info = NULL) {
  options_list <- list(errors = errors, warnings = warnings, info = info)
  options_list <- options_list[!unlist(lapply(options_list, is.null))]
  op <- getOption("lixoft_notificationOptions")
  for (i in seq_along(options_list)) {
    oname <- names(options_list)[[i]]
    oval <- options_list[[i]]
    if (!is.null(oval)) {
      if (!oval) {
        op[[oname]] <- 1
      } else {
        op[[oname]] <- 0
      }
    }
  }
  options(lixoft_notificationOptions = op)
  return(invisible(TRUE))
}

get_lixoft_options <- function() {
  op <- getOption("lixoft_notificationOptions")
  errors <- ifelse(op$errors == 1, FALSE, TRUE)
  warnings <- ifelse(op$warnings == 1, FALSE, TRUE)
  info <- ifelse(op$info == 1, FALSE, TRUE)
  return(list(errors = errors, warnings = warnings, info = info))
}

############################################################################

def.variable <- function(weight=NULL, prior=NULL, criterion=NULL, fix.param0=NULL, fix.param1=NULL) {
  go <- mlx.getObservationInformation()
  y.name <- names(go$mapping)
  y.type <- go$type[y.name]
  id <- n <- NULL
  for (yn in y.name) {
    id <- unique(c(id, go[[yn]]$id))
    n <- c(n, nrow(go[[yn]]))
  }
  N <- length(id)
  nc <- n[which(y.type=="continuous")]
  
  if (identical(criterion,"AIC")) {
    pen.coef <- rep(2, length(nc)+2)  
  } else if (identical(criterion,"BIC")) {
    pen.coef <- rep(log(N), length(nc)+2)
  } else if (identical(criterion,"BICc")) {
    pen.coef <- c(log(N), rep(log(sum(n)), length(nc)+1))
    #  pen.coef <- c(log(N), log(sum(n)), rep(log(sum(nc)), length(nc)))
    #    pen.coef <- c(log(N), log(sum(nc)), log(nc))
  } else {
    pen.coef <- rep(criterion, length(nc)+2)
  }
  
  if (is.null(weight$is.weight)) { 
    if (is.null(weight) & is.null(prior))
      is.weight <- F
    else
      is.weight <- T
  } else {
    is.weight <- weight$is.weight
  }
  
  w.cov <- weight$covariate
  if (!is.null(prior$covariate)) 
    w.cov <- -2*log(prior$covariate/(1-prior$covariate))/pen.coef[1]
  if (is.null(w.cov))
    w.cov <- 1
  if (length(w.cov)==1) {
    cov.model <- do.call(rbind, mlx.getIndividualParameterModel()$covariateModel)
    foo <- w.cov
    w.cov <- cov.model
    w.cov[] <- foo
  }
  
  w.cor <- weight$correlation
  if (!is.null(prior$correlation)) 
    w.cor <- -2*log(prior$correlation/(1-prior$correlation))/pen.coef[1]
  if (is.null(w.cor))
    w.cor <- 1
  if (length(w.cor)==1) {
    var.model <- mlx.getIndividualParameterModel()$variability$id
    n.param <- names(which(var.model))
    d.param <- length(n.param)
    w.cor <- matrix(w.cor, nrow=d.param, ncol=d.param, dimnames=list(n.param, n.param))
  }
  if (!is.null(w.cor))
    w.cor[fix.param0, ] <- w.cor[, fix.param0] <- 1e10
  w.cor <- w.cor*lower.tri(w.cor)
  w.cor[which(is.nan(w.cor))] <- 0
  
  if (!is.null(prior$variance)) 
    w.var <- -2*log(prior$variance/(1-prior$variance))/pen.coef[1]
  else
    w.var <- weight$variance 
  foo <- mlx.getIndividualParameterModel()$variability$id
  foo[] <- 1
  if (length(w.var)==1 & is.null(names(w.var)))
    foo[] <- w.var
  foo[names(w.var)] <- w.var
  w.var <- foo
  # w.var[fix.param0] <- 1e10
  # w.var[fix.param1] <- 0
  
  
  weight <- list(covariate=w.cov, correlation=w.cor, variance=w.var, is.weight=is.weight)
  
  return(list(n=n, N=N, weight=weight, pen.coef=pen.coef))
}

