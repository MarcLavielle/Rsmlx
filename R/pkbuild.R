#' Automatic PK model building
#'
#' Fit several structural PK models and select the best one based on a Bayesian Information Criterion.
#' Models to compare can be defined by rate constants and/or clearances and can include or not nonlinear elimination
#' models.
#' See http://rsmlx.webpopix.org/pkbuild/ for more details.
#' @param data a list with fields
#' \itemize{
#'   \item \code{dataFile}: path of a formatted data file
#'   \item \code{headerTypes}: a vector of strings
#'   \item \code{administration} ("iv", "bolus", "infusion", "oral", "ev"): route of administration 
#' }
#' @param project a Monolix project
#' @param stat ({FALSE}, TRUE): the statistical model is also built (using buildmlx)
#' @param param ({"clearance"}, "rate", "both): parameterization 
#' @param new.dir   name of the directory where the created files are stored 
#' (default is the current working directory) )
#' @param MM   ({FALSE}, TRUE): tested models include or not Michaelis Menten elimination models
#' @param linearization  {TRUE}/FALSE whether the computation of the likelihood is based on a linearization of the model (default=TRUE)
#' @param criterion  penalization criterion to optimize c("AIC", "BIC", {"BICc"}, gamma)
#' @param level an integer between 1 and 9 (used by setSettings)
#' @param settings.stat list of settings used by buildmlx (only if stat=TRUE)
#' 
#' @return A list of results 
#' @examples
#' \dontrun{
#' # Build a PK model for the warfarin PK data. 
#' # By default, only models using clearance (and inter compartmental clearances) are used
#' warf.pk1 <- pkbuild(data=warfarin)
#' 
#' # Models using elimination and transfer rate constants are used, 
#' # as well as nonlinear elimination models
#' warf.pk2 <- pkbuild(data=warfarin,  new.dir="warfarin", param="rate", MM=TRUE)
#' 
#' # Both models using clearances and rates are used. 
#' # Level is set to 7 in order to get accurate results.
#' warf.pk3 <- pkbuild(data=warfarin,  new.dir="warfarin", param="both", level=7)
#' }
#' @importFrom stats aggregate optim 
#' @export
pkbuild <- function(data=NULL, project=NULL, stat=FALSE, param="clearance", new.dir=".", 
                    MM=FALSE, linearization=T, criterion="BICc", level=NULL, settings.stat=NULL) {
  
  if (!initRsmlx()$status)
    return()
  
  RsmlxDemo1.project <- RsmlxDemo2.project <- warfarin.data  <- resMonolix <- NULL
  
  if (!is.null(project)) {
    r <- prcheck(project)
    data <- mlx.getData()[c('dataFile', 'headerTypes')]
    parameter <- mlx.getIndividualParameterModel()$name
    if ("ka" %in% parameter | "Tk0" %in% parameter)
      data$administration <- "oral"
    else
      data$administration <- "iv"
  }
  
  check.pkbuild(data, param, MM, level)
  
  if (param=="both") {
    r.cl <- pkbuild(data=data, param="clearance", new.dir=new.dir, MM=MM, level=level, linearization=linearization)
    r.k  <- pkbuild(data=data, param="rate", new.dir=new.dir, MM=MM, level=level, linearization=linearization)
    if (criterion == "AIC")
      test <- r.cl$res$aic[[1]] < r.k$res$aic[[1]]
    else if (criterion == "BIC")
      test <- r.cl$res$bic[[1]] < r.k$res$bic[[1]]
    else
      test <- r.cl$res$bicc[[1]] < r.k$res$bicc[[1]]
    
    if (test) 
      r.final <- r.cl
    else
      r.final <- r.k
    df.res <- rbind(r.cl$res, r.k$res)
    if (criterion == "AIC")
      r.final$res <- df.res[order(df.res$aic),]
    else if (criterion == "BIC")
      r.final$res <- df.res[order(df.res$bic),]
    else
      r.final$res <- df.res[order(df.res$bicc),]
    r.final$res <- r.final$res[!duplicated(r.final$res), ]
    row.names(r.final$bic) =NULL
    if (stat) 
      r.final <- pkbuild.stat(r.final, settings.stat)
    return(r.final)
  }
  
  if (!is.null(new.dir) && !dir.exists(new.dir))
    dir.create(new.dir)
  
  if (param=="rate") {
    p.lin1 <- c("V", "k")
    p.lin2 <- c("V", "k12", "k21", "k")
    p.lin3 <- c("V", "k12", "k21", "k13", "k31", "k")
    p.mm1 <- c("V", "Vm", "Km")
    p.mm2 <- c("V", "k12", "k21", "Vm", "Km")
    p.mm3 <- c("V", "k12", "k21", "k13", "k31", "Vm", "Km")
  } else {
    p.lin1 <- c("V", "Cl")
    p.lin2 <- c("V1", "V2", "Q", "Cl")
    p.lin3 <- c("V1", "V2", "V3", "Q2", "Q3", "Cl")
    p.mm1 <- c("V", "Vm", "Km")
    p.mm2 <- c("V1", "V2", "Q", "Vm", "Km")
    p.mm3 <- c("V1", "V2", "V3", "Q2", "Q3", "Vm", "Km")
  }
  
  r.model <- r.aic <- r.bic <- r.bicc <- NULL
  admin <- data$administration
  if (admin %in% c("oral", "ev")) {
    if (param=="rate") {
      p.abs1 <- c("ka", "V", "k")
      p.abs2 <- c("Tk0", "V", "k")
      p.abs3 <- c("Tlag", "ka", "V", "k")
      p.abs4 <- c("Tlag", "Tk0", "V", "k")
    } else {
      p.abs1 <- c("ka", "V", "Cl")
      p.abs2 <- c("Tk0", "V", "Cl")
      p.abs3 <- c("Tlag", "ka", "V", "Cl")
      p.abs4 <- c("Tlag", "Tk0", "V", "Cl")
    }
    r.abs1 <- compute.bic(parameter=p.abs1, data=data, new.dir=new.dir, level=level, linearization=linearization) 
    r.abs2 <- compute.bic(parameter=p.abs2, data=data, new.dir=new.dir, level=level, linearization=linearization) 
    par.ini <- c(Tlag=1, r.abs1$pop.est[1:3])
    names(par.ini)=gsub("_pop","",names(par.ini))
    r.abs3 <- compute.bic(parameter=p.abs3, data=data, new.dir=new.dir, level=level, linearization=linearization, par.ini=par.ini) 
    par.ini <- c(r.abs3$pop.est[1], r.abs2$pop.est[1:3])
    names(par.ini)=gsub("_pop","",names(par.ini))
    r.abs4 <- compute.bic(parameter=p.abs4, data=data, new.dir=new.dir, level=level, linearization=linearization, par.ini=par.ini) 
    r.model <- c(r.model, c(r.abs1$model, r.abs2$model, r.abs3$model, r.abs4$model))
    r.aic   <- c(r.aic, c(r.abs1$aic, r.abs2$aic, r.abs3$aic, r.abs4$aic))
    r.bic   <- c(r.bic, c(r.abs1$bic, r.abs2$bic, r.abs3$bic, r.abs4$bic))
    r.bicc   <- c(r.bicc, c(r.abs1$bicc, r.abs2$bicc, r.abs3$bicc, r.abs4$bicc))
    if (criterion == "AIC")
      oabs <- order(r.aic)
    else if (criterion == "BIC")
      oabs <- order(r.bic)
    else
      oabs <- order(r.bicc)
    eval(parse(text = paste0("p.absa <- p.abs",oabs[1])))
    eval(parse(text = paste0("p.absb <- p.abs",oabs[2])))
    p.absa <- p.absa[1:(length(p.absa)-2)]
    p.absb <- p.absb[1:(length(p.absb)-2)]
    p.lin1 <- c(p.absa, p.lin1)
    p.lin2a <- c(p.absa, p.lin2)
    p.lin2b <- c(p.absb, p.lin2)
    p.lin3a <- c(p.absa, p.lin3)
    p.lin3b <- c(p.absb, p.lin3)
    p.mm1 <- c(p.absa, p.mm1)
    p.mm2a <- c(p.absa, p.mm2)
    p.mm3a <- c(p.absa, p.mm3)
    eval(parse(text = paste0("r.lin1 <- r.abs",oabs[1])))
    
    par.ini <- r.lin1$pop.est[grep("_pop", names(r.lin1$pop.est))]
    names(par.ini)=gsub("_pop","",names(par.ini))
    if (param=="rate") {
      par.ini[['V']] <- par.ini[['V']]/2
      par.ini <- c(par.ini, k12=par.ini[["k"]], k21=par.ini[["k"]])
    } else {
      par.ini <- c(par.ini, V1=par.ini[["V"]]/2, V2=par.ini[["V"]]/2, Q=par.ini[["Cl"]]/2)
      par.ini <- par.ini[names(par.ini) != "V"]
    }
    #    browser()
    r.lin2a <- compute.bic(parameter=p.lin2a, data=data, new.dir=new.dir, level=level, linearization=linearization, par.ini=par.ini) 
    par.ini <- c(r.abs3$pop.est['Tlag_pop'], r.lin2a$pop.est[grep("_pop", names(r.lin2a$pop.est))])
    names(par.ini)=gsub("_pop","",names(par.ini))
    par.ini <- par.ini[p.lin2b]
    r.lin2b <- compute.bic(parameter=p.lin2b, data=data, new.dir=new.dir, level=level, linearization=linearization, par.ini=par.ini) 
    r.model <- c(r.model, r.lin2a$model, r.lin2b$model)
    r.aic <- c(r.aic, r.lin2a$aic, r.lin2b$aic)
    r.bic <- c(r.bic, r.lin2a$bic, r.lin2b$bic)
    r.bicc <- c(r.bicc, r.lin2a$bicc, r.lin2b$bicc)
    if (criterion == "AIC")
      test <- r.lin2a$aic < r.lin2b$aic
    else if (criterion == "BIC")
      test <- r.lin2a$bic < r.lin2b$bic
    else
      test <- r.lin2a$bicc < r.lin2b$bicc
    
    if (test) {
      r.lin2 <- r.lin2a
      p.lin2 <- p.lin2a
      p.lin3 <- p.lin3a
    } else {
      r.lin2 <- r.lin2b
      p.lin2 <- p.lin2b
      p.lin3 <- p.lin3b
    }
  } else {
    r.lin1 <- compute.bic(parameter=p.lin1, data=data, new.dir=new.dir, level=level, linearization=linearization) 
    par.ini <- r.lin1$pop.est[grep("_pop", names(r.lin1$pop.est))]
    names(par.ini)=gsub("_pop","",names(par.ini))
    if (param=="rate") {
      par.ini <- c(par.ini, k12=par.ini[["k"]]/2, k21=par.ini[["k"]]/2)
    } else {
      par.ini <- c(par.ini, V2=par.ini[["V"]]/2, Q=par.ini[["Cl"]]/2)
      par.ini <- par.ini[names(par.ini) != "V"]
    }
    r.lin2 <- compute.bic(parameter=p.lin2, data=data, new.dir=new.dir, level=level, linearization=linearization, par.ini=par.ini) 
    r.model <- c(r.model, r.lin1$model, r.lin2$model)
    r.aic <- c(r.aic, r.lin1$aic, r.lin2$aic)
    r.bic <- c(r.bic, r.lin1$bic, r.lin2$bic)
    r.bicc <- c(r.bicc, r.lin1$bicc, r.lin2$bicc)
  }
  if (criterion == "AIC")
    test <- r.lin2$aic < r.lin1$aic
  else if (criterion == "BIC")
    test <- r.lin2$bic < r.lin1$bic
  else
    test <- r.lin2$bicc < r.lin1$bicc
  if (test) {
    par.ini <- r.lin2$pop.est[grep("_pop", names(r.lin2$pop.est))]
    names(par.ini)=gsub("_pop","",names(par.ini))
    if (param=="rate") {
      par.ini <- c(par.ini, k13=par.ini[["k12"]]/2, k31=par.ini[["k21"]]/2)
    } else {
      par.ini <- c(par.ini, V3=par.ini[["V2"]]/2, Q3=par.ini[["Q"]]/2, Q2=par.ini[["Q"]])
      par.ini <- par.ini[names(par.ini) != "Q"]
    }
    r.lin3 <- compute.bic(parameter=p.lin3, data=data, new.dir=new.dir, level=level, linearization=linearization, par.ini=par.ini) 
    r.model <- c(r.model, r.lin3$model)
    r.aic <- c(r.aic, r.lin3$aic)
    r.bic <- c(r.bic, r.lin3$bic)
    r.bicc <- c(r.bicc, r.lin3$bicc)
    if (criterion == "AIC")
      test <- r.lin3$aic < r.lin2$aic
    else if (criterion == "BIC")
      test <- r.lin3$bic < r.lin2$bic
    else
      test <- r.lin3$bicc < r.lin2$bicc
    if (test) {
      cpt <- 3
      r.final <- r.lin3
    } else {
      cpt <- 2
      r.final <- r.lin2
    } 
  } else {
    cpt <- 1
    r.final <- r.lin1
  }
  
  elim <- "lin"
  if (MM) {
    if (cpt==1) {
      r.mm  <- compute.bic(parameter=p.mm1, data=data, new.dir=new.dir, level=level, linearization=linearization) 
    } else if (cpt==2) {
      r.mm  <- compute.bic(parameter=p.mm2, data=data, new.dir=new.dir, level=level, linearization=linearization) 
    } else {
      r.mm  <- compute.bic(parameter=p.mm3, data=data, new.dir=new.dir, level=level, linearization=linearization) 
    }
    r.model <- c(r.model, r.mm$model)
    r.aic <- c(r.aic, r.mm$aic)
    r.bic <- c(r.bic, r.mm$bic)
    r.bicc <- c(r.bicc, r.mm$bicc)
    if (criterion == "AIC")
      test <- r.mm$aic < r.final$aic
    else if (criterion == "BIC")
      test <- r.mm$bic < r.final$bic
    else
      test <- r.mm$bicc < r.final$bicc
    if (test) {
      r.final <- r.mm
    }
  }
  r.model <- unlist(lapply(as.character(r.model), function(x) basename(x)))
  df.res <- data.frame(model=r.model, AIC=r.aic, BIC=r.bic, BICc=r.bicc)
  if (criterion == "AIC")
    df.res <- df.res[order(df.res$AIC),]
  else if (criterion == "BIC")
    df.res <- df.res[order(df.res$BIC),]
  else
    df.res <- df.res[order(df.res$BICc),]
  r.final$res <- df.res
  row.names(r.final$res) =NULL
  r.final$pini <- NULL
  
  if (stat) 
    r.final <- pkbuild.stat(r.final, settings.stat)
  
  return(r.final)
}

pkbuild.stat <- function(r.final, settings.stat) {
  res.stat <- do.call("buildAll", c(list(project=r.final$project), settings.stat))
  mlx.loadProject(res.stat$project)
  r.final$project <- res.stat$project
  # r.final$covariate.model <- res.stat$covariate.model
  # r.final$correlation.model <- res.stat$correlation.model
  # r.final$error.model <- res.stat$error.model
  # r.final$pop.est <- mlx.getEstimatedPopulationParameters()
  return(r.final)
}

check.pkbuild <- function(data, param, MM, level) {
  if (is.null(data$dataFile)) 
    stop("A data file should be provided in data$dataFile", call.=FALSE)
  if (!file.exists(data$dataFile)) 
    stop(paste0(data$dataFile, " does not exists"), call.=FALSE)
  
  if (is.null(data$headerTypes)) 
    stop("Header types should be provided as a vector of strings in data$headerTypes", call.=FALSE)
  
  headerList = c('ID','TIME','AMT','ADM','RATE','TINF','Y','YTYPE', 'X','COV','CAT','OCC','MDV', 'EVID',
                 'ADDL','SS','II', 'DPT', 'AMOUNT', 'OBSERVATION', 'OBSID', 'CONTCOV', 'CATCOV', 'IGNORE',
                 'cens', 'limit', 'regressor', 'admid', 'date')
  
  if (!all(tolower(data$headerTypes) %in% tolower(headerList) ))
    stop("Valid header types should be provided", call.=FALSE)
  
  if (is.null(data$administration)) 
    stop("Route of administration should be provided in data (e.g. data$administation='oral')", call.=FALSE)
  
  admList <- c('IV', 'BOLUS', 'INFUSION', 'ORAL', 'EV')
  if (!(toupper(data$administration) %in% admList ))
    stop("A valid route of administration should be provided among \n{'iv', 'bolus', 'infusion', 'oral', 'ev'}", call.=FALSE)
  
  if (!is.logical(MM))
    stop("MM (Michaelis Menten) should be either TRUE or FALSE", call.=FALSE)
  
  if (!is.null(level) && !(level %in% 1:9))
    stop("level should be an integer between 1 and 9 (see the setSettings documentation)", call.=FALSE)
}
