stepAICmodified<-function (object, scope, scale = 0, direction = c("both", "backward", 
                                                           "forward"), trace = 1, keep = NULL, steps = 1000, use.start = FALSE, 
                   k = 2,prior=NULL) #Was +...
{
  mydeviance <- function(x) {#Was +...
    dev <- deviance(x)
    if (!is.null(dev)) 
      dev
    else extractAIC(x, k = 0)[2L]+extractAIC.getpenalisationprior(fit,prior)
  }
  cut.string <- function(string) {
    if (length(string) > 1L) 
      string[-1L] <- paste("\n", string[-1L], sep = "")
    string
  }
  re.arrange <- function(keep) {
    namr <- names(k1 <- keep[[1L]])
    namc <- names(keep)
    nc <- length(keep)
    nr <- length(k1)
    array(unlist(keep, recursive = FALSE), c(nr, nc), list(namr, 
                                                           namc))
  }
  step.results <- function(models, fit, object, usingCp = FALSE) {
    change <- sapply(models, "[[", "change")
    rd <- sapply(models, "[[", "deviance")
    dd <- c(NA, abs(diff(rd)))
    rdf <- sapply(models, "[[", "df.resid")
    ddf <- c(NA, abs(diff(rdf)))
    AIC <- sapply(models, "[[", "AIC")
    heading <- c("Stepwise Model Path \nAnalysis of Deviance Table", 
                 "\nInitial Model:", deparse(formula(object)), "\nFinal Model:", 
                 deparse(formula(fit)), "\n")
    aod <- if (usingCp) 
      data.frame(Step = change, Df = ddf, Deviance = dd, 
                 `Resid. Df` = rdf, `Resid. Dev` = rd, Cp = AIC, 
                 check.names = FALSE)
    else data.frame(Step = change, Df = ddf, Deviance = dd, 
                    `Resid. Df` = rdf, `Resid. Dev` = rd, AIC = AIC, 
                    check.names = FALSE)
    attr(aod, "heading") <- heading
    class(aod) <- c("Anova", "data.frame")
    fit$anova <- aod
    fit
  }
  Terms <- terms(object)
  object$formula <- Terms
  if (inherits(object, "lme")) 
    object$call$fixed <- Terms
  else if (inherits(object, "gls")) 
    object$call$model <- Terms
  else object$call$formula <- Terms
  if (use.start) 
    warning("'use.start' cannot be used with R's version of 'glm'")
  md <- missing(direction)
  direction <- match.arg(direction)
  backward <- direction == "both" | direction == "backward"
  forward <- direction == "both" | direction == "forward"
  if (missing(scope)) {
    fdrop <- numeric()
    fadd <- attr(Terms, "factors")
    if (md) 
      forward <- FALSE
  }
  else {
    if (is.list(scope)) {
      fdrop <- if (!is.null(fdrop <- scope$lower)) 
        attr(terms(update.formula(object, fdrop)), "factors")
      else numeric()
      fadd <- if (!is.null(fadd <- scope$upper)) 
        attr(terms(update.formula(object, fadd)), "factors")
    }
    else {
      fadd <- if (!is.null(fadd <- scope)) 
        attr(terms(update.formula(object, scope)), "factors")
      fdrop <- numeric()
    }
  }
  models <- vector("list", steps)
  if (!is.null(keep)) 
    keep.list <- vector("list", steps)
  n <- nobs(object, use.fallback = TRUE)
  fit <- object
  bAIC <- extractAIC(fit, scale, k = k)  #Was +...
  edf <- bAIC[1L]
  bAIC <- bAIC[2L]+extractAIC.getpenalisationprior(fit,prior)
  if (is.na(bAIC)) 
    stop("AIC is not defined for this model, so 'stepAIC' cannot proceed")
  if (bAIC == -Inf) 
    stop("AIC is -infinity for this model, so 'stepAIC' cannot proceed")
  nm <- 1
  Terms <- terms(fit)
  if (trace) {
    cat("Start:  AIC=", format(round(bAIC, 2)), "\n", cut.string(deparse(formula(fit))), 
        "\n\n", sep = "")
    utils::flush.console()
  }
  models[[nm]] <- list(deviance = mydeviance(fit), df.resid = n - 
                         edf, change = "", AIC = bAIC)
  if (!is.null(keep)) 
    keep.list[[nm]] <- keep(fit, bAIC)
  usingCp <- FALSE
  while (steps > 0) {
    steps <- steps - 1
    AIC <- bAIC
    ffac <- attr(Terms, "factors")
    if (!is.null(sp <- attr(Terms, "specials")) && !is.null(st <- sp$strata)) 
      ffac <- ffac[-st, ]
    scope <- factor.scope(ffac, list(add = fadd, drop = fdrop))
    aod <- NULL
    change <- NULL
    if (backward && length(scope$drop)) {
      aod <- dropterm(fit, scope$drop, scale = scale, trace = max(0, 
                                                                  trace - 1), k = k)#Was +...
      print("backward - no prior penality")
      print(aod)
      aod<-extractAIC.updateAICprior(fit,aod,prior,drop=TRUE)
      print("backward - with prior penality")
      print(aod)
      #browser()
      rn <- row.names(aod)
      row.names(aod) <- c(rn[1L], paste("-", rn[-1L], sep = " "))
      if (any(aod$Df == 0, na.rm = TRUE)) {
        zdf <- aod$Df == 0 & !is.na(aod$Df)
        nc <- match(c("Cp", "AIC"), names(aod))
        nc <- nc[!is.na(nc)][1L]
        ch <- abs(aod[zdf, nc] - aod[1, nc]) > 0.01
        if (any(is.finite(ch) & ch)) {
          warning("0 df terms are changing AIC")
          zdf <- zdf[!ch]
        }
        if (length(zdf) > 0L) 
          change <- rev(rownames(aod)[zdf])[1L]
      }
    }
    if (is.null(change)) {
      if (forward && length(scope$add)) {
        aodf <- addterm(fit, scope$add, scale = scale, 
                        trace = max(0, trace - 1), k = k)#Was +...
        print("forward - no prior penality")
        print(aodf)

        aodf<-extractAIC.updateAICprior(fit,aodf,prior,drop=FALSE)
        print("forward - with prior penality")
        print(aodf)
        #browser()
        
        rn <- row.names(aodf)
        row.names(aodf) <- c(rn[1L], paste("+", rn[-1L], 
                                           sep = " "))
        aod <- if (is.null(aod)) 
          aodf
        else rbind(aod, aodf[-1, , drop = FALSE])
      }
      attr(aod, "heading") <- NULL
      if (is.null(aod) || ncol(aod) == 0) 
        break
      nzdf <- if (!is.null(aod$Df)) 
        aod$Df != 0 | is.na(aod$Df)
      aod <- aod[nzdf, ]
      if (is.null(aod) || ncol(aod) == 0) 
        break
      nc <- match(c("Cp", "AIC"), names(aod))
      nc <- nc[!is.na(nc)][1L]
      o <- order(aod[, nc])
      if (trace) {
        print(aod[o, ])
        utils::flush.console()
      }
      if (o[1L] == 1) 
        break
      change <- rownames(aod)[o[1L]]
    }
    usingCp <- match("Cp", names(aod), 0) > 0
    fit <- update(fit, paste("~ .", change), evaluate = FALSE)
    fit <- eval.parent(fit)
    nnew <- nobs(fit, use.fallback = TRUE)
    if (all(is.finite(c(n, nnew))) && nnew != n) 
      stop("number of rows in use has changed: remove missing values?")
    Terms <- terms(fit)
    bAIC <- extractAIC(fit, scale, k = k)#Was +...
    edf <- bAIC[1L]
    bAIC <- bAIC[2L]+extractAIC.getpenalisationprior(fit,prior)
    if (trace) {
      cat("\nStep:  AIC=", format(round(bAIC, 2)), "\n", 
          cut.string(deparse(formula(fit))), "\n\n", sep = "")
      utils::flush.console()
    }
    if (bAIC >= AIC + 1e-07) 
      break
    nm <- nm + 1
    models[[nm]] <- list(deviance = mydeviance(fit), df.resid = n - 
                           edf, change = change, AIC = bAIC)
    if (!is.null(keep)) 
      keep.list[[nm]] <- keep(fit, bAIC)
  }
  if (!is.null(keep)) 
    fit$keep <- re.arrange(keep.list[seq(nm)])
  step.results(models = models[seq(nm)], fit, object, usingCp)
}

extractAIC.getpenalisationprior<-function(object,prior){
  variableyn<-prior
  variableyn[]<-0
  variableyn[attr(terms(object),"term.labels")]<-1
  penalprior<-+2*sum(variableyn*log((1-prior)/0.5)+(1-variableyn)*log(prior/0.5))
  return(penalprior)
}

extractAIC.updateAICprior<-function(object,objaod,prior,drop=TRUE){
  pen<-extractAIC.getpenalisationprior(object,prior)
  objaod["<none>","AIC"]<-objaod["<none>","AIC"]+pen
  for(namep in rownames(objaod)[2:dim(objaod)[1]]){
    if(drop){
      objaod[namep,"AIC"]<-objaod[namep,"AIC"]+pen+2*log(prior[namep]/0.5)-2*log((1-prior[namep])/0.5)
    }else{
      objaod[namep,"AIC"]<-objaod[namep,"AIC"]+pen-2*log(prior[namep]/0.5)+2*log((1-prior[namep])/0.5)
    }
  }
  return(objaod)
}


###TESTING
# set.seed(1001)
# n=250
# C1<-rnorm(n,0,1)
# C2<-rnorm(n,0,1)
# C3<-rnorm(n,0,1)
# C4<-rnorm(n,0,1)
# Y<-C1+0.001*C2+rnorm(n,0,0.5)
# plot(C1,Y)
# plot(C2,Y)
# plot(C3,Y)
# plot(C4,Y)
# 
# library(MASS)
# 
# ####################
# ### TESTING BACKWARD
# ####################
# #Reference
# resAV<-stepAIC(lm(Y~C1+C2+C3+C4),trace=-1,direction="backward")
# summary(resAV)
# #Result select C1 and C2
# 
# #Identical to no prior
# priorall<-matrix(c(0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5),ncol=4,nrow=2,byrow = TRUE)
# rownames(priorall)<-c("p1","p2")
# colnames(priorall)<-c("C1","C2","C3","C4")
# prior<-priorall["p1",]
# resNOPRIOR<-stepAICmodified(lm(Y~C1+C2+C3+C4),prior=prior,trace=-1,direction="backward")
# summary(resNOPRIOR)
# #Result select C1 and C2
# 
# #In favor of C1, not C2
# priorall<-matrix(c(0.9,0.1,0.5,0.6,0.5,0.5,0.5,0.5),ncol=4,nrow=2,byrow = TRUE)
# rownames(priorall)<-c("p1","p2")
# colnames(priorall)<-c("C1","C2","C3","C4")
# prior<-priorall["p1",]
# resNEWSTRONGC1SMALLC2<-stepAICmodified(lm(Y~C1+C2+C3+C4),prior=prior,trace=-1,direction="backward")
# summary(resNEWSTRONGC1SMALLC2)
# #Result select C1 only
# 
# #In favor of C2 VERY VERY not C1
# priorall<-matrix(c(1e-100,0.9,0.5,0.6,0.5,0.5,0.5,0.5),ncol=4,nrow=2,byrow = TRUE)
# rownames(priorall)<-c("p1","p2")
# colnames(priorall)<-c("C1","C2","C3","C4")
# prior<-priorall["p1",]
# resNEWSTRONGC2SMALLC1<-stepAICmodified(lm(Y~C1+C2+C3+C4),prior=prior,trace=-1,direction="backward")
# summary(resNEWSTRONGC2SMALLC1)
# #Result select C2 and C4 (not C1)
# 
# 
# #In favor of VERY VERY C4, C1, not C2
# priorall<-matrix(c(0.9,0.1,1-1e-10,0.5,0.5,0.5,0.5,0.5),ncol=4,nrow=2,byrow = TRUE)
# rownames(priorall)<-c("p1","p2")
# colnames(priorall)<-c("C1","C2","C3","C4")
# prior<-priorall["p1",]
# resNEWSTRONGC4C1SMALLC2<-stepAICmodified(lm(Y~C1+C2+C3+C4),prior=prior,trace=-1,direction="backward")
# summary(resNEWSTRONGC4C1SMALLC2)
# #Result select C1 and C3
# 
# 
# ####################
# ### TESTING FORWARD
# ####################
# #Reference
# resAV<-stepAIC(lm(Y~C3),trace=-1,scope = list(upper = ~C1+C2+C3+C4),direction="forward")
# summary(resAV)
# #Result select C1, C2 and C3
# 
# #Identical to no prior
# priorall<-matrix(c(0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5),ncol=4,nrow=2,byrow = TRUE)
# rownames(priorall)<-c("p1","p2")
# colnames(priorall)<-c("C1","C2","C3","C4")
# prior<-priorall["p1",]
# resNOPRIOR<-stepAICmodified(lm(Y~C3),prior=prior,trace=-1,scope = list(upper = ~C1+C2+C3+C4),direction="forward")
# summary(resNOPRIOR)
# #Result select C1, C2 and C3
# 
# #In favor of C1, not C2
# priorall<-matrix(c(0.9,0.1,0.5,0.6,0.5,0.5,0.5,0.5),ncol=4,nrow=2,byrow = TRUE)
# rownames(priorall)<-c("p1","p2")
# colnames(priorall)<-c("C1","C2","C3","C4")
# prior<-priorall["p1",]
# resNEWSTRONGC1SMALLC2<-stepAICmodified(lm(Y~C3),prior=prior,trace=-1,scope = list(upper = ~C1+C2+C3+C4),direction="forward")
# summary(resNEWSTRONGC1SMALLC2)
# #Result select C1 and C3
# 
# #In favor of C2 VERY VERY not C1
# priorall<-matrix(c(1e-100,0.9,0.5,0.6,0.5,0.5,0.5,0.5),ncol=4,nrow=2,byrow = TRUE)
# rownames(priorall)<-c("p1","p2")
# colnames(priorall)<-c("C1","C2","C3","C4")
# prior<-priorall["p1",]
# resNEWSTRONGC2SMALLC1<-stepAICmodified(lm(Y~C3),prior=prior,trace=-1,scope = list(upper = ~C1+C2+C3+C4),direction="forward")
# summary(resNEWSTRONGC2SMALLC1)
# #Result select C2, C3 and C4
# 
# 
# #In favor of VERY VERY C4, C1, not C2
# priorall<-matrix(c(0.9,0.1,1-1e-10,0.5,0.5,0.5,0.5,0.5),ncol=4,nrow=2,byrow = TRUE)
# rownames(priorall)<-c("p1","p2")
# colnames(priorall)<-c("C1","C2","C3","C4")
# prior<-priorall["p1",]
# resNEWSTRONGC4C1SMALLC2<-stepAICmodified(lm(Y~C3),prior=prior,trace=-1,scope = list(upper = ~C1+C2+C3+C4),direction="forward")
# summary(resNEWSTRONGC4C1SMALLC2)
# #Result select C1 and C3
# 
# 
# ####################
# ### TESTING BOTH
# ####################
# #Reference
# resAV<-stepAIC(lm(Y~C3),trace=-1,scope = list(upper = ~C1+C2+C3+C4, lower = ~1),direction="both")
# summary(resAV)
# #Result select C1, C2 
# 
# #Identical to no prior
# priorall<-matrix(c(0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5),ncol=4,nrow=2,byrow = TRUE)
# rownames(priorall)<-c("p1","p2")
# colnames(priorall)<-c("C1","C2","C3","C4")
# prior<-priorall["p1",]
# resNOPRIOR<-stepAICmodified(lm(Y~C3),prior=prior,trace=-1,scope = list(upper = ~C1+C2+C3+C4, lower = ~1),direction="both")
# summary(resNOPRIOR)
# #Result select C1, C2 
# 
# #In favor of C1, not C2
# priorall<-matrix(c(0.9,0.1,0.5,0.6,0.5,0.5,0.5,0.5),ncol=4,nrow=2,byrow = TRUE)
# rownames(priorall)<-c("p1","p2")
# colnames(priorall)<-c("C1","C2","C3","C4")
# prior<-priorall["p1",]
# resNEWSTRONGC1SMALLC2<-stepAICmodified(lm(Y~C3),prior=prior,trace=-1,scope = list(upper = ~C1+C2+C3+C4, lower = ~1),direction="both")
# summary(resNEWSTRONGC1SMALLC2)
# #Result select C1 
# 
# #In favor of C2 VERY VERY not C1
# priorall<-matrix(c(1e-100,0.9,0.5,0.6,0.5,0.5,0.5,0.5),ncol=4,nrow=2,byrow = TRUE)
# rownames(priorall)<-c("p1","p2")
# colnames(priorall)<-c("C1","C2","C3","C4")
# prior<-priorall["p1",]
# resNEWSTRONGC2SMALLC1<-stepAICmodified(lm(Y~C3),prior=prior,trace=-1,scope = list(upper = ~C1+C2+C3+C4, lower = ~1),direction="both")
# summary(resNEWSTRONGC2SMALLC1)
# #Result select C2 and C4
# 
# 
# #In favor of VERY VERY C4, C1, not C2
# priorall<-matrix(c(0.9,0.1,1-1e-10,0.5,0.5,0.5,0.5,0.5),ncol=4,nrow=2,byrow = TRUE)
# rownames(priorall)<-c("p1","p2")
# colnames(priorall)<-c("C1","C2","C3","C4")
# prior<-priorall["p1",]
# resNEWSTRONGC4C1SMALLC2<-stepAICmodified(lm(Y~C3),prior=prior,trace=-1,scope = list(upper = ~C1+C2+C3+C4, lower = ~1),direction="both")
# summary(resNEWSTRONGC4C1SMALLC2)
# #Result select C1 and C3



