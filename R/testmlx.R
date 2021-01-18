#' Statistical tests for model assessment
#'
#' Perform several statistical tests using the results of a Monolix run to assess the statistical
#' components of the model in use.
#' 
#' The tests used are:  
#' 1) F-tests (or, equivalently, correlation tests) to evaluate the effect of each covariate 
#' on each parameter ("covariate"), 
#' 2) correlation tests to assess the correlation structure of the random effects ("correlation"), 
#' 3) Shapiro-Wilk and Miao-Gel-Gastwirth tests to assess, respectively the normality and the
#' symmetry of the distribution of the random effects (""randomEffect"), 
#' 4) Shapiro-Wilk and Miao-Gel-Gastwirth tests to assess, respectively the normality and the
#' symmetry of the distribution of residual errors ("residual").
#' 
#' By default, the four tests are performed.
#' 
#' When several samples of the conditional distributions are used, two methods are proposed in order 
#' to take into the dependance of the samples for the Shapiro-Wilk and Miao-Gel-Gastwirth tests: 
#' "edf" computes an effective degrees of freedom, "BH" performs one test per replicates and adjust 
#' the smallest p-value using the Benjamini-Hochberg correction.
#' @param project a Monolix project
#' @param tests  a vector of strings: the list of tests to perform 
#' among c("covariate","randomEffect","correlation","residual")
#' @param plot  {FALSE}/TRUE  display some diagnostic plots associated to the tests (default=FALSE)
#' @param adjust method to take into account the dependency of MCMC sample  c({"edf"},"BH")
#' @param n.sample number of samples from the conditional distribution to be used (default = number of available samples in the project)
#' @return a list of data frames and ggplot objects if plot=TRUE
#' @examples
#' # RsmlxDemo2.mlxtran is a Monolix project for modelling the PK of warfarin using a PK model 
#' # with parameters ka, V, Cl.
#' 
#' #testmlx will perform statistical tests for the different component of the statistical model:
#' r1 <- testmlx(project="RsmlxDemo2.mlxtran")
#' 
#' #testmlx will perform statistical tests for the covariate model and the correlation model only.
#' r2 <- testmlx(project="RsmlxDemo2.mlxtran", tests=c("covariate","correlation"))
#' 
#' # See http://rsmlx.webpopix.org/userguide/testmlx/ for detailed examples of use of testmlx
#' # Download the demo examples here: http://rsmlx.webpopix.org/installation
#' 
#' 
#' @importFrom ggplot2 ggplot geom_point theme theme_set theme_bw aes geom_line xlab ylab facet_wrap facet_grid stat_ecdf aes_string
#'             geom_abline geom_boxplot geom_smooth
#' @importFrom gridExtra grid.arrange
#' @importFrom grDevices palette
#' @importFrom stats quantile anova binom.test cor cov dnorm lm logLik nlm p.adjust pchisq dchisq pnorm qnorm
#'             quantile rnorm runif sd shapiro.test spline t.test optimize median
#' @importFrom utils ls.str read.csv read.table write.table
#' @export

testmlx <- function(project, 
                    tests=c("covariate","randomEffect","correlation","residual"), 
                    plot=FALSE, adjust="edf", n.sample=NULL) 
{
  RsmlxDemo1.project <- RsmlxDemo2.project <- warfarin.data  <- resMonolix <- NULL
  
  r <- prcheck(project, f="test", tests=tests)
  if (r$demo)
    return(r$res)
  project <- r$project
  
  launched.tasks <- mlx.getLaunchedTasks()
  if (!launched.tasks[["populationParameterEstimation"]]) {
    cat("\nEstimation of the population parameters... \n")
    mlx.runPopulationParameterEstimation()
  }
  if (!launched.tasks[["conditionalDistributionSampling"]]) {
    cat("Sampling of the conditional distribution... \n")
    mlx.runConditionalDistributionSampling()
  }
  
  theme_set(theme_bw())
  
  #------------------------------
  if (!any(mlx.getIndividualParameterModel()$variability$id))
    stop("\nA least one parameter with random effects is required\n", call.=FALSE)
  
  method.adjust <- adjust
  res <- list()
  if ("covariate" %in% tests)
    res$covariate <- covariateTest(plot=plot, n.sample=n.sample)
  if ("residual" %in% tests)
    res$residual <- residualTest(plot=plot, method.adjust=method.adjust, n.sample=n.sample)
  if ("randomEffect" %in% tests)
    res$randomEffect <- randomEffectTest(plot=plot, method.adjust = method.adjust, n.sample=n.sample)
  if ("correlation" %in% tests)
    res$correlation <- correlationTest(plot=plot, n.sample=n.sample)
  return(res)
}

#----------------------------------------------------------
residualTest <- function(project=NULL, method.adjust="edf", n.sample=NULL, plot=FALSE) {
  
  if (!is.null(project)) 
    mlx.loadProject(project)
  
  if (is.null(mlx.getContinuousObservationModel()))  return(list())
  
  if (!is.null(n.sample) && n.sample=="mode") {
    residual <- getEstimatedResiduals()
    for (i.out in (1:length(residual))) {
      if (is.null(residual[[i.out]]$conditionalMode))
        stop("the conditional modes (EBE's) have not been estimated...", call.=FALSE)
      residual[[i.out]] <- data.frame(rep=1, residual=residual[[i.out]]$conditionalMode)
    }
    n.sample <- 1
  } else if (!is.null(n.sample) && n.sample=="mean") {
    residual <- getEstimatedResiduals()
    for (i.out in (1:length(residual))) {
      if (is.null(residual[[i.out]]$conditionalMean))
        stop("the conditional means have not been estimated...", call.=FALSE)
      residual[[i.out]] <- data.frame(rep=1, residual=residual[[i.out]]$conditionalMean)
    }
    n.sample <- 1
  } else {
    residual <- getSimulatedResiduals()
    for (i.out in (1:length(residual))) 
      nrep <- max(residual[[i.out]]$rep)
    if (is.null(n.sample))  n.sample <- nrep
    if (n.sample<1)
      stop("the number of replicates n.sample should be greater than or equal to 1", call. = FALSE)
    if (n.sample>nrep)
      stop(paste0("the number of replicates should be less than or equal to  ", nrep), call. = FALSE)
    
    residual[[i.out]] <- subset(residual[[i.out]], rep<=n.sample)
  }
  n.out <- length(residual)
  res.errorModel <- NULL
  for (i.out in (1:n.out)) {
    resi <- residual[[i.out]]
    nrep <- max(resi$rep)
    resi <- matrix(resi$residual, ncol=nrep)
    # adf1 <- adjust.df(resi)
    # ndf1 <- nrow(resi)*adf1
    ri1 <- sw.test(resi, method=method.adjust)
    ri2 <- mgg.test(resi, method=method.adjust)
    res.errorModel <- rbind(res.errorModel, c(normality=ri1, symmetry=ri2))
    #res.errorModel[[i.out]] <- list(normality=ri1, symmetry=ri2)
  }
  res.errorModel <- as.data.frame(res.errorModel)
  rownames(res.errorModel) <- names(residual)
  
  if (plot) {
    x <- seq(-3,3,length.out=100)
    dn <- data.frame(x,F=pnorm(x))
    pl <- list()
    for (i.out in (1:n.out)) {
      ri <- residual[[i.out]]$residual
      ni <- paste0("residual_",names(residual)[i.out])
      pl[[i.out]] <- ggplot() + stat_ecdf(aes_string(ri), geom = "step", size=1) + 
        geom_line(data=dn, aes_string(x,"F"), colour="red", size=1) + ylab("F") + 
        xlab(ni)
    }
    do.call(grid.arrange,c(pl,ncol=2))
    return(list(p.value=res.errorModel, plot=pl))
  } else
    return(list(p.value=res.errorModel))
}

#----------------------------------------------------------
randomEffectTest <- function(project=NULL, method.adjust="edf", n.sample=NULL, plot=FALSE) {
  
  if (!is.null(project)) 
    mlx.loadProject(project)
  
  if (!is.null(n.sample) && n.sample=="mode") {
    sim.randeff <- mlx.getEstimatedRandomEffects()$conditionalMode
    if (is.null(sim.randeff))
      stop("the conditional modes (EBE's) have not been estimated...", call.=FALSE)
    n.sample <- 1
  } else if (!is.null(n.sample) && n.sample=="mean") {
    sim.randeff <- mlx.getEstimatedRandomEffects()$conditionalMean
    if (is.null(sim.randeff))
      stop("the conditional means have not been estimated...", call.=FALSE)
    n.sample <- 1
  } else {
    sim.randeff <- mlx.getSimulatedRandomEffects()
  }
  if (is.null(sim.randeff$rep)) 
    sim.randeff$rep <- 1
  nrep <- max(sim.randeff$rep)
  if (is.null(n.sample))  n.sample <- nrep
  
  if (n.sample<1)
    stop("the number of replicates n.sample should be greater than or equal to 1", call. = FALSE)
  if (n.sample>nrep)
    stop(paste0("the number of replicates should be less than or equal to  ", nrep), call. = FALSE)
  sim.randeff <- subset(sim.randeff, rep<=n.sample)
  nrep <- n.sample
  
  pop <- mlx.getEstimatedPopulationParameters()
  
  var.param <- names(which(mlx.getIndividualParameterModel()$variability$id))
  var.randeff <- paste0("eta_",var.param)
  omega <- paste0("omega_",var.param)
  
  res.randeff <- NULL
  k <- 0
  for (nj in var.randeff) {
    k <- k+1
    refi <- matrix(sim.randeff[[nj]], ncol=nrep)/pop[omega[k]]
    
    # if (method.adjust=="edf") {
    #   adf1 <- adjust.df(refi)
    #   ndf1 <- nrow(refi)*adf1
    # } else {
    #   adf1 <- ndf1 <- NULL
    # }
    
    ri1 <- sw.test(refi, method=method.adjust)
    ri2 <- mgg.test(refi, method=method.adjust)
    res.randeff <- rbind(res.randeff, c(normality=ri1, symmetry=ri2))
  }
  res.randeff <- as.data.frame(res.randeff)
  rownames(res.randeff) <- var.randeff
  
  if (plot) {
    x <- seq(-3,3,length.out=100)
    dn <- data.frame(x,F=pnorm(x))
    pop.param <- mlx.getEstimatedPopulationParameters()
    names.pop.param <- gsub(" ","",names(pop.param))
    ind.param.omega <- grep("omega_",names.pop.param)
    iop.omega2 <- FALSE
    if (length(ind.param.omega)==0) {
      ind.param.omega <- grep("omega2_",names.pop.param)
      iop.omega2 <- TRUE
    }
    
    pl <- list()
    j <- 0
    for (nj in var.randeff) {
      j <- j+1
      nj0 <- grep(gsub("eta_","",nj),names.pop.param[ind.param.omega])
      if (iop.omega2)
        oj <- sqrt(pop.param[ind.param.omega[nj0]])
      else
        oj <- pop.param[ind.param.omega[nj0]]
      rj <- sim.randeff[[nj]]/oj
      
      pl[[j]] <- ggplot() + stat_ecdf(aes_string(rj), geom = "step", size=1) + 
        geom_line(data=dn, aes_string("x","F"), colour="red", size=1) + 
        ylab("F") + xlab(nj)
    }
    do.call(grid.arrange,c(pl,ncol=2))
    return(list(p.value=res.randeff, plot=pl))
  } else 
    return(list(p.value=res.randeff))
}

#----------------------------------------------------------
correlationTest <- function(project=NULL, n.sample=NULL, plot=FALSE) {
  
  if (!is.null(project)) 
    mlx.loadProject(project)
  
  if (!is.null(n.sample) && n.sample=="mode") {
    sim.randeff <- mlx.getEstimatedRandomEffects()$conditionalMode
    if (is.null(sim.randeff))
      stop("the conditional modes (EBE's) have not been estimated...", call.=FALSE)
    n.sample <- 1
  } else if (!is.null(n.sample) && n.sample=="mean") {
    sim.randeff <- mlx.getEstimatedRandomEffects()$conditionalMean
    if (is.null(sim.randeff))
      stop("the conditional means have not been estimated...", call.=FALSE)
    n.sample <- 1
  } else {
    sim.randeff <- mlx.getSimulatedRandomEffects()
  }
  if (is.null(sim.randeff$rep)) 
    sim.randeff$rep <- 1
  nrep <- max(sim.randeff$rep)
  if (is.null(n.sample))  n.sample <- nrep
  
  if (n.sample<1)
    stop("the number of replicates n.sample should be greater than or equal to 1", call. = FALSE)
  if (n.sample>nrep)
    stop(paste0("the number of replicates should be less than or equal to  ", nrep), call. = FALSE)
  sim.randeff <- subset(sim.randeff, rep<=n.sample)
  nrep <- n.sample
  
  
  var.param <- names(which(mlx.getIndividualParameterModel()$variability$id))
  var.randeff <- paste0("eta_",var.param)
  col.el <- which((names(sim.randeff) %in% var.randeff))
  nel <- length(col.el)
  if (nel>1) {
    nref <- names(sim.randeff)
    res.cor <- NULL
    for (i1 in (1:(nel-1))) {
      for (i2 in ((i1+1):nel)) {
        refi1 <- matrix(sim.randeff[,col.el[i1]],ncol=nrep)
        refi2 <- matrix(sim.randeff[,col.el[i2]],ncol=nrep)
        ri <- rowSums(refi1*refi2)
        ci <- cor(matrix(refi1,ncol=1),matrix(refi2,ncol=1))
        pv <- signif(t.test(ri)$p.value, 4)
        # pjc <- signif(anova(lm0, lmc)$`Pr(>F)`[2],4)
        # plrt <- signif(1-pchisq(2*c(logLik(lmc)-logLik(lm0)),1),4)
        # dnc <- data.frame(random.effect=ne,covariate=nc,p.value=pjc,p.ttest=pjc,p.lrt=plrt,in.model=g[[nj]][[nc]])
        if (is.null(res.cor))
          res.cor <- data.frame(randomEffect.1=nref[col.el[i1]], randomEffect.2=nref[col.el[i2]], correlation=ci, p.value=pv)
        else {
          resi <- data.frame(randomEffect.1=nref[col.el[i1]], randomEffect.2=nref[col.el[i2]], correlation=ci, p.value=pv)
          res.cor <- rbind(res.cor, resi)
        }
      }
    }
    
    if (plot) {
      pop.param <- mlx.getEstimatedPopulationParameters()
      names.pop.param <- gsub(" ","",names(pop.param))
      ind.param.omega <- grep("omega_",names.pop.param)
      iop.omega2 <- FALSE
      if (length(ind.param.omega)==0) {
        ind.param.omega <- grep("omega2_",names.pop.param)
        iop.omega2 <- TRUE
      }
      for (i in col.el) {
        ni <- names(sim.randeff)[i]
        ni0 <- grep(gsub("eta_","",ni),names.pop.param[ind.param.omega])
        if (iop.omega2)
          oi <- sqrt(pop.param[ind.param.omega[ni0]])
        else
          oi <- pop.param[ind.param.omega[ni0]]
        sim.randeff[,i] <- sim.randeff[,i]/oi
      }
      nel <- length(col.el)
      pl <- list()
      j <- 0
      
      for (i1 in (1:(nel-1))) {
        for (i2 in ((i1+1):nel)) {
          j <- j+1
          ni <- names(sim.randeff)[c(col.el[i1],col.el[i2])]
          mi <- as.matrix(sim.randeff[ni])
          vi <- eigen(t(mi)%*%mi)$vectors
          slopei <- vi[2]/vi[1]
          yb2 <- mean(sim.randeff[[ni[2]]])
          yb1 <- mean(sim.randeff[[ni[1]]])
          ini <- yb2 - slopei*yb1
          
          ci <- cov(mi)
          a1 <- ci[2]/ci[1]
          b1 <- yb2 - a1*yb1
          a2 <- ci[4]/ci[2]
          b2 <- yb2 - a2*yb1
          
          pl[[j]] <- ggplot(sim.randeff, aes_string(ni[1],ni[2])) + geom_point(color="blue") +
            #           geom_abline(intercept=ini, slope=slopei, color="green") + 
            geom_abline(intercept=b1, slope=a1, color="red") + 
            geom_abline(intercept=b2, slope=a2, color="cyan") + 
            ylab(ni[2]) + xlab(ni[1])
        }
      }
      do.call(grid.arrange,c(pl,ncol=ceiling(sqrt(length(pl)))))
      return(list(p.value=res.cor, plot=pl))
    } else 
      return(list(p.value=res.cor))
  }
}

#----------------------------------------------------------
covariateTest <- function(project=NULL, n.sample=NULL, plot=FALSE) {
  
  if (!is.null(project)) 
    mlx.loadProject(project)
  
  cov.info <- mlx.getCovariateInformation()
  if (is.null(cov.info)) return(list())
  cov.names <- cov.info$name
  cov.types <- cov.info$type
  j.cat <- grep("categorical",cov.types)
  j.cont <- grep("continuous",cov.types)
  cont.cov <- cov.names[j.cont]
  cat.cov <- cov.names[j.cat]
  covariates <- cov.info$covariate
  covariates[cat.cov] <-  lapply(covariates[cat.cov],as.factor)
  #covariates["id"] <- NULL
  covariates <- covariates[order(covariates$id),]
  
  # sim.randeff <- mlx.getSimulatedRandomEffects()
  # if (is.null(sim.randeff$rep)) 
  #   sim.randeff$rep <- 1
  # nrep <- max(sim.randeff$rep)
  
  var.param <- names(which(mlx.getIndividualParameterModel()$variability$id))
  var.randeff <- paste0("eta_",var.param)
  ind.dist <- mlx.getIndividualParameterModel()$distribution[var.param]
  n.param <- length(var.param)
  # m.indparam <- mlx.getEstimatedIndividualParameters()$conditionalMean[var.param]
  # m.randeff <- mlx.getEstimatedRandomEffects()$conditionalMean[var.randeff]
  
  if (!is.null(n.sample) && n.sample=="mode") {
    m.indparam <- mlx.getEstimatedIndividualParameters()$conditionalMode[c("id",var.param)]
    m.randeff <- mlx.getEstimatedRandomEffects()$conditionalMode[c("id",var.randeff)]
    if (is.null(m.indparam))
      stop("the conditional modes (EBE's) have not been estimated...", call.=FALSE)
    n.sample <- 1
  } else if (!is.null(n.sample) && n.sample=="mean") {
    m.indparam <- mlx.getEstimatedIndividualParameters()$conditionalMean[c("id",var.param)]
    m.randeff <- mlx.getEstimatedRandomEffects()$conditionalMean[c("id",var.randeff)]
    if (is.null(m.indparam))
      stop("the conditional means have not been estimated...", call.=FALSE)
    n.sample <- 1
  } else {
    m.indparam <- mlx.getSimulatedIndividualParameters()
    m.randeff <- mlx.getSimulatedRandomEffects()
  }
  if (is.null(m.indparam$rep)) {
    m.indparam$rep <- 1
    m.randeff$rep <- 1
  }
  nrep <- max(m.indparam$rep)
  if (is.null(n.sample))  n.sample <- nrep
  
  if (n.sample<1)
    stop("the number of replicates n.sample should be greater than or equal to 1", call. = FALSE)
  if (n.sample>nrep)
    stop(paste0("the number of replicates should be less than or equal to  ", nrep), call. = FALSE)
  m.indparam <- subset(m.indparam, rep<=n.sample)[c("id",var.param)]
  m.randeff <- subset(m.randeff, rep<=n.sample)[c("id",var.randeff)]
  
  g=mlx.getIndividualParameterModel()$covariateModel
  lnj <- NULL
  for (nj in var.param) {
    dj <- tolower(ind.dist[nj])
    yj <- m.indparam[[nj]]
    n.yj <- nj
    if (dj == "lognormal") {
      yj <- log(yj)
      n.yj <- paste0("log.",nj)
    } else if (dj == "logitnormal") {
      yj <- log(yj/(1-yj))
      n.yj <- paste0("logit.",nj)
    } else if (dj == "probitnormal") {
      yj <- qnorm(yj)
      n.yj <- paste0("probit.",nj)
    } 
    m.indparam[[nj]] <- yj
    lnj <- c(lnj, n.yj)
  }
  m.indparam <- aggregate(m.indparam[var.param],list(m.indparam$id),mean)
  m.indparam <- m.indparam[order(m.indparam[,1]),]
  
  d1 <- NULL  
  j <- 0
  for (nj in var.param) {
    j <- j+1
    yj <- m.indparam[[nj]]
    lm0 <- lm(yj ~1)
    for (nc in cov.names) {
      lmc <- lm(yj ~ covariates[[nc]])
      pjc <- signif(anova(lm0, lmc)$`Pr(>F)`[2],4)
      plrt <- signif(1-pchisq(2*c(logLik(lmc)-logLik(lm0)),1),4)
      dnc <- data.frame(parameter=lnj[j],covariate=nc,p.value=pjc,p.ttest=pjc,p.lrt=plrt,in.model=g[[nj]][[nc]])
      d1 <- rbind(d1,dnc)
    }
  }
  foo <- d1[order(d1$parameter, d1$p.value),]
  d1 <- foo[order(foo$in.model, decreasing=TRUE),]
  
  m.randeff <- aggregate(m.randeff[var.randeff],list(m.randeff$id),mean)
  m.randeff <- m.randeff[order(m.randeff[,1]),]
  
  d2 <- NULL
  for (ne in var.randeff) {
    nj <- gsub("eta_","",ne)
    yj <- m.randeff[[ne]]
    lm0 <- lm(yj ~1)
    for (nc in cov.names) {
      lmc <- lm(yj ~ covariates[[nc]])
      pjc <- signif(anova(lm0, lmc)$`Pr(>F)`[2],4)
      plrt <- signif(1-pchisq(2*c(logLik(lmc)-logLik(lm0)),1),4)
      dnc <- data.frame(random.effect=ne,covariate=nc,p.value=pjc,p.ttest=pjc,p.lrt=plrt,in.model=g[[nj]][[nc]])
      d2 <- rbind(d2,dnc)
    }
  }
  d2 <- d2[order(d2$in.model, decreasing=TRUE),]
  
  if (plot) {
    sim.param <- mlx.getSimulatedIndividualParameters()
    sim.param$rep <- NULL
    tr.param <- NULL
    for (nj in var.param) {
      dj <- tolower(ind.dist[nj])
      yj <- sim.param[[nj]]
      n.yj <- nj
      if (dj == "lognormal") {
        yj <- log(yj)
        n.yj <- paste0("log.",nj)
      } else if (dj == "logitnormal") {
        yj <- log(yj/(1-yj))
        n.yj <- paste0("logit.",nj)
      } else if (dj == "probitnormal") {
        yj <- qnorm(yj)
        n.yj <- paste0("probit.",nj)
      } 
      sim.param[[nj]] <- yj
      names(sim.param)[names(sim.param)==nj] <- n.yj
      tr.param <- c(tr.param, n.yj)
    }
    r=merge(sim.param,covariates)
    pl <- list()
    j <- 0
    for (nj in tr.param) {
      yj <- m.indparam[[nj]]
      for (nc in cov.names) {
        j <- j+1
        if (nc %in% cont.cov) {
          pl[[j]] <- ggplot(r, aes_string(nc,nj)) + geom_point(color="blue") +
            geom_smooth(method="lm", color="red", se=FALSE) + xlab(nc) + ylab(nj)
        } else {
          pl[[j]] <- ggplot(r, aes_string(nc,nj)) + geom_boxplot() 
        }
      }
    }
    do.call(grid.arrange,c(pl,nrow=length(tr.param)))
    return(list(p.value.parameters=d1, p.value.randomEffects=d2, plot=pl))
  } else 
    return(list(p.value.parameters=d1, p.value.randomEffects=d2))
}

#----------------------------------------------------------
mgg.u <- function(x, en=NULL) {
  if (is.null(en))
    en <- length(x)
  dmx <- mean(x)-median(x)
  T <- sqrt(en)/0.9468922*dmx/mean(abs(x-median(x)))
  pT <- pnorm(T)
  return(2*min(pT, 1-pT))
}

sw.u <- function(x, en=NULL) {
  n <- length(x)
  if (is.null(en))
    en <- n
  if (n <4)
    stop("n should be greater than 3")
  mt <- qnorm(((1:n)- 3/8)/(n+0.25))
  c  <- 1/sqrt(sum(mt^2))*mt
  xn <- 1/sqrt(n)
  at <- vector(length=n)
  at[n] <- c[n] + 0.221157*xn - 0.147981*xn^2 - 2.071190*xn^3 + 4.434685*xn^4 - 2.706056*xn^5
  at[n-1] <- c[n-1] + 0.042981*xn - 0.293762*xn^2 - 1.752461*xn^3 + 5.682633*xn^4 - 3.582663*xn^5
  
  if (n<=5) {
    phi <- (sum(mt^2) - 2*(mt[n]^2))/(1 - 2*(at[n]^2))
    at[2:(n-1)] <- mt[2:(n-1)]/sqrt(phi)
    at[1] <- -at[n]
  } else {
    phi <- (sum(mt^2) - 2*(mt[n]^2) - 2*(mt[n-1]^2))/(1 - 2*(at[n]^2) - 2*(at[n-1]^2))
    at[3:(n-2)] <- mt[3:(n-2)]/sqrt(phi)
    at[1:2] <- -at[c(n, n-1)]
  }
  xs <- sort(x)
  W <- sum(at*xs)^2/sum((x-mean(x))^2)
  
  if (en <= 11) {
    mu <-  0.5440 - 0.39978*en + 0.025054*en^2 - 0.0006714*en^3
    sigma <- exp( 1.3822 - 0.77857*en + 0.062767*en^2 - 0.0020322*en^3 ) 
    gamma <- -2.273 + 0.459*en
    w <- -log(gamma - log(1-W))
  } else {
    ln <- log(en)
    mu <- - 1.5861 - 0.31082*ln - 0.083751*(ln^2) + 0.0038915*(ln^3)
    sigma <- exp( -0.4803 - 0.082676*ln + 0.0030302*(ln^2)) 
    w <- log(1-W)
  }
  return(1-pnorm(w,mu,sigma))
}


mgg.test <- function (x,  method="edf") {
  my.test(x=x, method=method, FUN=mgg.u)
}

sw.test <- function (x,  method="edf") {
  my.test(x=x, method=method, FUN=sw.u)
}

my.test <- function (x,  method="edf", FUN=NULL) {
  if (is.vector(x)) 
    x <- matrix(x,ncol=1)
  x <- (x -mean(x))/sd(x)
  K <- ncol(x)
  if (K==1)
    return(FUN(x))
  
  if (method=="ortho") {
    rho <- mean((rowSums(x)^2-rowSums(x^2)))/(K*(K-1))
    Q <- matrix(rho,K,K)  + diag(rep(1-rho, K))
    z <- x%*%solve(chol(Q))
    return(FUN(z))
  }
  if (method=="cov") {
    z <- x%*%solve(chol(cov(x)))
    return(FUN(z))
  }
  if (method=="BH") {
    p <- NULL
    for (k in (1:K)) 
      p <- c(p, FUN(x[,k]))
    p <- p.adjust(sort(p), method="BH")[1]
    return(p)
  }
  if (method=="edf") {
    en <- nrow(x)*adjust.df(x)
    return(FUN(x, en))
  }
} 


# ----------------------------------------
adjust.df <- function(e) {
  e1.df <- function(x,y,K) {
    n <- length(y)
    M <- K/x
    ofv <- n*log(M) - sum(dchisq(x=y/M, df=x, log=TRUE))
    return(ofv)
  }
  K <- ncol(e)
  if (K==1) {
    nu <- 1
  } else {
    s <- rowSums(e^2)
    s <- s[s>0]
    nu <- optimize(e1.df, c(1,K), y=s, K=K)$minimum
  }
  return(nu)
}



