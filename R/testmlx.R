#' Statistical tests for model assessment
#'
#' Perform several statistical tests using the results of a Monolix run to assess the statistical
#' components of the model in use.
#' 
#' The tests used are:  1) F-tests (or, equivalently, correlation tests) to evaluate the effect of 
#' each covariate on each parameter ("covariate"), 2) Shapiro-Wilk and symmetry tests to assess 
#' the distribution of the random effects (""randomEffect"), 3) correlation tests to assess the 
#' correlation structure of the random effects ("correlation"), 4) Shapiro-Wilk and symmetry 
#' tests to assess the distribution of the residual errors ("residual").
#' 
#' By default, the four tests are performed
#' @param project a Monolix project
#' @param tests  a vector of strings: the list of tests to perform 
#' among c("covariate","randomEffect","correlation","residual")
#' @param plot  {FALSE}/TRUE  display some diagnostic plots associated to the tests (default=FALSE)
#' @return a list of data frames and ggplot objects if plot=TRUE
#' @examples
#' \dontrun{
#' r = testmlx(project="PBPKproject.mlxtran")
#' r = testmlx(project="PBPKproject.mlxtran", tests=c("covariate","correlation"), plot=TRUE)
#' }
#' @importFrom ggplot2 ggplot geom_point theme theme_set theme_bw aes geom_line xlab ylab facet_wrap facet_grid stat_ecdf aes_string
#'             geom_abline geom_boxplot geom_smooth
#' @importFrom gridExtra grid.arrange
#' @importFrom grDevices palette
#' @importFrom stats quantile anova binom.test cor cov dnorm lm logLik nlm p.adjust pchisq pnorm qnorm
#'             quantile rnorm runif sd shapiro.test spline t.test
#' @importFrom utils ls.str read.csv read.table write.table
#' @export

testmlx <- function(project, 
                    tests=c("covariate","randomEffect","correlation","residual"), 
                    plot=FALSE) 
{
  if(!file.exists(project)){
    message(paste0("ERROR: project '", project, "' does not exists"))
    return(invisible(FALSE))}
  
  loadProject(project)   
  
  launched.tasks <- getLaunchedTasks()
  if (!launched.tasks[["populationParameterEstimation"]]) {
    cat("\nEstimation of the population parameters... \n")
    runPopulationParameterEstimation()
  }
  if (!launched.tasks[["conditionalDistributionSampling"]]) {
    cat("Sampling of the conditional distribution... \n")
    runConditionalDistributionSampling()
  }
  
  #------------------------------
  if (!any(getIndividualParameterModel()$variability$id))
    stop("\nA least one parameter with random effects is required\n", call.=FALSE)
  
  if (plot)
    ggplot2::theme_set(theme_bw())
  
  res <- list()
  if ("covariate" %in% tests)
  res$covariate <- covariateTest(plot=plot)
  if ("residual" %in% tests)
    res$residual <- residualTest(plot=plot)
  if ("randomEffect" %in% tests)
    res$randomEffect <- randomEffectTest(plot=plot)
  if ("correlation" %in% tests)
    res$correlation <- correlationTest(plot=plot)
  return(res)
}

#----------------------------------------------------------
residualTest <- function(project=NULL, plot=FALSE) {
  
  if (!is.null(project)) 
    loadProject(project)
 
  residual <- getSimulatedResiduals()
  n.out <- length(residual)
  res.errorModel <- list()
  for (i.out in (1:n.out)) {
    resi <- residual[[i.out]]
    nrep <- max(resi$rep)
    resi <- matrix(resi$residual, ncol=nrep)
    ri1 <- swmlx.test(resi)
    ri2 <- vdwmlx.test(resi)
    res.errorModel[[i.out]] <- list(normality=ri1, symmetry=ri2)
  }
  
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
randomEffectTest <- function(project=NULL, plot=FALSE) {
  
  if (!is.null(project)) 
    loadProject(project)
  
  sim.randeff <- getSimulatedRandomEffects()
  if (is.null(sim.randeff$rep)) 
    sim.randeff$rep <- 1
  nrep <- max(sim.randeff$rep)
  
  var.param <- names(which(getIndividualParameterModel()$variability$id))
  var.randeff <- paste0("eta_",var.param)
  res.randeff <- NULL
  for (nj in var.randeff) {
    refi <- matrix(sim.randeff[[nj]],ncol=nrep)
    ri1 <- swmlx.test(refi)
    ri2 <- vdwmlx.test(refi)
    res.randeff <- rbind(res.randeff, 
                         data.frame(randomEffect=nj,normality=ri1, symmetry=ri2))
  }
  
  if (plot) {
    x <- seq(-3,3,length.out=100)
    dn <- data.frame(x,F=pnorm(x))
    pop.param <- getEstimatedPopulationParameters()
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
correlationTest <- function(project=NULL, plot=FALSE) {
  
  if (!is.null(project)) 
    loadProject(project)
  
  sim.randeff <- getSimulatedRandomEffects()
  if (is.null(sim.randeff$rep)) 
    sim.randeff$rep <- 1
  nrep <- max(sim.randeff$rep)
  
  var.param <- names(which(getIndividualParameterModel()$variability$id))
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
        if (is.null(res.cor))
          res.cor <- data.frame(randomEffect.1=nref[col.el[i1]], randomEffect.2=nref[col.el[i2]], correlation=ci, p.value=pv)
        else {
          resi <- data.frame(randomEffect.1=nref[col.el[i1]], randomEffect.2=nref[col.el[i2]], correlation=ci, p.value=pv)
          res.cor <- rbind(res.cor, resi)
        }
      }
    }
    
    if (plot) {
      pop.param <- getEstimatedPopulationParameters()
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
covariateTest <- function(project=NULL, plot=FALSE) {
  
  if (!is.null(project)) 
    loadProject(project)
  
  sim.randeff <- getSimulatedRandomEffects()
  if (is.null(sim.randeff$rep)) 
    sim.randeff$rep <- 1
  nrep <- max(sim.randeff$rep)
  
  var.param <- names(which(getIndividualParameterModel()$variability$id))
  var.randeff <- paste0("eta_",var.param)
  ind.dist <- getIndividualParameterModel()$distribution[var.param]
  n.param <- length(var.param)
  m.indparam <- getEstimatedIndividualParameters()$conditionalMean[var.param]
  m.randeff <- getEstimatedRandomEffects()$conditionalMean[var.randeff]
  
  cov.info <- getCovariateInformation()
  cov.names <- cov.info$name
  cov.types <- cov.info$type
  j.cat <- grep("categorical",cov.types)
  j.cont <- grep("continuous",cov.types)
  cont.cov <- cov.names[j.cont]
  cat.cov <- cov.names[j.cat]
  covariates <- cov.info$covariate
  covariates[cat.cov] <-  lapply(covariates[cat.cov],as.factor)
  #covariates["id"] <- NULL
  
  d1 <- NULL
  g=getIndividualParameterModel()$covariateModel
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
    names(m.indparam)[names(m.indparam)==nj] <- n.yj
    lm0 <- lm(yj ~1)
    for (nc in cov.names) {
      lmc <- lm(yj ~ covariates[[nc]])
      pjc <- signif(anova(lm0, lmc)$`Pr(>F)`[2],4)
      dnc <- data.frame(parameter=n.yj,covariate=nc,p.value=pjc,in.model=g[[nj]][[nc]])
      d1 <- rbind(d1,dnc)
    }
  }
  d1 <- d1[order(d1$in.model, decreasing=TRUE),]
  
  d2 <- NULL
  for (ne in var.randeff) {
    nj <- gsub("eta_","",ne)
    yj <- m.randeff[[ne]]
    lm0 <- lm(yj ~1)
    for (nc in cov.names) {
      lmc <- lm(yj ~ covariates[[nc]])
      pjc <- signif(anova(lm0, lmc)$`Pr(>F)`[2],4)
      dnc <- data.frame(random.effect=ne,covariate=nc,p.value=pjc,in.model=g[[nj]][[nc]])
      d2 <- rbind(d2,dnc)
    }
  }
  d2 <- d2[order(d2$in.model, decreasing=TRUE),]
  
  if (plot) {
    sim.param <- getSimulatedIndividualParameters()
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
vdwmlx.test <- function(x) {
  if (is.matrix(x))
    x <- apply(x, MARGIN = 1, mean)
  i0 <- which(x==0)
  x[i0] <- rnorm(length(i0))*sd(x)*0.000001
  mx <- 0
  x1 <- x[x>mx] - mx
  x2 <- mx-x[x<mx]
  n1 <- length(x1)
  n2 <- length(x2)
  n <- length(x)
  
  r <- rank(c(x1,x2))  
  r1 <- r[1:n1]        
  r2 <- r[(n1+1):n]    
  
  m1 <- mean(qnorm(r1/(n+1)))
  m2 <- mean(qnorm(r2/(n+1)))
  s2 <- sum(qnorm((1:n)/(n+1))^2)/(n-1)
  T <- (n1*m1^2 + n2*m2^2)/s2  # test statistic ~ Chi2(1)
  p1 <- 1 - pchisq(T,1)         # p-value
  p2 <- binom.test(sum(x>0),length(x))$p.value
  p <- min(1,2*min(p1, p2))
  
  return(p)
}

swmlx.test <- function(x) {
  if (is.vector(x)) {
    p <- shapiro.test(x)$p.value
  } else {
    K <- ncol(x)
    #   cu <- cov(x)
    #   du <- mean(diag(cu))
    #   cu <- cu + diag(mean(cu)-diag(cu))
    #   mu <- mean(cu)
    #   cu <- matrix(mu,nrow=K,ncol=K)
    #   cu <- cu + diag(du-diag(cu))
    #   v=x%*%solve(chol(cu))
    #   p <- shapiro.test(v)$p.value
    # 
    p <- NULL
    for (k in (1:K))
      p <- c(p, shapiro.test(x[,k])$p.value)
    p <- p.adjust(sort(p), method="BH")[1]
    p <- shapiro.test(rowMeans(x))$p.value
  }
  return(p)
}


