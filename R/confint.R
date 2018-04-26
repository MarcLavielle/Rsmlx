confint.mlx <- function(project=NULL, level=0.95, linearization=FALSE)
{
  if (!is.null(project))
    loadProject(project)
  
  if (level<=0 | level>=1)
    stop("Level of the confidence interval should be strictly between 0 and 1", call.=FALSE)
  
  launched.tasks <- getLaunchedTasks()
  if (!linearization) {
    if (!("stochasticApproximation" %in% launched.tasks[["standardErrorEstimation"]]) ) {
      runStandardErrorEstimation(linearization=FALSE)
      cat("Estimation of the standard errors ... \n")
    }    
    se <- as.numeric(unlist(getEstimatedStandardErrors()$stochasticApproximation))
    names(se) <- names(getEstimatedStandardErrors()$stochasticApproximation)
  } else {
    if (!("linearization" %in% launched.tasks[["standardErrorEstimation"]])) {
      cat("Estimation of the standard errors ... \n")
      runStandardErrorEstimation(linearization=TRUE)
    }
    se <- getEstimatedStandardErrors()$linearization
  }
  
  param <- getEstimatedPopulationParameters()
  pname <- names(param)
  
  io <- match(pname, names(se))
  se <- se[io]
  
  indParam.name <- getIndividualParameterModel()$name
  indParam.dist <- getIndividualParameterModel()$distribution
  p2.name <- sub("_pop","",pname)
  i.pop <- match(indParam.name,p2.name)
  i.log <- i.pop[grep("logNormal", indParam.dist)]
  i.logit <- i.pop[grep("logitNormal", indParam.dist)]
  i.probit <- i.pop[grep("probitNormal", indParam.dist)]
  
  i.omega <- c(grep("^omega_",pname),grep("^omega2_",pname))
  i.corr <- grep("^corr_",pname)
  
  tr <- rep("N", length(param))
  tr[i.log] <- "L"
  tr[i.logit] <- "G"
  tr[i.probit] <- "P"
  tr[i.omega] <- "L"
  tr[i.corr] <- "R"
  
  set <- se
  mut <- param
  iL <- which(tr=="L")
  set[iL] <- se[iL]/param[iL]
  mut[iL] <- log(param[iL])
  iG <- which(tr=="G")
  set[iG] <- se[iG]/(param[iG]*(1-param[iG]))
  mut[iG] <- log(param[iG]/(1-param[iG]))
  iR <- which(tr=="R")
  set[iR] <- se[iR]*2/(1 - param[iR]^2)
  mut[iR] <- log((param[iR]+1)/(1-param[iR]))
  iP <- which(tr=="P")
  set[iP] <- se[iP]/dnorm(qnorm(param[iP]))
  mut[iP] <- qnorm(param[iP])
  
  cit <- cbind(mut + qnorm((1-level)/2)*set, mut + qnorm((1+level)/2)*set)
  ci <- cit
  ci[iL,] <- exp(ci[iL,])
  ci[iG,] <- 1/(1+exp(-ci[iG,]))
  ci[iP,] <- pnorm(ci[iP,])
  ci[iR,] <- (1-exp(-ci[iR,]))/(1+exp(-ci[iR,]))
  ci <- as.data.frame(ci)
  names(ci) <- c("lower", "upper")
  return(list(confint=ci, level=level, estimates=param))
}

