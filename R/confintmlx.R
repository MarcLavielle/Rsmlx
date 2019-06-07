#' Confidence intervals for population parameters 
#'
#' Compute confidence intervals for the population parameters estimated by Monolix.
#' 
#' The method used for computing the confidence intervals can be either based on the 
#' standard errors derived from an estimation of the Fisher Information Matrix ("fim"),
#' on the profile likelihood ("proflike") or on nonparametric bootstrap estimate ("bootstrap").
#' \code{method="fim"} is used by default.
#' 
#' When method="fim", the FIM can be either estimated using a linearization of the model 
#' or a stochastic approximation. When method="proflike", the observed likelihood can be 
#' either estimated using a linearization of the model or an importance sampling Monte Carlo 
#' procedure. When method="bootstrap", the bootstrap estimates are obtained using the bootmlx
#' function
#' 
#' @param project a Monolix project
#' @param method  method c({"fim"}, "proflike", "bootstrap")
#' @param parameters list of parameters for which confidence intervals are computed (default="all")
#' @param level  confidence level, a real number between 0 and 1 (default=0.90)
#' @param linearization  {TRUE}/FALSE  whether the calculation of the standard errors 
#' or the profile likelihood  is based on a linearization of the model (default=TRUE) 
#' @param nboot number of bootstrat replicates (default=100, used when method="bootstrap")
#' @param settings a list of settings for the profile likelihood method:
#' \itemize{
#' \item \code{max.iter} maximum number of iterations to find the solution (default=10)
#' \item \code{tol.LL} absolute tolerance for -2LL  (default=0.001)
#' \item \code{tol.param} relative tolerance for the parameter (default=0.01)
#' \item \code{print} {TRUE}/FALSE display the results (default=TRUE)
#' }
#' @return a list with the computed confidence intervals, the method used and the level.
#' @examples
#' # RsmlxDemo2.mlxtran is a Monolix project for modelling the PK of warfarin using a PK model 
#' # with parameters ka, V, Cl.
#' 
#' # confintmlx will compute a 90% confidence interval for all the population parameters 
#' # using the population estimates obtained by Monolix and the Fisher Information Matrix 
#' # estimated by linearization
#' r1 <- confintmlx(project="RsmlxDemo2.mlxtran") 
#' 
#' # 95% confidence intervals are now computed, using the FIM estimated by Monolix using a 
#' # stochastic approximation algorithm:
#' r2 <- confintmlx(project="RsmlxDemo2.mlxtran", linearization=FALSE, level=0.95) 
#' 
#' # Confidence intervals are computed for ka_pop and omega_ka only, 
#' # using the profile likelihood method:
#' r <- confintmlx(project    = "RsmlxDemo2.mlxtran", 
#'                 method     = "proflike", 
#'                 parameters = c("ka_pop","omega_ka")) 
#' 
#' # Confidence intervals are computed using 200 bootstrap samples:
#' r3 <- confintmlx(project="RsmlxDemo2.mlxtran", method="bootstrap", nboot=200)
#' 
#' # See http://rsmlx.webpopix.org/userguide/confintmlx/ for detailed examples of use of confintmlx
#' # Download the demo examples here: http://rsmlx.webpopix.org/installation
#'
#' 
#' @importFrom stats qchisq
#' @export
confintmlx <- function(project, parameters="all", method="fim", level=0.90, 
                       linearization=TRUE, nboot=100, settings=NULL)
{
  
  r <- prcheck(project, f="conf", level=level, method=method )
  if (r$demo)
    return(r$res)
  project <- r$project
  
  if (level<=0 | level>=1)
    stop("Level of the confidence interval should be strictly between 0 and 1", call.=FALSE)
  
  launched.tasks <- mlx.getLaunchedTasks()
  if (!launched.tasks[["populationParameterEstimation"]]) {
    cat("\nEstimation of the population parameters... \n")
    mlx.runPopulationParameterEstimation()
  }
  
  
  parameters <- unlist(parameters)
  if (method=="proflike") {
    if (identical(parameters, "all"))
      parameters=NULL
    settings$dLLthreshold <- qchisq(level,1)
    r <- llp(project, parameters, settings)
    r$level <- level
    r$method <- "proflike"
    return(r)
  }
  
  if (method=="bootstrap") {
    r.boot <- bootmlx(project, nboot=nboot, settings=list(plot=FALSE, level=level))
    c.inf <- apply(r.boot,MARGIN=2, quantile,(1-level)/2)
    c.sup <- apply(r.boot,MARGIN=2, quantile,(1+level)/2)
    lp <- mlx.loadProject(project) 
    if (!lp) return()
    c.est <- mlx.getEstimatedPopulationParameters()
    ci <- data.frame(estimate=mlx.getEstimatedPopulationParameters(), lower=c.inf, upper=c.sup)
    return(list(confint=ci, level=level, method="bootstrap" ))
  }
  
  
  if (!linearization) {
    if (!("stochasticApproximation" %in% launched.tasks[["standardErrorEstimation"]]) ) {
      mlx.runStandardErrorEstimation(linearization=FALSE)
      cat("Estimation of the standard errors ... \n")
    }    
    se <- as.numeric(unlist(mlx.getEstimatedStandardErrors()$stochasticApproximation))
    names(se) <- names(mlx.getEstimatedStandardErrors()$stochasticApproximation)
  } else {
    if (!("linearization" %in% launched.tasks[["standardErrorEstimation"]])) {
      cat("Estimation of the standard errors ... \n")
      mlx.runStandardErrorEstimation(linearization=TRUE)
    }
    se <- mlx.getEstimatedStandardErrors()$linearization
  }
  
  param <- mlx.getEstimatedPopulationParameters()
  pname <- names(param)
  
  io <- match(pname, names(se))
  se <- se[io]
  
  indParam.name <- mlx.getIndividualParameterModel()$name
  indParam.dist <- mlx.getIndividualParameterModel()$distribution
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
  ci <- cbind(param, ci)
  ci <- as.data.frame(ci)
  names(ci) <- c("estimate","lower", "upper")
  if (!identical(parameters, "all"))
    ci <- ci[parameters,]
  return(list(confint=ci, level=level, method="fim"))
}

