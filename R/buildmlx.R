#' Automatic model building
#'
#' Iterative procedure to accelerate and optimize the process of model building by identifying 
#' at each step how best to improve some of the model components. This method allows to find 
#' the optimal statistical model which minimizes some information criteria in very few steps.
#' @param project a string: the initial Monolix project
#' @param final.project  a string: the final Monolix project (default adds "_built" to the original project)
#' @param model  component of the model to optimize, default=c("residualError", "covariate", "correlation")
#' @param penalization  penalization criteria to optimize c({"BIC"}, "AIC")
#' @param max.iter maximum number of iterations (default=20)
#' @param covToTransform  list of covariates to be log-transformed (default="none")
#' @param paramToUse  list of parameters which could be function of covariates (default="all")
#' @param linearization  {TRUE}/FALSE whether the computation of the likelihood is based on a linearization of the model (default=TRUE)
#' @param seqcc {TRUE}/FALSE whether the covariate model is built before the correlation model (default=TRUE) 
#' @param nb.model number of models to display at each iteration (default=1)
#' @param print {TRUE}/FALSE display the results (default=TRUE)
#' @param direction method for covariate search c({"full"}, "both", "backward", "forward"), (default="full")
#' @param steps number of iteration for stepAIC/BIC (default=1000)
#' @param p.min minimum p-value (for the correlation test) for keeping a covariate in a model  (default=0.2)
#' @return a new Monolix project with a new statistical model.
#' @importFrom MASS stepAIC 
#' @importFrom stats coef 
#' @examples
#' \dontrun{
#' r = buildmlx("warfPK.mlxtran")
#' }
#' @export
buildmlx <- function(project, final.project=NULL, model=c("residualError", "covariate", "correlation"), 
                     penalization="BIC", max.iter=20, covToTransform="none", paramToUse="all", linearization=TRUE,
                     seqcc=TRUE, nb.model=1, print=TRUE, direction='full', steps=1000, p.min=0.2)
                     #lambda='cv', settings=NULL)
{
  
# #' @param lambda penalization coefficient for "lasso" (default="cv")
# #' @param settings when penalization="lasso", settings used for package glmnet 
  
  
  if(!file.exists(project)){
    message(paste0("ERROR: project '", project, "' does not exists"))
    return(invisible(FALSE))}
  
  #  initializeMlxConnectors(software = "monolix")
  
  if (length(penalization)==1) {
    if (penalization=="lasso") {
      penalization <- c(penalization, "BIC", "BIC")
    } else {
      penalization <- rep(penalization, 3)
    }
  } else if (length(penalization)==2) 
    penalization <- c(penalization, penalization[2])
  
  penalization <- penalization[1:3]
  
  
  loadProject(project)   
  
  launched.tasks <- getLaunchedTasks()
  if (!launched.tasks[["populationParameterEstimation"]]) {
    cat("\nEstimation of the population parameters using the initial model ... \n")
    runPopulationParameterEstimation()
  }
  if (!launched.tasks[["conditionalDistributionSampling"]]) {
    cat("Sampling of the conditional distribution using the initial model ... \n")
    runConditionalDistributionSampling()
  }
  
  #------------------------------
  if (!any(getIndividualParameterModel()$variability$id))
    stop("\nA least one parameter with random effects is required\n", call.=FALSE)
  iop.error <- "residualError" %in% model
  if (is.null(getContinuousObservationModel()))
    iop.error <- FALSE
  iop.covariate <- "covariate" %in% model
  if (is.null(getCovariateInformation()))
    iop.covariate <- FALSE
  iop.correlation <- "correlation" %in% model
  if (sum(getIndividualParameterModel()$variability$id)==1)
    iop.correlation <- FALSE
  if (!any(c(iop.error, iop.covariate, iop.correlation)))
    stop("\nThere is no statistical model to build...\n", call.=FALSE)
  
  
  #----  LL & linearization
  # iop.ll <- as.vector(getScenario()$tasks['logLikelihoodEstimation'])
  # lin.ll <- getScenario()$linearization
  iop.ll <- 1
  lin.ll <- linearization
  method.ll <- ifelse (lin.ll,"linearization","importanceSampling")
  if (iop.ll) {
    if (!(method.ll %in% launched.tasks[["logLikelihoodEstimation"]]))  {
      if (lin.ll & !launched.tasks[["conditionalModeEstimation"]])
        runConditionalModeEstimation()
      runLogLikelihoodEstimation(linearization = lin.ll)
      cat("Estimation of the log-likelihood of the initial model ... \n")
    }
    ll <- getEstimatedLogLikelihood()
    ll.ini <- ll[[method.ll]][penalization[2]]
  }
  
  p.ini <- getPopulationParameterInformation()
  ind.omega <- grep("omega_",p.ini[['name']])
  omega.ini <- p.ini[ind.omega,]
  
  error.model <- getContinuousObservationModel()$errorModel
  covariate.model <- getIndividualParameterModel()$covariateModel
  correlation.model <- getIndividualParameterModel()$correlationBlocks$id
  if (is.null(correlation.model))
    correlation.model <- list()
  
  if (print) {
    cat("____________________________________________\n")
    cat("Initialization:\n")
    if (iop.error) {
      cat("\nResidual error model:\n")
      print(error.model)
    }
    if (iop.covariate) {
      cat("\nCovariate model:\n")
      print(covariate.model)
    }
    if (iop.correlation) {
      cat("\nCorrelation model:\n")
      print(correlation.model)
    }
    if (iop.ll) {
      cat("\nEstimated log-likelihood:\n")
      print(lapply(ll,round,2))
    }
  }
  
  #-----------------------------------------
  
  if (is.null(final.project))
    final.project <- paste0(sub(pattern = "(.*)\\..*$", replacement = "\\1", project),"_built.mlxtran")
  final.dir <- sub(pattern = "(.*)\\..*$", replacement = "\\1", final.project)
  
  res.dir <- getProjectSettings()$directory
  if (dir.exists(res.dir)) {
    list_of_files <- c(list.files(res.dir) , ".Internals")
    if (dir.exists(final.dir))
      unlink(final.dir, recursive=TRUE)
  }
  dir.create(final.dir)
  foo <- file.copy(file.path(res.dir,list_of_files), final.dir, recursive=TRUE)
  foo <- file.rename(file.path(final.dir,".Internals",basename(project)),
                     file.path(final.dir,".Internals",basename(final.project)))
  
  setProjectSettings(directory=final.dir)
  saveProject(final.project)
  loadProject(final.project)
  
  #--------------
  
  stop.test <- FALSE
  if (seqcc==TRUE)
    corr.test <- FALSE
  else
    corr.test <- TRUE
  iter <- 0
  
  cov.names0 <- cov.names <- NULL
  
  ll.prev <- Inf
  ll.new <- ll.ini
  while (!stop.test) {
    iter <- iter + 1
    
    error.model0 <- error.model
    covariate.model0 <- covariate.model
    correlation.model0 <- correlation.model
    if (iop.ll)
      ll0 <- ll
    
    if (iop.error) {
      res.error <- errorModelSelection(penalization=penalization[3], nb.model=nb.model)
      error.model <- getContinuousObservationModel()$errorModel
      for (k in (1:length(error.model)))
        error.model[[k]] <- as.character(res.error[[k]]$error.model[1])
    }
    
    if (iop.covariate) {
      res.covariate <- covariateModelSelection(penalization=penalization[1], nb.model=nb.model,
                                               covToTransform=covToTransform, direction=direction, 
                                               steps=steps, p.min=p.min, paramToUse=paramToUse)
                                              # lambda=lambda,  settings=settings)
      covariate.model <- res.covariate$model
      e <- res.covariate$residuals
      cov.names <- lapply(covariate.model, function(x) {sort(names(which(x)))})
      cov.names0 <- lapply(covariate.model0, function(x) {sort(names(which(x)))})
    } else {
      e <- getSimulatedRandomEffects()
      e$id <- e$rep <- NULL
    }
    
    if (iop.correlation & !corr.test) {
      if (isTRUE(all.equal(cov.names0,cov.names)) &
          isTRUE(all.equal(error.model0,error.model))) {
        corr.test <- TRUE
        cat("Start building correlation model too\n")
      }
    }
    
    if (iop.correlation & corr.test) {
      res.correlation <- correlationModelSelection(e=e, penalization=penalization[2],
                                                   nb.model=nb.model, corr0=correlation.model0)
      correlation.model <- res.correlation$blocks[[1]]
      if (is.null(correlation.model))  correlation.model <- list()
    } else {
      res.correlation <- getIndividualParameterModel()$correlationBlocks$id
    }
    #-------------------------------
    
    
    if (!iop.correlation | corr.test) {
      if (isTRUE(all.equal(cov.names0,cov.names)) &
          isTRUE(all.equal(error.model0,error.model)) &
          isTRUE(all.equal(correlation.model0,correlation.model))) {
        stop.test <- TRUE
        cat("No difference between two successive iterations\n")
      }
    }
    
    if (!stop.test) {
      
      if (print) {
        cat("____________________________________________\n")
        cat(paste0("Iteration ",iter,":\n"))
        if (iop.covariate) {
          cat("\nCovariate model:\n")
          print(res.covariate$res)
        }
        if (iop.correlation) {
          cat("\nCorrelation model:\n")
          print(res.correlation)
        }
        if (iop.error) {
          cat("\nResidual error model:\n")
          print(res.error)
        }
      }
      
      setInitialEstimatesToLastEstimates()
      p.ini <- getPopulationParameterInformation()
      p.ini[ind.omega,] <- omega.ini
      setPopulationParameterInformation(p.ini)
      
      setPopulationParameterEstimationSettings(simulatedAnnealing=FALSE)
      
      if (iop.error)
        setErrorModel(error.model)
      
      if (iop.covariate) {
        if (length(res.covariate$add.covariate) >0) {
          for (k in 1:length(res.covariate$add.covariate))
            eval(parse(text=res.covariate$add.covariate[[k]]))
        }
        setCovariateModel(covariate.model)
      }
      
      if (iop.correlation & corr.test)
        setCorrelationBlocks(correlation.model)
      
      #-------------------------------
      
      if (dir.exists(final.dir))
        unlink(final.dir, recursive=TRUE)
      
      if (ll.new > ll.prev) {
        g <- getGeneralSettings()
        g$autochains <- FALSE
        g$nbchains <- g$nbchains+1
        setGeneralSettings(g)
      }
      ll.prev <- ll.new
      
      saveProject(final.project)
      cat(paste0("Run scenario for model ",iter," ... \n"))
      #runScenario()
      cat("Estimation of the population parameters... \n")
      runPopulationParameterEstimation()
      cat("Sampling from the conditional distribution... \n")
      runConditionalDistributionSampling()
      if (iop.ll) {
        cat("Estimation of the log-likelihood... \n")
        if (lin.ll)
          runConditionalModeEstimation()
        runLogLikelihoodEstimation(linearization = lin.ll)
      }
      ll.new <- ll[[method.ll]][penalization[2]]
      
      loadProject(final.project)
      
      if (stop.test)
        ll <- ll0
      else
        ll <- getEstimatedLogLikelihood()
      
      if (print & iop.ll) {
        cat("\nEstimated log-likelihood:\n")
        print(lapply(ll,round,2))
      }
      if (iter >= max.iter) {
        stop.test <- TRUE
        if (print)
          cat("Maximum number of iterations reached\n")
      }
    }
  }
  
  if (print) {
    cat("____________________________________________\n")
    cat("Final model:\n")
    if (iop.covariate) {
      cat("\nCovariate model:\n")
      print(covariate.model)
    }
    if (iop.correlation) {
      cat("\nCorrelation model:\n")
      print(correlation.model)
    }
    if (iop.error) {
      cat("\nResidual error model:\n")
      print(error.model)
    }
    if (iop.ll) {
      cat("\nEstimated log-likelihood:\n")
      print(lapply(ll,round,2))
    }
  }
  
  if (iop.covariate) {
    foo <- lapply(res.covariate$model,function(x) {which(x)})
    cov.model <- unique(unlist(lapply(foo,function(x) {names(x)})))
    cov.type <- getCovariateInformation()$type[cov.model]
    cov.cont <- names(cov.type[cov.type=="continuous"])
    for (ck in cov.cont) {
      cck <- paste0("c",ck)
      covk <- getCovariateInformation()$covariate[[ck]]
      tr.str <- paste0(cck,' = "',ck,"-",signif(mean(covk),digits=2),'"')
      tr.str <- paste0("addContinuousTransformedCovariate(",tr.str,")")
      eval(parse(text=tr.str))
      g=getIndividualParameterModel()$covariateModel
      cg <- lapply(g, function(x) {foo <- x[ck]; x[ck]<-x[cck]; x[cck]=foo; return(x)})
      setCovariateModel(cg)
    }
  }
  saveProject(final.project)
  res <- list(project=final.project)
  if (iop.covariate)
    res <- c(res, list(covariate.model=covariate.model))
  if (iop.correlation)
    res <- c(res, list(correlation.model=correlation.model))
  if (iop.error)
    res <- c(res, list(error.model=error.model))
  return(res)
}
