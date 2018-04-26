#' Build a model
#'
#' Get the individual predictions obtained with the simulated individual parameters :
#' @return a list of data frames (one data frame per output).
#' @examples
#' \dontrun{
#' r = getSimulatedPredictions()
#' names(r)
#'     "Cc"  "E"
#' }
#' @export
buildmlx <- function(project=NULL, final.project=NULL, max.iter=100,
                     model=c("residualError", "covariate", "correlation"), penalization="BIC",
                     nb.model=1, seq=TRUE, trans.cov="all", param.cov="all", print.out=TRUE)
{
  
  initializeMlxConnectors(software = "monolix")
  
  if (is.null(project)) {
    dp <- getProjectSettings()$directory
    if (!is.null(dp))
      project <- paste0(basename(dp),".mlxtran")
    else
      project <- "temp"
    saveProject(project)
  } else {
    loadProject(project)
  }
  
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
  iop.ll <- as.vector(getScenario()$tasks['logLikelihoodEstimation'])
  lin.ll <- getScenario()$linearization
  method.ll <- ifelse (lin.ll,"linearization","importanceSampling")
  if (iop.ll) {
    if (!(method.ll %in% launched.tasks[["logLikelihoodEstimation"]]))  {
      runLogLikelihoodEstimation(linearization = lin.ll)
      cat("Estimation of the log-likelihood of the initial model ... \n")
    }
    ll <- getEstimatedLogLikelihood()
  }
  
  error.model <- getContinuousObservationModel()$errorModel
  covariate.model <- getIndividualParameterModel()$covariateModel
  correlation.model <- getIndividualParameterModel()$correlationBlocks$id
  if (is.null(correlation.model))
    correlation.model <- list()
  
  if (print.out) {
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
    final.project <- paste0(sub(pattern = "(.*)\\..*$", replacement = "\\1", project),"_builtmlx.mlxtran")
  final.dir <- sub(pattern = "(.*)\\..*$", replacement = "\\1", final.project)
  
  res.dir <- getProjectSettings()$directory
  if (dir.exists(res.dir)) {
    list_of_files <- c(list.files(res.dir) , ".Internals")
    if (dir.exists(final.dir))
      unlink(final.dir, recursive=TRUE)
  }
  dir.create(final.dir)
  #browser()
  # foo = tryCatch( {
  #   file.copy(file.path(res.dir,list_of_files), final.dir, recursive=TRUE)}
  #   , error=function(e) {
  #     return(list())        }
  # )
  # if (length(foo)==0)
  foo <- file.copy(file.path(res.dir,list_of_files), final.dir, recursive=TRUE)
  foo <- file.rename(file.path(final.dir,".Internals",project),file.path(final.dir,".Internals",final.project))
  
  setProjectSettings(directory=final.dir)
  saveProject(final.project)
  loadProject(final.project)
  
  #--------------
  
  stop.test <- FALSE
  if (seq==TRUE)
    corr.test <- FALSE
  else
    corr.test <- TRUE
  iter <- 0
  
  cov.names0 <- cov.names <- NULL
  
  while (!stop.test) {
    iter <- iter + 1
    
    error.model0 <- error.model
    covariate.model0 <- covariate.model
    correlation.model0 <- correlation.model
    if (iop.ll)
      ll0 <- ll
    
    if (iop.error) {
      res.error <- errorModelSelection(penalization=penalization, nb.model=nb.model)
      error.model <- getContinuousObservationModel()$errorModel
      for (k in (1:length(error.model)))
        error.model[[k]] <- as.character(res.error[[k]]$error.model[1])
    }
    
    if (iop.covariate) {
      res.covariate <- covariateModelSelection(penalization=penalization, nb.model=nb.model,
                                               trans.cov=trans.cov, param.cov=param.cov)
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
      res.correlation <- correlationModelSelection(e=e, penalization=penalization,
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
      
      if (print.out) {
        cat("____________________________________________\n")
        cat(paste0("Iteration ",iter,":\n"))
        if (iop.error) {
          cat("\nResidual error model:\n")
          print(res.error)
        }
        if (iop.covariate) {
          cat("\nCovariate model:\n")
          print(res.covariate$res)
        }
        if (iop.correlation) {
          cat("\nCorrelation model:\n")
          print(res.correlation)
        }
      }
      
      setInitialEstimatesToLastEstimates()
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
      
      saveProject(final.project)
      cat(paste0("Run scenario for model ",iter," ... \n"))
      #runScenario()
      cat("Estimation of the population parameters... \n")
      runPopulationParameterEstimation()
      cat("Sampling from the conditional distribution... \n")
      runConditionalDistributionSampling()
      if (iop.ll) {
        cat("Estimation of the log-likelihood... \n")
        runLogLikelihoodEstimation(linearization = lin.ll)
      }
      
      loadProject(final.project)
      
      if (print.out) {
        if (iop.ll) {
          cat("\nEstimated log-likelihood:\n")
          if (stop.test)
            ll <- ll0
          else
            ll <- getEstimatedLogLikelihood()
          print(lapply(ll,round,2))
        }
      }
      if (iter >= max.iter) {
        stop.test <- TRUE
        cat("Maximum number of iterations reached\n")
      }
    }
  }
  
  if (print.out) {
    cat("____________________________________________\n")
    cat("Final model:\n")
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
  loadProject(project)
  res <- list(project=final.project)
  if (iop.error)
    res <- c(res, list(res.error=res.error))
  if (iop.covariate)
    res <- c(res, list(res.covariate=res.covariate$res))
  if (iop.correlation)
    res <- c(res, list(res.correlation=res.correlation))
  return(res)
}
