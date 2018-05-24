#' Automatic model building
#'
#' buildmlx uses SAMBA (Stochastic Approximation for Model Building Algorithm), an iterative procedure to accelerate and optimize the process of model building by identifying 
#' at each step how best to improve some of the model components. This method allows to find 
#' the optimal statistical model which minimizes some information criterion in very few steps.
#' 
#' Penalization criterion can be either a custom penalization of the form gamma*(number of parameters),
#' AIC (gamma=2) or BIC (gamma=log(N)).
#' 
#' Several strategies can be used for building the covariate model at each iteration of the algorithm: 
#' \code{direction="full"} means that all the possible models are compared (default when the number of covariates
#' is less than 10). Othrwise, \code{direction} is the mode of stepwise search of \code{stepAIC {MASS}}, 
#' can be one of "both", "backward", or "forward", with a default of "both" when there are at least 10 covariates.
#' 
#' See http://rsmlx.webpopix.org for more details.
#' @param project a string: the initial Monolix project
#' @param final.project  a string: the final Monolix project (default adds "_built" to the original project)
#' @param model  components of the model to optimize c("residualError", "covariate", "correlation"), (default="all")
#' @param paramToUse  list of parameters possibly function of covariates (default="all")
#' @param covToTest  components of the covariate model that can be modified   (default="all")
#' @param covToTransform  list of (continuous) covariates to be log-transformed (default="none")
#' @param criterion  penalization criterion to optimize c({"BIC"}, "AIC", gamma)
#' @param direction method for covariate search c({"full"}, "both", "backward", "forward"), (default="full" or "both")
#' @param max.iter maximum number of iterations (default=20)
#' @param exp.iter  number of iterations during the exploratory phase (default=1)
#' @param print {TRUE}/FALSE display the results (default=TRUE)
#' @param nb.model number of models to display at each iteration (default=1)
#' @param linearization  TRUE/{FALSE} whether the computation of the likelihood is based on a linearization of the model (default=FALSE)
#' @param seqcc TRUE/{FALSE} whether the covariate model is built before the correlation model (default=FALSE) 
#' @param p.max maximum p-value (for the correlation test) for keeping a covariate in a model  (default=1)
#' @param steps maximum number of iteration for stepAIC (default=1000)
#' @return a new Monolix project with a new statistical model.
#' @importFrom MASS stepAIC 
#' @importFrom stats coef as.formula model.matrix
#' @examples
#' \dontrun{
#' # the 3 covariates available ("wt","age", "sex") are considered:
#' r = buildmlx("warfPK.mlxtran")
#'   
#' # log-transformation of the continuous covariates are also considered
#' # the covariate model is built for V and Cl only 
#' r = buildmlx("warfPK.mlxtran", covToTransform="all", paramToUse=c("V", "Cl"))  
#' }
#' @export
buildmlx <- function(project, final.project=NULL, model="all", 
                     paramToUse="all", covToTest="all", covToTransform="none", 
                     criterion="BIC", direction=NULL,
                     max.iter=20, print=TRUE, nb.model=1, linearization=FALSE, 
                     seqcc=FALSE, p.max=1, steps=1000, exp.iter=2)
  #lambda='cv', glmnet.settings=NULL)
{
  
  if (!grepl("\\.",project))
    project <- paste0(project,".mlxtran")
  if(!file.exists(project)){
    message(paste0("ERROR: project '", project, "' does not exists"))
    return(invisible(FALSE))}
  lp <- loadProject(project) 
  if (!lp) return(invisible(FALSE))
  
  if (length(getIndividualParameterModel()$variability)>1)
    stop("Multiple levels of variability are not supported in this version of buildmlx", call.=FALSE)
  #  initializeMlxConnectors(software = "monolix")
  
  # if (length(criterion)==1) {
  #   if (criterion=="lasso") {
  #     criterion <- c(criterion, "BIC", "BIC")
  #   } else {
  #     criterion <- rep(criterion, 3)
  #   }
  # } else if (length(criterion)==2) 
  #   criterion <- c(criterion, criterion[2])
  # 
  # criterion <- criterion[1:3]
  # 
  
  r <- check(project, final.project, covToTransform, paramToUse, covToTest, model, direction)
  covToTransform <- r$covToTransform
  paramToUse <- r$paramToUse
  final.project <- r$final.project
  covToTest <- r$covToTest
  model <- r$model
  direction <- r$direction
  
  final.dir <- sub(pattern = "(.*)\\..*$", replacement = "\\1", final.project)
  if (dir.exists(final.dir)) {
    unlink(final.dir, recursive=TRUE)
  }
  
  ptm <- proc.time()
  Sys.sleep(0.1)
  project.dir <- getProjectSettings()$directory
  if (!dir.exists(project.dir))
    dir.create(project.dir)
  buildmlx.dir <- file.path(getProjectSettings()$directory,"buildmlx")
  Sys.sleep(0.1)
  if (dir.exists(buildmlx.dir))
    unlink(buildmlx.dir, recursive=TRUE)
  
  Sys.sleep(0.1)
  dir.create(buildmlx.dir)
  summary.file = file.path(buildmlx.dir,"summary.txt")
  
  launched.tasks <- getLaunchedTasks()
  Sys.sleep(0.1)
  dir.create(final.dir)
  
  if (!launched.tasks[["populationParameterEstimation"]]) {
    lineDisplay <- "\nEstimation of the population parameters using the initial model ... \n"
    if (print) cat(lineDisplay)
    runPopulationParameterEstimation()
  }
  if (!launched.tasks[["conditionalDistributionSampling"]]) {
    lineDisplay <- "Sampling of the conditional distribution using the initial model ... \n"
    if (print) cat(lineDisplay)
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
  iop.ll <- 1
  lin.ll <- linearization
  method.ll <- ifelse (lin.ll,"linearization","importanceSampling")
  if (iop.ll) {
    if (!(method.ll %in% launched.tasks[["logLikelihoodEstimation"]]))  {
      if (lin.ll & !launched.tasks[["conditionalModeEstimation"]])
        runConditionalModeEstimation()
      lineDisplay <- ("Estimation of the log-likelihood of the initial model ... \n")
      if (print) cat(lineDisplay)
      runLogLikelihoodEstimation(linearization = lin.ll)
    }
    ll.ini <- computecriterion(criterion, method.ll)
    ll <- getEstimatedLogLikelihood()[[method.ll]]
    names(ll)[which(names(ll)=="standardError")] <- "s.e."
    if (is.numeric(criterion))
      ll['criterion'] <- ll.ini
  }
  
  p.ini <- getPopulationParameterInformation()
  rownames(p.ini) <- p.ini$name
  ind.omega <- grep("omega_",p.ini[['name']])
  omega <- p.ini$name[ind.omega]
  omega.ini <- p.ini[ind.omega,]
  
  error.model <- getContinuousObservationModel()$errorModel
  obs.dist <- getContinuousObservationModel()$distribution
  covariate.model <- getIndividualParameterModel()$covariateModel
  cov.ini <- names(covariate.model[[1]])
  correlation.model <- getIndividualParameterModel()$correlationBlocks$id
  if (is.null(correlation.model))
    correlation.model <- list()
  
  #-----------------------------------------
  
  res.dir <- getProjectSettings()$directory
  if (dir.exists(res.dir)) 
    list_of_files <- c(list.files(res.dir) , ".Internals")
  foo <- file.copy(file.path(res.dir,list_of_files), final.dir, recursive=TRUE)
  foo <- file.rename(file.path(final.dir,".Internals",basename(project)),
                     file.path(final.dir,".Internals",basename(final.project)))
  
  setProjectSettings(directory=final.dir)
  saveProject(final.project)
  loadProject(final.project)
  
  if (!is.null(r$idir)) {
    if (print)
      cat(paste0('\n\ndirection: = "',direction, '" will be used for the covariate search\n\n'))
    sink(summary.file)
    cat(paste0('\n\ndirection: = "',direction, '" will be used for the covariate search\n\n'))
    sink()
  }
  
  
  
  if (print) {
    cat("____________________________________________\nInitialization:")
    if (iop.covariate) {
      cat("\nCovariate model:\n")
      print(formatCovariateModel(covariate.model, cov.ini))
    }
    if (iop.correlation) {
      cat("\nCorrelation model:\n")
      print(correlation.model)
    }
    if (iop.error) {
      cat("\nResidual error model:\n")
      print(formatErrorModel(error.model))
    }
    if (iop.ll) {
      cat(paste0("\nEstimated criteria (",method.ll,"):\n"))
      print(round(ll,2))
    }
  }
  if (file.exists(summary.file)) sink(summary.file, append=TRUE)
  else sink(summary.file)
  cat("____________________________________________\nInitialization:")
  if (iop.covariate) {
    cat("\nCovariate model:\n")
    print(formatCovariateModel(covariate.model, cov.ini))
  }
  if (iop.correlation) {
    cat("\nCorrelation model:\n")
    print(correlation.model)
  }
  if (iop.error) {
    cat("\nResidual error model:\n")
    print(formatErrorModel(error.model))
  }
  if (iop.ll) {
    cat(paste0("\nEstimated criteria (",method.ll,"):\n"))
    print(round(ll,2))
  }
  sink()
  
  
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
  sp0 <- NULL
  while (!stop.test) {
    iter <- iter + 1
    
    obs.dist0 <- obs.dist
    error.model0 <- error.model
    covariate.model0 <- covariate.model
    correlation.model0 <- correlation.model
    if (iop.ll)
      ll0 <- ll
    
    if (iop.error) {
      res.error <- errorModelSelection(criterion=criterion, nb.model=nb.model)
      if (nb.model==1)
        error.model <- res.error
      else {
        error.model <- getContinuousObservationModel()$errorModel
        for (k in (1:length(error.model)))
          error.model[[k]] <- as.character(res.error[[k]]$error.model[1])
      }
    }
    
    if (iop.covariate) {
      res.covariate <- covariateModelSelection(criterion=criterion, nb.model=nb.model,
                                               covToTransform=covToTransform, covToTest=covToTest, direction=direction, 
                                               steps=steps, p.max=p.max, paramToUse=paramToUse, sp0=sp0)
      res.covariate$res <- sortCov(res.covariate$res, cov.ini)
      if (iter>exp.iter) sp0 <- res.covariate$sp
      covToTransform <- setdiff(covToTransform, res.covariate$tr0)
      covariate.model <- res.covariate$model
      e <- res.covariate$residuals
      cov.names <- lapply(covariate.model, function(x) {sort(names(which(x)))})
      cov.names0 <- lapply(covariate.model0, function(x) {sort(names(which(x)))})
    } else {
      e <- getSimulatedRandomEffects()
      #     e$id <- e$rep <- NULL
    }
    
    if (iop.correlation & !corr.test) {
      if (isTRUE(all.equal(cov.names0,cov.names)) &
          isTRUE(all.equal(error.model0,error.model))) {
        corr.test <- TRUE
        lineDisplay <- "Start building correlation model too ... \n"
        sink(summary.file, append=TRUE); cat(lineDisplay); sink(); if (print) cat(lineDisplay)
        
      }
    }
    
    if (iop.correlation & corr.test) {
      res.correlation <- correlationModelSelection(e=e, criterion=criterion, nb.model=nb.model, 
                                                   corr0=correlation.model0, seqcc=seqcc)
      if (nb.model==1) correlation.model <- res.correlation
      else  correlation.model <- res.correlation[[1]]$blocks
      if (is.null(correlation.model))  correlation.model <- list()
    } else {
      res.correlation <- getIndividualParameterModel()$correlationBlocks$id
    }
    #-------------------------------
    
    if (max.iter>0) {
      if (!iop.correlation | corr.test) {
        if (isTRUE(all.equal(cov.names0,cov.names)) &
            isTRUE(all.equal(error.model0,error.model)) &
            isTRUE(all.equal(obs.dist0,obs.dist)) &
            isTRUE(all.equal(correlation.model0,correlation.model))) {
          stop.test <- TRUE
          lineDisplay <- "No difference between two successive iterations\n"
          sink(summary.file, append=TRUE); cat(lineDisplay); sink(); if (print) cat(lineDisplay)
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
        sink(summary.file, append=TRUE)
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
        sink()
      }
    } else if (!stop.test | (nb.model>1 & max.iter==0)) {
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
      sink(summary.file, append=TRUE)
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
      sink()
    }
    
    if (!stop.test) {
      setInitialEstimatesToLastEstimates()
      p.ini <- getPopulationParameterInformation()
      rownames(p.ini) <- p.ini$name
      p.ini[omega,] <- omega.ini
      setPopulationParameterInformation(p.ini)
      
      #  setPopulationParameterEstimationSettings(simulatedAnnealing=FALSE)
      
      if (iop.error) {
        emodel <- error.model
        odist <- getContinuousObservationModel()$distribution
        for (k in (1:length(emodel))) {
          if (identical(emodel[[k]],"exponential")) {
            emodel[[k]] <- "constant"
            odist[[k]] <- "lognormal"
          } else {
            odist[[k]] <- "normal"
          }
        }
        setErrorModel(emodel)
        setObservationDistribution(odist)
      }
      
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
      
      if (max.iter>0) {
        if (ll.new > ll.prev) {
          g <- getGeneralSettings()
          g$autochains <- FALSE
          g$nbchains <- g$nbchains+1
          setGeneralSettings(g)
        }
        ll.prev <- ll.new
        g=getConditionalDistributionSamplingSettings()
        g$nbminiterations <- max(100, g$nbminiterations)
        setConditionalDistributionSamplingSettings(g)
        
        
        buildmlx.dir.iter <- file.path(buildmlx.dir,paste0("iteration",iter))
        buildmlx.project.iter <- paste0(buildmlx.dir.iter,".mlxtran")
        saveProject(buildmlx.project.iter)
      }
      
      saveProject(final.project)
      if (dir.exists(final.dir))
        unlink(final.dir, recursive=TRUE)
      
      
      if (max.iter>0) {
        
        lineDisplay <- paste0("Run scenario for model ",iter," ... \nEstimation of the population parameters... \n")
        sink(summary.file, append=TRUE); cat(lineDisplay); sink(); if (print) cat(lineDisplay)
        
        runPopulationParameterEstimation()
        lineDisplay <- "Sampling from the conditional distribution... \n"
        sink(summary.file, append=TRUE); cat(lineDisplay); sink(); if (print) cat(lineDisplay)
        
        runConditionalDistributionSampling()
        if (iop.ll) {
          lineDisplay <- "Estimation of the log-likelihood... \n"
          sink(summary.file, append=TRUE); cat(lineDisplay); sink(); if (print) cat(lineDisplay)
          if (lin.ll)
            runConditionalModeEstimation()
          runLogLikelihoodEstimation(linearization = lin.ll)
        }
        ll.new <- computecriterion(criterion, method.ll)
        
        if (!dir.exists(buildmlx.dir.iter))
          dir.create(buildmlx.dir.iter)
        
        if (dir.exists(final.dir)) 
          list_of_files <- setdiff(c(list.files(final.dir) , ".Internals"), "buildmlx")
        foo <- file.copy(file.path(final.dir,list_of_files), buildmlx.dir.iter, recursive=TRUE)
        foo <- file.rename(file.path(buildmlx.dir.iter,".Internals",basename(final.project)),
                           file.path(buildmlx.dir.iter,".Internals",basename(buildmlx.project.iter)))
        
        #loadProject(final.project)
        
        if (stop.test)
          ll <- ll0
        else {
          ll <- getEstimatedLogLikelihood()[[method.ll]]
          names(ll)[which(names(ll)=="standardError")] <- "s.e."
          if (is.numeric(criterion))
            ll['criterion'] <- ll.new
        }
        
        if (iop.ll) {
          if (print) {
            cat(paste0("\nEstimated criteria (",method.ll,"):\n"))
            print(round(ll,2))
          }
          sink(summary.file, append=TRUE)
          cat(paste0("\nEstimated criteria (",method.ll,"):\n"))
          print(round(ll,2))
          sink()
        }
        if (iter >= max.iter) {
          stop.test <- TRUE
          lineDisplay <- "Maximum number of iterations reached\n"
          sink(summary.file, append=TRUE); cat(lineDisplay); sink(); if (print) cat(lineDisplay)
        }
      }
    }
    if (max.iter==0) stop.test <- TRUE
    saveProject(final.project)
  }
  
  if (print) {
    cat("____________________________________________\n")
    cat("Final model:\n")
    if (iop.covariate) {
      cat("\nCovariate model:\n")
      print(formatCovariateModel(covariate.model))
    }
    if (iop.correlation) {
      cat("\nCorrelation model:\n")
      print(correlation.model)
    }
    if (iop.error) {
      cat("\nResidual error model:\n")
      print(formatErrorModel(error.model))
    }
    if (iop.ll & max.iter>0) {
      cat(paste0("\nEstimated criteria (",method.ll,"):\n"))
      print(round(ll,2))
    }
  }
  sink(summary.file, append=TRUE)
  cat("____________________________________________\n")
  cat("Final model:\n")
  if (iop.covariate) {
    cat("\nCovariate model:\n")
    print(formatCovariateModel(covariate.model))
  }
  if (iop.correlation) {
    cat("\nCorrelation model:\n")
    print(correlation.model)
  }
  if (iop.error) {
    cat("\nResidual error model:\n")
    print(formatErrorModel(error.model))
  }
  if (iop.ll & max.iter>0) {
    cat(paste0("\ncriterion (",method.ll,"):\n"))
    print(round(ll,2))
  }
  sink()
  test.del <- FALSE
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
      test.del <- TRUE
    }
  }
  if (test.del & dir.exists(final.dir)) {
    unlink(final.dir, recursive=TRUE)
  }
  
  saveProject(final.project)
  # con        = file(summary.file, open = "r")
  # lines      = readLines(con, warn=FALSE)
  # close(con)
  # summary.file = file.path(final.dir,"buildmlx.txt")
  # write(lines, file=summary.file)
  
  res <- list(project=final.project)
  if (iop.covariate)
    res <- c(res, list(covariate.model=covariate.model))
  if (iop.correlation)
    res <- c(res, list(correlation.model=correlation.model))
  if (iop.error)
    res <- c(res, list(error.model=error.model))
  
  dt <- proc.time() - ptm
  sink(summary.file, append=TRUE)
  cat("____________________________________________\n")
  cat(paste0("total time: ", round(dt["elapsed"], digits=1),"s\n"))
  cat("____________________________________________\n")
  sink()
  if (print) {
    cat("____________________________________________\n")
    cat(paste0("total time: ", round(dt["elapsed"], digits=1),"s\n"))
    cat("____________________________________________\n")
  }
  return(res)
}


check <- function(project, final.project, covToTransform, paramToUse, covToTest, model, direction) {
  
  if (is.null(final.project)) 
    final.project <- paste0(sub(pattern = "(.*)\\..*$", replacement = "\\1", project),"_built.mlxtran")
  if (!grepl("\\.",final.project))
    final.project <- paste0(final.project,".mlxtran")
  if (!grepl("\\.mlxtran",final.project))
    stop(paste0(final.project, " is not a valid name for a Monolix project (use the .mlxtran extension)"), call.=FALSE)
  
  
  cov.info <- getCovariateInformation()
  cov.names <- cov.info$name
  cov.types <- cov.info$type
  j.trans <- grep("transformed",cov.types)
  j.cont <- grep("continuous",cov.types)
  cont.cov <- cov.names[j.cont]
  if (identical(covToTransform,"all"))
    covToTransform = cov.names[j.cont]
  else if (identical(covToTransform,"none"))
    covToTransform=NULL
  if (!is.null(covToTransform)) {
    ncat0 <- names(which(cov.types[covToTransform]=="categorical"))
    if (length(ncat0)>0) {
      warning(paste0(ncat0, " is a categorical covariate and will not be transformed..."), call.=FALSE)
      covToTransform <- setdiff(covToTransform, ncat0)
    }
    ncov0 <- covToTransform[(!(covToTransform %in% cov.names))]
    if (length(ncov0)>0) {
      warning(paste0(ncov0, " is not a valid covariate and will be ignored"), call.=FALSE)
      covToTransform <- setdiff(covToTransform, ncov0)
    }
  }
  if (identical(covToTest,"all"))
    covToTest = cov.names
  else if (identical(covToTest,"none"))
    covToTest=NULL
  if (!is.null(covToTest)) {
    ncov0 <- covToTest[(!(covToTest %in% cov.names))]
    if (length(ncov0)>0)  stop(paste0(ncov0, " is not a valid covariate"), call.=FALSE)
  }
  
  ind.dist <- getIndividualParameterModel()$distribution
  param.names <- names(ind.dist)
  if (identical(paramToUse,"all"))
    paramToUse <- param.names
  nparam0 <- paramToUse[(!(paramToUse %in% param.names))]
  if (length(nparam0)>0) 
    stop(paste0(nparam0, " is not a valid parameter"), call.=FALSE)
  
  model.names <- c("residualError", "covariate", "correlation")
  if (identical(model,"all")) model <- model.names
  mod0 <- model[(!(model %in% model.names))]
  if (length(mod0)>0) 
    stop(paste0(mod0, " is not a valid model component"), call.=FALSE)
  
  idir <- NULL
  if ("covariate" %in% model) {
    dir.names <- c("full", "both", "backward", "forward")
    if (is.null(direction)) {
      nbcov <- length(getCovariateInformation()$name)
      direction <- ifelse(nbcov<=10,"full","both")
      idir <- direction
    }
    dir0 <- direction[(!(direction %in% dir.names))]
    if (length(dir0)>0) 
      stop(paste0(dir0, " is not a valid direction name"), call.=FALSE)
  }
  
  return(list(covToTransform=covToTransform, paramToUse=paramToUse, covToTest=covToTest, 
              final.project=final.project, model=model, direction=direction, idir=idir))
}

#------------------------

formatCovariateModel <- function(m, cov.ini=NULL) {
  i0 <- which(unlist(lapply(m,function(x) {identical(x,"none")})))
  if (length(i0)>0)
    m <- m[-i0]
  param.names <- names(m)
  cov.names <- names(m[[1]])
  if (length(cov.names)>0) {
    if (is.vector(m[[1]]) | (!is.null(nrow(m[[1]])) && nrow(m[[1]])==1)) {
      mr <- matrix(unlist(m),ncol=length(cov.names),byrow=TRUE)*1
      row.names(mr) <- param.names
      colnames(mr) <- cov.names
      mr <- as.data.frame(mr[,order(colnames(mr))])
    } else {
      mr <- list()
      nr <- nrow(m[[1]])
      for (k in 1:nr) {
        mk <- lapply(m, function(x) x[k,])
        mrk <- matrix(unlist(mk),ncol=length(cov.names),byrow=TRUE)*1
        row.names(mrk) <- param.names
        colnames(mrk) <- cov.names
        mrk <- mrk[,order(colnames(mrk))]
        mr[[k]] <- as.data.frame(mrk)
      }
    }
    if (!is.null(cov.ini))
      mr <- sortCov(mr, cov.ini)
    return(mr)
  }
}

sortCov <- function(r, cov.ini) {
  n1 <- cov.ini
  n3 <- c("ll", "df", "criterion")
  if (is.data.frame(r)) {
    n2 <- setdiff(names(r), c(n1,n3))
    rs <- r[,c(n1,n2)]
    #    rs <- r[c(cov.ini, setdiff(names(r),cov.ini))]
  } else {
    n2 <- setdiff(names(r[[1]]), c(n1,n3))
    rs <- lapply(r, function(x) x[c(n1,n2,n3)])
  }
  return(rs)
}

formatErrorModel <- function(m) {
  out.names <- names(m)
  return(m)
}

computecriterion <- function(criterion, method.ll) {
  ll <- getEstimatedLogLikelihood()[[method.ll]]
  if (identical(criterion,"AIC")) cr <- ll[["AIC"]]
  else if (identical(criterion,"BIC")) cr <- ll[["BIC"]]
  else {
    d <- (ll[["AIC"]] - ll[["-2LL"]])/2
    cr <- ll[["-2LL"]]+d*criterion
  }
  return(cr)
}

