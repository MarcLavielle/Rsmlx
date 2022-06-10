#' Automatic statistical model building
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
#lambda='cv', glmnet.settings=NULL)

#' See http://rsmlx.webpopix.org for more details.
#' @param project a string: the initial Monolix project
#' @param final.project  a string: the final Monolix project (default adds "_built" to the original project)
#' @param model  components of the model to optimize c("residualError", "covariate", "correlation"), (default="all")
#' @param prior  list of prior probabilities for each component of the model (default=NULL)
#' @param weight list of penalty weights for each component of the model (default=NULL)
#' @param coef.w1 multiplicative weight coefficient used for the first iteration only (default=0.5)
#' @param paramToUse  list of parameters possibly function of covariates (default="all")
#' @param covToTest  components of the covariate model that can be modified   (default="all")
#' @param covToTransform  list of (continuous) covariates to be log-transformed (default="none")
#' @param center.covariate TRUE/{FALSE} center the covariates of the final model (default=FALSE) 
#' @param criterion  penalization criterion to optimize c("AIC", "BIC", {"BICc"}, gamma)
#' @param test  {TRUE}/FALSE  perform additional statistical tests for building the model (default=TRUE)
#' @param ll  {TRUE}/FALSE  compute the observe likelihood and the criterion to optimize at each iteration
#' @param linearization  TRUE/{FALSE} whether the computation of the likelihood is based on a linearization of the model (default=FALSE)
#' @param fError.min  minimum fraction of residual variance for combined error model (default = 1e-3)
#' @param seq.cov TRUE/{FALSE} whether the covariate model is built before the correlation model  
#' @param seq.cov.iter number of iterations before building the correlation model (only when seq.cov=F, default=0) 
#' @param seq.corr {TRUE}/FALSE whether the correlation model is built iteratively (default=TRUE) 
#' @param p.max  maximum p-value used for removing non significant relationships between covariates and individual parameters (default=0.1)
#' @param p.min vector of 3 minimum p-values used for testing the components of a new model (default=c(0.075, 0.05, 0.1))
#' @param direction method for covariate search c({"full"}, "both", "backward", "forward"), (default="full" or "both")
#' @param steps maximum number of iteration for stepAIC (default=1000)
#' @param n.full  maximum number of covariates for an exhaustive comparison of all possible covariate models (default=10)
#' @param max.iter maximum number of iterations (default=20)
#' @param explor.iter  number of iterations during the exploratory phase (default=2)
#' @param print {TRUE}/FALSE display the results (default=TRUE)
#' @param nb.model number of models to display at each iteration (default=1)
#' @return a new Monolix project with a new statistical model.
#' @examples
#' # RsmlxDemo1.mlxtran is a Monolix project for modelling the pharmacokinetics (PK) of warfarin 
#' # using a PK model with parameters ka, V, Cl.
#' 
#' # By default, buildmlx will compute the best statistical model in term of BIC, i.e , 
#' # the best covariate model, the best correlation model for the three random effects and the best 
#' # residual error model in terms of BIC. 
#' # In this example, three covariates (wt, age, sex) are available with the data and will be used 
#' # for building the covariate model for the three PK parameters:
#' r1 <- buildmlx(project="RsmlxDemo1.mlxtran")
#'   
#' # Here, the covariate model will be built for V and Cl only and log-transformation of all 
#' # continuous covariates will also be considered:
#' r2 <- buildmlx(project="RsmlxDemo1.mlxtran", paramToUse=c("V", "Cl"), covToTransform="all") 
#' 
#' # Only the covariate model will be  built, using AIC instead of BIC:
#' r3 <- buildmlx(project="RsmlxDemo1.mlxtran", model="covariate", criterion="AIC") 
#' 
#' # See http://rsmlx.webpopix.org/userguide/buildmlx/ for detailed examples of use of buildmlx
#' # Download the demo examples here: http://rsmlx.webpopix.org/installation
#' 
#' 
#' @importFrom MASS addterm dropterm 
#' @importFrom stats coef as.formula model.matrix deviance formula extractAIC factor.scope nobs terms update update.formula
#' @importFrom utils data write.csv packageVersion
#' @importFrom dplyr filter select rename arrange bind_rows mutate
#' @importFrom dplyr %>%
#' @export
buildmlx <- function(project=NULL, final.project=NULL, model="all", prior=NULL, weight=NULL, coef.w1=0.5,
                     paramToUse="all", covToTest="all", covToTransform="none", center.covariate=FALSE, 
                     criterion="BICc", linearization=FALSE, ll=T, test=T, 
                     direction=NULL, steps=1000, n.full=10,
                     max.iter=20, explor.iter=2, fError.min=1e-3, 
                     seq.cov=FALSE, seq.cov.iter=0, seq.corr=TRUE, 
                     p.max=0.1, p.min=c(0.075, 0.05, 0.1),
                     print=TRUE, nb.model=1)
{
  
  ptm <- proc.time()
  
  dashed.line <- "--------------------------------------------------\n"
  plain.line <-  "__________________________________________________\n"
  dashed.short <- "-----------------------\n"
  plain.short  <- "_______________________\n"
  
  op.original <- options()
  op.new <- options()
  op.new$lixoft_notificationOptions$warnings <- 1   #hide the warning messages
  options(op.new)
  
  
  RsmlxDemo1.project <- RsmlxDemo2.project <- warfarin.data  <- resMonolix <- NULL
  pi <- 4*atan(1)
  
  if (!is.null(project)) {
    r <- prcheck(project, f="build", paramToUse=paramToUse, model=model)
    if (r$demo)
      return(r$res)
    project <- r$project
  } else {
    project <- mlx.getProjectSettings()$project
  }
  
  method.ll <- iop.ll <- pen.coef <- NULL
  r <- buildmlx.check(project, final.project, model, paramToUse, covToTest, covToTransform, center.covariate, 
                      criterion, linearization, ll, test, direction, steps, max.iter, explor.iter, 
                      seq.cov, seq.cov.iter, seq.corr, p.max, p.min, print, nb.model, prior, weight, n.full)
  if (!is.null(r$change)) return(list(change=F))
  for (j in 1:length(r))
    eval(parse(text=paste0(names(r)[j],"= r[[j]]")))
  
  r <- def.variable(weight=weight, prior=prior, criterion=criterion)
  for (j in 1:length(r))
    eval(parse(text=paste0(names(r)[j],"= r[[j]]")))
  
  is.weight <- weight$is.weight
  is.prior <- NULL
  
  final.dir <- sub(pattern = "(.*)\\..*$", replacement = "\\1", final.project)
  if (dir.exists(final.dir)) 
    unlink(final.dir, recursive=TRUE)
  
  project.dir <- mlx.getProjectSettings()$directory
  if (!dir.exists(project.dir))
    dir.create(project.dir)
  buildmlx.dir <- file.path(mlx.getProjectSettings()$directory,"buildmlx")
  Sys.sleep(0.1)
  if (dir.exists(buildmlx.dir))
    unlink(buildmlx.dir, recursive=TRUE)
  
  Sys.sleep(0.1)
  dir.create(buildmlx.dir)
  summary.file = file.path(buildmlx.dir,"summary.txt")
  
  Sys.sleep(0.1)
  if (!dir.exists(final.dir))
    dir.create(final.dir, recursive=T)
  
  
  
  #-------------------------------------------------
  
  to.cat <- paste0("\n", dashed.line, "\nBuilding:\n")
  if (model$covariate)  to.cat <- c(to.cat, "  -  The covariate model\n")
  if (model$correlation)  to.cat <- c(to.cat, "  -  The correlation model\n")
  if (model$residualError)  to.cat <- c(to.cat, "  -  The residual error model\n")
  to.cat <- c(to.cat, "\n")
  print.result(print, summary.file, to.cat=to.cat, to.print=NULL) 
  
  print.line <- F
  launched.tasks <- mlx.getLaunchedTasks()
  # if (!launched.tasks[["populationParameterEstimation"]]) {
  #   to.cat <- paste0(plain.line,"\nEstimation of the population parameters using the initial model ... \n")
  #   print.result(print, summary.file, to.cat=to.cat, to.print=NULL) 
  #   print.line <- T
  #   mlx.runPopulationParameterEstimation()
  # }
  
  #-------------------------------------------------
  
  p.ini <- mlx.getPopulationParameterInformation()
  rownames(p.ini) <- p.ini$name
  ind.omega <- grep("omega_",p.ini[['name']])
  omega <- p.ini$name[ind.omega]
  omega.ini <- p.ini[ind.omega,]
  
  error.model <- mlx.getContinuousObservationModel()$errorModel
  obs.dist <- mlx.getContinuousObservationModel()$distribution
  covariate.model <- mlx.getIndividualParameterModel()$covariateModel
  cov.ini <- names(covariate.model[[1]])
  correlation.model <- lapply(mlx.getIndividualParameterModel()$correlationBlocks$id, sort)
  if (length(correlation.model)==0)
    correlation.model <- NULL
  
  error.model.ini <- error.model
  covariate.model.ini <- covariate.model
  correlation.model.ini <- correlation.model
  
  to.cat <- ("- - - Initialization - - -\n")
  if (!print.line)
    to.cat <- paste0(plain.line, to.cat)
  print.result(print, summary.file, to.cat=to.cat, to.print=NULL) 
  
  if (model$covariate) {
    to.cat <- "\nCovariate model:\n"
    to.print <- formatCovariateModel(covariate.model, cov.ini)
    print.result(print, summary.file, to.cat=to.cat, to.print=to.print) 
  }
  if (model$correlation) {
    to.cat <- "\nCorrelation model:\n"
    to.print <- ifelse(!is.null(correlation.model), correlation.model, "NULL")
    print.result(print, summary.file, to.cat=to.cat, to.print=to.print) 
  }
  if (model$residualError) {
    to.cat <- "\nResidual error model:\n"
    to.print <- formatErrorModel(error.model)
    print.result(print, summary.file, to.cat=to.cat, to.print=to.print) 
  }
  
  
  if (!launched.tasks[["populationParameterEstimation"]]) {
    to.cat <- "\nEstimation of the population parameters using the initial model ... \n"
    print.result(print, summary.file, to.cat=to.cat, to.print=NULL) 
    mlx.runPopulationParameterEstimation()
  }
  if (!launched.tasks[["conditionalDistributionSampling"]]) {
    to.cat <- "Sampling of the conditional distribution using the initial model ... \n"
    print.result(print, summary.file, to.cat=to.cat, to.print=NULL) 
    mlx.runConditionalDistributionSampling()
  }
  
  gi <- mlx.getSimulatedIndividualParameters()
  gi <- gi %>% filter(rep==gi$rep[nrow(gi)]) %>% select(-rep)
  lin.ll <- method.ll=="linearization"
  if (iop.ll) {
    if (!(method.ll %in% launched.tasks[["logLikelihoodEstimation"]]))  {
      if (lin.ll & !launched.tasks[["conditionalModeEstimation"]])
        mlx.runConditionalModeEstimation()
      
      to.cat <- "Estimation of the log-likelihood of the initial model ... \n"
      print.result(print, summary.file, to.cat=to.cat, to.print=NULL) 
      
      mlx.runLogLikelihoodEstimation(linearization = lin.ll)
    }
    
    ll.ini <- compute.criterion(criterion, method.ll, weight, pen.coef)
    ll <- formatLL(mlx.getEstimatedLogLikelihood()[[method.ll]], criterion, ll.ini, is.weight, is.prior)
    list.criterion <- ll.ini
    
    to.cat <- paste0("\nEstimated criteria (",method.ll,"):\n")
    to.print <- round(ll,2)
    print.result(print, summary.file, to.cat=to.cat, to.print=to.print) 
    #    print.result(print, summary.file, to.cat=plain.line, to.print=NULL) 
    
  }
  
  mlx.saveProject(final.project)
  
  #-----------------------------------------
  if (!is.null(covToTransform)) {
    covariate <- mlx.getCovariateInformation()$covariate
    for (cov.name in covToTransform) {
      covkj <- covariate[[cov.name]]
      lckj <- paste0("logt",toupper(substr(cov.name, 1, 1)), substr(cov.name, 2, nchar(cov.name)))
      if (!(lckj %in% mlx.getCovariateInformation()$name)) {
        tr.str <- paste0(lckj,' = "log(',cov.name,"/",signif(mean(covkj),digits=2),')"')
        trs <- paste0("lixoftConnectors::addContinuousTransformedCovariate(",tr.str,")")
        eval(parse(text=trs))
        #        colnames(weight$covariate) <- gsub(cov.name, lckj, colnames(weight$covariate))
        foo <- colnames(weight$covariate)
        weight$covariate <- cbind(weight$covariate, weight$covariate[,cov.name])
        colnames(weight$covariate) <- c(foo, lckj)
      }
      covToTest <- c(covToTest, lckj)
      
    }
    covToTest <- unique(setdiff(covToTest, covToTransform))
    covToTransform <- NULL
    mlx.saveProject(final.project)
    
    if (!mlx.getLaunchedTasks()[["populationParameterEstimation"]]) {
      to.cat <- "\nEstimation of the population parameters using the transformed covariates ... \n"
      print.result(print, summary.file, to.cat=to.cat, to.print=NULL) 
      mlx.runPopulationParameterEstimation()
    }
    if (!mlx.getLaunchedTasks()[["conditionalDistributionSampling"]]) {
      to.cat <- "Sampling of the conditional distribution using the the transformed covariates ... \n"
      print.result(print, summary.file, to.cat=to.cat, to.print=NULL) 
      mlx.runConditionalDistributionSampling()
    }
  }
  
  #--------------
  
  stop.test <- FALSE
  corr.test <- FALSE
  iter <- 0
  
  cov.names0 <- cov.names <- NULL
  
  if (identical(covToTest,"all"))
    covFix = NULL
  else
    covFix <- setdiff(mlx.getCovariateInformation()$name, covToTest)
  
  
  if (iop.ll) {
    ll.prev <- Inf
    ll.new <- ll.ini
  }
  sp0 <- NULL
  cov.test <- NULL
  e <- NULL
  while (!stop.test) {
    iter <- iter + 1
    
    if (iter==1){
      weight.covariate <- weight$covariate*coef.w1
      weight.correlation <- weight$correlation*coef.w1
    } else {
      weight.covariate <- weight$covariate
      weight.correlation <- weight$correlation
    }
    
    to.cat <- paste0(plain.line,"- - - Iteration ",iter," - - -\n")
    print.result(print, summary.file, to.cat=to.cat, to.print=NULL) 
    
    obs.dist0 <- obs.dist
    error.model0 <- error.model
    covariate.model0 <- covariate.model
    correlation.model0 <- correlation.model
    if (iop.ll)
      ll0 <- ll
    
    if (model$residualError) {
      res.error <- errorModelSelection(pen.coef=pen.coef[-c(1, 2)], nb.model=nb.model, f.min=fError.min)
      if (nb.model==1)
        error.model <- res.error
      else {
        error.model <- mlx.getContinuousObservationModel()$errorModel
        for (k in (1:length(error.model)))
          error.model[[k]] <- as.character(res.error[[k]]$error.model[1])
      }
    }
    
    pmax.cov <-  p.max
    if (model$covariate) {
      #      pmax.cov <- ifelse(iter <= 1, 1, p.max) 
      res.covariate <- covariateModelSelection(pen.coef=pen.coef[1], nb.model=nb.model, weight=weight.covariate,
                                               covToTransform=covToTransform, covFix=covFix, direction=direction, 
                                               steps=steps, p.max=pmax.cov, paramToUse=paramToUse, sp0=sp0, iter=iter,
                                               correlation.model = correlation.model, n.full=n.full, eta=e)
      res.covariate$res <- sortCov(res.covariate$res, cov.ini)
      if (iter>explor.iter) sp0 <- res.covariate$sp
      covToTransform <- setdiff(covToTransform, res.covariate$tr0)
      covariate.model <- res.covariate$model
      e <- res.covariate$residuals
      
      if (nb.model==1)
        cov.select <- rownames(res.covariate$res)
      else
        cov.select <- rownames(res.covariate$res[[1]])
      
      cov.names <- lapply(covariate.model[cov.select], function(x) {sort(names(which(x)))})
      cov.names0 <- lapply(covariate.model0[cov.select], function(x) {sort(names(which(x)))})
      
      cov.test <- covariate.test(cov.test, covToTest, covToTransform, paramToUse)
      
    } else {
      e <- mlx.getSimulatedRandomEffects()
    }
    
    if (model$correlation & !corr.test) {
      if (isTRUE(all.equal(cov.names0,cov.names))) # & isTRUE(all.equal(error.model0,error.model))) 
        corr.test <- TRUE
      if (!seq.cov & iter>seq.cov.iter)
        corr.test <- TRUE
      
      if (corr.test & (seq.cov==T | seq.cov.iter>0)) {
        to.cat <- "Start building correlation model too ... \n"
        print.result(print, summary.file, to.cat=to.cat, to.print=NULL) 
      }
    }
    if (model$correlation & corr.test) {
      pen.corr <- ifelse(iter <= 1, pen.coef[1], pen.coef[1]) 
      res.correlation <- correlationModelSelection(e0=e, pen.coef=pen.corr, nb.model=nb.model, 
                                                   corr0=correlation.model0, seqmod=seq.corr, weight=weight.correlation)
      if (nb.model==1) 
        correlation.model <- res.correlation
      else  
        correlation.model <- res.correlation[[1]]$block
      if (length(correlation.model)==0)
        correlation.model <- NULL
    } else {
      res.correlation <- lapply(mlx.getIndividualParameterModel()$correlationBlocks$id, sort)
    }
    if (length(res.correlation)==0)
      res.correlation <- NULL
    
    
    #-------------------------------
    
    if (max.iter>0 | nb.model>1) {
      eq.cov <- isTRUE(all.equal(cov.names0,cov.names)) & (pmax.cov == p.max)
      eq.err <-  isTRUE(all.equal(error.model0,error.model))
      eq.corr <- isTRUE(all.equal(correlation.model0,correlation.model))
      eq.dist <- isTRUE(all.equal(obs.dist0,obs.dist))
      if (!model$correlation | corr.test) {
        if ( eq.cov & eq.err & eq.dist & eq.corr ) 
          stop.test <- TRUE
        if ( model$covariate & eq.cov & eq.dist & eq.corr ) 
          stop.test <- TRUE
      }
      if (stop.test) {
        to.cat <- "\nNo difference between two successive iterations\n"
        print.result(print, summary.file, to.cat=to.cat, to.print=NULL) 
      }
      
      if (!stop.test | nb.model>1) {
        if (model$covariate) {
          to.cat <- "\nCovariate model:\n"
          to.print <- res.covariate$res
          print.result(print, summary.file, to.cat=to.cat, to.print=to.print) 
          
        }
        if (model$correlation) {
          to.cat <- "\nCorrelation model:\n"
          if (!is.null(res.correlation)) 
            to.print <- res.correlation
          else
            to.print <- "NULL"
          print.result(print, summary.file, to.cat=to.cat, to.print=to.print) 
        }
        if (model$residualError) {
          to.cat <- "\nResidual error model:\n"
          to.print <- res.error
          print.result(print, summary.file, to.cat=to.cat, to.print=to.print) 
        }
        
      }
    } 
    
    if (!stop.test) {
      p.est <- mlx.getEstimatedPopulationParameters()
      mlx.setInitialEstimatesToLastEstimates(fixedEffectsOnly = T)
      p.ini <- mlx.getPopulationParameterInformation()
      rownames(p.ini) <- p.ini$name
      i.omega <- which(grepl("omega_", p.ini$name) & (!identical(p.ini$method,"FIXED")))
      p.ini$initialValue[i.omega] <- p.est[p.ini$name[i.omega]]*3
      #          p.ini[omega,] <- omega.ini
      jcor <- grep("corr_",p.ini$name)
      if (length(jcor)>0)  p.ini <- p.ini[-jcor,]
      mlx.setPopulationParameterInformation(p.ini)
      
      if (model$residualError) {
        emodel <- error.model
        odist <- mlx.getContinuousObservationModel()$distribution
        for (k in (1:length(emodel))) {
          if (identical(emodel[[k]],"exponential")) {
            emodel[[k]] <- "constant"
            odist[[k]] <- "lognormal"
          } else {
            odist[[k]] <- "normal"
          }
        }
        mlx.setErrorModel(emodel)
        mlx.setObservationDistribution(odist)
      }
      
      if (model$covariate) {
        if (length(res.covariate$add.covariate) >0) {
          for (k in 1:length(res.covariate$add.covariate))
            eval(parse(text=res.covariate$add.covariate[[k]]))
        }
        mlx.setCovariateModel (covariate.model)
      }
      
      if (model$correlation & corr.test)
        mlx.setCorrelationBlocks(correlation.model)
      
      #-------------------------------
      
      if (max.iter>0) {
        if (iop.ll) { 
          if (ll.new > ll.prev) {
            g <- mlx.getGeneralSettings()
            g$autochains <- FALSE
            g$nbchains <- g$nbchains+1
            mlx.setGeneralSettings(g)
          }
          ll.prev <- ll.new
        }
        g=mlx.getConditionalDistributionSamplingSettings()
        g$nbminiterations <- max(100, g$nbminiterations)
        
        #    if (!is.null(g$enablemaxiterations)) {
        # g$enablemaxiterations <- T
        # g$nbmaxiterations <- 500
        #     g$nbminiterations <- 150
        #     g$ratio=0.03
        # }
        
        mlx.setConditionalDistributionSamplingSettings(g)
        
        # g <- getPopulationParameterEstimationSettings()
        # g$nbexploratoryiterations <-  400
        # g$nbsmoothingiterations <- 200
        # g$smoothingautostop <- g$exploratoryautostop <- F
        # g$simulatedannealing <- T
        # setPopulationParameterEstimationSettings(g)
        
      }
      
      mlx.saveProject(final.project)
      if (dir.exists(final.dir))
        unlink(final.dir, recursive=TRUE)
      mlx.loadProject(final.project)
      
      if (max.iter>0) {
        
        to.cat <- paste0("\nRun scenario for model ",iter," ... \nEstimation of the population parameters... \n")
        print.result(print, summary.file, to.cat=to.cat, to.print=NULL) 
        
        #        mlx.runPopulationParameterEstimation(parameter=gi)
        mlx.runPopulationParameterEstimation()
        to.cat <- "Sampling from the conditional distribution... \n"
        print.result(print, summary.file, to.cat=to.cat, to.print=NULL) 
        
        mlx.runConditionalDistributionSampling()
        gi <- mlx.getSimulatedIndividualParameters()
        gi <- gi %>% filter(rep==gi$rep[nrow(gi)]) %>% select(-rep)
        if (iop.ll) {
          to.cat <- "Estimation of the log-likelihood... \n"
          print.result(print, summary.file, to.cat=to.cat, to.print=NULL) 
          if (lin.ll & !launched.tasks[["conditionalModeEstimation"]])
            mlx.runConditionalModeEstimation()
          mlx.runLogLikelihoodEstimation(linearization = lin.ll)
          ll.new <- compute.criterion(criterion, method.ll, weight, pen.coef)
          list.criterion <- c(list.criterion, ll.new)
        }
        
        buildmlx.project.iter <- file.path(buildmlx.dir,paste0("iteration",iter,".mlxtran"))
        mlx.saveProject(buildmlx.project.iter)
        
        if (iop.ll) {
          if (stop.test)
            ll <- ll0
          else {
            ll <- formatLL(mlx.getEstimatedLogLikelihood()[[method.ll]], criterion, ll.new, is.weight)
          }
          
          to.cat <- paste0("\nEstimated criteria (",method.ll,"):\n")
          to.print <- round(ll,2)
          print.result(print, summary.file, to.cat=to.cat, to.print=to.print) 
        }
        if (iter >= max.iter) {
          stop.test <- TRUE
          to.cat <- "Maximum number of iterations reached\n"
          print.result(print, summary.file, to.cat=to.cat, to.print=NULL) 
        }
      }
    }
    if (max.iter==0) 
      stop.test <- TRUE
    mlx.saveProject(final.project)
  }
  
  change.error.model <- NULL
  if (iop.ll) {
    ll.min <- min(list.criterion)
    iter.opt <- which.min(list.criterion)
    if (iter.opt==1)
      buildmlx.project.iter <- project
    else {
      buildmlx.project.iter <- file.path(buildmlx.dir,paste0("iteration",iter.opt-1,".mlxtran"))
    }
    mlx.loadProject(buildmlx.project.iter)
    mlx.saveProject(final.project)
  } else {
    iter.opt <- iter
  }
  
  if (test) {
    
    if (!iop.ll) {
      g <- as.list(mlx.getLaunchedTasks())$logLikelihoodEstimation
      if (!linearization & !("importanceSampling" %in% g)) 
        mlx.runLogLikelihoodEstimation(linearization=F)
      if (linearization & !("linearization" %in% g)) 
        mlx.runLogLikelihoodEstimation(linearization=T)
      ll.min <- compute.criterion(criterion, method.ll, weight, pen.coef)
    }
    
    to.cat <- paste0(plain.line,"- - - Further tests - - -\n")
    print.result(print, summary.file, to.cat=to.cat, to.print=NULL) 
    if (model$covariate) {
      g0 <- mlx.getIndividualParameterModel()
      covariate <- random.effect <- p.ttest <- p.lrt <- in.model <- p.value <- NULL
      r.test <- covariate.test(cov.test, covToTest, covToTransform, paramToUse)
      r.test <- r.test %>% filter(!in.model) %>% select(-in.model) 
      if (is.weight) {
        w.cov <- weight$covariate[cbind(r.test[['parameter']], r.test[['covariate']])]
        r.test <- r.test %>%  mutate(p.value = p.weight(p.value, w.cov, pen.coef[1]))
      }
      r.cov0 <- res.covariate$r.cov0
      for (j in 1:nrow(r.test)) {
        pj <- r.test$parameter[j]
        if (!is.null(r.cov0[[pj]])) {
          if (r.test$covariate[j] %in% r.cov0[[pj]])
            r.test$p.value[j] <- 1
        }
      }
      i.min <- which(as.numeric(r.test$p.value) < p.min[1])
      if (length(i.min)>0) {
        
        to.cat <- paste0(plain.short,"Add parameters/covariates relationships:\n")
        to.print <- r.test[i.min,]
        print.result(print, summary.file, to.cat=to.cat, to.print=to.print) 
        g <- mlx.getIndividualParameterModel()
        stop.test <- F
        for (i in i.min) {
          param.i <- r.test$parameter[i]
          cov.i <- r.test$covariate[i]
          g$covariateModel[[param.i]][cov.i] <- T
        }
        mlx.setIndividualParameterModel(g)
        iter <- iter+1
        to.cat <- paste0("\nRun scenario for model ",iter," ... \nEstimation of the population parameters... \n")
        print.result(print, summary.file, to.cat=to.cat, to.print=NULL) 
        buildmlx.project.iter <- file.path(buildmlx.dir,paste0("iteration",iter,".mlxtran"))
        mlx.saveProject(buildmlx.project.iter)
        mlx.runPopulationParameterEstimation()
      } else {
        stop.test <- F
      }
      if (any(mlx.getObservationInformation()$type != "continuous"))
        mlx.runStandardErrorEstimation(linearization=F)
      else
        mlx.runStandardErrorEstimation(linearization=T)
      #mlx.runStandardErrorEstimation(linearization = lin.ll)
      r.test <- mlx.getTests()$wald
      
      # if (is.weight | is.prior) {
      #   w.cov <- weight$covariate[cbind(r.test[['parameter']], r.test[['covariate']])]
      #   r.test <- r.test %>%  mutate(p.value = p.weight(p.value, w.cov, pen.coef[1]))
      # }
      g <- mlx.getIndividualParameterModel()
      n.param <- g$name
      n.cov <- names(g$covariateModel[[1]])
      
      pv <- as.numeric(gsub("<", "", r.test$p.value))
      pv[which(is.nan(pv))] <- 0.99999
      list.ipc <- NULL
      for (np in n.param) {
        gp <- g$covariateModel[[np]]
        ngp <- names(which(gp))
        if (length(ngp) > 0) {
          for (nc in ngp) {
            g$covariateModel[[np]][nc] <- F
            ipc <- grep(paste0("beta_",np,"_",nc), r.test$parameter)
            pv[ipc] <- p.weight(pv[ipc], weight$covariate[np, nc], pen.coef[1])
            if (min(pv[ipc]) < p.min[2])
              g$covariateModel[[np]][nc] <- T
            else
              list.ipc <- c(list.ipc, ipc)
          }
        }
      }
      
      if (identical(g$covariateModel, g0$covariateModel)) 
        stop.test <- T
      
      if (length(list.ipc) >0) {
        to.cat <- paste0(plain.short,"Remove parameters/covariates relationships:\n")
        method <- statistics <- parameter <- NULL
        to.print <- (r.test %>% select(-c(method, statistics)) %>%
                       rename(coefficient=parameter))[list.ipc,]
        print.result(print, summary.file, to.cat=to.cat, to.print=to.print) 
      }
      
      if (!stop.test) {
        if (length(list.ipc) >0) {
          mlx.setIndividualParameterModel(g)
          iter <- iter+1
          to.cat <- paste0("\nRun scenario for model ",iter," ... \nEstimation of the population parameters... \n")
          print.result(print, summary.file, to.cat=to.cat, to.print=NULL) 
          buildmlx.project.iter <- file.path(buildmlx.dir,paste0("iteration",iter,".mlxtran"))
          mlx.saveProject(buildmlx.project.iter)
          mlx.runPopulationParameterEstimation()
        }
        if (lin.ll) {
          if(!launched.tasks[["conditionalModeEstimation"]]) 
            mlx.runConditionalModeEstimation()
        } else {
          to.cat <- "Sampling from the conditional distribution... \n"
          print.result(print, summary.file, to.cat=to.cat, to.print=NULL) 
          mlx.runConditionalDistributionSampling()
        }
        to.cat <- "Estimation of the log-likelihood... \n"
        print.result(print, summary.file, to.cat=to.cat, to.print=NULL) 
        mlx.runLogLikelihoodEstimation(linearization = lin.ll)
        ll.new <- compute.criterion(criterion, method.ll, weight, pen.coef)
        ll <- formatLL(mlx.getEstimatedLogLikelihood()[[method.ll]], criterion, ll.new, is.weight, is.prior)
        to.cat <- paste0("\nEstimated criteria (",method.ll,"):\n")
        to.print <- round(ll,2)
        print.result(print, summary.file, to.cat=to.cat, to.print=to.print) 
        
        if (ll.new < ll.min) {
          ll.min <- ll.new
          mlx.saveProject(final.project)
        } else {
          mlx.loadProject(final.project)
        }
      }
      
      # ---   Wald tests
      g <- as.list(mlx.getLaunchedTasks())$standardErrorEstimation
      if (!linearization & !("stochasticApproximation" %in% g))
        mlx.runStandardErrorEstimation(linearization=F)
      if (!("linearization" %in% g))
        mlx.runStandardErrorEstimation(linearization=T)
      
      r.test <- mlx.getTests()$wald
      g <- mlx.getIndividualParameterModel()
      n.param <- g$name
      n.cov <- names(g$covariateModel[[1]])
      
      pv <- as.numeric(gsub("<", "", r.test$p.value))
      pv[which(is.nan(pv))] <- 0
      
      list.ipc <- NULL
      for (np in n.param) {
        gp <- g$covariateModel[[np]]
        ngp <- names(which(gp))
        if (length(ngp) > 0) {
          for (nc in ngp) {
            ipc <- grep(paste0("beta_",np,"_",nc), r.test$parameter)
            pv[ipc] <- p.weight(pv[ipc], weight$covariate[np, nc], pen.coef[1])
            if (max(pv[ipc]) > p.min[2]) {
              g$covariateModel[[np]][nc] <- F
              list.ipc <- c(list.ipc, ipc)
            }
          }
        }
      }
      if (length(list.ipc) >0) {
        to.cat <- paste0(plain.short,"Remove parameters/covariates relationships:\n")
        method <- statistics <- parameter <- NULL
        to.print <- (r.test %>% select(-c(method, statistics)) %>%
                       rename(coefficient=parameter))[list.ipc,]
        print.result(print, summary.file, to.cat=to.cat, to.print=to.print) 
        
        mlx.setIndividualParameterModel(g)
        iter <- iter+1
        to.cat <- paste0("\nRun scenario for model ",iter," ... \nEstimation of the population parameters... \n")
        print.result(print, summary.file, to.cat=to.cat, to.print=NULL) 
        buildmlx.project.iter <- file.path(buildmlx.dir,paste0("iteration",iter,".mlxtran"))
        mlx.saveProject(buildmlx.project.iter)
        mlx.runPopulationParameterEstimation()
        if (lin.ll) {
          if(!launched.tasks[["conditionalModeEstimation"]]) 
            mlx.runConditionalModeEstimation()
        } else {
          to.cat <- "Sampling from the conditional distribution... \n"
          print.result(print, summary.file, to.cat=to.cat, to.print=NULL) 
          mlx.runConditionalDistributionSampling()
        }
        to.cat <- "Estimation of the log-likelihood... \n"
        print.result(print, summary.file, to.cat=to.cat, to.print=NULL) 
        mlx.runLogLikelihoodEstimation(linearization = lin.ll)
        ll.new <- compute.criterion(criterion, method.ll, weight, pen.coef)
        ll <- formatLL(mlx.getEstimatedLogLikelihood()[[method.ll]], criterion, ll.new, is.weight, is.prior)
        to.cat <- paste0("\nEstimated criteria (",method.ll,"):\n")
        to.print <- round(ll,2)
        print.result(print, summary.file, to.cat=to.cat, to.print=to.print) 
        
        if (ll.new < ll.min) {
          ll.min <- ll.new
          mlx.saveProject(final.project)
        }
      }
      
      # ---  end Wald tests
      
    }
    if (model$correlation) {
      test.cor <- T
      while (test.cor) {
        mlx.loadProject(final.project)
        p.cortest <- NULL
        r.test <- correlationTest()$p.value %>% filter(!in.model) %>% rename(p.value=p.cortest)
        param1 <- gsub("eta_","",r.test$randomEffect.1)
        param2 <- gsub("eta_","",r.test$randomEffect.2)
        w.cor <- weight$correlation[cbind(param1, param2)]+weight$correlation[cbind(param2, param1)]
        r.test <- r.test %>% mutate(p.value = p.weight(p.value, w.cor, pen.coef[1]))
        
        i.min <- which(as.numeric(r.test$p.value) < p.min[3])
        g <- mlx.getIndividualParameterModel()
        param.list <- unlist(g$correlationBlocks$id)
        if (length(i.min)>0) {
          i.min <- i.min[which.min(r.test$p.value[i.min])]
          param1 <- gsub("eta_","",r.test$randomEffect.1[i.min])
          param2 <- gsub("eta_","",r.test$randomEffect.2[i.min])
          test.cor <- F
          if ( !(param1 %in% param.list) & !(param2 %in% param.list) ) {
            l.block <- length(g$correlationBlocks$id)+1
            g$correlationBlocks$id[[l.block]] <- c(param1, param2)
            test.cor <- T
          }
          if ( !(param1 %in% param.list) & (param2 %in% param.list) ) {
            l.block <- grep(param2, g$correlationBlocks$id)
            g$correlationBlocks$id[[l.block]] <- c(g$correlationBlocks$id[[l.block]], param1)
            test.cor <- T
          }
          if ( (param1 %in% param.list) & !(param2 %in% param.list) ) {
            l.block <- grep(param1, g$correlationBlocks$id)
            g$correlationBlocks$id[[l.block]] <- c(g$correlationBlocks$id[[l.block]], param2)
            test.cor <- T
          }
          if (test.cor) {
            
            to.cat <- paste0(plain.short, "Add correlation:\n")
            to.print <- r.test[i.min,]
            print.result(print, summary.file, to.cat=to.cat, to.print=to.print) 
            to.print <- g$correlationBlocks$id
            print.result(print, summary.file, to.cat=NULL, to.print=to.print) 
            gi <- mlx.getPopulationParameterInformation()
            gi$initialValue[which(gi$name==paste0("omega_",param1))] <- 3*gi$initialValue[which(gi$name==paste0("omega_",param1))]
            gi$initialValue[which(gi$name==paste0("omega_",param2))] <- 3*gi$initialValue[which(gi$name==paste0("omega_",param2))]
            mlx.setPopulationParameterInformation(gi)
            mlx.setIndividualParameterModel(g)
            
            ###
            iter <- iter+1
            buildmlx.project.iter <- file.path(buildmlx.dir,paste0("iteration",iter,".mlxtran"))
            mlx.saveProject(buildmlx.project.iter)
            to.cat <- paste0("\nRun scenario for model ",iter,"  ... \nEstimation of the population parameters... \n")
            print.result(print, summary.file, to.cat=to.cat, to.print=NULL) 
            mlx.runPopulationParameterEstimation()
            if (lin.ll) {
              if(!launched.tasks[["conditionalModeEstimation"]]) 
                mlx.runConditionalModeEstimation()
            } else {
              to.cat <- "Sampling from the conditional distribution... \n"
              print.result(print, summary.file, to.cat=to.cat, to.print=NULL) 
              mlx.runConditionalDistributionSampling()
            }
            to.cat <- "Estimation of the log-likelihood... \n"
            print.result(print, summary.file, to.cat=to.cat, to.print=NULL) 
            mlx.runLogLikelihoodEstimation(linearization = lin.ll)
            ll.new <- compute.criterion(criterion, method.ll, weight, pen.coef)
            ll.disp <- formatLL(mlx.getEstimatedLogLikelihood()[[method.ll]], criterion, ll.new, is.weight, is.prior)
            
            to.cat <- paste0("\nEstimated criteria (",method.ll,"):\n")
            to.print <- round(ll.disp,2)
            print.result(print, summary.file, to.cat=to.cat, to.print=to.print) 
            
            if (ll.new < ll.min) {
              ll.min <- ll.new
              mlx.saveProject(final.project)
            } else {
              test.cor <- F
            }
            
          } else {
            test.cor <- F
          }
        } else {
          test.cor <- F
        }
      }
      
      mlx.loadProject(final.project)
      p <- mlx.getEstimatedPopulationParameters()
      if (any(mlx.getObservationInformation()$type != "continuous")) {
        mlx.runStandardErrorEstimation(linearization=F)
        se <- mlx.getEstimatedStandardErrors()$stochasticApproximation$se
        names(se) <- mlx.getEstimatedStandardErrors()$stochasticApproximation$parameter
      } else {
        mlx.runStandardErrorEstimation(linearization=T)
        se <- mlx.getEstimatedStandardErrors()$linearization$se
        names(se) <- mlx.getEstimatedStandardErrors()$linearization$parameter
      }
      z <- as.numeric(p)/as.numeric(se[names(p)])
      names(z) <- names(p)
      pv <- pnorm(-abs(z))*2
      pv.corr <- pv[grep("corr_", names(pv))]
      if (length(which(pv.corr>p.min[2])) > 0) {
        mlx.saveProject(buildmlx.project.iter)
        ind <- mlx.getEstimatedIndividualParameters()$saem
        pv.block <- strsplit(gsub("corr_", "", names(pv.corr)),"_")
        ind.mod <- mlx.getIndividualParameterModel()
        cb <- ind.mod$correlationBlocks$id
        cl <- unlist(cb)
        pv.tot <- rep(0,length(cl))
        names(pv.tot) <- cl
        for (k in 1: length(pv.corr))
          pv.tot[pv.block[[k]]] <- pv.tot[pv.block[[k]]] + log(pv.corr[k])
        param0 <- names(which.max(pv.tot))
        cb <- lapply(cb, function(x) setdiff(x, param0))
        i.cb <- which(lapply(cb, length)==1)
        if (length(i.cb)>0)
          cb[[i.cb]] <- NULL
        ind.mod$correlationBlocks$id <- cb
        to.cat <- paste0(plain.short,"Test correlation model:\n")
        to.print <- cb
        print.result(print, summary.file, to.cat=to.cat, to.print=to.print) 
        
        mlx.setInitialEstimatesToLastEstimates(fixedEffectsOnly = TRUE)
        mlx.setIndividualParameterModel(ind.mod)
        iter <- iter+1
        buildmlx.project.iter <- file.path(buildmlx.dir,paste0("iteration",iter,".mlxtran"))
        mlx.saveProject(buildmlx.project.iter)
        to.cat <- paste0("Run scenario for model ",iter,"  ... \nEstimation of the population parameters... \n")
        print.result(print, summary.file, to.cat=to.cat, to.print=NULL) 
        mlx.runPopulationParameterEstimation(parameters=ind)
        if (lin.ll) {
          if(!launched.tasks[["conditionalModeEstimation"]]) 
            mlx.runConditionalModeEstimation()
        } else {
          to.cat <- "Sampling from the conditional distribution... \n"
          print.result(print, summary.file, to.cat=to.cat, to.print=NULL) 
          mlx.runConditionalDistributionSampling()
        }
        to.cat <- "Estimation of the log-likelihood... \n"
        print.result(print, summary.file, to.cat=to.cat, to.print=NULL) 
        mlx.runLogLikelihoodEstimation(linearization = lin.ll)
        ll.new <- compute.criterion(criterion, method.ll, weight, pen.coef)
        ll <- formatLL(mlx.getEstimatedLogLikelihood()[[method.ll]], criterion, ll.new, is.weight, is.prior)
        to.cat <- paste0("\nEstimated criteria (",method.ll,"):\n")
        to.print <- round(ll,2)
        print.result(print, summary.file, to.cat=to.cat, to.print=to.print) 
        if (ll.new < ll.min) {
          ll.min <- ll.new
          mlx.saveProject(final.project)
        }
      }
    }
    
  } 
  
  if (model$covariate & nb.model>1)
    res.covariate$res <- sortCov(res.covariate$res[[1]], cov.ini)
  
  mlx.loadProject(final.project)
  
  if (model$covariate)
    covariate.model.print <- formatCovariateModel(mlx.getIndividualParameterModel()$covariateModel)
  if (model$correlation) {
    correlation.model.print <- lapply(mlx.getIndividualParameterModel()$correlationBlocks$id, sort)
    if (length(correlation.model.print)==0)
      correlation.model.print <- NULL
  }
  if (model$residualError)
    error.model.print <- formatErrorModel(mlx.getContinuousObservationModel()$errorModel)
  if (iop.ll) {
    ll.final <- compute.criterion(criterion, method.ll, weight, pen.coef)
    ll <- formatLL(mlx.getEstimatedLogLikelihood()[[method.ll]], criterion, ll.final, is.weight, is.prior)
  }
  
  to.cat <- paste0(plain.line,"\nFinal statistical model:\n")
  print.result(print, summary.file, to.cat=to.cat, to.print=NULL) 
  
  if (model$covariate) {
    to.cat <- "\nCovariate model:\n"
    to.print <- covariate.model.print
    print.result(print, summary.file, to.cat=to.cat, to.print=to.print) 
  }
  if (model$correlation) {
    to.cat <- "\nCorrelation model:\n"
    if (!is.null(correlation.model.print))
      to.print <- correlation.model.print
    else
      to.print <- "NULL"
    print.result(print, summary.file, to.cat=to.cat, to.print=to.print) 
  }
  if (model$residualError) {
    to.cat <- "\nResidual error model:\n"
    to.print <- error.model.print
    print.result(print, summary.file, to.cat=to.cat, to.print=to.print) 
  }
  if (iop.ll & max.iter>0) {
    to.cat <- paste0("\nEstimated criteria (",method.ll,"):\n")
    to.print <- round(ll,2)
    print.result(print, summary.file, to.cat=to.cat, to.print=to.print) 
  }
  
  test.del <- FALSE
  if (model$covariate & center.covariate) {
    foo <- lapply(res.covariate$model,function(x) {which(x)})
    cov.model <- unique(unlist(lapply(foo,function(x) {names(x)})))
    cov.type <- mlx.getCovariateInformation()$type[cov.model]
    cov.cont <- names(cov.type[cov.type=="continuous"])
    for (ck in cov.cont) {
      cck <- paste0("c",ck)
      covk <- mlx.getCovariateInformation()$covariate[[ck]]
      tr.str <- paste0(cck,' = "',ck,"-",signif(mean(covk),digits=2),'"')
      tr.str <- paste0("lixoftConnectors::addContinuousTransformedCovariate(",tr.str,")")
      eval(parse(text=tr.str))
      g=mlx.getIndividualParameterModel()$covariateModel
      cg <- lapply(g, function(x) {foo <- x[ck]; x[ck]<-x[cck]; x[cck]=foo; return(x)})
      mlx.setCovariateModel (cg)
      test.del <- TRUE
      covariate.model <- cg
    }
  }
  if (test.del & dir.exists(final.dir)) {
    unlink(final.dir, recursive=TRUE)
  }
  
  g=mlx.getScenario()
  g$tasks[[2]]=TRUE
  mlx.setScenario(g)
  mlx.saveProject(final.project)
  mlx.loadProject(final.project)
  
  error.model <- mlx.getContinuousObservationModel()$errorModel
  covariate.model <- mlx.getIndividualParameterModel()$covariateModel
  correlation.model <- lapply(mlx.getIndividualParameterModel()$correlationBlocks$id, sort)
  if (length(correlation.model)==0)
    correlation.model <- NULL
  
  dt <- proc.time() - ptm
  res <- list(project=final.project, niter=iter.opt, time=dt["elapsed"])
  if (model$covariate)
    res <- c(res, list(covariate.model=covariate.model))
  if (model$correlation)
    res <- c(res, list(correlation.model=correlation.model))
  if (model$residualError)
    res <- c(res, list(error.model=error.model))
  
  to.cat <- paste0("\ntotal time: ", round(dt["elapsed"], digits=1),"s\n", plain.line)
  print.result(print, summary.file, to.cat=to.cat, to.print=NULL) 
  
  
  res$change <- !(identical(error.model,error.model.ini) & 
                    identical(covariate.model,covariate.model.ini) &
                    identical(correlation.model, correlation.model.ini))
  res$change.error.model <- change.error.model
  res$weight <- weight
  options(op.original)
  return(res)
}

buildmlx.check <- function(project, final.project, model, paramToUse, covToTest, covToTransform, center.covariate, 
                           criterion, linearization, ll, test, direction, steps, max.iter, explor.iter, 
                           seq.cov, seq.cov.iter, seq.corr, p.max, p.min, print, nb.model, prior, weight, n.full) {
  
  if (length(mlx.getIndividualParameterModel()$variability)>1)
    stop("Multiple levels of variability are not supported in this version of buildmlx", call.=FALSE)
  
  if (is.null(final.project)) 
    final.project <- paste0(sub(pattern = "(.*)\\..*$", replacement = "\\1", project),"_built.mlxtran")
  if (!grepl("\\.",final.project))
    final.project <- paste0(final.project,".mlxtran")
  if (!grepl("\\.mlxtran",final.project))
    stop(paste0(final.project, " is not a valid name for a Monolix project (use the .mlxtran extension)"), call.=FALSE)
  
  if (!is.logical(test))
    stop(" 'test' should be boolean", call.=FALSE)
  if (!is.logical(center.covariate))
    stop(" 'center.covariate' should be boolean", call.=FALSE)
  if (!is.logical(linearization))
    stop(" 'linearization' should be boolean", call.=FALSE)
  if (!is.logical(ll))
    stop(" 'll' should be boolean", call.=FALSE)
  if (!is.logical(print))
    stop(" 'print' should be boolean", call.=FALSE)
  if (!is.numeric(criterion) && !(criterion %in% c("AIC", "BIC", "BICc")))
    stop(" 'criterion' should be in {'AIC', 'BIC', 'BICc'} or be numerical > 0", call.=FALSE)
  if (is.numeric(criterion) && criterion<=0)
    stop(" 'criterion' should be in {'AIC', 'BIC', 'BICc'} or be numerical > 0", call.=FALSE)
  if (!is.numeric(steps) | steps<=0)
    stop(" 'steps' should be numerical > 0", call.=FALSE)
  if ((round(max.iter)!=max.iter) | max.iter<=0)
    stop(" 'max.iter' should be an integer > 0", call.=FALSE)
  if ((round(explor.iter)!=explor.iter) | explor.iter<=0)
    stop(" 'explor.iter' should be an integer > 0", call.=FALSE)
  if ((round(n.full)!=n.full) | n.full<=0)
    stop(" 'n.full' should be an integer > 0", call.=FALSE)
  if (!is.logical(seq.cov))
    stop(" 'seq.cov' should be boolean", call.=FALSE)
  if ((round(seq.cov.iter)!=seq.cov.iter) | seq.cov.iter<0)
    stop(" 'seq.cov.iter' should be an integer >= 0", call.=FALSE)
  if (!is.logical(seq.corr))
    stop(" 'seq.corr' should be boolean", call.=FALSE)
  if (!is.numeric(p.max) | p.max<=0 | p.max>1)
    stop(" 'p.max' should be probability > 0", call.=FALSE)
  if (!is.numeric(p.min) | (length(p.min) != 3))
    stop(" 'p.min' should be a 3-vector of probabilities", call.=FALSE)
  if (min(p.min)<=0 | max(p.min) >1)
    stop(" 'p.min' should be a 3-vector of probabilities", call.=FALSE)
  if ((round(nb.model)!=nb.model) | nb.model<=0)
    stop(" 'nb.model' should be an integer > 0", call.=FALSE)
  if (!is.logical(print))
    stop(" 'print' should be boolean", call.=FALSE)
  
  
  cov.info <- mlx.getCovariateInformation()
  cov.names <- cov.info$name
  cov.types <- cov.info$type
  j.trans <- grep("transformed",cov.types)
  j.cont <- grep("continuous",cov.types)
  cont.cov <- cov.names[j.cont]
  # j.strat <- grep("stratification",cov.types)
  # strat.cov <- cov.names[j.strat]
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
  
  #covToTest <- setdiff(covToTest, strat.cov)
  
  covToTransform <- intersect(covToTest, covToTransform)
  if (length(covToTransform)==0) 
    covToTransform <- NULL
  
  ind.dist <- mlx.getIndividualParameterModel()$distribution
  param.names <- names(ind.dist)
  if (identical(paramToUse,"all"))
    paramToUse <- param.names
  nparam0 <- paramToUse[(!(paramToUse %in% param.names))]
  if (length(nparam0)>0) 
    stop(paste0(nparam0, " is not a valid parameter"), call.=FALSE)
  
  model.names <- c("residualError", "covariate", "correlation")
  if (identical(model,"all")) model <- model.names
  # if (is.null(cov.names))  model <- setdiff(model, "covariate")
  # if (sum(mlx.getIndividualParameterModel()$variability$id) < 2) model <- setdiff(model, "correlation")
  mod0 <- model[(!(model %in% c(model.names, "variance")))]
  if (length(mod0)>0) 
    stop(paste0(mod0, " is not a valid model component"), call.=FALSE)
  foo <- model
  model <- list(F, F, F)
  names(model) <- model.names
  model[foo] <- T
  
  #------------------------------
  if (!any(mlx.getIndividualParameterModel()$variability$id))
    stop("\nA least one parameter with random effects is required\n", call.=FALSE)
  if (is.null(mlx.getContinuousObservationModel()))
    model$residualError <- FALSE
  if (is.null(mlx.getCovariateInformation()))
    model$covariate <- FALSE
  if (sum(mlx.getIndividualParameterModel()$variability$id) < 2)
    model$correlation <- FALSE
  if (!any(unlist(model))) {
    return(list(change=F))
    warning("\nThere is no statistical model to build...\n", call.=FALSE)
  }
  
  idir <- NULL
  if (model$covariate) {
    dir.names <- c("full", "both", "backward", "forward")
    if (is.null(direction)) {
      nbcov <- length(mlx.getCovariateInformation()$name)
      direction <- ifelse(nbcov<=n.full,"full","both")
      idir <- direction
    }
    dir0 <- direction[(!(direction %in% dir.names))]
    if (length(dir0)>0) 
      stop(paste0(dir0, " is not a valid direction name"), call.=FALSE)
  }
  
  if (ll=="none") {
    warning(" Use 'll=FALSE' instead of 'll='none'' ", call.=FALSE)
    ll = FALSE
  }
  
  if (ll=="importanceSampling") {
    warning(" Use 'linearization=FALSE' instead of 'll='importanceSampling'' ", call.=FALSE)
    ll = TRUE
    linearization=FALSE
  }
  
  if (ll=="linearization") {
    warning(" Use 'linearization=TRUE' instead of 'll='linearization'' ", call.=FALSE)
    ll = TRUE
    linearization=TRUE
  }
  
  if (linearization)
    method.ll<-"linearization"
  else
    method.ll <- "importanceSampling"
  iop.ll <- ll
  
  
  if (model$covariate) {
    p.cov <- prior$covariate
    w.cov <- weight$covariate
    cov.model <- do.call(rbind, mlx.getIndividualParameterModel()$covariateModel)
    if (!is.null(p.cov) & !is.null(w.cov)) {
      warning("Covariate model: only 'weight' or 'prior' can be defined, not both. 'weight' will be used and prior will be ignored", call.=FALSE)
      p.cov <- NULL
    }
    #if (!is.null(p.cov) & is.null(w.cov)) { 
    # if (length(p.cov)==1) {
    #   foo <- p.cov
    #   p.cov <- cov.model
    #   p.cov[is.logical(cov.model)] <- foo
    if (length(p.cov) > 1) {
      if (!identical(sort(colnames(cov.model)), sort(colnames(p.cov))) | !identical(sort(rownames(cov.model)), sort(rownames(p.cov))))
        stop("prior$covariate should be a matrix whose column names are the names of the covariates and whose row names are the names of the parameters", call.=FALSE)
      else
        prior$covariate <- p.cov[,colnames(cov.model)]
    }
    # if (is.null(p.cov) & is.null(w.cov)) 
    #   w.cov <- 1
    # if (!is.null(w.cov)) {
    #   if (length(w.cov)==1) {
    #     foo <- w.cov
    #     w.cov <- cov.model
    #     w.cov[is.logical(cov.model)] <- foo
    #   } else {
    if (length(w.cov) > 1) {
      if (!identical(sort(colnames(cov.model)), sort(colnames(w.cov))) | !identical(sort(rownames(cov.model)), sort(rownames(w.cov))))
        stop("weight$covariate should be a matrix whose column names are the names of the covariates and whose row names are the names of the parameters", call.=FALSE)
      else 
        weight$covariate <- w.cov[,colnames(cov.model)]
    }
    
  } else {
    seq.cov <- F
  }
  
  if (model$correlation) {
    p.cor <- prior$correlation
    w.cor <- weight$correlation
    if (!is.null(p.cor) & !is.null(w.cor)) {
      warning("Correlation model: only 'weight' or 'prior' can be defined, not both. 'weight' will be used and prior will be ignored", call.=FALSE)
      p.cor <- NULL
    }
    var.model <- mlx.getIndividualParameterModel()$variability$id
    n.param <- names(which(var.model))
    d.param <- length(n.param)
    if (!is.null(p.cor)) { 
      if (length(p.cor) > 1) 
        if (!all(n.param %in% colnames(p.cor)) | !all(n.param %in% rownames(p.cor)) |
            (!identical(p.cor, t(p.cor)) & !identical(p.cor, lower.tri(p.cor)*p.cor))) {
          print(p.cor)
          stop("prior$correlation should be a symetrical or triangular inferior square matrix whose column names and row names are the names of the parameters with variability", call.=FALSE)
        } 
    }
    
    if (!is.null(w.cor)) {
      if (length(w.cor)==1) 
        if (length(p.cor) > 1) 
          if (!all(n.param %in% colnames(w.cor)) | !all(n.param %in% rownames(w.cor)) |
              (!identical(w.cor, t(w.cor)) & !identical(w.cor, lower.tri(w.cor)*w.cor))) {
            print(w.cor)
            stop("weight$correlation should be a symetrical or triangular inferior square matrix whose column names and row names are the names of the parameters with variability", call.=FALSE)
          }
    }
    
    weight$correlation <- w.cor
    prior$correlation <- p.cor
  }
  
  return(list(covToTransform=covToTransform, paramToUse=paramToUse, covToTest=covToTest, 
              final.project=final.project, model=model, direction=direction, idir=idir,
              seq.cov=seq.cov, seq.corr=seq.corr, iop.ll=iop.ll, method.ll=method.ll, weight=weight, prior=prior))
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
      if (length(cov.names)==1)
        mr <- as.data.frame(mr)
      else {
        if (length(param.names)>1)
          mr <- as.data.frame(mr[,order(colnames(mr))])
        else {
          mr <- data.frame(mr)
        }
      }
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
  } else {
    n2 <- setdiff(names(r[[1]]), c(n1,n3))
    rs <- lapply(r, function(x) x[c(n1,n2,n3)])
  }
  return(rs)
}

formatErrorModel <- function(m) {
  i.ln <- which(mlx.getContinuousObservationModel()$distribution=="logNormal")
  if (length(i.ln) > 0)
    m[i.ln] <- "exponential"
  return(m)
}

formatLL <- function(ll, criterion, cr, is.weight, is.prior=F) {
  is.prior <- F
  llr <- ll[c('AIC', 'BIC', 'BICc')]
  if (is.prior) {
    llr['criterion'] <- cr
  } else {
    if (is.numeric(criterion)) {
      if (is.weight)
        llr['w.criterion'] <- cr
      else
        llr['criterion'] <- cr
    } else if (is.weight){
      llr[paste0("w",criterion)] <- cr
    }
  }
  llr["s.e."] <- ll["standardError"]
  return(llr)
}

compute.criterion <- function(criterion, method.ll, weight=NULL, pen.coef=NULL) {
  ofv <- mlx.getEstimatedLogLikelihood()[[method.ll]][["OFV"]]
  
  ind.model <- mlx.getIndividualParameterModel()
  i1 <- names(which(ind.model$variability$id))
  i0 <- names(which(!ind.model$variability$id))
  
  cov.model <- do.call(rbind, ind.model$covariateModel)
  if (ncol(cov.model) > 0)
    pen.covariate <- sum((cov.model*weight$covariate)[i1,])*pen.coef[1] + sum((cov.model*weight$covariate)[i0,])*pen.coef[2]
  else
    pen.covariate <- 0
  
  cB <- ind.model$correlationBlocks$id
  cor.model <- weight$correlation*0
  for (k in seq_along(cB))
    cor.model[cB[[k]], cB[[k]]] <- 1
  pen.correlation <- sum((lower.tri(cor.model)*cor.model)*weight$correlation)*pen.coef[1]
  
  v <- ind.model$variability$id
  pen.variance <- sum(v*weight$variance[names(v)])*pen.coef[1]
  
  pen.pop <- length(v)*pen.coef[2]
  
  error.model <- mlx.getContinuousObservationModel()$errorModel
  pen.error <- 0
  for (k in seq_along(error.model))
    pen.error <- pen.error + pen.coef[2+k]*(1 + grepl("combined",error.model[[k]]))
  
  cr <- ofv + pen.pop + pen.covariate + pen.correlation + pen.variance + pen.error
  return(cr)
}

compute.criterion.old <- function(criterion, method.ll, weight=NULL, pen.coef=NULL) {
  ll <- mlx.getEstimatedLogLikelihood()[[method.ll]]
  if (identical(criterion,"AIC"))       cr <- ll[["AIC"]]
  else if (identical(criterion,"BIC"))  cr <- ll[["BIC"]]
  else if (identical(criterion,"BICc")) cr <- ll[["BICc"]]
  else {
    d <- (ll[["AIC"]] - ll[["OFV"]])/2
    cr <- ll[["OFV"]]+d*criterion
  }
  cov.model <- do.call(rbind, mlx.getIndividualParameterModel()$covariateModel)
  w.covariate <- sum(cov.model*(weight$covariate - 1))*pen.coef[1]
  
  cB <- mlx.getIndividualParameterModel()$correlationBlocks$id
  cor.model <- weight$correlation*0
  for (k in seq_along(cB))
    cor.model[cB[[k]], cB[[k]]] <- 1
  w.correlation <- sum((lower.tri(cor.model)*cor.model)*(weight$correlation-1))*pen.coef[1]
  
  v <- mlx.getIndividualParameterModel()$variability$id
  w.variance <- sum(v*(weight$variance[names(v)]-1))*pen.coef[1]
  
  return(cr + w.covariate + w.correlation + w.variance)
}

print.result <- function(print, summary.file, to.cat=NULL,to.print=NULL) {
  
  if (file.exists(summary.file)) 
    sink(summary.file, append=TRUE)
  else  
    sink(summary.file)
  if (!is.null(to.cat))
    cat(to.cat)
  if (!is.null(to.print))
    print(to.print)
  sink()
  
  if (print) {
    if (!is.null(to.cat))
      cat(to.cat)
    if (!is.null(to.print))
      print(to.print)
  }
}

covariate.test <- function(cov.test, covToTest, covToTransform, paramToUse) {
  p.ttest <- random.effect <- covariate <- parameter <- id <- NULL
  r.test <- covariateTest()$p.value.randomEffects  %>% 
    rename(p.value=p.ttest) %>%
    mutate(parameter=gsub("eta_","", random.effect)) %>%
    select(-random.effect)
  r.test <- r.test %>% 
    filter(covariate %in% covToTest & parameter %in% paramToUse) 
  if (!is.null(cov.test))
    r.test$p.value <- pmin(cov.test$p.value, r.test$p.value)
  r.test <- r.test[,c("parameter", "covariate", "p.value", "in.model")]
  return(r.test)
}


p.weight <- function(p, pw, coef) {
  p <- pmax(p, 2.2e-16)
  p <- pmin(p, 0.99999)
  A <- pmax(p/(1-p)*(exp(pw*coef/2)-1)/(exp(coef/2)-1)  , 0)
  p <- A/(1+A)
  # r <- exp(coef*(pw-1)/2)
  # p <- p*r/(1 - p + p*r)
  return(p)
}

