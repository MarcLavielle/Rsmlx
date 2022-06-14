#' Automatic complete statistical model building
#'
#' buildAll builds the complete statistical model by iteratively calling functions
#'  buildmlx and buildVar 
#'  
#' Penalization criterion can be either a custom penalization of the form gamma*(number of parameters),
#' AIC (gamma=2) or BIC (gamma=log(N)).
#' 
#' See http://rsmlx.webpopix.org for more details.
#' @param project a string: the initial Monolix project
#' @param final.project  a string: the final Monolix project (default adds "_buildAll" to the original project)
#' @param model  components of the model to optimize c("residualError", "covariate", "correlation"), (default="all")
#' @param prior  list of prior probabilities for each component of the model (default=NULL)
#' @param weight list of penalty weights for each component of the model (default=NULL)
#' @param coef.w1 multiplicative weight coefficient used for the first iteration only (default=0.5)
#' @param cv.min  value of the coefficient of variation below which an individual parameter is considered fixed (default=0.001)
#' @param fError.min  minimum fraction of residual variance for combined error model (default = 1e-3)
#' @param paramToUse  list of parameters possibly function of covariates (default="all")
#' @param fix.param1  parameters with variability that cannot be removed (default=NULL)
#' @param fix.param0  parameters without variability that cannot be added (default=NULL)
#' @param covToTest  components of the covariate model that can be modified   (default="all")
#' @param covToTransform  list of (continuous) covariates to be log-transformed (default="none")
#' @param center.covariate TRUE/{FALSE} center the covariates of the final model (default=FALSE) 
#' @param criterion  penalization criterion to optimize c("AIC", "BIC", {"BICc"}, gamma)
#' @param linearization  TRUE/{FALSE} whether the computation of the likelihood is based on a linearization of the model (default=FALSE, deprecated)
#' @param test  {TRUE}/FALSE  perform additional statistical tests for building the model (default=TRUE)
#' @param ll  {TRUE}/FALSE  compute the observe likelihood and the criterion to optimize at each iteration
#' @param seq.cov TRUE/{FALSE} whether the covariate model is built before the correlation model  
#' @param seq.cov.iter number of iterations before building the correlation model (only when seq.cov=F, default=0) 
#' @param seq.corr {TRUE}/FALSE whether the correlation model is built iteratively (default=TRUE) 
#' @param p.max  maximum p-value used for removing non significant relationships between covariates and individual parameters (default=0.1)
#' @param p.min minimum p-values used for testing the components of a new model (default=c(0.075, 0.05, 0.1))
#' @param remove  try to remove random effects (default=T)
#' @param add  try to add random effects (default=T)
#' @param delta  maximum difference in criteria for testing a new model (default=c(30,10,5))
#' @param omega.set settings to define how a variance varies during iterations of SAEM
#' @param pop.set1  Monolix settings 1
#' @param pop.set2  Monolix settings 2
#' @param direction method for covariate search c({"full"}, "both", "backward", "forward"), (default="full" or "both")
#' @param steps maximum number of iteration for stepAIC (default=1000)
#' @param max.iter maximum number of iterations (default=20)
#' @param explor.iter  number of iterations during the exploratory phase (default=2)
#' @param print {TRUE}/FALSE display the results (default=TRUE)
#' @param nb.model number of models to display at each iteration (default=1)
#' @return a new Monolix project with a new statistical model.
#' @examples
#' \dontrun{
#' # Build the complete statistical model using the default settings
#' r1 <- buildAll(project="warfarinPK_project.mlxtran")
#' 
#' # Force parameter Tlag to be fixed (no variability) and parameter Cl to vary
#' r2 <- buildAll(project="warfarinPK_project.mlxtran", fix.param0="Tlag", fix.param1="Cl")
#' 
#' # Estimate the log-likelihood by linearization of the model (faster)
#' r3 <- buildAll(project="warfarinPK_project.mlxtran", linearization=T)
#' 
#' }
#' 
#' # See http://rsmlx.webpopix.org/userguide/buildmlx/ for detailed examples of use of buildmlx
#' # Download the demo examples here: http://rsmlx.webpopix.org/installation
#' 
#' @importFrom MASS stepAIC 
#' @importFrom stats coef as.formula model.matrix
#' @importFrom utils data write.csv packageVersion
#' @importFrom dplyr filter select rename arrange bind_rows rename
#' @importFrom dplyr %>%
#' @export
buildAll <- function(project=NULL, final.project=NULL, model="all", prior=NULL, weight=NULL, coef.w1=0.5, cv.min=0.001, fError.min=1e-3,
                     paramToUse="all", covToTest="all", covToTransform="none", center.covariate=FALSE, 
                     criterion="BICc", linearization=FALSE, ll=T, test=T, direction=NULL, steps=1000,
                     max.iter=20, explor.iter=2, seq.cov=FALSE, seq.corr=TRUE, seq.cov.iter=0, 
                     p.max=0.1, p.min=c(0.075, 0.05, 0.1), print=TRUE, nb.model=1,
                     fix.param1=NULL, fix.param0=NULL, remove=T, add=T, delta=c(30,10,5), 
                     omega.set=NULL, pop.set1=NULL, pop.set2=NULL) {
  
  ptm <- proc.time()
  
  dashed.line <- "--------------------------------------------------\n"
  plain.line <-  "__________________________________________________\n"
  dashed.short <- "-----------------------\n"
  plain.short  <- "_______________________\n"
  
  op.original <- options()
  op.new <- options()
  op.new$lixoft_notificationOptions$warnings <- 1   #hide the warning messages
  options(op.new)
  
  # is.weight <- !is.null(weight)
  # is.prior <- !is.null(prior)
  
  if (!is.null(project)) {
    project <- prcheck(project)$project
    mlx.loadProject(project)
  }  else 
    project <- mlx.getProjectSettings()$project
  
  in.model <- p.ttest <- random.effect <- covariate <- param <- pen.coef <- NULL
  
  if (is.null(final.project))
    final.project <- gsub(".mlxtran", "_buildAll.mlxtran", project)
  
  dir.project <- gsub(".mlxtran","", project)
  dir.built <- file.path(dir.project, "buildAll")
  if (dir.exists(dir.built))
    unlink(dir.built, recursive=TRUE)
  dir.create(dir.built)
  
  #  start with "full variance" model
  g <- mlx.getIndividualParameterModel()
  if (any(!g$variability$id[1:length(g$variability$id)])) {
    g$variability$id[1:length(g$variability$id)] <- T
    mlx.setIndividualParameterModel(g)
    project.ini <- file.path(dir.built, "project_ini.mlxtran")
    mlx.saveProject(project.ini)
  } else {
    project.ini <- project
  }
  g <- mlx.getScenario()
  g$tasks['plots'] <- F
  mlx.setScenario(g)
  
  # ---------------------
  # pset1 <- list(nbexploratoryiterations=75, nbsmoothingiterations=75, simulatedannealing=F, smoothingautostop=F, exploratoryautostop=F)
  # if (!is.null(pop.set1))
  #   pset1 <- modifyList(pset1, pop.set1[intersect(names(pop.set1), names(pset1))])
  # pop.set1 <- mlx.getPopulationParameterEstimationSettings()
  # pop.set1 <- modifyList(pop.set1, pset1[intersect(names(pset1), names(pop.set1))])
  # 
  # pset2 <- list(nbsmoothingiterations=25)
  # if (!is.null(pop.set2))
  #   pset2 <- modifyList(pset2, pop.set2[intersect(names(pop.set2), names(pset2))])
  # pop.set2 <- modifyList(pop.set1, pset2[intersect(names(pset2), names(pop.set1))])
  # ---------------------
  
  r <- def.variable(weight=weight, prior=prior, criterion=criterion, fix.param0=fix.param0, fix.param1=fix.param1)
  for (j in 1:length(r))
    eval(parse(text=paste0(names(r)[j],"= r[[j]]")))
  weight0 <- weight
  #weight0$variance <- 1
  
  method.ll <- ifelse(linearization,"linearization","importanceSampling")
  
  #  iterate buildmlx (SAMBA) and buildvar until convergence
  change <- T
  iter <- 0
  relToTest <- NULL
  #  if (1>2) {
  while (change) {
    iter <- iter+1
    if (iter==1 || r.buildvar$change) {
      if (iter==1) {
        project.ini.build <- project.ini
        loadProject(project.ini)
      } else {
        coef.w1=1
        project.ini.build <- project.final.builtvar
        seq.corr <- seq.cov <- F
        seq.cov.iter=0
        covToTransform <- NULL
        covToTest <- r.build$covToTest
      }
      project.final.built <- file.path(dir.built, paste0("project_built",iter,".mlxtran"))
      mlx.loadProject(project.ini.build)
      cov1 <- mlx.getIndividualParameterModel()$covariateModel
      cor1 <- rep(F, length(cov1))
      names(cor1) <- names(cov1)
      cor1[unlist(mlx.getIndividualParameterModel()$correlationBlocks$id)] <- T
      r.build <- buildmlx(project=project.ini.build, final.project=project.final.built, covToTest=covToTest, prior=NULL, weight=weight, coef.w1=coef.w1,
                          covToTransform=covToTransform, seq.corr=seq.corr, seq.cov=seq.cov, seq.cov.iter=seq.cov.iter, p.max=p.max, p.min=p.min,
                          model=model, paramToUse=paramToUse, center.covariate=center.covariate, criterion=criterion, 
                          linearization=linearization, ll=ll, test=test, direction=direction, steps=steps, fError.min=fError.min,
                          max.iter=max.iter, explor.iter=explor.iter, nb.model=nb.model, print=print)
      if (test) {
      if (!mlx.getLaunchedTasks()$conditionalDistributionSampling)
        mlx.runConditionalDistributionSampling()
      pv.re <- covariateTest()$p.value.randomEffects
      if (!is.null(pv.re))
        relToTest <- rbind(relToTest, covariateTest()$p.value.randomEffects %>% filter(in.model==F & p.ttest<p.min[1]) %>% select(c(random.effect, covariate, p.ttest)))
      }
      cov2 <- mlx.getIndividualParameterModel()$covariateModel
      cor2 <- rep(F, length(cov2))
      names(cor2) <- names(cov2)
      cor2[unlist(mlx.getIndividualParameterModel()$correlationBlocks$id)] <- T
      var2 <- unlist(mlx.getIndividualParameterModel()$variability$id)
      weight <- r.build$weight
    } else {
      r.build$change <- F
    }
    if ((iter==1 || r.build$change) & ("variance" %in% model | identical(model,"all"))) {
      fix0 <- fix.param0
      fix1 <- fix.param1
      if (iter >= 2) {
        list.n <- setdiff(names(cov1), c(fix.param0, fix.param1))
        for (nc in list.n) {
          if (!var2[nc] & identical(cov1[nc], cov2[nc]))
            fix0 <- c(fix0, nc)
          if (var2[nc] & identical(cov1[nc], cov2[nc])) {
            if (!cor1[nc] | identical(cor1[nc], cor2[nc]))
              fix1 <- c(fix1, nc)
          }
        }
      }
      if (length(setdiff(names(cov1), c(fix0, fix1)))>0){
        project.ini.buildvar <- project.final.built
        project.final.builtvar <- file.path(dir.built, paste0("project_builtvar",iter,".mlxtran"))
        r.buildvar <- buildVar(project.ini.buildvar, final.project=project.final.builtvar, linearization=linearization, cv.min=cv.min,
                               fix.param1=fix1, fix.param0=fix0, criterion=criterion, print=print, prior=NULL, weight=weight,
                               remove=remove, add=add, delta=delta, omega.set=omega.set, pop.set1=pop.set1, pop.set2=pop.set2)
      } else {
        r.buildvar <- list(change = F)
      }
    } else {
      r.buildvar <- list(change=F)
    }
    change <- r.build$change | r.buildvar$change 
  }
  mlx.saveProject(final.project)
  # }
  # relToTest <- rbind(relToTest, covariateTest()$p.value.randomEffects %>% filter(in.model==F & p.ttest<p.min[1]) %>% select(c(random.effect, covariate)))
  # mlx.loadProject(final.project)
  
  if (!identical(model, "all") & !("covariate" %in% model))
    relToTest <- NULL
  
  add.test <- NULL
  test.all <- F
  if (!is.null(relToTest) & test.all==T) {
    param0 <- names(which(!mlx.getIndividualParameterModel()$variability$id))
    relToTest <- unique(relToTest) %>% mutate(param=gsub("eta_","", random.effect)) %>% filter(param %in% param0) %>% select(-random.effect)
    if (nrow(relToTest)>0) {
      if (print) {
        cat(paste0("\n",dashed.line))
        cat("\nTesting additional relationships between covariates and parameters\n")
      }
      linearization.iter <- T
      if (any(mlx.getObservationInformation()$type != "continuous")) 
        linearization <- linearization.iter <-  F
      method.ll <- ifelse(linearization, "linearization", "importanceSampling")
      method.ll.iter <- ifelse(linearization.iter, "linearization", "importanceSampling")
      
      g <- as.list(mlx.getLaunchedTasks())
      if (linearization | linearization.iter) {
        if (!("linearization" %in% g[['logLikelihoodEstimation']])) {
          if (!g$conditionalModeEstimation) {
            mlx.runConditionalModeEstimation()
          }
          mlx.runLogLikelihoodEstimation(linearization=TRUE)
        }
      }
      
      project.relToTest <- file.path(dir.built, paste0("project_relToTest.mlxtran"))
      mlx.saveProject(project.relToTest)
      
      ll.ini <- compute.criterion(criterion, method.ll, weight, pen.coef)
      ll.min <- compute.criterion(criterion, method.ll.iter, weight, pen.coef)
      
      ind.param <- mlx.getEstimatedIndividualParameters()$saem
      mlx.setInitialEstimatesToLastEstimates()
      pop.set.ini <- mlx.getPopulationParameterEstimationSettings()
      scenario.ini <- mlx.getScenario()
      mlx.setPopulationParameterEstimationSettings(pop.set1)
      g.min <- g.ini <- mlx.getIndividualParameterModel()
       for (k in 1:nrow(relToTest)) {
        if (print)  print(relToTest[k,])
        g <- g.min
        g$covariateModel[[relToTest$param[k]]][[relToTest$covariate[k]]] <- T
        mlx.setIndividualParameterModel(g)
        mlx.runPopulationParameterEstimation(parameters=ind.param)
        if (linearization.iter) {
          mlx.runConditionalModeEstimation()
          mlx.runLogLikelihoodEstimation(linearization=TRUE)
        } else {
          mlx.runConditionalDistributionSampling()
          mlx.runLogLikelihoodEstimation(linearization=FALSE)
        }
        
        ll.iter <- compute.criterion(criterion, method.ll.iter, weight, pen.coef)
        print(c(ll.min, ll.iter))
        if (ll.iter < ll.min) {
          ll.min <- ll.iter
          g.min <- g
          ind.param <- mlx.getEstimatedIndividualParameters()$saem
          mlx.setInitialEstimatesToLastEstimates()
        }
      }
      if (!identical(g.ini, g.min)) {
        add.test <- list( relToTest=relToTest, g.ini=g.ini, g.min=g.min)
        mlx.loadProject(final.project)
        mlx.saveProject(project.relToTest)
        mlx.setPopulationParameterEstimationSettings(pop.set.ini)
        mlx.setIndividualParameterModel(g.min)
        if (print)
          cat("Fitting the model... \n")
        mlx.runPopulationParameterEstimation()
        if (linearization) {
          mlx.runConditionalModeEstimation()
          mlx.runLogLikelihoodEstimation(linearization=TRUE)
        } else {
          mlx.runConditionalDistributionSampling()
          mlx.runLogLikelihoodEstimation(linearization=FALSE)
        }
        r <- def.variable(weight=weight, prior=prior, criterion=criterion)
        for (j in 1:length(r))
          eval(parse(text=paste0(names(r)[j],"= r[[j]]")))
        
        ll.new <- compute.criterion(criterion, method.ll, weight, pen.coef)
        ll.disp <- formatLL(mlx.getEstimatedLogLikelihood()[[method.ll]], criterion, ll.new, weight$is.weight)
        if (is.numeric(criterion))
          ll.disp['criterion'] <- ll.new
        if (print) {
          cat(paste0("\nEstimated criteria (",method.ll,"):\n"))
          print(round(ll.disp,2))
        }
        if (ll.new < ll.ini) {
          mlx.saveProject(final.project)
          # if (print)
          r.build <- buildmlx(project=final.project, covToTest=covToTest, prior=NULL, weight=weight,
                              covToTransform=covToTransform, seq.corr=F, seq.cov=F, seq.cov.iter=0, p.max=p.max, p.min=p.min,
                              model=model, paramToUse=paramToUse, center.covariate=center.covariate, criterion=criterion, 
                              linearization=linearization, ll=ll, test=test, direction=direction, steps=steps, fError.min=fError.min,
                              max.iter=max.iter, explor.iter=explor.iter, nb.model=nb.model, print=print)
          mlx.saveProject(final.project)
          unlink(r.build$project)
          unlink(gsub("mlxtran", "mlxproperties", r.build$project))
          unlink(gsub(".mlxtran", "", r.build$project), recursive=T)
        }
      }
    }
    mlx.loadProject(final.project)
  }
  
  if (print) {
    param0 <- names(which(!mlx.getIndividualParameterModel()$variability$id))
    param1 <- names(which(mlx.getIndividualParameterModel()$variability$id))
    covariate.model.print <- formatCovariateModel(mlx.getIndividualParameterModel()$covariateModel)
    correlation.model.print <- lapply(mlx.getIndividualParameterModel()$correlationBlocks$id, sort)
    error.model.print <- formatErrorModel(mlx.getContinuousObservationModel()$errorModel)
    
    cat(paste0("\n",dashed.line,"\nFinal complete model:\n"))
    cat("\nVariance model: \n")
    cat("Parameters without variability:", param0, "\n")
    cat("Parameters with variability   :", param1, "\n")
    cat("\nCovariate model:\n")
    print(covariate.model.print)
    cat("\nCorrelation model:\n")
    if (length(correlation.model.print)>0)
      print(correlation.model.print)
    else
      print("NULL")
    if (length(error.model.print)>0) {
      cat("\nResidual error model:\n")
      print(error.model.print)
    }
    if (ll) {
      ll.final <- compute.criterion(criterion, method.ll, r.build$weight, pen.coef)
      ll <- formatLL(mlx.getEstimatedLogLikelihood()[[method.ll]], criterion, ll.final, weight$is.weight)
      cat(paste0("\nEstimated criteria (",method.ll,"):\n"))
      print(round(ll,2)) 
    }
    dt <- proc.time() - ptm
    cat(paste0("\ntotal time: ", round(dt["elapsed"], digits=1),"s\n", dashed.line, "\n"))
    
  }
  
#  if (!is.null(add.test))  browser()
#  return(list(project=final.project, add.test=add.test))
  return(list(project=final.project))
}




