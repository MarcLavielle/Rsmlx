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
#' @param final.project  a string: the final Monolix project (default adds "_built" to the original project)
#' @param model  components of the model to optimize c("residualError", "covariate", "correlation"), (default="all")
#' @param paramToUse  list of parameters possibly function of covariates (default="all")
#' @param fix.param1  parameters with variability that cannot be removed (default=NULL)
#' @param fix.param0  parameters without variability that cannot be added (default=NULL)
#' @param covToTest  components of the covariate model that can be modified   (default="all")
#' @param covToTransform  list of (continuous) covariates to be log-transformed (default="none")
#' @param center.covariate TRUE/{FALSE} center the covariates of the final model (default=FALSE) 
#' @param criterion  penalization criterion to optimize c("AIC", "BIC", {"BICc"}, gamma)
#' @param linearization  TRUE/{FALSE} whether the computation of the likelihood is based on a linearization of the model (default=FALSE, deprecated)
#' @param ll  {TRUE}/FALSE  compute the observe likelihood and the criterion to optimize at each iteration
#' @param seq.cov TRUE/{FALSE} whether the covariate model is built before the correlation model  
#' @param seq.cov.iter number of iterations before building the correlation model (only when seq.cov=F, default=0) 
#' @param seq.corr {TRUE}/FALSE whether the correlation model is built iteratively (default=TRUE) 
#' @param p.max  maximum p-value used for removing non significant relationships between covariates and individual parameters (default=0.1)
#' @param p.min minimum p-values used for testing the components of a new model (default=c(0.075, 0.05, 0.1))
#' @param remove  try to remove random effects (default=T)
#' @param add  try to add random effects (default=T)
#' @param delta  maximum difference in criteria for testing a new model (default=c(30,5))
#' @param omega.set settings to define how a variance varies during iterations of SAEM
#' @param pop.set1  Monolix settings 1
#' @param pop.set2  Monolix settings 2
#' @param pen.cov multiplicative penalty term for the covariate model (default=1)
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
buildAll <- function(project, final.project=NULL, model="all", 
                     paramToUse="all", covToTest="all", covToTransform="none", center.covariate=FALSE, 
                     criterion="BICc", linearization=FALSE, ll=T, pen.cov=1, direction=NULL, steps=1000,
                     max.iter=20, explor.iter=2, seq.cov=FALSE, seq.corr=TRUE, seq.cov.iter=0, 
                     p.max=0.1, p.min=c(0.075, 0.05, 0.1), print=TRUE, nb.model=1,
                     fix.param1=NULL, fix.param0=NULL, remove=T, add=T, delta=c(30,5), 
                     omega.set=NULL, pop.set1=NULL, pop.set2=NULL) {
  
  in.model <- p.ttest <- random.effect <- covariate <- param <- NULL
  
  if (is.null(final.project))
    final.project <- gsub(".mlxtran", "_buildAll.mlxtran", project)
  
  dir.project <- gsub(".mlxtran","", project)
  dir.built <- file.path(dir.project, "buildAll")
  if (dir.exists(dir.built))
    unlink(dir.built, recursive=TRUE)
  dir.create(dir.built)
  
  #  start with "full variance" model
  mlx.loadProject(project)
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
  pset1 <- list(nbexploratoryiterations=75, nbsmoothingiterations=75, simulatedannealing=F, smoothingautostop=F, exploratoryautostop=F)
  if (!is.null(pop.set1))
    pset1 <- modifyList(pset1, pop.set1[intersect(names(pop.set1), names(pset1))])
  pop.set1 <- mlx.getPopulationParameterEstimationSettings()
  pop.set1 <- modifyList(pop.set1, pset1[intersect(names(pset1), names(pop.set1))])
  
  pset2 <- list(nbsmoothingiterations=25)
  if (!is.null(pop.set2))
    pset2 <- modifyList(pset2, pop.set2[intersect(names(pop.set2), names(pset2))])
  pop.set2 <- modifyList(pop.set1, pset2[intersect(names(pset2), names(pop.set1))])
  # ---------------------
  
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
      } else {
        project.ini.build <- project.final.builtvar
        seq.corr <- seq.cov <- F
        seq.cov.iter=0
      }
      project.final.built <- file.path(dir.built, paste0("project_built",iter,".mlxtran"))
      mlx.loadProject(project.ini.build)
      cov1 <- mlx.getIndividualParameterModel()$covariateModel
      cor1 <- rep(F, length(cov1))
      names(cor1) <- names(cov1)
      cor1[unlist(mlx.getIndividualParameterModel()$correlationBlocks$id)] <- T
      r.build <- buildmlx(project=project.ini.build, final.project=project.final.built, covToTest=covToTest,
                          covToTransform=covToTransform, seq.corr=seq.corr, seq.cov=seq.cov, seq.cov.iter=seq.cov.iter, p.max=p.max, p.min=p.min,
                          model=model, paramToUse=paramToUse, center.covariate=center.covariate, criterion=criterion, 
                          linearization=linearization, ll=ll, pen.cov=pen.cov, direction=direction, steps=steps,
                          max.iter=max.iter, explor.iter=explor.iter, nb.model=nb.model, print=print)
      relToTest <- rbind(relToTest, covariateTest()$p.value.randomEffects %>% filter(in.model==F & p.ttest<p.min[1]) %>% select(c(random.effect, covariate, p.ttest)))
      cov2 <- mlx.getIndividualParameterModel()$covariateModel
      cor2 <- rep(F, length(cov2))
      names(cor2) <- names(cov2)
      cor2[unlist(mlx.getIndividualParameterModel()$correlationBlocks$id)] <- T
      var2 <- unlist(mlx.getIndividualParameterModel()$variability$id)
    } else {
      r.build$change <- F
    }
    if (iter==1 || r.build$change) {
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
        r.buildvar <- buildVar(project.ini.buildvar, final.project=project.final.builtvar, linearization=linearization,
                               fix.param1=fix1, fix.param0=fix0, criterion=criterion, print=print,
                               remove=remove, add=add, delta=delta, omega.set=omega.set, pop.set1=pop.set1, pop.set2=pop.set2)
      } else {
        r.buildvar$change <- F
      }
    } else {
      r.buildvar$change <- F
    }
    change <- r.build$change | r.buildvar$change 
  }
  mlx.saveProject(final.project)
  # }
  # relToTest <- rbind(relToTest, covariateTest()$p.value.randomEffects %>% filter(in.model==F & p.ttest<p.min[1]) %>% select(c(random.effect, covariate)))
  # mlx.loadProject(final.project)
  
  
  if (!is.null(relToTest)) {
    param0 <- names(which(!mlx.getIndividualParameterModel()$variability$id))
    relToTest <- unique(relToTest) %>% mutate(param=gsub("eta_","", random.effect)) %>% filter(param %in% param0) %>% select(-random.effect)
    if (nrow(relToTest)>0) {
      
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
      
      ll.ini <- computecriterion(criterion, method.ll)
      ll.min <- computecriterion(criterion, method.ll.iter)
      
      ind.param <- mlx.getEstimatedIndividualParameters()$saem
      mlx.setInitialEstimatesToLastEstimates()
      pop.set.ini <- mlx.getPopulationParameterEstimationSettings()
      scenario.ini <- mlx.getScenario()
      mlx.setPopulationParameterEstimationSettings(pop.set1)
      g.min <- g.ini <- mlx.getIndividualParameterModel()
      for (k in 1:nrow(relToTest)) {
        print(relToTest[k,])
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
        ll.iter <- computecriterion(criterion, method.ll.iter)
        print(c(ll.min, ll.iter))
        if (ll.iter < ll.min) {
          ll.min <- ll.iter
          g.min <- g
          ind.param <- mlx.getEstimatedIndividualParameters()$saem
          mlx.setInitialEstimatesToLastEstimates()
        }
      }
      if (!identical(g.ini, g.min)) {
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
        ll.new <- computecriterion(criterion, method.ll)
        ll.disp <- formatLL(mlx.getEstimatedLogLikelihood()[[method.ll]])
        if (is.numeric(criterion))
          ll.disp['criterion'] <- ll.new
        if (print) {
          cat(paste0("\nEstimated criteria (",method.ll,"):\n"))
          print(round(ll.disp,2))
        }
        if (ll.new < ll.ini) {
          mlx.saveProject(final.project)
          r.build <- buildmlx(project=final.project, covToTest=covToTest,
                              covToTransform=covToTransform, seq.corr=F, seq.cov=F, seq.cov.iter=0, p.max=p.max, p.min=p.min,
                              model=model, paramToUse=paramToUse, center.covariate=center.covariate, criterion=criterion, 
                              linearization=linearization, ll=ll, pen.cov=pen.cov, direction=direction, steps=steps,
                              max.iter=max.iter, explor.iter=explor.iter, nb.model=nb.model, print=print)
          mlx.saveProject(final.project)
        }
      }
    }
    mlx.loadProject(final.project)
  }
  
  return(list(project=final.project))
}




