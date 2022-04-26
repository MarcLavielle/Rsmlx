#' Automatic model variance building
#'
#' buildVar is designed to build the best variance model for the random effects by selecting
#' which individual parameters vary and which ones are fixed.
#' 
#' Penalization criterion can be either a custom penalization of the form gamma*(number of parameters),
#' AIC (gamma=2) or BIC (gamma=log(N)).
#' 
#' See http://rsmlx.webpopix.org for more details.
#' @param project a string: the initial Monolix project
#' @param final.project  a string: the final Monolix project (default adds "_var" to the original project)
#' @param prior  named vector of prior probabilities (default=NULL)
#' @param weight  named vector of weights (default=NULL)
#' @param cv.min  value of the coefficient of variation below which an individual parameter is considered fixed (default=0.001)
#' @param fix.param1  parameters with variability that cannot be removed (default=NULL)
#' @param fix.param0  parameters without variability that cannot be added (default=NULL)
#' @param criterion  penalization criterion to optimize c("AIC", "BIC", {"BICc"}, gamma)
#' @param linearization  TRUE/{FALSE} whether the computation of the likelihood is based on a linearization of the model (default=FALSE)
#' @param remove  {TRUE}/FALSE try to remove random effects (default=TRUE)
#' @param add  {TRUE}/FALSE try to add random effects (default=TRUE)
#' @param delta  maximum difference in criteria for testing a new model (default=c(30,5))
#' @param omega.set settings to define how a variance varies during iterations of SAEM
#' @param pop.set1  Monolix settings 1
#' @param pop.set2  Monolix settings 2
#' @param print {TRUE}/FALSE display the results (default=TRUE)
#' @return a new Monolix project with a new inter individual variability model.
#' @examples
#' \dontrun{
#' # Build the variability model using the default settings
#' r1 <- buildVar(project="warfarinPK_project.mlxtran")
#' 
#' # Force parameter Tlag to be fixed (no variability) and parameter Cl to vary
#' r2 <- buildVar(project="warfarinPK_project.mlxtran", fix.param0="Tlag", fix.param1="Cl")
#' 
#' # Estimate the log-likelihood by linearization of the model (faster)
#' r3 <- buildVar(project="warfarinPK_project.mlxtran", linearization=T)
#' 
#' }
#' 
#' # See http://rsmlx.webpopix.org/userguide/buildvar/ for detailed examples of use of buildvar
#' # Download the demo examples here: http://rsmlx.webpopix.org/installation
#' 
#' 
#' @importFrom stats coef as.formula model.matrix integrate
#' @importFrom utils modifyList 
#' @export

buildVar <- function(project=NULL, weight=NULL, prior=NULL, cv.min=0.001, final.project=NULL, 
                     fix.param1=NULL, fix.param0=NULL, criterion="BICc", linearization=F, remove=T, add=T,
                     delta=c(30,5), omega.set=NULL, pop.set1=NULL, pop.set2=NULL, print=TRUE) {
  
  
  ptm <- proc.time()
  
  dashed.line <- "--------------------------------------------------\n"
  plain.line <-  "__________________________________________________\n"
  dashed.short <- "-----------------------\n"
  plain.short  <- "_______________________\n"
  
  op.original <- options()
  op.new <- options()
  op.new$lixoft_notificationOptions$warnings <- 1   #hide the warning messages
  options(op.new)
  
  is.weight <- !is.null(weight)
  is.prior <- !is.null(prior)
  
  if (!is.null(project)) 
    project <- prcheck(project)$project
  else 
    project <- mlx.getProjectSettings()$project
  
  oset <- project.built <- project.temp <- pen.coef <- NULL
  r <- buildVar.check(project=project, final.project=final.project, fix.param1=fix.param1, fix.param0=fix.param0, 
                      criterion=criterion, linearization=linearization, remove=remove, add=add, delta=delta, 
                      omega.set=omega.set, pop.set1=pop.set1, pop.set2=pop.set2, print=print, prior=prior, weight=weight) 
  
  for (j in 1:length(r))
    eval(parse(text=paste0(names(r)[j],"= r[[j]]")))
  
  r <- def.variable(weight, prior, criterion)
  for (j in 1:length(r))
    eval(parse(text=paste0(names(r)[j],"= r[[j]]")))
  w.var <- weight$variance
  
  if (print) cat(paste0("\n",dashed.line,"\nBuilding the variance model\n"))
  
  mlx.saveProject(project.built)
  
  final.dir <- sub(pattern = "(.*)\\..*$", replacement = "\\1", final.project)
  if (dir.exists(final.dir)) 
    unlink(final.dir, recursive=TRUE)
  Sys.sleep(0.1)
  dir.create(final.dir, recursive=T)
  
  project.dir <- mlx.getProjectSettings()$directory
  if (!dir.exists(project.dir))
    dir.create(project.dir)
  
  # buildVar.dir <- file.path(mlx.getProjectSettings()$directory,"buildVar")
  # Sys.sleep(0.1)
  # if (dir.exists(buildVar.dir))
  #   unlink(buildVar.dir, recursive=TRUE)
  # 
  # Sys.sleep(0.1)
  # dir.create(buildVar.dir)
  # summary.file = file.path(buildVar.dir,"summary.txt")
  
  linearization.iter <- T
  if (any(mlx.getObservationInformation()$type != "continuous")) {
    linearization <- linearization.iter <-  F
  }
  
  method.ll <- ifelse(linearization, "linearization", "importanceSampling")
  method.ll.iter <- ifelse(linearization.iter, "linearization", "importanceSampling")
  pop.set0 <- mlx.getPopulationParameterEstimationSettings()
  scenario0 <- mlx.getScenario()
  
  g <- mlx.getIndividualParameterModel()
  gv <- gv0 <- g$variability$id
  change.ini <- F
  if (!is.null(fix.param1))
    gv[fix.param1] <- T
  if (!is.null(fix.param0))
    gv[fix.param0] <- F
  if (!identical(gv,gv0)) {
    g$variability$id <- gv
    mlx.setIndividualParameterModel(g)
    change.ini <- T
  }
  
  p.param1 <- sapply(mlx.getIndividualParameterModel()$name,function(x) NULL)
  lpar <- max(nchar(mlx.getIndividualParameterModel()$name))
  
  #---------------------------------------------
  print.line <- F
  if (remove) {
    g <- mlx.getIndividualParameterModel()$variability$id
    param0 <- names(which(!g))
    param1 <- names(which(g))
    
    g <- as.list(mlx.getLaunchedTasks())
    if (print) {
      cat(paste0("\n", plain.line))
      print.line <- T
    }
    if (!g[['populationParameterEstimation']]) {
      if (print)  cat("\nEstimating the population parameters\n\n")
      print.line <- F
      mlx.runPopulationParameterEstimation()
    }
    pop.built <- mlx.getEstimatedPopulationParameters()
    
    list.param1 <- names(which(!unlist(lapply(p.param1[param1], function(x) all(param1 %in% x)))))
    list.param1 <- setdiff(list.param1, fix.param1)
    if (length(list.param1)>0) {
      r.param1 <- compute.cv(pop.built, list.param1)
      j0 <- which(r.param1$cv < cv.min)
      if (length(j0) > 0) {
        if (print) {
          if (!print.line) cat(paste0("\n", plain.line))
          cat("\nParameters without variability:", param0, "\n")
          cat("Parameters with variability   :", param1, "\n\n")
          print(r.param1$o)
          print(r.param1$cv)
          cat(paste0("\nRemove variability on ",list.param1[j0], "\n"))
          print.line <- F
        }
        param0 <- c(param0, list.param1[j0])
        param1 <- setdiff(param1, list.param1[j0])
        update.project(project, project.built, param0, param1, NULL, pop.set1)
        list.param1 <- names(which(!unlist(lapply(p.param1[param1], function(x) all(param1 %in% x)))))
        list.param1 <- setdiff(list.param1, fix.param1)
      }
      list.omega1 <- paste0("omega_", list.param1)
      list.omega1 <- list.omega1[order(pop.built[list.omega1])]
    }
  }
  #---------------------------------------------
  
  g <- as.list(mlx.getLaunchedTasks())
  if (!g[['populationParameterEstimation']]) {
    if (print) {
      if (!print.line) cat(plain.line)
      cat("\nFitting the initial model\n")
      print.line <- F
    }
    mlx.runPopulationParameterEstimation()
  }
  
  p0.pop <- mlx.getEstimatedPopulationParameters()
  omega.value <- p0.pop[grep("omega", names(p0.pop))]
  p0.ind <- mlx.getEstimatedIndividualParameters()$saem
  N <- nrow(p0.ind)
  
  g <- as.list(mlx.getLaunchedTasks())
  if (linearization | linearization.iter) {
    if (!("linearization" %in% g[['logLikelihoodEstimation']])) {
      if (!g$conditionalModeEstimation) {
        mlx.runConditionalModeEstimation()
      }
      mlx.runLogLikelihoodEstimation(linearization=TRUE)
    }
  }
  
  if (!linearization | !linearization.iter) {
    if (!("importanceSampling" %in% g[['logLikelihoodEstimation']])) {
      mlx.runConditionalDistributionSampling()
      mlx.runLogLikelihoodEstimation(linearization=FALSE)
    }
  }
  BICc.built <- compute.criterion(criterion, method.ll, weight, pen.coef)
  BICc.built.iter <- compute.criterion(criterion, method.ll.iter, weight, pen.coef)
  
  mlx.setInitialEstimatesToLastEstimates(fixedEffectsOnly = FALSE)
  mlx.saveProject(project.built)
  
  #BICc.built <- ll0$BICc
  pop.built <- p0.pop
  ind.built <- p0.ind
  g <- mlx.getIndividualParameterModel()$variability$id
  param0.built <- names(which(!g))
  param1.built <-names(which(g))
  
  beginIter <- oset$beginIter 
  nbIter <- oset$nbIter
  nbEst <- oset$nbEst
  endIter <- oset$endIter
  r.omega <- oset$r.omega
  
  pset1 <- list(nbexploratoryiterations=75, nbsmoothingiterations=75, simulatedannealing=F, smoothingautostop=F, exploratoryautostop=F)
  if (!is.null(pop.set1))
    pset1 <- modifyList(pset1, pop.set1[intersect(names(pop.set1), names(pset1))])
  pop.set1 <- mlx.getPopulationParameterEstimationSettings()
  pop.set1 <- modifyList(pop.set1, pset1[intersect(names(pset1), names(pop.set1))])
  
  pset2 <- list(nbsmoothingiterations=25)
  if (!is.null(pop.set2))
    pset2 <- modifyList(pset2, pop.set2[intersect(names(pop.set2), names(pset2))])
  pop.set2 <- modifyList(pop.set1, pset2[intersect(names(pset2), names(pop.set1))])
  pop.set2$nbexploratoryiterations <- beginIter + nbIter + endIter
  pop.set3 <- pop.set2
  pop.set3$nbburningiterations <- pop.set3$nbexploratoryiterations <- pop.set3$nbsmoothingiterations <- 0
  
  stop <- F 
  change.any <- change.ini
  change.remove <- change.add <- F
  p.last.remove <- p.last.add <- NULL
  k <- 0
  p.last.add <- NULL
  
  while (!stop) {
    k <- k+1
    if (print) {
      if (!print.line)  cat(plain.line)
      print.line <- F
      cat("Iteration ", k, "\n\n")
    }
    mlx.loadProject(project.built)
    if (remove) {
      if (print) {
        #       cat(plain.short)
        cat("removing variability...\n")
      }
      g <- mlx.getIndividualParameterModel()$variability$id
      param0 <- names(which(!g))
      param1 <- names(which(g))
      
      change <- F
      n1 <- NULL
      stop.remove <- F
      m <- 0
      while (!stop.remove) {
        list.param1 <- names(which(!unlist(lapply(p.param1[param1], function(x) all(param1 %in% x)))))
        if (length(setdiff(list.param1, fix.param1))>0) {
          list.param1 <- setdiff(list.param1, fix.param1)
          m <- m+1
          if (print) {
            cat("\n")
            cat(dashed.short)
            cat("Step ", m, "\n")
            cat("Parameters without variability:", param0, "\n")
            cat("Parameters with variability   :", param1, "\n\n")
            if (!linearization & linearization.iter)
              cat("Criterion (linearization): ",round(BICc.built.iter,1), "\n")
            else
              cat("Criterion: ",round(BICc.built.iter,1), "\n")
          }
          list.omega1 <- paste0("omega_", list.param1)
          list.omega1 <- list.omega1[order(pop.built[list.omega1])]
          #   var.mod <- mlx.getIndividualParameterModel()$variability$id
          op.min <- BICc.min <-  pop.min <- ind.min <- NULL
          n.remove <- 0
          for (j in 1:length(list.omega1)) {
            op <- list.omega1[j]
            pj <- gsub("omega_", "", op)
            if (!identical(pj, p.last.add) | m>1) {
              if (print)
                cat("trying to remove", sprintf(paste0("%-",lpar+6,"s"), op), ": ")
              
              mlx.loadProject(project.built)
              mlx.setPopulationParameterEstimationSettings(pop.set2)
              mlx.saveProject(project.built)
              re <-   runPopEstOmega(omega=op, o.ini=pop.built[op], o.final=pop.built[op]*r.omega, 
                                     K.ini=beginIter, K.trans=nbIter, K.final=endIter, K.est=0, 
                                     p.ind=ind.built, SAEM=F)
              
              #            print(mlx.getEstimatedPopulationParameters())
              if (linearization.iter) {
                mlx.runConditionalModeEstimation()
                mlx.runLogLikelihoodEstimation(linearization=TRUE)
              } else {
                # mlx.setInitialEstimatesToLastEstimates()
                # mlx.setPopulationParameterEstimationSettings(pop.set3)
                # mlx.runPopulationParameterEstimation()
                mlx.runConditionalDistributionSampling()
                mlx.runLogLikelihoodEstimation(linearization=FALSE)
              }
              BICc1.iter <- compute.criterion(criterion, method.ll.iter, weight, pen.coef) - w.var[pj]*pen.coef[1]
              
              #  browser()
              if (print)
                cat(sprintf("%.1f", BICc1.iter), "\n")
              if (BICc1.iter - BICc.built.iter > delta[1] ) {
                p.param1[[pj]] <- param1
              } else if (BICc1.iter - BICc.built.iter < delta[2] ) {
                n.remove <- n.remove+1
                BICc.min[n.remove] <- BICc1.iter
                op.min[[n.remove]] <- op
                pmin <- mlx.getEstimatedPopulationParameters()
                pop.min[[n.remove]] <- pmin[names(pmin) != op]
                ind.min[[n.remove]] <- mlx.getEstimatedIndividualParameters()$saem
              }
            }
          }
          if (n.remove>0) {
            if (print)
              cat("\nCriterion: ",round(BICc.built,1))
            order.min <- order(BICc.min)
            j.min <- 0
            test.min <- T
            while (test.min) {
              j.min <- j.min+1
              i.min <- order.min[j.min]
              
              p.min <- gsub("omega_","",op.min[[i.min]])
              # r <- c(BICc.built, BICc.min[i.min])
              # if (length(param0) == 0)
              #   names(r) <- c("none", paste0("test_",p.min))
              # else
              #   names(r) <- c(paste0(param0, collapse='_'), paste0("test_",p.min))
              param0 <- c(param0, p.min)
              param1 <- setdiff(param1, p.min)
              update.project(project.built, project.temp, param0, param1, pop.min[[i.min]], pop.set1)
              
              if (print)
                cat("\nfitting the model with no variability on ", param0, ": ")
              
              if (BICc.min[i.min] > -Inf) {
                mlx.runPopulationParameterEstimation(parameters=ind.min[[i.min]])
                if (linearization) {
                  mlx.runConditionalModeEstimation()
                  mlx.runLogLikelihoodEstimation(linearization=TRUE)
                } else {
                  mlx.runConditionalDistributionSampling()
                  mlx.runLogLikelihoodEstimation(linearization=FALSE)
                }
                #        print(mlx.getEstimatedPopulationParameters())
                BICc0 <- compute.criterion(criterion, method.ll, weight, pen.coef)
              } else {
                BICc0 <- -Inf
              }
              if (print)
                cat(round(BICc0, 1))
              # r <- round(c(r, BICc0), 1)
              # names(r)[length(r)] <- paste0(param0, collapse='_')
              #     print(r)
              if (BICc0 < BICc.built) {
                if (print)
                  cat("\nvariability on", p.min, "removed\n")
                change <- T
                test.min <- F
                p.last.remove <- p.min
                BICc.built <- BICc0
                
                if (linearization.iter & !linearization) {
                  mlx.runConditionalModeEstimation()
                  mlx.runLogLikelihoodEstimation(linearization=TRUE)
                } else if (!linearization & linearization.iter) {
                  mlx.runConditionalDistributionSampling()
                  mlx.runLogLikelihoodEstimation(linearization=FALSE)
                }
                BICc.built.iter <- compute.criterion(criterion, method.ll.iter, weight, pen.coef)
                
                pop.built <- mlx.getEstimatedPopulationParameters()
                ind.built <- mlx.getEstimatedIndividualParameters()$saem
                # print(head(ind.built))
                # browser()
                param0.built <- param0
                param1.built <- param1
                mlx.setInitialEstimatesToLastEstimates(fixedEffectsOnly = FALSE)
                mlx.setPopulationParameterEstimationSettings(pop.set2)
                mlx.saveProject(project.built)
              } else {
                param1 <- c(param1, p.min)
                param0 <- setdiff(param0, p.min)
                
              }
              if (j.min == n.remove & test.min == T) {
                if (print)
                  cat("\nno more variability can be removed\n")
                test.min  <- F
                stop.remove=T
              }
              
            }
            
          } else {
            if (print)
              cat("\nno more variability can be removed\n")
            stop.remove=T
          }
          
        } else {
          if (print)
            cat("\nno more variability can be removed\n")
          stop.remove=T
        }
      } 
      update.project(project.built, project.built, param0.built, param1.built, pop.built, pop.set1)
      change.remove <- change
      if (change.remove)
        change.any <- T
    }
    
    if (add & (change.remove | k==1)) {
      if (print) {
        cat(plain.short)
        cat("\nadding variability...\n")
      }
      
      mlx.loadProject(project.built)
      change <- F
      stop.add <- F
      m <- 0
      p.last.add <- NULL
      while (!stop.add) {
        mod.ind <- mlx.getIndividualParameterModel()$variability$id
        param0 <-  names(mod.ind)[which(!mod.ind)]
        param1 <-  names(mod.ind)[which(mod.ind)]
        if (length(setdiff(param0, fix.param0))>0) {
          m <- m+1
          if (print) {
            cat("\n")
            cat(dashed.short)
            cat("Step ", m, "\n")
            cat("Parameters without variability:", param0, "\n")
            cat("Parameters with variability   :", param1, "\n\n")
            if (!linearization & linearization.iter)
              cat("Criterion (linearization): ",round(BICc.built.iter,1), "\n")
            else
              cat("Criterion: ",round(BICc.built.iter,1), "\n")
          }
          param0 <- setdiff(param0, fix.param0)
          
          ind.mod <- mlx.getIndividualParameterModel()
          op.min <- BICc.min <-  pop.min <- ind.min <- NULL
          n.add <- 0
          for (j in 1: length(param0)) {
            pj <- param0[j]
            if (!identical(pj, p.last.remove) | m>1) {
              op <- paste0("omega_",pj)
              if (print)
                cat("trying to add", sprintf(paste0("%-",lpar+6,"s"), op), ": ")
              mlx.loadProject(project.built)
              ind.modj <- ind.mod
              ind.modj$variability$id[pj] <- T
              
              mlx.setIndividualParameterModel(ind.modj)
              mlx.setPopulationParameterEstimationSettings(pop.set2)
              if (op %in% names(omega.value)) 
                re <-  runPopEstOmega(omega=op, r.o=1/r.omega, o.final=omega.value[op], p.ind=ind.built, 
                                      K.ini=beginIter, K.trans=nbIter, K.final=0, K.est=nbEst,
                                      SAEM=F)
              else
                re <-  runPopEstOmega(omega=op, r.o=1/r.omega, p.ind=ind.built, 
                                      K.ini=beginIter, K.trans=nbIter, K.final=0, K.est=nbEst,
                                      SAEM=F)
              popj <- mlx.getEstimatedPopulationParameters()
              omega.value[op] <- popj[op]
              
              if (linearization.iter) {
                mlx.runConditionalModeEstimation()
                mlx.runLogLikelihoodEstimation(linearization=TRUE)
              } else {
                # mlx.setInitialEstimatesToLastEstimates()
                # mlx.setPopulationParameterEstimationSettings(pop.set3)
                # mlx.runPopulationParameterEstimation()
                mlx.runConditionalDistributionSampling()
                mlx.runLogLikelihoodEstimation(linearization=FALSE)
              }
              BICc1.iter <- compute.criterion(criterion, method.ll.iter, weight, pen.coef) 
              cvj <- compute.cv(popj, pj)$cv
              if (cvj < cv.min) {
                if (print)
                  cat(sprintf("%.1f", BICc1.iter), " ; cv = ", cvj, "\n")
              } else {
                if (print)
                  cat(sprintf("%.1f", BICc1.iter), "\n")
                if (BICc1.iter - BICc.built.iter < delta[2]) {
                  n.add <- n.add+1
                  BICc.min[n.add] <- BICc1.iter
                  op.min[[n.add]] <- op
                  pmin <- mlx.getEstimatedPopulationParameters()
                  pop.min[[n.add]] <- pmin
                  ind.min[[n.add]] <- mlx.getEstimatedIndividualParameters()$saem
                }
              }
            }
          }
          if (n.add>0) {
            if (print)
              cat("\nCriterion: ",round(BICc.built,1))
            order.min <- order(BICc.min)
            j.min <- 0
            test.min <- T
            while (test.min) {
              j.min <- j.min+1
              i.min <- order.min[j.min]
              
              p.min <- gsub("omega_","",op.min[[i.min]])
              # r <- c(BICc.built, BICc.min[i.min])
              # names(r) <- c(paste0(c(fix.param0,param0), collapse='_'), paste0("test_",p.min))
              
              param1 <- c(param1, p.min)
              param0 <- setdiff(param0, p.min)
              fparam0 <- c(fix.param0,param0)
              update.project(project.built, project.temp, param0, param1, pop.min[[i.min]], pop.set1)
              
              if (print) {
                if (length(fparam0)>0)
                  cat("\nfitting the model with no variability on ", fparam0, ": ")
                else
                  cat("\nfitting the model with variability on all parameters: ")
              }
              mlx.runPopulationParameterEstimation(parameters=ind.min[[i.min]])
              if (linearization) {
                mlx.runConditionalModeEstimation()
                mlx.runLogLikelihoodEstimation(linearization=TRUE)
              } else {
                mlx.runConditionalDistributionSampling()
                mlx.runLogLikelihoodEstimation(linearization=FALSE)
              }
              
              BICc0 <- compute.criterion(criterion, method.ll, weight, pen.coef)
              if (print)
                cat(round(BICc0,1))
              # r <- round(c(r, BICc0), 1)
              # names(r)[length(r)] <- ifelse(length(fparam0) == 0, "none", paste0(fparam0, collapse='_'))
              # print(r)
              
              if (BICc0 < BICc.built) {
                if (print)
                  cat("\nvariability on", p.min, "added\n")
                p.last.add <- p.min
                change <- T
                test.min <- F
                BICc.built <- BICc0
                
                if (linearization.iter & !linearization) {
                  mlx.runConditionalModeEstimation()
                  mlx.runLogLikelihoodEstimation(linearization=TRUE)
                } else if (!linearization & linearization.iter) {
                  mlx.runConditionalDistributionSampling()
                  mlx.runLogLikelihoodEstimation(linearization=FALSE)
                }
                BICc.built.iter <- compute.criterion(criterion, method.ll.iter, weight, pen.coef)
                
                pop.built <- mlx.getEstimatedPopulationParameters()
                ind.built <- mlx.getEstimatedIndividualParameters()$saem
                param0.built <- param0
                param1.built <- param1
                io <- grep("omega_", names(pop.built))
                omega.value[names(pop.built)[io]] <- pop.built[io]
                mlx.setInitialEstimatesToLastEstimates(fixedEffectsOnly = FALSE)
                mlx.setPopulationParameterEstimationSettings(pop.set2)
                mlx.saveProject(project.built)
              } else {
                param0 <- c(param0, p.min)
                param1 <- setdiff(param1, p.min)
              }
              if (j.min == n.add & test.min == T) {
                if (print)
                  cat("\nno more variability can be added\n")
                test.min  <- F
                stop.add=T
              }
            }
          } else {
            if (print)
              cat("\nno more variability can be added\n")
            stop.add <- T
          }
        } else {
          if (print)
            cat("\nno more variability can be added\n")
          stop.add <- T
        }
      }
      #      browser()
      update.project(project.built, project.built, param0.built, param1.built, pop.built, pop.set1)
      stop <- !change
      if (change)
        change.any <- T
    } else {
      stop <- T
    }
    
    
  }
  
  g <- mlx.getIndividualParameterModel()$variability$id
  param0 <- names(which(!g))
  param1 <- names(which(g))
  
  if (change.any) {
    mlx.setPopulationParameterEstimationSettings(pop.set0)
    mlx.setScenario(scenario0)
    if (remove) {
      p.param1 <- sapply(mlx.getIndividualParameterModel()$name,function(x) NULL)
      
      list.param1 <- names(which(!unlist(lapply(p.param1[param1], function(x) all(param1 %in% x)))))
      list.param1 <- setdiff(list.param1, fix.param1)
      if (length(list.param1)>0) {
        #        list.omega1 <- paste0("omega_", list.param1)
        # p.param1 <- pop.built[paste0(list.param1,"_pop")]
        # o.param1 <- pop.built[paste0("omega_", list.param1)]
        # d.param1 <- mlx.getIndividualParameterModel()$distribution[list.param1]
        r.param1 <- compute.cv(pop.built, list.param1)
        j0 <- which(r.param1$cv < cv.min)
        if (length(j0) > 0) {
          if (print) {
            # cat("Parameters without variability:", param0, "\n")
            # cat("Parameters with variability   :", param1, "\n")
            cat(paste0("\n",plain.line,"\n"))
            print(r.param1$o)
            print(r.param1$cv)
            cat("remove variability on ",list.param1[j0], "\n\n")
          }
          
          param0 <- c(param0, list.param1[j0])
          param1 <- setdiff(param1, list.param1[j0])
          update.project(project.built, project.built, param0, param1, NULL, pop.set0)
          # mlx.saveProject(project.built)
          # mlx.runPopulationParameterEstimation(parameters=ind.built)
        }
      }
    }
    
    if (print) {
      cat(paste0("\n",plain.line))
      cat("\nFinal variance model: \n")
      cat("\nParameters without variability:", param0, "\n")
      cat("Parameters with variability   :", param1, "\n")
      cat("\nFitting the final model using the original settings... \n")
    }
    
    mlx.saveProject(project.built)
    mlx.runPopulationParameterEstimation(parameters=ind.built)
    
    
    if (linearization) {
      mlx.runConditionalModeEstimation()
      mlx.runLogLikelihoodEstimation(linearization=TRUE)
    } else {
      mlx.runConditionalDistributionSampling()
      mlx.runLogLikelihoodEstimation(linearization=FALSE)
    }
    
    if (print) {
      ll.final <- compute.criterion(criterion, method.ll, weight, pen.coef)
      ll <- formatLL(mlx.getEstimatedLogLikelihood()[[method.ll]], criterion, ll.final, weight$is.weight)
      cat(paste0("\nEstimated criteria (",method.ll,"):\n"))
      print(round(ll,2))
    }
    
    # ll <- formatLL(mlx.getEstimatedLogLikelihood()[[method.ll]], criterion, ll.final, is.weight, is.prior)
    # to.cat <- paste0("\nEstimated criteria (",method.ll,"):\n")
    # to.print <- round(ll,2)
    # print.result(print, NULL, to.cat=to.cat, to.print=to.print) 
    
  } else {
    if (print) {
      cat(paste0("\n",plain.line))
      cat("\nFinal variance model: \n")
      cat("\nParameters without variability:", param0, "\n")
      cat("Parameters with variability   :", param1, "\n")
    }
    mlx.loadProject(project)
    mlx.saveProject(project.built)
  }
  dir.temp <- gsub(".mlxtran","",project.temp)
  prop.temp <- gsub(".mlxtran",".mlxproperties",project.temp)
  if (dir.exists(dir.temp))
    unlink(dir.temp, recursive=TRUE)
  if (file.exists(project.temp))
    foo <- file.remove(project.temp)
  if (file.exists(prop.temp))
    foo <- file.remove(prop.temp)
  
  
  dt <- proc.time() - ptm
  if (print)
    cat(paste0("\ntotal time: ", round(dt["elapsed"], digits=1),"s\n\n"))
  
  options(op.original)
  return(list(project=project.built, niter=k, change=change.any, 
              variability.model=mlx.getIndividualParameterModel()$variability$id, time=dt["elapsed"]))
}

runPopEstOmega <- function(omega=NULL, o.ini=NULL, o.final=NULL, r.o=0.001, K.ini=0, 
                           K.trans=30, K.final=0, K.est=0, p.ind=NULL, SAEM=F) {
  #  browser()
  if (!is.null(o.ini) & is.null(o.final))
    o.final <- o.ini*r.o
  if (is.null(o.ini) & !is.null(o.final))
    o.ini <- o.final/r.o
  estim.ini <- ifelse(is.null(o.ini), T, F) 
  
  pop.set <- mlx.getPopulationParameterEstimationSettings()
  pop.set$simulatedannealing <- F
  estimates <- NULL
  
  if (estim.ini) {
    for (j in 1:length(omega))
      .hiddenCall(paste0("lixoftConnectors::setPopulationParameterInformation(",omega[j],"= list(initialValue = 0.3, method = 'MLE'))"))
    pop.set$nbexploratoryiterations <- 20
    mlx.setPopulationParameterEstimationSettings(pop.set)
    if (is.null(p.ind)) 
      mlx.runPopulationParameterEstimation()
    else
      mlx.runPopulationParameterEstimation(parameters=p.ind)
    o.final <- mlx.getEstimatedPopulationParameters()[omega]
    o.ini <- o.final/r.o
    if (SAEM)
      estimates <- rbind(estimates, mlx.getSAEMiterations()$estimates)
  }
  
  test.o <- T
  while (test.o) {
    for (j in 1:length(omega))
      .hiddenCall(paste0("lixoftConnectors::setPopulationParameterInformation(",omega[j],"= list(initialValue = o.ini[j], method = 'transition', 
                                                  beginIter = K.ini, nbIter = K.trans, target = o.final[j]))"))
    pop.set <- mlx.getPopulationParameterEstimationSettings()
    pop.set$nbexploratoryiterations <- K.ini + K.trans + K.final
    mlx.setPopulationParameterEstimationSettings(pop.set)
    if (is.null(p.ind)) 
      mlx.runPopulationParameterEstimation()
    else
      mlx.runPopulationParameterEstimation(parameters=p.ind)
    
    o.est <- mlx.getEstimatedPopulationParameters()[omega]
    j.o <- which( abs(o.est - o.final)/o.est > 0.01)
    if (length(j.o) > 0) {
      o.final[j.o] <- exp(log((2*o.final[j.o] + o.ini[j.o])/3))
    } else {
      test.o <- F
    }
  }
  
  if (SAEM)
    estimates <- rbind(estimates, mlx.getSAEMiterations()$estimates)
  
  if (K.est > 0) {
    proj <- mlx.getProjectSettings()$directory
    p.ind <- read.csv(file.path(proj, "IndividualParameters", "estimatedIndividualParameters.txt"))
    p.ind <- p.ind[, c(1, grep("_SAEM", names(p.ind)))]
    names(p.ind) <- gsub("_SAEM", "", names(p.ind))
    p.pop <- read.csv(file.path(proj, "populationParameters.txt"))
    
    for (j in 1:length(omega))
      .hiddenCall(paste0("lixoftConnectors::setPopulationParameterInformation(",omega[j],"= list(method = 'MLE'))"))
    j.fix <- which(o.final<o.ini)
    if (length(omega)>1 & length(j.fix)>0)
      for (j in j.fix)
        .hiddenCall(paste0("lixoftConnectors::setPopulationParameterInformation(",omega[j],"= list(method = 'FIXED'))"))
    g <- mlx.getPopulationParameterInformation()
    i <- match(p.pop$parameter, g$name)
    g[i, 'initialValue'] <- p.pop$value
    mlx.setPopulationParameterInformation(g)
    pop.set$nbexploratoryiterations <- K.est
    mlx.setPopulationParameterEstimationSettings(pop.set)
    mlx.runPopulationParameterEstimation(parameters=p.ind)
    if (SAEM)
      estimates <- bind_rows(estimates, mlx.getSAEMiterations()$estimates)
  }
  if (SAEM)
    return(estimates)
}

update.project <- function(project, project.temp, param0, param1, pop=NULL, pop.set) {
  mlx.loadProject(project)
  ind.mod <- mlx.getIndividualParameterModel()
  cb <- ind.mod$correlationBlocks$id
  if (!is.null(cb)) {
    cb <- lapply(cb, function(x) setdiff(x, param0))
    i.cb <- which(lapply(cb, length)==1)
    if (length(i.cb)>0)
      cb[[i.cb]] <- NULL
    ind.mod$correlationBlocks$id <- cb
    mlx.setIndividualParameterModel(ind.mod)
  }
  ind.mod$variability$id[param0] <- F
  ind.mod$variability$id[param1] <- T
  mlx.setIndividualParameterModel(ind.mod)
  if (!is.null(pop)) {
    pop.inf <- mlx.getPopulationParameterInformation()
    i <- match(names(pop), pop.inf$name)
    j <- which(!is.na(i))
    pop.inf[i[j], 'initialValue'] <- pop[j]
    mlx.setPopulationParameterInformation(pop.inf)
  }
  mlx.setPopulationParameterEstimationSettings(pop.set)
  mlx.saveProject(project.temp)
}


buildVar.check <- function(project, final.project, fix.param1, fix.param0, 
                           criterion, linearization, remove, add, delta, 
                           omega.set, pop.set1, pop.set2, print, weight, prior) {
  
  if (length(mlx.getIndividualParameterModel()$variability)>1)
    stop("Multiple levels of variability are not supported in this version of buildVar", call.=FALSE)
  
  if (is.null(final.project)) 
    final.project <- paste0(sub(pattern = "(.*)\\..*$", replacement = "\\1", project),"_var.mlxtran")
  if (!grepl("\\.",final.project))
    final.project <- paste0(final.project,".mlxtran")
  if (!grepl("\\.mlxtran",final.project))
    stop(paste0(final.project, " is not a valid name for a Monolix project (use the .mlxtran extension)"), call.=FALSE)
  
  
  project.temp <- gsub(".mlxtran" ,"_temp.mlxtran", project)
  project.built <- final.project
  
  oset <- list(beginIter=0, endIter=25, nbIter=25, nbEst=50, r.omega=0.001)
  if (!is.null(omega.set)) {
    if (is.null(names(omega.set)) | !(names(omega.set) %in% names(oset)) )
      stop(" 'omega.set' is not a valid list of settings", call.=FALSE)
    oset <- modifyList(oset, omega.set[intersect(names(omega.set), names(oset))])
  }
  
  p <- mlx.getIndividualParameterModel()$name
  if (!is.null(fix.param1) && !(fix.param1 %in% p))
    stop(" 'fix.param1' is not a valid list of parameter names", call.=FALSE)
  if (!is.null(fix.param0) && !(fix.param0 %in% p))
    stop(" 'fix.param0' is not a valid list of parameter names", call.=FALSE)
  if (!is.logical(linearization))
    stop(" 'linearization' should be boolean", call.=FALSE)
  if (!is.logical(print))
    stop(" 'print' should be boolean", call.=FALSE)
  if (!is.logical(remove))
    stop(" 'remove' should be boolean", call.=FALSE)
  if (!is.logical(add))
    stop(" 'add' should be boolean", call.=FALSE)
  if (!is.numeric(criterion) && !(criterion %in% c("AIC", "BIC", "BICc")))
    stop(" 'criterion' should be in {'AIC', 'BIC', 'BICc'} or be numerical > 0", call.=FALSE)
  if (is.numeric(criterion) && criterion<=0)
    stop(" 'criterion' should be in {'AIC', 'BIC', 'BICc'} or be numerical > 0", call.=FALSE)
  if (!is.numeric(delta) | (length(delta) != 2))
    stop(" 'delta' should be a 2-vector of thresholds >= 0", call.=FALSE)
  if (min(delta)< 0 )
    stop(" 'delta' should be a 2-vector of thresholds >= 0", call.=FALSE)
  
  
  if (!is.null(prior) & !is.list(prior))
    prior <- list(variance=prior)
  if (!is.null(weight) & !is.list(weight))
    weight <- list(variance=weight)
  
  if (!is.null(prior$variance) & !is.null( weight$variance)) {
    warning("Variance model: only 'weight' or 'prior' can be defined, not both. 'weight' will be used and prior will be ignored", call.=FALSE)
    prior$variance <- NULL
  }
  p.var <- prior$variance
  w.var <- weight$variance
  if (!is.null(p.var) & !is.null(w.var)) {
    warning("Variance model: only 'weight' or 'prior' can be defined, not both. 'weight' will be used and prior will be ignored", call.=FALSE)
    prior$variance <- NULL
  }
  var.model <- mlx.getIndividualParameterModel()$variability$id
  if (!is.null(p.var) & !identical(names(var.model), names(p.var))) 
    stop("prior$variance should be a named vector whose names are the names of the parameters", call.=FALSE)
  if (length(w.var)>1 & is.null(names(w.var)))
    stop("weight$variance should be a scalar or a named vector whose names are parameter names", call.=FALSE)
  if (!all(names(w.var) %in% names(var.model)))
    stop("weight$variance should be a scalar or a named vector whose names are parameter names", call.=FALSE)
  
  # if (!is.null(p.var) & is.null(w.var)) { 
  #   if (length(p.var)==1) {
  #     foo <- p.var
  #     p.var <- var.model
  #     p.var[] <- foo
  #   } else if (!identical(names(var.model), names(p.var)))
  #     stop("prior$variance should be a named vector whose names are the names of the parameters", call.=FALSE)
  # }
  # prior$variance <- p.var
  # 
  # if (is.null(p.var)) {
  #   foo <- var.model
  #   foo[] <- 1
  #   if (length(w.var)==1 & is.null(names(w.var))) 
  #     foo[] <- w.var
  #   if (length(w.var)>1 & is.null(names(w.var))) 
  #     stop("weight$variance should be a scalar or a named vector whose names are parameter names", call.=FALSE)
  #   if (!all(names(w.var) %in% names(var.model)))
  #     stop("weight$variance should be a scalar or a named vector whose names are parameter names", call.=FALSE)
  #   foo[names(w.var)] <- w.var
  #   weight$variance <- foo
  # }
  # 
  return(list(final.project=final.project, oset=oset, project.temp=project.temp, 
              project.built=project.built, weight=weight, prior=prior))
}


compute.cv <- function(vp, lp) {
  f0 <- function(x, mu=0, sigma=1) {
    f <- 1/sqrt(2*pi)/sigma*exp(-((log(x/(1-x)) - mu)^2)/(2*sigma^2))/(x*(1-x)) }
  f1 <- function(x, mu=0, sigma=1) {
    f <- x*f0(x, mu, sigma) }
  f2 <- function(x, mu=0, sigma=1) {
    f <- (x^2)*f0(x, mu, sigma) }
  
  p <- vp[paste0(lp,"_pop")]
  o <- vp[paste0("omega_", lp)]
  d <- mlx.getIndividualParameterModel()$distribution[lp]
  
  cv <- o/abs(p)
  j1 <- which(d=="logNormal")
  if (length(j1) > 0)
    cv[j1] <- sqrt(exp(o[j1]^2) - 1)
  j2 <- which(d=="logitNormal")
  for (k in seq_along(j2)) {
    lk <- mlx.getIndividualParameterModel()$limits[[lp[j2[k]]]]
    ok <- o[j2[k]]
    muk <- log((p[j2[k]] - lk[1])/(lk[2]-p[j2[k]]))
    mk <- integrate(f1, lower = 0, upper = 1, mu=muk, sigma=ok)$value
    vk <- integrate(f2, lower = 0, upper = 1, mu=muk, sigma=ok)$value - mk^2
    cv[j2[k]] <- (lk[2]-lk[1])*sqrt(vk)/(mk*(lk[2]-lk[1])+lk[1])
  }
  names(cv) <- gsub("omega", "cv", names(o))
  return(list(cv=cv, o=o, p=p))
}


# if (swap) {
#   stop.swap <- F
#   m <- 0
#   while (!stop.swap) {
#     mod.ind <- mlx.getIndividualParameterModel()$variability$id
#     param0 <- names(mod.ind)[which(!mod.ind)]
#     param1 <- names(mod.ind)[which(mod.ind)]
#     if (length(param0)>0) {
#       m <- m+1
#       if (m>0)
#         stop.swap <- T
#       if (print)
#         print(c(k, m, "C", round(BICc.built)))
#       BICc.min <- Inf
#       ind.mod <- mlx.getIndividualParameterModel()
#       for (j0 in 1: length(param0)) {
#         pj0 <- param0[j0]
#         op0 <- paste0("omega_",pj0)
#         if (!(op0 %in% names(omega.value)))
#           omega.value[op0] <- 0.3
#         for (j1 in 1: length(param1)) {
#           pj1 <- param1[j1]
#           op1 <- paste0("omega_",pj1)
#           op <-c(op0, op1)
#           if (print)
#             print(c(pj0, pj1))
#           mlx.loadProject(project.built)
#           ind.modj <- ind.mod
#           ind.modj$variability$id[pj0] <- T
#           
#           # var.modj <- ind.modj$variability$id
#           # var.modj[pj1] <- F
#           # num.mod <- BinToDec(var.modj*1)
#           # i.mod1 <- which(num.mod == list.mod1)
#           # #if (length(i.mod1)>0)
#           # list.mod1 <- c(list.mod1, num.mod)
#           # list.algo1 <- c(list.algo1, 3)
#           
#           mlx.setIndividualParameterModel(ind.modj)
#           mlx.setPopulationParameterEstimationSettings(pop.set2)
#           o.final <- omega.value[op]*c(1, r.omega)
#           o.ini <- omega.value[op]*c( r.omega, 1)
#           runPopEstOmega(omega=op, o.ini=o.ini, o.final=o.final, p.ind=ind.built, 
#                          K.ini=beginIter, K.trans=nbIter, K.final=0, K.est=nbEst)
#           
#           if (linearization) {
#             mlx.runConditionalModeEstimation()
#             mlx.runLogLikelihoodEstimation(linearization=TRUE)
#             ll1 <- mlx.getEstimatedLogLikelihood()[['linearization']]
#           } else {
#             # mlx.setInitialEstimatesToLastEstimates()
#             # mlx.setPopulationParameterEstimationSettings(pop.set3)
#             # mlx.runPopulationParameterEstimation()
#             mlx.runConditionalDistributionSampling()
#             mlx.runLogLikelihoodEstimation(linearization=FALSE)
#             ll1 <- mlx.getEstimatedLogLikelihood()[['importanceSampling']]
#           }
#           BICc1 <- ll1$BICc 
#           #if (length(i.mod1)>0)
#           #list.BICc1 <- c(list.BICc1, BICc1)
#           if (print)
#             print(BICc1)
#           if (BICc1 < BICc.min) {
#             BICc.min <- BICc1
#             op.min <- op
#             pop.min <- mlx.getEstimatedPopulationParameters()
#             pop.min <- pop.min[names(pop.min) != op[2]]
#             ind.min <- mlx.getEstimatedIndividualParameters()$saem
#           }
#         }
#       }
#       p.min <- gsub("omega_","",op.min)
#       r <- c(BICc.built, BICc.min)
#       names(r) <- c(paste0(param0, collapse='_'), paste0("test_",paste0(p.min, collapse='_')))
#       
#       param1 <- c(setdiff(param1, p.min[2]), p.min[1])
#       param0 <- c(setdiff(param0, p.min[1]), p.min[2])
#       update.project(project.built, project.built, param0, param1, pop.min, pop.set1)
#       
#       # num.mod <- BinToDec(mlx.getIndividualParameterModel()$variability$id*1)
#       # i.mod0 <- which(num.mod == list.mod0)
#       # #if (length(i.mod0)>0)
#       # list.mod0 <- c(list.mod0, i.mod0)
#       # list.algo0 <- c(list.algo0, 3)
#       
#       mlx.runPopulationParameterEstimation(parameters=ind.min)
#       if (linearization) {
#         mlx.runConditionalModeEstimation()
#         mlx.runLogLikelihoodEstimation(linearization=TRUE)
#         ll0 <- mlx.getEstimatedLogLikelihood()[['linearization']]
#       } else {
#         mlx.runConditionalDistributionSampling()
#         mlx.runLogLikelihoodEstimation(linearization=FALSE)
#         ll0 <- mlx.getEstimatedLogLikelihood()[['importanceSampling']]
#       }
#       
#       BICc0 <- ll0$BICc
#       #if (length(i.mod0)>0)
#       #list.BICc0 <- c(list.BICc0, BICc0)
#       r <- c(r, BICc0)
#       names(r)[length(r)] <- ifelse(length(param0) == 0, "none", paste0(param0, collapse='_'))
#       if (print)
#         print(r)
#       
#       if (BICc0 < BICc.built) {
#         change <<- T
#         BICc.built <<- BICc0
#         pop.built <<- mlx.getEstimatedPopulationParameters()
#         ind.built <<- mlx.getEstimatedIndividualParameters()$saem
#         param0.built <<- param0
#         param1.built <<- param1
#         mlx.setInitialEstimatesToLastEstimates(fixedEffectsOnly = FALSE)
#         mlx.setPopulationParameterEstimationSettings(pop.set2)
#         mlx.saveProject(project.built)
#       } else {
#         stop.swap <- T
#       }
#     } else {
#       stop.swap <- T
#     }
#   }
#   update.project(project.built, project.built, param0.built, param1.built, pop.built, pop.set1)
#   change.swap <- change
#   if (change.swap)
#     change.any <- T
# } 
# 
# 
