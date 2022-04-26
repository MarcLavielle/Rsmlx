errorModelSelection <- function(project=NULL, pen.coef=NULL, nb.model=1) {
  
  if (!is.null(project))
    mlx.loadProject(project)
  
  go <- mlx.getObservationInformation()
  go$name <- names(go$mapping)
  go$prediction <- mlx.getContinuousObservationModel()$prediction
  go$type <- go$type[1:length(go$name)]
  ic <- which(go$type=="continuous")
  y.name <- go$name[ic]
  y.obs <- go[y.name]
  y.pred <- go$prediction[y.name]
  n.out <- length(ic)
  pred <- getSimulatedPredictions()
  
  # while( (sum(apply(pred[[1]],2,is.nan))+sum(apply(pred[[1]],2,is.na)) + sum(apply(pred[[1]],2,is.infinite))) > 0) {
  #   browser()
  #   runConditionalDistributionSampling()
  #   pred <- getSimulatedPredictions()
  # }
  
  # obs.model <- mlx.getContinuousObservationModel()
  # obs.names <- mlx.getData()$observationNames
  # fit.info <- mlx.getFit()
  # for (noj in seq_along(obs.names)) {
  #   nojk <- grep(obs.names[noj],fit.info$data)
  #   if (length(nojk)>0)
  #     obs.names[noj] <- fit.info$model[nojk]
  # }
  # i.contObs <- which(mlx.getData()$observationTypes=="continuous") 
  # i.contModel <- which(names(obs.model$prediction) %in% obs.names[i.contObs])
  # n.out <- length(i.contModel)
  #d <- mlx.getObservationInformation()
  
  
  if (nb.model==1) {
    res.errorModel <- NULL
    for (i.out in (1:n.out)) {

      # name.predi <- obs.model$prediction[[i.contModel[i.out]]]
      # name.obsi <- names(obs.model$prediction[i.contModel[i.out]])
      # y.obsi <- d[[name.obsi]][[name.obsi]]
      # y.predi <- pred[[name.predi]][[name.predi]]
      y.obsi <- y.obs[[y.name[i.out]]][[y.name[i.out]]]
      y.predi <- pred[[y.pred[i.out]]][[y.pred[i.out]]]
      test <- T
      kt <- 0
      while(test | kt<10) {
        kt <- kt+1
        resi <- tryCatch(computeBIC(y.obs=y.obsi, y.pred=y.predi, pen.coef=pen.coef[i.out], nb.model=nb.model),
                         error = function(e) e)
        if (inherits(resi, "error")) {
          print("error !!!")
          mlx.saveProject("foo.mlxtran")
          mlx.loadProject("foo.mlxtran")
          pred <- getSimulatedPredictions()
          y.predi <- pred[[y.pred[i.out]]][[y.pred[i.out]]]
        } else {
          test <- F
        }
      }
      res.errorModel <- c(res.errorModel, as.character(resi[['error.model']]))
      names(res.errorModel)[i.out] <- y.name[i.out]
    }
  } else {
    res.errorModel <- list()
    for (i.out in (1:n.out)) {
      # name.predi <- obs.model$prediction[[i.contModel[i.out]]]
      # name.obsi <- names(obs.model$prediction[i.contModel[i.out]])
      # y.obsi <- d[[name.obsi]][[name.obsi]]
      # y.predi <- pred[[name.predi]][[name.predi]]
      y.obsi <- y.obs[[y.name[i.out]]][[y.name[i.out]]]
      y.predi <- pred[[y.pred[i.out]]][[y.pred[i.out]]]
      resi <- tryCatch(computeBIC(y.obs=y.obsi, y.pred=y.predi, pen.coef=pen.coef[i.out], nb.model=nb.model),
                       error = function(e) e)
      if (inherits(resi, "error")) {
        print("error !!!")
        mlx.saveProject("foo.mlxtran")
        mlx.loadProject("foo.mlxtran")
        pred <- getSimulatedPredictions()
        y.predi <- pred[[y.pred[i.out]]][[y.pred[i.out]]]
        resi <- computeBIC(y.obs=y.obsi, y.pred=y.predi, pen.coef=pen.coef[i.out], nb.model=nb.model)
      }
      res.errorModel[[i.out]] <- resi
      names(res.errorModel)[i.out] <- y.name[i.out]
    }
  }
  return(res.errorModel)
}

e.min1 <- function(x,y.pred,y.obs) {
  sigma2 <- (x[1]+x[2]*y.pred)^2
  e <- sum((y.obs-y.pred)^2/sigma2) + sum(log(sigma2))
  return(e)
}

e.min2 <- function(x,y.pred,y.obs) {
  sigma2 <- x[1]^2+(x[2]*y.pred)^2
  e <- sum((y.obs-y.pred)^2/sigma2) + sum(log(sigma2))
  return(e)
}

computeBIC <- function(y.obs, y.pred, pen.coef, nb.model) {
  
  y.pos <- (min(y.obs) >= 0) 
  nrep <- length(y.pred)/length(y.obs)
  y.obs <- rep(y.obs, nrep)
  
  if (y.pos) {
    iy <- (y.pred>0)
    y.obs <- y.obs[iy]
    y.pred <- y.pred[iy]
  }
  
  # if (sum(is.na(y.pred))+sum(is.nan(y.pred))>0)
  #   browser()
  a.cons <- sqrt(mean((y.obs-y.pred)^2, na.rm=T))
  
  x.min <- nlm(e.min2,c(a.cons/10,0.2),y.pred,y.obs)
  a.comb2 <- abs(x.min$estimate[1])
  b.comb2 <- abs(x.min$estimate[2])
  
  x.min <- nlm(e.min1,c(a.comb2,b.comb2),y.pred,y.obs)
  a.comb1 <- x.min$estimate[1]
  b.comb1 <- x.min$estimate[2]

  x.min <- nlm(e.min2,c(a.comb1,b.comb1),y.pred,y.obs)
  a.comb2 <- abs(x.min$estimate[1])
  b.comb2 <- abs(x.min$estimate[2])

    if (y.pos == TRUE ) {
    b.prop <- sqrt(mean(((y.obs/y.pred-1)^2)))
    
    a.expo <- sqrt(mean((log(y.obs)-log(y.pred))^2))
    
    error.model=c("constant", "proportional", "combined1","combined2","exponential")
    sigma2 <- cbind(a.cons^2, (b.prop*y.pred)^2, (a.comb1+b.comb1*y.pred)^2, 
                    a.comb2^2+(b.comb2*y.pred)^2, a.expo^2)
    sigma2 <- cbind(a.cons^2, (b.comb2*y.pred)^2, (a.comb1+b.comb1*y.pred)^2, 
                    a.comb2^2+(b.comb2*y.pred)^2, a.expo^2)
    df <- c(1, 1, 2, 2, 1)
  } else {
    error.model=c("constant","combined1","combined2")
    sigma2 <- cbind(a.cons^2, (a.comb1+b.comb1*y.pred)^2, a.comb2^2+(b.comb2*y.pred)^2)
    df <- c(1, 2, 2)
  }
  
  ll <- pen <- bic <- NULL
    for (k in 1:length(error.model)) {
    if (error.model[k] == "exponential")
      ll <- c(ll, - 0.5*sum( (log(y.obs)-log(y.pred))^2/(a.expo^2) + log(2*pi*(a.expo^2)) + 2*log(y.obs) )/nrep)
    else
      ll <- c(ll, - 0.5*sum( (y.obs-y.pred)^2/sigma2[,k] + log(2*pi*sigma2[,k]) )/nrep)
  }
  bic <- -2*ll + pen.coef*df
  
  E <- data.frame(error.model=error.model,
                  ll=ll, df=df, criterion= bic)
  nb.model <- min(nb.model, length(bic))
  E <- E[order(bic)[1:nb.model],]
  row.names(E) <- 1:nrow(E)
  return(E)
}


