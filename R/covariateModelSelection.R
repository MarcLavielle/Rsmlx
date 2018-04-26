covariateModelSelection <- function(penalization="BIC", nb.model=1, trans.cov="all", param.cov="all") {
  
  project.folder <- getProjectSettings()$directory
  sp.file <- file.path(project.folder,"IndividualParameters","simulatedIndividualParameters.txt")
  sp.df <- read.csv(sp.file)
  if (is.null(sp.df$rep))
    sp.df$rep <- 1
  nrep <- max(sp.df$rep)
  
  ind.dist <- getIndividualParameterModel()$distribution
  param.names <- names(ind.dist)
  if (identical(param.cov,"all"))
    param.cov <- param.names
  n.param <- length(param.names)
  #sim.parameters <- sp.df[c("rep","id",param.names)]
  sim.parameters <- getSimulatedIndividualParameters()
  
  cov.info <- getCovariateInformation()
  cov.names <- cov.info$name
  cov.types <- cov.info$type
  j.trans <- grep("transformed",cov.types)
  j.cont <- grep("continuous",cov.types)
  cont.cov <- cov.names[j.cont]
  if (identical(trans.cov,"all"))
    trans.cov = cov.names[j.cont]
  else if (identical(trans.cov,"none"))
    trans.cov=NULL
  tcov.names <- unique(c(cov.names[j.trans], setdiff(cont.cov,trans.cov)))
  cov.types <- gsub("transformed","",cov.types)
  covariates <- sp.df[c("rep","id",cov.names)]
  cov.cat <- cov.names[cov.types == "categorical"]
  covariates[cov.cat] <-  lapply(covariates[cov.cat],as.factor)
  
  indvar <- getIndividualParameterModel()$variability$id
  indvar[setdiff(param.names, param.cov)] <- FALSE
  
  r <- res <- list()
  for (j in (1:n.param)) {
    dj <- ind.dist[j]
    nj <- names(dj)
    if (indvar[j]) {
      #    print(nj)
      yj <- sp.df[nj]
      if (tolower(dj) == "lognormal") {
        yj <- log(yj)
        names(yj) <- paste0("log.",nj)
      } else if (tolower(dj) == "logitnormal") {
        yj <- log(yj/(1-yj))
        names(yj) <- paste0("logit.",nj)
      } else if (tolower(dj) == "probitnormal") {
        yj[[1]] <- qnorm(yj[[1]])
        names(yj) <- paste0("probit.",nj)
      } 
      
      r[[j]] <- lm.all(yj,covariates,tcov.names,penalization=penalization,nb.model=nb.model)
      res[[j]] <- r[[j]]$res
    } else {
      r[[j]] <- list(model="fixed")
      res[[j]] <- "none"
    }
    r[[j]]$p.name <- nj
  }
  names(res) <- param.names
  e <- as.data.frame(lapply(r[indvar], function(x) {x$model$residuals}))
  names(e) <- paste0("e_",lapply(r[indvar], function(x) {x$p.name}))
  
  covariate.model <- getIndividualParameterModel()$covariateModel
  covariate <- getCovariateInformation()$covariate
  js <- 0
  trs <- list()
  for (k in (1:n.param)) {
    covariate.model[[k]][1:length(covariate.model[[k]])] <- FALSE
    if (indvar[k]) {
      ck <- attr(r[[k]]$model$terms,"term.labels")
      if (length(ck)>0) {
        for (j in (1:length(ck))) {
          ckj <- ck[j]
          if (identical(substr(ckj,1,4),"log.")) {
            js <- js+1
            ckj.name <- sub("log.","",ckj)
            covkj <- covariate[[ckj.name]]
            lckj <- paste0("l",ckj.name)
            tr.str <- paste0(lckj,' = "log(',ckj.name,"/",signif(mean(covkj),digits=2),')"')
            trs[[js]] <- paste0("addContinuousTransformedCovariate(",tr.str,")")	
            #eval(parse(text=tr.str))
            covariate.model[[k]][lckj] <- TRUE
          } else {
            covariate.model[[k]][ckj] <- TRUE
          }
        }
      }
    }
  }
  
  return(list(model=covariate.model, residuals=e, res=res, add.covariate=trs))
}

#-----------------------------------

lm.all <- function(y, x, tr.names=NULL, penalization=penalization, nb.model=nb.model) {
  if (!is.null(x$id)) {
    N <- length(unique(x$id))
    nrep <- nrow(x)/N
  } else {
    N <-  nrow(x)
    nrep <- 1
  }
  x$id <- x$rep <- NULL
  nx <- ncol(x)
  l <- x
  s <- rep("0:1", nx)
  j.num <- which(!sapply(x, is.factor))
  if (length(j.num)>0) {
    if (length(j.num)==1)
      j0 <- which(min(x[,j.num])>0)
    else
      j0 <- which(sapply(x[,j.num],min)>0)
    j0.num <- j.num[j0]
    s[j0.num] <- "0:2"
    l[,j0.num] <- log(x[,j0.num])
    names(l)[j0.num] <- paste0("log.",names(x)[j0.num])
  }
  x[,j.num] <- scale(x[j.num], scale=FALSE)
  l[,j.num] <- scale(l[j.num], scale=FALSE)
  
  if (length(tr.names)>0){
    j.newc <- which((names(x) %in% tr.names))
    if (length(j.newc)>0)
      s[j.newc] <- "0:1"
  }
  
  s <- paste(s,collapse=",")
  s <- paste0("G <- expand.grid(", s, ")")
  eval(parse(text=s))
  
  if (length(tr.names)>0){
    jc <- which(!(names(x) %in% tr.names))
    c.names <- names(x)[jc]
    for (k in 1:length(j.newc)) {
      jk <- which(paste0("l",c.names)==tr.names[k])
      if (length(jk)>0) {
        j <- which(names(x) == c.names[jk])
        tj <- which(names(x) == tr.names[k])
        i0 <- which(G[,j]>0 & G[,tj]>0)
        G <- G[-i0,]
      }
    }
  }
  ng <- nrow(G)
  d  <- ncol(G)
  names(G) <- names(x)
  
  ll <- df <- bic <- NULL
  if (penalization=="BIC")
    pen.bic <- log(N)
  else if (penalization=="AIC")
    pen.bic <- 2
  else 
    pen.bic <- penalization
  
  for (k in 1:ng) {
    xk <- data.frame(y=y[[1]])
    Gk <- G[k,]
    j1 <- which(Gk==1)
    if (length(j1)>0)
      xk[names(x)[j1]] <- x[j1]
    j2 <- which(Gk==2)
    if (length(j2)>0)
      xk[names(l)[j2]] <- l[j2]
    llk= tryCatch( {
      lmk <- lm(y ~ ., data=xk)
      logLik(lmk)[1]/nrep }
      , error=function(e) {
        return(-Inf)        }      
    )    
    dfk <- sum(Gk>0)
    bick <- -2*llk + pen.bic*dfk #+ any(Gk==2)*0.001
    ll <- c(ll , llk)
    df <- c(df, dfk)
    bic <- c(bic , bick)
  }
  
  bic <- round(bic, digits=3)
  i0 <- rep(1,ng)
  mG <- ncol(G)
  for (k in 1:(ng-1)) {
    if (i0[k]==1) {
      ik <- which(bic[(k):ng]==bic[k]) + k-1
      sk <- .rowSums(G[ik,]==2, n=length(ik), m=mG)
      ik0 <- ik[which(sk==0)]
      if (length(ik0)==0)
        ik0 <- ik[order(sk)[1]]
      i0[ik] <- 0
      i0[ik0] <- 1
      i0[k] <- 1
    }
  }
  res <- data.frame(ll=round(ll,digits=3), df=df, criteria=bic)
  res <- res[i0==1,]
  G <- G[i0==1,]
  bic <- bic[i0==1]
  
  #ill <- which(!duplicated(ll))
  eval(parse(text=paste0(names(y)," <- y[[1]]")))
  obic <- order(bic)
  k.min <- obic[1]
  if (mG==1) 
    Gkmin <- G[k.min]
  else
    Gkmin <- G[k.min,]
  
  j1 <- which(Gkmin==1)
  j2 <- which(Gkmin==2)
  if (length(j1)>0) {
    for (k in (1:length(j1)))
      eval(parse(text=paste0(names(x)[j1[k]]," <- x[[j1[k]]]")))
  }
  if (length(j2)>0) {
    for (k in (1:length(j2)))
      eval(parse(text=paste0(names(l)[j2[k]]," <- l[[j2[k]]]")))
  }
  list.x <- c("1",names(x)[j1],names(l)[j2])
  form1 <- paste0(names(y), "~",  paste(list.x, collapse = "+")) 
  eval(parse(text=paste0("lm.min <- lm(",form1,")")))
  lm.min$covsel=Gkmin
  
  res <- cbind(G, res)
  nb.model <- min(nb.model, length(bic))
  res <- res[obic[1:nb.model],]
  row.names(res) <- 1:nrow(res)
  
  return(list(model=lm.min, res=res))
}
