covariateModelSelection <- function(pen.coef=NULL, nb.model=1, covToTransform=NULL, covFix=NULL, 
                                    direction="both", paramToUse="all", steps=1000, p.max=1, 
                                    sp0=NULL, iter=1, correlation.model=NULL) {
  
  # correlation.model=NULL
  project.folder <- mlx.getProjectSettings()$directory
  sp.file <- file.path(project.folder,"IndividualParameters","simulatedIndividualParameters.txt")
  
  sp.df <- read.res(sp.file)
  if (is.null(sp.df$rep))
    sp.df$rep <- 1
  
  if (!is.null(sp0)) {
    nrep0 <- max(sp0$rep) 
    sp.df$rep <- sp.df$rep + nrep0
    dn <- setdiff(names(sp.df), names(sp0))
    sp0[dn] <- sp.df[dn]
    sp.df <- rbind(sp0, sp.df)
  }
  nrep <- max(sp.df$rep)
  
  ind.dist <- mlx.getIndividualParameterModel()$distribution
  param.names <- names(ind.dist)
  if (identical(paramToUse,"all"))
    paramToUse <- param.names
  n.param <- length(param.names)
  #sim.parameters <- sp.df[c("rep","id",param.names)]
  sim.parameters <- mlx.getSimulatedIndividualParameters()
  
  cov.info <- mlx.getCovariateInformation()
  cov.names <- cov.info$name
  cov.types <- cov.info$type
  j.strat <- grep("stratification",cov.types)
  if (length(j.strat)>0) {
    cov.names <- cov.names[-j.strat]
    cov.types <- cov.types[-j.strat]
  }
  j.trans <- grep("transformed",cov.types)
  #j.trans <- NULL
  j.cont <- grep("continuous",cov.types)
  cont.cov <- cov.names[j.cont]
  if (identical(covToTransform,"all"))
    covToTransform = cov.names[j.cont]
  tcov.names <- unique(c(cov.names[j.trans], setdiff(cont.cov,covToTransform)))
  cov.types <- gsub("transformed","",cov.types)
  covariates <- sp.df[c("rep","id",cov.names)]
  cov.cat <- cov.names[cov.types == "categorical"]
  covariates[cov.cat] <-  lapply(covariates[cov.cat],as.factor)
  
  indvar <- mlx.getIndividualParameterModel()$variability$id
  indvar[setdiff(param.names, paramToUse)] <- FALSE
  
  cov.model <- mlx.getIndividualParameterModel()$covariateModel
  r <- res <- list()
  eps <- 1e-15
  for (j in (1:n.param)) {
    dj <- ind.dist[j]
    nj <- names(dj)
    if (indvar[j]) {
      yj <- sp.df[nj]
      if (tolower(dj) == "lognormal") {
        yj <- log(yj+eps)
        names(yj) <- paste0("log.",nj)
      } else if (tolower(dj) == "logitnormal") {
        yj <- log((yj+eps)/(1-yj+eps))
        names(yj) <- paste0("logit.",nj)
      } else if (tolower(dj) == "probitnormal") {
        yj[[1]] <- qnorm(yj[[1]])
        names(yj) <- paste0("probit.",nj)
      } 
      
      if (length(covFix)>0) {
        cmj <- cov.model[[nj]][covFix]
        cov0 <- names(which(!cmj))
        cov1 <- names(which(cmj))
      } else {
        cov0 <- cov1 <- NULL
      }
      
      r[[j]] <- lm.all(yj, covariates, tcov.names, pen.coef=pen.coef, nb.model=nb.model, 
                       direction=direction, steps=steps, p.max=p.max, cov0=cov0, cov1=cov1, iter=iter)
      res[[j]] <- r[[j]]$res
      #   res[[j]][,c("ll","df","criterion")] <- NULL
      
      names(res[[j]]) <- gsub("log[.]","l",names(res[[j]]))
    } else {
      r[[j]] <- list(model="fixed")
      res[[j]] <- "none"
    }
    r[[j]]$p.name <- nj
  }
  
  names(res) <- param.names
  e <- as.data.frame(lapply(r[indvar], function(x) {x$model$residuals}))
  names(e) <- paste0("eta_",lapply(r[indvar], function(x) {x$p.name}))
  if (!is.null(sp.df["id"]))
    e <- cbind(sp.df["id"], e)
  if (!is.null(sp.df["rep"]))
    e <- cbind(sp.df["rep"], e)
  
  if (length(correlation.model) >0 ) {
    #    en <- as.data.frame(scale(e))
    for (ic in (1:length(correlation.model))) {
      pic <- correlation.model[[ic]]
      ipic <- match(pic, gsub("eta_","", names(e)))
      epic <- e[,ipic]
      gic <- solve(cov(epic))
      jic <- match(pic, names(ind.dist))
      #   print(head(epic))
      for (itk in 1:2) {
        jk <- 0
        for (j in jic) {
          jk <- jk+1
          dj <- ind.dist[j]
          nj <- names(dj)
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
          
          if (length(covFix)>0) {
            cmj <- cov.model[[nj]][covFix]
            cov0 <- names(which(!cmj))
            cov1 <- names(which(cmj))
          } else {
            cov0 <- cov1 <- NULL
          }
          ejc <- as.matrix(epic[,-jk])%*%matrix(gic[jk, -jk], ncol=1)/gic[jk, jk]
          yjc <- yj + ejc
          
          r[[j]] <- lm.all(yjc, covariates, tcov.names, pen.coef=pen.coef, nb.model=nb.model, 
                           direction=direction, steps=steps, p.max=p.max, cov0=cov0, cov1=cov1, iter=iter)
          res[[j]] <- r[[j]]$res
          #   res[[j]][,c("ll","df","criterion")] <- NULL
          
          names(res[[j]]) <- gsub("log[.]","l",names(res[[j]]))
          r[[j]]$p.name <- nj
          e[paste0("eta_",nj)] <- epic[,jk] <- r[[j]]$model$residuals - ejc
        }
        # print(head(epic))
      }
    }
  }
  
  
  covariate.model <- mlx.getIndividualParameterModel()$covariateModel
  covariate <- mlx.getCovariateInformation()$covariate
  js <- 0
  trs <- list()
  tr0 <- NULL
  for (k in (1:n.param)) {
    if (!identical(res[[k]],"none")) {
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
              trs[[js]] <- paste0("lixoftConnectors::addContinuousTransformedCovariate(",tr.str,")")
              tr0 <- unique(c(tr0,ckj.name))
              #eval(parse(text=tr.str))
              covariate.model[[k]][lckj] <- TRUE
            } else {
              covariate.model[[k]][ckj] <- TRUE
            }
          }
        }
      }
    }
  }
  #res <- formatCovariateModel(getIndividualParameterModel()$covariateModel)
  res <- formatCovariateModel(res)
  return(list(model=covariate.model, residuals=e, res=res, add.covariate=trs, sp=sp.df, tr0=tr0))
}

#-----------------------------------

lm.all <- function(y, x, tr.names=NULL, pen.coef=NULL, nb.model=NULL,
                   direction='both',steps = 1000, p.max=1, cov0=NULL, cov1=NULL, iter=1) {
  if (!is.null(x$id)) {
    N <- length(unique(x$id))
    nrep <- nrow(x)/N
  } else {
    x$id <- 1
    N <-  nrow(x)
    nrep <- 1
  }
  if (is.null(x$rep)) x$rep <- 1
  
  if (p.max<=1) {
    nx <- setdiff(names(x), c("id","rep"))
    yc <- y[[1]][order(x$rep, x$id)]
    yc <- rowMeans(matrix(yc, ncol=nrep))
    xc <- subset(x,rep==1)
    xc <- xc[order(xc$id),,drop=FALSE]
    #  xc <- x[seq(1,nrow(x),by=nrep),,drop=FALSE]
    lm0 <- lm(yc ~1)
    nxc <- setdiff(nx, cov0)
    pjc <- NULL
    for (nc in nxc) {
      lmc <- lm(yc ~ xc[[nc]])
      pjc <- c(pjc, signif(anova(lm0, lmc)$`Pr(>F)`[2],4))
    }
    # if (length(which(pjc<p.max)) == 0)
    #   list.c <- which(pjc<max(pjc))
    # else
    #   pjc <- p.adjust(pjc, method="BH")
    list.c <- which(pjc>p.max)
    # browser()
    
    cov0 <- c(cov0, nxc[list.c])
    direction <- ifelse(length(setdiff(nx, cov0))<=10,"full",direction)
    # x <- x[, -which(names(x) %in% nx[list.c])]
  } else
    list.c <- NULL
  
  x$id <- x$rep <- NULL
  nx <- ncol(x)
  l <- x
  s <- rep("0:1", nx)
  names(s) <- names(x)
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
    # x[,j.num] <- scale(x[j.num], scale=TRUE)
    # l[,j.num] <- scale(l[j.num], scale=TRUE)
    l[,j.num] <- scale(l[j.num], scale=FALSE)
  } else
    j0.num <- vector(length=0)
  
  #print(direction)
  #############% stepAIC/stepBIC
  if (direction != 'full') {
    # create a data frame with all possible covariates (all transformations)
    l_data = data.frame(y=y[[1]])
    # not use non transformed covariates
    j0.num = j0.num[!names(j0.num) %in% tr.names]
    
    # not drop data frame to a vector!
    
    l_data = cbind.data.frame(l_data,x,l[,j0.num,drop=FALSE])
    
    # apply stepAIC
    llk = tryCatch( {
      # stepAIC forward, backward, both
      f.sature <- "y ~ ."
      if (length(cov0)>0)
        f.sature <- paste0(f.sature, "-", paste(cov0,collapse="-"))
      f.sature <- as.formula(f.sature)
      model.sature=lm(f.sature,l_data)
      f.cst <- "y ~ 1"
      if (length(cov1)>0)
        f.cst <- paste0(f.cst, "+", paste(cov1,collapse="+"))
      f.cst <- as.formula(f.cst)
      model.cst=lm(f.cst,l_data)
      #  model.sature=lm(y~.,l_data)
      # model.cst=lm(y~1,data=l_data) 
      
      
      if(direction=='backward'){
        lm.sat=stepAIC(model.sature, direction='backward', trace = FALSE, k=nrep*pen.coef,
                       scope=list(upper=model.sature, lower=model.cst), steps=steps)
      }
      else {
        lm.cst=stepAIC(model.cst, direction=direction, trace = FALSE, k=nrep*pen.coef,
                       scope=list(upper=model.sature, lower=model.cst), steps=steps)
      }
    }
    , error=function(e) {
      print('Error in stepAIC')
      return(-Inf)        }      
    )
    ## after tryCatch
    # Which covariates are used (written in 1-row data frame G)
    Gnames = names(x)
    G <- data.frame(matrix(ncol = length(Gnames), nrow = 1))
    colnames(G) <- Gnames
    
    # exact used covariates (with original names of cats), fill in G
    usedcovariates = names(llk$model)[-1] #names(coef(llk)) 
    for (i in names(G)){
      if (i %in% usedcovariates ) G[1,i]=1 else G[1,i]=0
      if (paste0('log.',i) %in% usedcovariates ) G[1,i]=G[1,i]+2  
      # = 3 if both log and non-log are used
    }
    
    ll = logLik(llk)/nrep
    df = length(coef(llk))-1 # except Intercept
    criterion = -2*ll + df*pen.coef
    
    res <- data.frame(ll=round(ll,digits=3), df=df,  criterion=round(criterion,digits=3))
    j0.num <- j0.num[!(names(j0.num) %in% tr.names)]
    if (length(j0.num)>0){
      G[names(l)[j0.num]] <- 0
      i2 <- (G[names(x)[j0.num]]==2)
      G[names(l)[j0.num]][i2] <- 1
      G[names(x)[j0.num]][i2] <- 0
    }
    res <- cbind(G==1, res)
    #  res <- cbind(G, res)
    
    #  browser()
    return(list(model=llk, res=res))
  } 
  
  if (length(tr.names)>0){
    j.newc <- which((names(x) %in% tr.names))
    if (length(j.newc)>0)
      s[j.newc] <- "0:1"
  }
  
  if (!is.null(cov0))
    s[cov0] <- "0"
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
  names(G) <- names(x)
  
  if (length(cov0)>0) {
    i0 <- which(rowSums(G[cov0])==0)
    G <- G[i0,]
  }
  
  if (length(cov1)>0) {
    i1 <- which(rowSums(G[cov1]==1)==length(cov1))
    G <- G[i1,]
  }
  
  ng <- nrow(G)
  d  <- ncol(G)
  #  print(c(ng,d))
  
  ll <- df <- bic <- bic.cor <- NULL
  corb <- log(iter^2/(iter^2+3))
  
  #iop.mean <- ((iter ==2))
    iop.mean <- F
  if (iop.mean) {
    yg <- colMeans(matrix(y[[1]], nrow=nrep))
    xg <- x[seq(1,nrow(x),by=nrep),]
    if (nrow(l)>0)
      lg <- l[seq(1,nrow(l),by=nrep),]
    nrepg <- 1
    #  pen.coef <- pen.coef + log(nrep)
  } else {
    yg <- y[[1]]
    xg <- x
    lg <- l
    nrepg <- nrep
  }
  
  for (k in 1:ng) {
    xk <- data.frame(y=yg)
    Gk <- G[k,,drop=FALSE]
    j1 <- which(Gk==1)
    if (length(j1)>0)
      xk[names(x)[j1]] <- xg[j1]
    j2 <- which(Gk==2)
    if (length(j2)>0)
      xk[names(l)[j2]] <- lg[j2]
    llk= tryCatch( {
      lmk <- lm(y ~ ., data=xk)
      logLik(lmk)[1]/nrepg }
      , error=function(e) {
        return(-Inf)        }      
    )    
    #dfk <- sum(Gk>0)
    dfk <- lmk$rank
    bick <- -2*llk + pen.coef*dfk + any(Gk==2)*0.001
    ll <- c(ll , llk)
    df <- c(df, dfk)
    bic <- c(bic , bick)
    #  bic.cor <- c(bic.cor , bick-corb*dfk)
    bic.cor <- c(bic.cor , bick)
  }
  
  bic <- round(bic.cor, digits=3)
  i0 <- rep(1,ng)
  mG <- ncol(G)
  for (k in seq_len(ng-1)) {
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
  res <- data.frame(ll=round(ll,digits=3), df=df, criterion=bic)
  # print(res %>% group_by(df) %>% summarize(max(ll)))
  # browser()
  res <- res[i0==1,]
  G <- G[i0==1,,drop=FALSE]
  bic <- bic[i0==1]
  
  #ill <- which(!duplicated(ll))
  eval(parse(text=paste0(names(y)," <- y[[1]]")))
  obic <- order(bic)
  k.min <- obic[1]
  # if (mG==1) 
  #   Gkmin <- G[k.min]
  # else
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
  #  if (names(y)=="log.Cl") {print(cbind(G,res)); browser()}
  nb.model0 <- min(nb.model, length(bic))
  res <- res[obic[1:nb.model0],]
  G <- G[obic[1:nb.model],,drop=FALSE]
  j0.num <- j0.num[!(names(j0.num) %in% tr.names)]
  if (length(j0.num)>0){
    G[names(l)[j0.num]] <- 0
    i2 <- (G[names(x)[j0.num]]==2)
    G[names(l)[j0.num]][i2] <- 1
    G[names(x)[j0.num]][i2] <- 0
  }
  if (nb.model>nb.model0)
    res[(nb.model0+1):nb.model,c("ll","df","criterion")] <- NA
  res <- cbind(G==1, res)
  if (nb.model==1)
    res[,c("ll","df","criterion")] <- NULL
  row.names(res) <- 1:nrow(res)
  
  return(list(model=lm.min, res=res))
}
