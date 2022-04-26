correlationModelSelection <- function(e0=NULL, pen.coef=NULL, nb.model=1, corr0=NULL, 
                                      seqmod=TRUE, prior=NULL, cor.list=NULL, weight=NULL) {
  
  # if (criterion=="BICc")  criterion="BIC"
  
  p.name <- mlx.getIndividualParameterModel()$name
  id <- NULL
  e <- e0 %>% arrange(rep, id)
  if (is.null(e)) {
    project.folder <- mlx.getProjectSettings()$directory
    sp.file <- file.path(project.folder,"IndividualParameters","simulatedRandomEffects.txt")
    e <- read.res(sp.file)
    # 
    # browser()
    if (is.null(e$rep)) 
      e$rep <- 1
    eta.names <- paste0("eta_",p.name)
    e <- e[c("rep","id",eta.names)]
  } 
  # else {
  #   e$rep <- mlx.getSimulatedRandomEffects()$rep
  # }
  nrep <- max(e$rep)
  N <- nrow(e)/nrep
  e$rep <- e$id <- NULL
  e.name <- names(which(mlx.getIndividualParameterModel()$variability$id))
  e.var <- paste0("eta_",e.name)
  e.var <- e.var[e.var %in% names(e)]
  e <- e[e.var]
  e <- scale(e)
  e <- as.data.frame(e)
  n.param <- ncol(e)
  alpha.cor <- 0.01
  if (n.param>1) {
    weight <- weight[e.name, e.name]
    C <- cov(e)
    C <- (C + alpha.cor*diag(diag(C)))/(1+alpha.cor)
    C0 <- diag(diag(C))
    ll <- -log(det(C0))*N/2
    df <- 0
    CC <- list(C0)
    r <- cor.list(v=1:n.param)
    for (rk in r) {
      Ck <- C0
      for (j in 1:length(rk)) {
        rjk <- rk[[j]] 
        Ck[rjk, rjk] <- C[rjk, rjk]
      }
      df <- c(df, sum(weight*(Ck !=0)))
      ll <- c(ll, -log(det(Ck))*N/2)
      CC <- c(CC, list(Ck))
    }
    
    bic <- -2*ll + pen.coef*df
    # pvl <- 1
    # for (k in (2:length(PP))) 
    #   pvl <- c(pvl, min(PP[[k]][which(PP[[k-1]]==0 & PP[[k]]>0)]))
    
    if (seqmod) {
      if (length(corr0)==0)
        bl=1
      else
        bl <- unlist(lapply(corr0,length))
      bl[which.max(bl)] <- bl[which.max(bl)]+1
      bl <- bl-1
      bm <- sum(bl*(bl+1)/2)
      bic[df>bm] <- Inf
      #     pvl[df>bm] <- 1
    }
    
    obic <- order(bic)
    #  if (obic[1] < length(bic)) {
    #   if (pvl[obic[1]+1] < p.min) {
    #     ib <- obic[1]+1
    #     obic <- c(ib, setdiff(obic, ib))
    #   }
    # }
    
    nb.model <- min(nb.model, length(bic))
    #   E <- data.frame(ll=ll, df=df, criterion=bic, pv=pvl)
    E <- data.frame(ll=ll, df=df, criterion=bic)
    # print(E)
    E <- E[obic[1:nb.model],]
    row.names(E) <- 1:nrow(E)
    # rct <- cortest(C,e,pen.bic,n.param,nrep)
    
    correlation.model <- vector(mode = "list", length = nb.model)
    for (j in 1:nb.model) {
      Cj <- CC[[obic[j]]]
      u <- rep(1,n.param)
      cm <- list()
      m <- 0
      for (k in (1:(n.param-1)))
        if (u[k]==1) {
          ik <- which(Cj[k,k:n.param] != 0) + k-1
          if (length(ik)>1) {
            m <- m+1
            cm[[m]] <- ik
            u[ik] <- 0
          }
        }
      if (length(cm)>0) {
        for (k in 1:length(cm))
          cm[[k]] <- e.name[cm[[k]]]
        correlation.model[[j]] <- lapply(cm, sort)
      }
    }
    if (nb.model==1) {
      return(correlation.model[[1]])
    } else {
      r <- list()
      for (j in (1:nb.model)) {
        r[[j]] <- list(block=correlation.model[[j]], res=E[j,])
      }
      return(r)
    }
  } else {
    return(NULL)
  }
}

#-------------------------------------------

my.dmvnorm <- function(x, mean = rep(0, p), sigma = diag(p), log = FALSE) 
{
  if (is.vector(x)) 
    x <- matrix(x, ncol = length(x))
  p <- ncol(x)
  if (!missing(mean)) {
    if (!is.null(dim(mean))) 
      dim(mean) <- NULL
    if (length(mean) != p) 
      stop("mean and sigma have non-conforming size")
  }
  if (!missing(sigma)) {
    if (p != ncol(sigma)) 
      stop("x and sigma have non-conforming size")
    if (!isSymmetric(sigma, tol = sqrt(.Machine$double.eps), 
                     check.attributes = FALSE)) 
      stop("sigma must be a symmetric matrix")
  }
  dec <- tryCatch(chol(sigma), error = function(e) e)
  if (inherits(dec, "error")) {
    x.is.mu <- colSums(t(x) != mean) == 0
    logretval <- rep.int(-Inf, nrow(x))
    logretval[x.is.mu] <- Inf
  }
  else {
    tmp <- backsolve(dec, t(x) - mean, transpose = TRUE)
    rss <- colSums(tmp^2)
    logretval <- -sum(log(diag(dec))) - 0.5 * p * log(2 * pi) - 0.5 * rss
  }
  names(logretval) <- rownames(x)
  if (log) 
    logretval
  else exp(logretval)
}


#-------------------------------------------

cor.list <- function(l=NULL, v=NULL, l0=list(), r=list()){
  if (length(l)>1)  
    r[[length(r)+1]] <- c(l0,list(l))
  nv <- length(v)
  lk <- l
  if (nv>0) {
    for (k in 1:(nv)) {
      lk <- c(l, v[k])
      if (k<nv)
        vk <- v[(k+1):nv]
      else
        vk <- NULL
      r <- cor.list(l=lk, v=vk, l0=l0, r=r)
      if (length(lk)>1)
        r <- cor.list(l=NULL, v=vk, l0=c(l0, list(lk)), r=r)
    }
  }  
  return(r)
}
