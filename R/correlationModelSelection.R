correlationModelSelection <- function(e=NULL, criterion="BIC", nb.model=1, corr0=NULL, seqcc=TRUE) {
  
  if (criterion=="BICc")  criterion="BIC"
  
  p.name <- getIndividualParameterModel()$name
  if (is.null(e)) {
    project.folder <- getProjectSettings()$directory
    sp.file <- file.path(project.folder,"IndividualParameters","simulatedRandomEffects.txt")
    e <- read.res(sp.file)
    if (is.null(e$rep)) 
      e$rep <- 1
    eta.names <- paste0("eta_",p.name)
    e <- e[c("rep","id",eta.names)]
  } 
  # else {
  #   e$rep <- getSimulatedRandomEffects()$rep
  # }
  
  nrep <- max(e$rep)
  N <- nrow(e)/nrep
  e$rep <- e$id <- NULL
  e.name <- names(which(getIndividualParameterModel()$variability$id))
  e.var <- paste0("eta_",e.name)
  e.var <- e.var[e.var %in% names(e)]
  e <- e[e.var]
  e <- as.data.frame(scale(e))
  n.param <- ncol(e)
  if (n.param>1) {
    
    if (criterion=="BIC")
      pen.bic <- log(N)
    else if (criterion=="AIC")
      pen.bic <- 2
    else 
      pen.bic <- criterion
    
    C <- cov(e)
    
    ll <- sum(my.dmvnorm(e, sigma=diag(rep(1,n.param)), log=T))/nrep
    # ll <- -0.5*(sum(e^2)/nrep + N*n.param*log(2*pi))
    df <- 0
    Ck <- diag(rep(1,n.param))
    row.names(Ck) <- names(e)
    colnames(Ck) <- names(e)
    
    A <- abs(C)
    A[lower.tri(A, diag=T)] <- 0
    b <- 1:n.param
    test <- FALSE
    CC <- list(Ck)
    k <- 1
    while (test==F){
      k <- k+1
      m <- which(A == max(A), arr.ind = TRUE)
      i1 <- which(b==b[m[1]])
      b[i1] <- b[m[2]]
      #print(b)
      ib <- which(b==b[m[1]])
      A[ib,ib] <- 0
      if (length(unique(b))==1)
        test <- TRUE
      for (j in (1:n.param)) {
        ij <- which(b==b[j])
        Ck[ij,ij] <- C[ij,ij]
      }
      llk <- sum(my.dmvnorm(e, sigma=Ck, log=T))/nrep
      dfk <- (length(which(Ck !=0))-n.param)/2
      df <- c(df, dfk)
      ll <- c(ll , llk)
      CC[[k]] <- Ck
    }
    bic <- -2*ll + pen.bic*df
    
    if (seqcc==TRUE) {
      if (length(corr0)==0)
        bl=1
      else
        bl <- unlist(lapply(corr0,length))
      bl[which.max(bl)] <- bl[which.max(bl)]+1
      bl <- bl-1
      bm <- sum(bl*(bl+1)/2)
      bic[df>bm] <- Inf
    }
    obic <- order(bic)
    nb.model <- min(nb.model, length(bic))
    E <- data.frame(ll=ll, df=df, criterion=bic)
    E <- E[order(bic)[1:nb.model],]
    row.names(E) <- 1:nrow(E)
    
    # rct <- cortest(C,e,pen.bic,n.param,nrep)
    
    correlation.model <- list()
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
      }
      correlation.model[[j]] <- cm
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