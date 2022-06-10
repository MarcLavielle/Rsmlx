#' Compute initial population PK parameters
#'
#' Use the pooled PK data to derive population PK parameters for a "standard" PK model
#' (i.e. a model of the Monolix PK library). 
#' The structural model is automatically defined using the names of the PK parameters.
#' Allowed names are: 'Tlag', 'Mtt', 'Ktr', 'ka', 'Tk0', 'V', 'V1', 'V2', 'V3', 'Q', 'Q2', 'Q3', 
#' 'Cl', 'k', 'k12', 'k21', 'k13', 'k31', 'Vm', 'Km'.
#' 
#' A Monolix project is then automatically created using these values as initial population parameters.
#' 
#' See http://rsmlx.webpopix.org/pkpopini/ for more details.
#' @param data a list with fields
#' \itemize{
#'   \item \code{dataFile}: path to a formatted data file
#'   \item \code{headerTypes}: a vector of strings
#' }
#' @param project a Monolix project
#' @param parameter a vector of strings (names of the PK parameters)
#' @param new.project name of the new Monolix project  (a default name is created if not provided)
#' @param new.dir   name of the directory where the created files are stored 
#' (default is the current working directory) )
#' @param par.ini   a vector of PK parameter values 
#' 
#' @return A list of results 
#' @examples
#' \dontrun{
#' # Create in the working directory a Monolix project for a 1 cpt model with 
#' # lag time, 0 order absorption and linear elimination
#' warf.ini1 <- pkpopini(data=warfarin, param=c("Tlag", "Tk0", "V", "Cl")) 
#' 
#' # Create in directory 'warfarin' a Monolix project called 'warfPK2.mlxtran' 
#' # for a 2 cpt model with 1st order absorption and nonlinear elimination
#' warf.ini3 <- pkpopini(data=warfarin, param=c("ka", "V", "k12", "k21", "Vm", "Km"), 
#'                       new.dir="warfarin", new.project="warfPK2.mlxtran") 
#' }
#' @export
pkpopini <- function(data=NULL, project=NULL, parameter=NULL, new.project=NULL, new.dir=NULL, par.ini=NULL) {

  if (!initRsmlx()$status)
    return()
  
  if (!is.null(project)) {
    r <- prcheck(project)
    data <- mlx.getData()[c('dataFile', 'headerTypes')]
    if (is.null(parameter))
      parameter <- mlx.getIndividualParameterModel()$name
    y.name <- mlx.getObservationInformation()$name
  } else {
    y.name <- "y"
  }
  
  if (is.null(data$administration)) {
    if ("ka" %in% parameter | "Tk0" %in% parameter)
      data$administration <- "oral"
    else
      data$administration <- "iv"
  }
  check.popini(data, parameter, par.ini)
  
  if (!is.null(new.dir) && !dir.exists(new.dir))
    dir.create(new.dir)
  
  if (length(which(c("Q","Q2") %in% parameter))>0)
    param <- "clearance"
  else
    param <- "rate"
  
  # param <- "rate"
  
  if (param=="clearance"  & is.null(par.ini)) {
    param.rate <- parameter
    param.rate <- replace(param.rate, param.rate=="Cl", "k")
    param.rate <- replace(param.rate, param.rate=="V1", "V")
    param.rate <- replace(param.rate, param.rate=="Q", "k12")
    param.rate <- replace(param.rate, param.rate=="Q2", "k12")
    param.rate <- replace(param.rate, param.rate=="V2", "k21")
    param.rate <- replace(param.rate, param.rate=="Q3", "k13")
    param.rate <- replace(param.rate, param.rate=="V3", "k31")
    p.rate <- pkpopini(parameter=param.rate, new.dir=new.dir, data=data, new.project=new.project)
    fr <- file.remove(p.rate$project)
    p1 <- p.rate$pop.ini
    par.ini <- p1[names(p1) %in% parameter]
    par.ini <- c(par.ini, V1=p1[['V']], V2=p1[['k12']]/p1[['k21']]*p1[['V']])
    if ("Cl" %in% parameter)  par.ini <- c(par.ini, Cl=p1[['k']]*p1[['V']])
    if ("Q"  %in% parameter)  par.ini <- c(par.ini, Q=p1[['k12']]*p1[['V']])
    if ("Q2" %in% parameter)  par.ini <- c(par.ini, Q2=p1[['k12']]*p1[['V']], 
                                           Q3=p1[['k13']]*p1[['V']], V3=p1[['k13']]/p1[['k31']]*p1[['V']])
    par.ini <- par.ini[parameter]
    r <- pkpopini(parameter=parameter, data=data, new.dir=new.dir, par.ini=par.ini, new.project=new.project)
    return(r)
  }
  
  if (is.null(new.project))
    new.project <- paste0("pk_",paste0(parameter, collapse=""),'.mlxtran')
  if (is.null(new.dir))
    new.dir <- "."
  new.project <- file.path(new.dir,new.project)
  
  new.model = tryCatch( {
    whichPKmodel(parameter, lib=TRUE)  }
    , error=function(e) {
      new.model <- paste0("pk_",paste0(parameter, collapse=""),'.txt')
      new.model <- file.path(new.dir,new.model)
      mlx.generatePKmodel(parameter=parameter, model=new.model)
      return(new.model)
    }
  )
  
  #mlx.newProject(data = c(data, list(observationTypes = list(y1="continuous"), nbSSDoses = 10)), modelFile = new.model)
  mlx.newProject(data = data, modelFile = new.model)
  g <- mlx.getContinuousObservationModel()
  eval(parse(text=paste0('lixoftConnectors::setErrorModel(',names(g$errorModel),'= "combined2")')))
  
  r <- readDataRsmlx(data=data)
  #r <- readDatamlx(data=data)
  pini <- compute.ini(r, parameter)
  if (param=="rate" & is.null(par.ini)) 
    popt <- pop.opt(pini)
  else
    popt <- par.ini
  
  pop.ini <- mlx.getPopulationParameterInformation()
  j.pop <- which(pop.ini$name %in% paste0(parameter,"_pop"))
  j.pop <- match(paste0(parameter,"_pop"), pop.ini$name[j.pop])
  pop.ini$initialValue[j.pop] <- popt
  mlx.setPopulationParameterInformation(pop.ini)
  
  mlx.saveProject(projectFile = new.project)
  
  return(list(pop.ini=popt, project=new.project, model=new.model, data=data))
}

check.popini <- function(data, parameter=NULL, par.ini=NULL) {
  if (is.null(data$dataFile)) 
    stop("A data file should be provided in data$dataFile", call.=FALSE)
  if (!file.exists(data$dataFile)) 
    stop(paste0(data$dataFile, " does not exists"), call.=FALSE)
  
  if (is.null(data$headerTypes)) 
    stop("Header types should be provided as a vector of strings in data$headerTypes", call.=FALSE)
  
  headerList = c('ID','TIME','AMT','ADM','RATE','TINF','Y','YTYPE', 'X','COV','CAT','OCC','MDV', 'EVID',
                 'ADDL','SS','II', 'DPT', 'AMOUNT', 'OBSERVATION', 'OBSID', 'CONTCOV', 'CATCOV', 'IGNORE',
                 'cens', 'limit', 'regressor', 'admid', 'date')
  
  if (!all(tolower(data$headerTypes) %in% tolower(headerList) ))
    stop("Valid header types should be provided", call.=FALSE)
  
  paramList <- c('Tlag', 'Mtt', 'Ktr', 'ka', 'Tk0', 'V', 'V1', 'V2', 'V3', 'Q', 'Q2', 'Q3', 'Cl', 'k', 'k12',
                 'k21', 'k13', 'k31', 'Vm', 'Km')
  
  if (!all(parameter %in% paramList ))
    stop("Valid PK parameter names should be provided", call.=FALSE)
  
  if (!is.null(par.ini) && !all(names(par.ini) %in% parameter ))
    stop("Valid initial parameters should be provided", call.=FALSE)
  
}


readDataRsmlx  <- function(data = NULL, project=NULL, out.data=FALSE, nbSSDoses=10, obs.rows=FALSE,
                           datafile=NULL, header=NULL, infoProject=NULL, addl.ss=NULL){
  id <- NULL
  observationName <- NULL
  datas=NULL
  time <- NULL
  
  if (!is.null(addl.ss)) 
    warning("addl.ss is deprecated. Use nbSSDoses instead.", call. = FALSE)
  if (!is.null(nbSSDoses))
    addl.ss <- nbSSDoses 
  
  
  
  header <- data$headerTypes
  datafile <- data$dataFile
  
  headerList      = c('ID','TIME','AMT','ADM','RATE','TINF','Y','YTYPE',
                      'X','COV','CAT','OCC','MDV','EVID','ADDL','SS','II',
                      'CENS', 'LIMIT', 'DATE')
  newList         = tolower(headerList)
  newList[3:4] <- c('amount','type')
  newHeader       = vector(length=length(header))
  ixdose = NULL
  
  header=unlist(strsplit(header, ",")) 
  header <- toupper(header)
  header[header=="DPT"]="ADM"
  header[header=="ADMID"]="ADM"
  header[header=="AMOUNT"]="AMT"
  header[header=="OBSERVATION"]="Y"
  header[header=="OBSID"]="YTYPE"
  header[header=="CONTCOV"]="COV"
  header[header=="CATCOV"]="CAT"
  header[header=="REGRESSOR"]="X"
  nlabel = length(header)
  
  icov <- icat <- iid <- iamt <- iy <- iytype <- ix <- iocc <- imdv <- NULL
  ievid <- iaddl <- iii <- iss <- iadm <- irate <- itinf <- ilimit <- icens <- NULL
  
  for (i in 1:length(headerList)) {
    hi=headerList[[i]]
    eval(parse(text= paste0("i",hi," <- NULL")))
    ih <- which(header %in% hi)
    if (length(ih)==0)
      ih=NULL
    ii <- paste0( 'i', tolower(hi) )
    eval(parse(text= paste0(ii,"=ih")))
    if (!is.null(ih))
      newHeader[ih]=newList[i]      
  }
  
  if (length(iocc)>1) 
    stop("Multiple levels of occasions are not supported", call.=FALSE)
  
  # iss <- ists
  
  ##************************************************************************
  #       Gestion du format du fichier de donnees
  #*************************************************************************  
  fileFormat=NULL
  if (is.list(infoProject)) {
    if (!is.null(infoProject$dataformat))
      fileFormat = infoProject$dataformat      
  }
  if (is.null(fileFormat)) {
    tmp=unlist(strsplit(datafile, "\\."))
    e= tmp[[length(tmp)]]
    if (e=='csv')
      fileFormat="csv"
    else
      fileFormat="space"
  }
  if (tolower(fileFormat)=="csv") {
    h1 = Rsmlx.read.table(file = datafile, sep = ",", nrows = 1)
    h2 = Rsmlx.read.table(file = datafile, sep = ";", nrows = 1)
    if (length(h1)>length(h2)) 
      delimiter=','
    else 
      delimiter=';'
    
  } else if (tolower(fileFormat)=="space") {
    delimiter=""
  } else if (tolower(fileFormat)==" ") {
    delimiter=""
  } else if (tolower(fileFormat)=="tab") {
    delimiter='\t'
  } else if (tolower(fileFormat)==";"||tolower(fileFormat)=="semicolumn"||tolower(fileFormat)=="semicolon") {
    delimiter=';'
  } else if (tolower(fileFormat)=="\t") {
    delimiter='\t'
  } else
    delimiter=','
  
  catNames<-NULL
  headerTest =Rsmlx.read.table(file = datafile, comment.char = "", sep = delimiter, nrows = 1, stringsAsFactors = FALSE)
  if(headerTest[1,1]=="#") {
    headerToUse<-headerTest[,-1]
    dataNoHeader    =  tryCatch(
      Rsmlx.read.table(file = datafile,comment.char = "#", sep=delimiter,stringsAsFactors=FALSE)
      , error=function(e) {
        error<-  geterrmessage()
        message(paste0("WARNING: reading data using delimiter '",delimiter,"' failed: ", geterrmessage()))
        return( Rsmlx.read.table(file = datafile,comment.char = "#",stringsAsFactors=FALSE))
      }      
    )    
    data<- dataNoHeader
    names(data)<- headerToUse
    
  } else {
    data = tryCatch(
      (if(!is.null(icat)) {
        colCatType<-rep(NA,length(headerTest))
        catNames <-headerTest[icat]
        names(catNames)<-NULL
        catNames<-trimws(unlist(catNames))
        colCatType[icat]<-rep("character",length(icat))
        Rsmlx.read.table(file = datafile, comment.char="", header = TRUE, sep=delimiter,colClasses = colCatType)
      } else {
        Rsmlx.read.table(file = datafile, comment.char="", header = TRUE, sep=delimiter)
      }), error=function(e) {
        error<-  geterrmessage()
        message(paste0("WARNING: reading data using delimiter '",delimiter,"' failed: ", geterrmessage()))
        return( Rsmlx.read.table(file = datafile, comment.char="", header = TRUE))
      }      
    )
  }
  
  if (out.data) {
    infoProject$delimiter <- delimiter
    infoProject$dataheader <- header
    return(list(data=data, infoProject=infoProject))
  }
  
  #---remove rows containing NA-------
  if (!is.null(iid)) {
    narowsData <- which(is.na(data[iid])) # removed in ID column only
    if(length(narowsData)>0)
      data <- data[-narowsData,]
  }
  
  #-------------------------------
  if (!is.null(iii))
    levels(data[[iii]])[levels(data[[iii]])=="."]=0
  
  #-------------------------------
  if (!is.null(iocc)) {
    names(data)[iocc] <- "occ"
    data$OCC <- as.factor(data[[iocc]])
    icat <- c(icat,which(names(data)=="OCC"))
  }
  
  if (!is.null(iid)) names(data)[iid] <- "id"
  if (!is.null(itime)) names(data)[itime] <- "time"
  
  S       = data
  S0      = names(data)
  i.new <- c(icov,icat,ix,iocc)
  newHeader[i.new] = S0[i.new]
  
  if (!is.null(iid)) {
    iobs1   = which(S[[iy]] !='.')
    if (!is.null(imdv)) {
      for (kmdv in (1: length(imdv)))
        iobs1 <- iobs1[S[iobs1,imdv[kmdv]]!=1]
    }
    if (!is.null(ievid)) {
      levels(S[,ievid])[levels(S[,ievid])=="."] <- "0"  
      S[,ievid] <- as.numeric(as.character(S[,ievid]))
      iobs1 <- iobs1[S[iobs1,ievid]==0]
    }
    i0 <- c(grep(' .',S[iobs1,iy],fixed=TRUE),grep('. ',S[iobs1,iy],fixed=TRUE))
    if (length(i0)>0)
      iobs1 <- iobs1[-i0]
    #keep  the data only for  ids which have observation
    if (!is.null(iytype)) {
      idObs <-NULL
      ytype <- factor(S[iobs1,iytype])
      l.ytype <- levels(ytype)
      if (is.null(observationName))
        observationName <- paste0("y",l.ytype)
      n.y <- length(observationName)
      for (in.y in (1:n.y)) 
        idObs <-c(idObs,as.character(S[iobs1[which(ytype==l.ytype[in.y])],iid]))
      idObs<-unique(idObs)
    } else {
      idObs<-unique(S[[iid]][iobs1])
    }
    idObsRows<-which(S[[iid]]%in%idObs)
    S <- S[idObsRows,]
    
    iduf    = unique(S[[iid]])
    iuf    = match(iduf, S[[iid]])
    ia    = sort( iuf)
    ib    = match( ia, iuf)
    
    
    foo = sort(ib)
    ic  = match(foo, ib)
    ids = iduf[ib]
    
    iu    = ia
    # idnum = as.factor(ic[idnumf])
    idnum = as.factor(S[[iid]])
  } else {
    iduf <- 1
    idnum <- as.factor(1)
    idObsRows <- nrow(S)
  }
  
  iduf = as.factor(iduf)
  N     = length(iduf)
  
  if (is.null(itime)) 
    itime=ix[1]
  
  if (!is.null(itime)) {
    t=S[[itime]]
  } else {
    # no time,  no regressor
    t = 1:nrow(S) 
  }
  nx=length(ix)
  
  if (!is.null(icat)) {
    for (j in (1:length(icat))) {
      Scatj <- S[[icat[j]]]  
      Scatj <- gsub(" ", "", Scatj, fixed = TRUE)
      S[[icat[j]]] <- as.factor(Scatj)  
    }
  }
  
  if (!is.null(iocc)) {
    socc <- S[,c(iid,itime,iocc)]
    socc['time'] <- socc[names(S)[itime]]
    socc1 <- socc[with(socc, order(id, time,occ)), 1:3 ]
    socc2 <- socc[with(socc, order(id,occ, time)), 1:3 ]
    if (!identical(socc1,socc2))
      stop("Overlapping occasions are not handled by Simulx", call.=FALSE)
    #stop("Only occasions within a same period of time are supported", call.=FALSE)
  }
  ##************************************************************************
  #       TREATMENT FIELD
  #**************************************************************************
  #  i1.evid <- NULL
  u.evid <- NULL
  if (!is.null(iamt)) {
    i1 = which(S[[iamt]] != '.')
    if (!is.null(ievid)) {
      i1 <- i1[S[i1,ievid]!=0 & S[i1,ievid]!=2]
      #      i1.evid <- i1[S[i1,ievid]==4]
    }
    i0 <- c(grep(' .',S[i1,iamt],fixed=TRUE),grep('. ',S[i1,iamt],fixed=TRUE))
    if (length(i0)>0)
      i1 <- i1[-i0]
    si1 <- S[i1,]
    si1[[iamt]] <- as.numeric(as.character(si1[[iamt]]))
    if (!is.null(iadm)) si1[[iadm]] <- as.numeric(as.character(si1[[iadm]]))
    if (!is.null(ievid)) si1[[ievid]] <- as.numeric(as.character(si1[[ievid]]))
    ixdose <- c(iamt, irate, itinf, iadm, ievid)
    #si1dose <- data.frame(sapply(si1[,ixdose], function(x) as.numeric(as.character(x))))
    si1dose <- si1[,ixdose]
    if (length(ixdose)==1)
      u=data.frame(idnum[i1], t[i1], si1dose)
    else
      u=cbind(list(idnum[i1], t[i1]), si1dose)
    names(u) = c('id',newHeader[[itime]],newHeader[ixdose])
    #
    u.addl <- NULL
    u.ss <- NULL
    if (!is.null(iaddl)) {
      levels(S[,iaddl])[levels(S[,iaddl])=="."] <- "0"  
      addl <- as.numeric(as.character(S[i1,iaddl]))
      levels(S[,iii])[levels(S[,iii])=="."] <- "0"  
      ii <- as.numeric(as.character(S[i1,iii]))
      j.addl <- which(addl>0)
      if (length(j.addl)>0) {
        for (j in (1:length(j.addl))) {
          k <- j.addl[j]
          adk <- addl[k]
          uk <- u[rep(k, adk),]
          uk$time <- u$time[k] + ii[k]*seq(1:adk)
          u.addl <- rbind(u.addl,uk)
        }
      }
      j.addl <- which(addl<0)
      if (length(j.addl)>0) {
        for (j in (1:length(j.addl))) {
          k <- j.addl[j]
          adk <- -addl[k]
          uk <- u[rep(k, adk),]
          uk$time <- u$time[k] - ii[k]*seq(1:adk)
          u.addl <- rbind(u.addl,uk)
        }
      }
    }
    i1.ss <- NULL
    if (!is.null(iss)) {
      ss <- S[i1,iss]
      ii <- as.numeric(as.character(S[i1,iii]))
      j.ss <- which(ss==1)
      if (length(j.ss)>0) {
        i1.ss <- i1[j.ss]
        for (j in (1:length(j.ss))) {
          k <- j.ss[j]
          uk <- u[rep(k, addl.ss),]
          uk$time <- u$time[k] - ii[k]*seq(1:addl.ss)
          u.ss <- rbind(u.ss,uk)
          #          u.evid <- rbind(u.evid,uk[addl.ss,c(1,2)])
        }
      }
    }
    #    if (!is.null(ievid)) {   }
    u <- rbind(u,u.addl)
    u <- rbind(u,u.ss)
    # u <- u[order(u$id,u$time),]
    
    if ("evid" %in% names(u) && any(u$evid==4))
      stop("Washout (EVID=4) is not supported", call.=FALSE)
    
    if (!is.null(u$rate) )
      u$rate <- as.numeric(as.character(u$rate))
    if (!is.null(u$tinf) )
      u$tinf <- as.numeric(as.character(u$tinf))
    if(nrow(u)) {
      datas   = c(datas,list(treatment = u))
    }
  }
  
  ##************************************************************************
  #       OBSERVATION FIELD
  #**************************************************************************
  iobs1   = which(S[[iy]] != '.')
  if (!is.null(imdv)) {
    for (kmdv in (1: length(imdv)))
      iobs1 <- iobs1[S[iobs1,imdv[kmdv]]!=1]
  }
  if (!is.null(ievid))
    iobs1 <- iobs1[S[iobs1,ievid]==0]
  i0 <- c(grep(' .',S[iobs1,iy],fixed=TRUE),grep('. ',S[iobs1,iy],fixed=TRUE))
  if (length(i0)>0)
    iobs1 <- iobs1[-i0]
  
  yvalues = data.frame(id=idnum[iobs1], time=t[iobs1], y=as.numeric(as.character(S[iobs1,iy])))
  if (!is.null(icens))
    yvalues['cens'] <- S[iobs1,icens]
  if (!is.null(ilimit))
    yvalues['limit'] <- S[iobs1,ilimit]
  if (!is.null(iytype)){ 
    ytype <- factor(S[iobs1,iytype])
    l.ytype <- levels(ytype)
    if (is.null(observationName))
      observationName <- paste0("y",l.ytype)
    n.y <- length(observationName)
    # n.y <- length(l.ytype)
    # if (length(observationName)<n.y)
    #   observationName <- paste0("y",l.ytype)
    y<- obsRows <- list()
    for (in.y in (1:n.y)) {
      idObs.i <- which(ytype==l.ytype[in.y])
      y[[in.y]] <- yvalues[idObs.i,]
      names(y[[in.y]])[3] <- observationName[in.y]
      attr(y[[in.y]],'type') <- "longitudinal"
      obsRows[[in.y]] <-  idObsRows[iobs1[idObs.i]]
    }
    names(obsRows) <- observationName
    #     yvalues$ytype <- observationName[S[[iytype]][iobs1]]
  } else {
    y <- yvalues
    if (is.null(observationName))
      observationName <- "y"
    names(y)[3] <- observationName
    attr(y,'type') <- "longitudinal"
    y <- list(y)
    obsRows <- list(idObsRows[iobs1])
    names(obsRows) <- observationName
    
  }
  names(y) <- observationName
  datas <- c(datas,y)
  # datas$observation = y
  
  ##************************************************************************
  #       REGRESSOR FIELD
  #**************************************************************************
  
  if (nx>0) {
    Sx <- S[ix]
    Dx <- data.frame(id=idnum, time=t, Sx)
    ix.num <- which(!sapply(Sx,is.numeric))
    if (!is.null(ix.num)) {
      jx <- NULL
      for (k in (ix.num)) {
        jx <- c(jx, which(Sx[[ix.num[k]]] == '.'))
      }
      # if (!is.null(jx))
      #   Dx <- Dx[-jx,]
      # for (k in (ix.num))
      #   Dx[[k+2]] <- as.numeric(as.character(Dx[[k+2]]))        
    }
    datas$regressor <- subset(Dx, !duplicated(cbind(id,time)))
  }
  
  ##************************************************************************
  #       OCCASION FIELD
  #**************************************************************************
  
  if (!is.null(iocc)) {
    S[[iocc]] <- as.numeric(as.character(S[[iocc]]))
    ov <-data.frame(id=idnum, time=t, occ=S[iocc])
    oo=ov[-2]
    u <- unique(oo)
    #io=match(data.frame(t(u)), data.frame(t(oo)))
    io=match(data.frame(t(as.numeric(rownames(u)))), data.frame(t(as.numeric(rownames(oo)))))
    ov.io <- ov[io,]
    datas$occasion <- ov.io
  } else {
    ov.io <- NULL
    if (!is.null(ievid)) {
      io <- sort(unique(c(iu,which(S[[ievid]]==4))))
      ov.io <- S[io,c(iid,itime)]
    }
    if (!is.null(u.evid)) 
      ov.io <- u.evid
    if (!is.null(ov.io)) {
      ov.id <- ov.io[,1]
      occ <- rep(1,length(ov.id))
      if (length(ov.id)>=2) {
        for (i in (2:length(ov.id))) {
          if (ov.id[i]==ov.id[i-1])
            occ[i] <- occ[i-1] + 1
          else
            occ[i] <- 1
        }     
      }
      ov.io$occ <- occ
      datas$occasion <- ov.io
      iocc <- length(names(data))+1
    }
  }
  
  if (!is.null(datas$occasion)) {
    if (length(unique(datas$occasion$occ)) == 1)
      datas$occasion <- iocc <- ov.io$occ <- NULL
  }
  
  ##************************************************************************
  #       COVARIATE FIELD
  #*************************************************************************
  
  cdf <- cdf.iov <- NULL
  nc = length(icov) + length(icat)
  if (nc>0) {
    ic <- c(icov,icat)
    cdf <- data.frame(id=iduf)
    if (!is.null(iocc)) 
      cdf.iov <- ov.io[,c(1,2)]
    k1 <- 1
    k2 <- 2
    for (k in (1:nc)) {
      if (is.null(iocc) | dim(unique(S[,c(iid,ic[k])]))[1]==N) {
        k1 <- k1+1
        cdf[[k1]] <- S[[ic[k]]][iu]
        names(cdf)[k1]=names(S)[ic[k]]    
      } else {
        k2 <- k2+1
        cdf.iov[[k2]] <- S[[ic[k]]][io]
        names(cdf.iov)[k2]=names(S)[ic[k]]    
      }
    }
    if (dim(cdf)[2]>1 & k1>1)
      datas[["covariate"]] = cdf
    if (!is.null(cdf.iov) & k2>2)
      datas[["covariate.iov"]] = cdf.iov
  }
  
  for (k in (1:length(datas))) {
    dk <- datas[[k]]
    if (is.data.frame(dk)) {
      ik <- match(dk$id,iduf)
      if (!is.null(dk$time))
        ok <- order(ik,dk$time)
      else
        ok <- order(ik)
      datas[[k]] <- dk[ok,]
      if (is.null(iid))
        datas[[k]]$id <- NULL
      imk <- match(names(datas)[k], names(obsRows))
      if (!is.na(imk))
        obsRows[[imk]] <- obsRows[[imk]][ok]
    }
  }
  if (obs.rows)
    datas$obsRows <- obsRows
  if (!is.null(iid)) {
    datas$id <- iduf  
    datas$N <- N
  }
  foo <- catNames[catNames %in% names(cdf)]
  if (length(foo) >0) 
    datas$catNames <- foo
  foo <-catNames[catNames %in% names(cdf.iov)]
  if (length(foo) >0) 
    datas$catNames.iov <- foo
  if (!is.null(datas$covariate.iov)) {
    datas[["covariate.iiv"]] <- datas[["covariate"]]
    datas[["catNames.iiv"]] <- datas[["catNames"]]
    datas[["catNames"]] <- datas[["covariate"]] <- NULL
  }
  return(datas)
}

Rsmlx.read.table <- function(file, sep = "", ...){
  
  if (!is.null(file) && !file.exists(file))
    return()
  
  
  testedDelimiters <- if (sep != "") c(sep) else c()
  testedDelimiters <- c(testedDelimiters, "\t", ",", ";", " ")
  if (is.null(file))
    testedDelimiters <- c(testedDelimiters, "")
  
  out <- table(0,0)
  count <- 1
  
  while(ncol(out) == 1 && count <= length(testedDelimiters)){
    
    if (!is.null(file))
      try(out <- read.table(file = file, sep = testedDelimiters[count], ...), silent = TRUE)
    else
      try(out <- read.table(sep = testedDelimiters[count], ...), silent = TRUE)
    
    count <- count + 1
    
  }
  
  return(if (ncol(out) == 1) NULL else out)
  
}


