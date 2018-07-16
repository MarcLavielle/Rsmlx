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
  
  if (!is.null(project)) {
    r <- prcheck(project)
    data <- getData()[c('dataFile', 'headerTypes')]
    if (is.null(parameter))
      parameter <- getIndividualParameterModel()$name
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
    whichPKmodel(parameter)  }
    , error=function(e) {
      new.model <- paste0("pk_",paste0(parameter, collapse=""),'.txt')
      new.model <- file.path(new.dir,new.model)
      generatePKmodel(parameter=parameter, model=new.model)
      return(new.model)
    }
  )
  
  
  
  
  newProject(data = c(data, list(observationTypes = "continuous", nbSSDoses = 10)), modelFile = new.model)
  g <- getContinuousObservationModel()
  eval(parse(text=paste0('setErrorModel(',names(g$errorModel),'= "combined2")')))
  
  r <- readDatamlx(data=data)
  pini <- compute.ini(r, parameter)
  if (param=="rate") 
    popt <- pop.opt(pini)
  else
    popt <- par.ini
  pop.ini <- getPopulationParameterInformation()
  j.pop <- which(pop.ini$name %in% paste0(parameter,"_pop"))
  j.pop <- match(paste0(parameter,"_pop"), pop.ini$name[j.pop])
  pop.ini$initialValue[j.pop] <- popt
  setPopulationParameterInformation(pop.ini)
  
  saveProject(projectFile = new.project)
  
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

