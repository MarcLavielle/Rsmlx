#' Initialize Rsmlx library
#' 
#' @param errors (\emph{optional}) (\emph{boolean}) Display/Hide lixoft connectors errors (default TRUE)
#' @param warnings (\emph{optional}) (\emph{boolean}) Display/Hide lixoft connectors warnings (default TRUE)
#' @param info (\emph{optional}) (\emph{boolean}) Display/Hide lixoft connectors info (default FALSE)
#' @param path  (\emph{optional}) path to the installation directory of the Lixoft suite
#' If no path is given, the one written in the <user home>/lixoft/lixoft.ini file is used
#' @return A list:
#' \itemize{
#'   \item \code{software}: the software that is used (should be monolix with Rsmlx)
#'   \item \code{path}: the path to MonolixSuite
#'   \item \code{version}: the version of MonolixSuite that is used
#'   \item \code{status}: boolean equaling TRUE if the initialization has been successful.
#' }
#' @examples
#' \dontrun{
#' initRsmlx()  # print the info about Monolix and lixoftConnectors
#' initRsmlx(path="C:/ProgramData/Lixoft/MonolixSuite2019R1")  # use MonolixSuite 2019R1
#' }
#' @export
initRsmlx <- function(path=NULL, errors = TRUE, warnings = TRUE, info = FALSE){
  packinfo <- utils::installed.packages()
  # if (!is.element("lixoftConnectors", packinfo[,1]))
  #   stop("You need to install the lixoftConnectors package in order to use Rsmlx", call. = FALSE)
  
  if (any(as.logical(c(errors, warnings, info)))) {
    op = options()
    if (!errors) op$lixoft_notificationOptions$errors = 1
    if (!warnings) op$lixoft_notificationOptions$warnings = 1
    if (!info) op$lixoft_notificationOptions$info = 1
    options(op)
  }
  
  lixoftConnectorsState <- mlx.getLixoftConnectorsState(quietly = TRUE)
  
  if (!is.null(lixoftConnectorsState)){
    
    if (lixoftConnectorsState$software == "monolix"  && is.null(path)) {
      status=TRUE
    } else {
      status = mlx.initializeLixoftConnectors(path=path)
    }
    
  } else {
    status = mlx.initializeLixoftConnectors(path=path)
  }
  lixoftConnectorsState <- mlx.getLixoftConnectorsState(quietly = TRUE)
  lixoftConnectorsState$status <- status
  if (is.null(path))
    return(lixoftConnectorsState)
  else
    return(invisible(lixoftConnectorsState))
  
}

################################################################################
# Check project, initialize connectors, load project
################################################################################
.loadProject <- function(project, f=NULL, settings=NULL, model=NULL, paramToUse=NULL,
                         parameters=NULL, level=NULL, tests=NULL, nboot=NULL, method=NULL) {
  if (identical(substr(project,1,9),"RsmlxDemo")) {
    RsmlxDemo1.project <- RsmlxDemo2.project <- warfarin.data  <- resMonolix <- NULL
    rm(RsmlxDemo1.project, RsmlxDemo2.project, warfarin.data, resMonolix)
    eval(parse(text="data(RsmlxDemo)"))
    tmp.dir <- tempdir()
    write(RsmlxDemo1.project, file=file.path(tmp.dir,"RsmlxDemo1.mlxtran"))
    write(RsmlxDemo2.project, file=file.path(tmp.dir,"RsmlxDemo2.mlxtran"))
    write.csv(warfarin.data, file=file.path(tmp.dir,"warfarin_data.csv"), quote=FALSE, row.names = FALSE)
    project <- file.path(tmp.dir,project)
    demo <- TRUE
  } else {
    demo <- FALSE
  }
  if (demo & !is.null(f)) {
    if (f=="boot") {
      if (is.null(settings))
        res <- resMonolix$r1.boot
      else if (!is.null(settings$N) & is.null(settings$covStrat))
        res <- resMonolix$r2.boot
      else
        res <- resMonolix$r3.boot
    } else if (f=="build") {
      if (identical(model,"all") & identical(paramToUse,"all")) 
        res <- resMonolix$r1.build
      else if (identical(model,"all")) 
        res <- resMonolix$r2.build
      else 
        res <- resMonolix$r3.build
    } else if (f=="conf") {
      if (method == "fim" & level==0.90)
        res <- resMonolix$r1.conf
      else if (method == "fim" & level==0.95)
        res <- resMonolix$r2.conf
      else if (method == "proflike")
        res <- resMonolix$r3.conf
      else
        res <- resMonolix$r4.conf
    } else if (f=="cov") {
      if (identical(method,"COSSAC") & identical(paramToUse,"all")) 
        res <- resMonolix$r1.cov
      else if (identical(method,"SCM")) 
        res <- resMonolix$r2.cov
      else 
        res <- resMonolix$r3.cov
    } else if (f=="test") {
      if (length(tests)==4) 
        res <- resMonolix$r1.test
      else 
        res <- resMonolix$r2.test
    } else if (f=="set")
      res="foo"
  } else {
    
    if (!initRsmlx()$status)
      return()
    
    if (!grepl("\\.",project))
      project <- paste0(project,".mlxtran")
    
    if(!file.exists(project))
      stop(paste0("Project '", project, "' does not exist"), call.=FALSE)
    
    lp <- mlx.loadProject(project)
    if (!lp) 
      stop(paste0("Could not load project '", project, "'"), call.=FALSE)
    res <- NULL
  }
  
  return(list(project=project, demo=demo, res=res))
  #  return(project)
}
