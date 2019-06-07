#' Easy tuning of the settings of a Monolix project
#'
#' Use a single accuracy level, between 1 and 9, to automatically tune all the settings  
#' of a Monolix project.
#' When the accuray level is equal to 1, the algorithms are very fast but the results may be
#' not precise.  When the accuray level is equal to 9, the algorithms are slow but the results 
#' are accurate. Default Monolix settings are obtained with level=5.
#' @param project a string: a Monolix project (the loaded project if NULL)
#' @param new.project  a string: the new created Monolix project (default is the original project)
#' @param level  an integer between 1 and 9 (default=5)
#' @examples
#' \dontrun{
#' # RsmlxDemo1.mlxtran is a Monolix project for modelling the PK of warfarin.

#' # All settings of the project are set so that algorithms used by Monolix converge as 
#' # quickly as possible possible:
#' setSettings(project="RsmlxDemo1.mlxtran", level=1)
#'
#' # A new project will be created with settings set in order to obtain the most 
#' # precise results possible:
#' new.project= file.path(tempdir(),"RsmlxDemoNew.mlxtran")
#' setSettings(project="RsmlxDemo1.mlxtran", new.project=new.project, level=9)
#' 
#' # See http://rsmlx.webpopix.org/userguide/setSettings/ for detailed examples of use of setSettings
#' # Download the demo examples here: http://rsmlx.webpopix.org/installation

#' }
#' @export

setSettings  <- function(project=NULL, new.project=NULL, level=5) {
  
  if (!initRsmlx())
    return()
  
  if (!is.null(project)){
    r <- prcheck(project, f="set")
    if (r$demo)
      return(invisible(FALSE))
    project <- r$project
  }
  
  if (!is.null(project) & is.null(new.project)) 
    new.project <- project
  n <- 9
  
  level <- round(level)
  if (level<1 | level>n)
    stop(paste0("level should be an integer beween 1 and ",n), call. = FALSE)
  
  tr.str <- paste0("lixoftConnectors::setPopulationParameterEstimationSettings(exploratoryautostop=TRUE, smoothingautostop=TRUE)")
  eval(parse(text=tr.str))
  tr.str <- paste0("lixoftConnectors::setGeneralSettings(autochains=TRUE)")
  eval(parse(text=tr.str))
  
  #  mlx.setGeneralSettings(autochains=TRUE)
  #  mlx.setPopulationParameterEstimationSettings(exploratoryautostop=TRUE, smoothingautostop=TRUE)
  
  set.list <- c(
    "minindivforchains", 			
    "nbexploratoryiterations",   	
    "nbsmoothingiterations",    	
    "exploratoryinterval", 		
    "smoothinginterval",			
    "tauomega",				
    "tauerrormodel",
    "nboptimizationiterations",
    "optimizationtolerance",
    "nbfixediterations",
    "nbminiterations",
    "ratio",
    "nbsimulatedparameters",
    "miniterations",
    "maxiterations",
    "nboptimizationiterationsmode",
    "optimizationtolerancemode")
  
  p <- matrix(c(
    20,  50,  200,
    500, 500, 1000,
    50, 200, 1000,
    50, 150, 300,
    20, 50, 200,
    0.9, 0.95, 0.98,
    0.9, 0.95, 0.98,
    10, 20, 40,
    0.0005, 0.0001, 0.00005,
    5000, 10000, 50000,
    20, 50, 200,
    0.1, 0.05, 0.01,
    5, 10, 25,
    50, 50, 200,
    200, 200, 400,
    100, 200, 300,
    1e-5, 1e-6, 1e-7
  ), ncol=3, byrow=TRUE)
  
  
  np <- dim(p)[1]
  ind.int <- rep(TRUE, np)  
  ind.int[c(6, 7, 9, 12, 17)] <- FALSE
  
  x <- c(1, round((n+1)/2), n)
  
  g1 <- mlx.getGeneralSettings()
  g2 <- mlx.getPopulationParameterEstimationSettings()
  g3 <- mlx.getLogLikelihoodEstimationSettings()
  g4 <- mlx.getConditionalDistributionSamplingSettings()
  g5 <- mlx.getStandardErrorEstimationSettings()
  g6 <- mlx.getConditionalModeEstimationSettings()
  
  for (j in 1:np) {
    pj <- spline(x,p[j,],n,method="hyman")$y[level]
    if (ind.int[j])
      pj <- round(pj)
    if (set.list[j] %in% names(g1))
      g1[set.list[j]] <- pj
    else if (set.list[j] %in% names(g2))
      g2[set.list[j]] <- pj
    else if (set.list[j] %in% names(g3))
      g3[set.list[j]] <- pj
    else if (set.list[j] %in% names(g4))
      g4[set.list[j]] <- pj
    else if (set.list[j] %in% names(g5))
      g5[set.list[j]] <- pj
    else if (set.list[j] %in% names(g6))
      g6[set.list[j]] <- pj
  }
  mlx.setGeneralSettings(g1)
  g2$simulatedannealingiterations <- NULL
  mlx.setPopulationParameterEstimationSettings(g2)
  mlx.setLogLikelihoodEstimationSettings(g3)
  mlx.setConditionalDistributionSamplingSettings(g4)
  mlx.setStandardErrorEstimationSettings(g5)
  mlx.setConditionalModeEstimationSettings(g6)
  
  if (!is.null(new.project)) 
    mlx.saveProject(projectFile = new.project)
  else
    eval(parse(text=paste0("lixoftConnectors::saveProject()")))
  
  
}
