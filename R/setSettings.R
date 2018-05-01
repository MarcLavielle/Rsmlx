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
#' setSettings(project="theophylline.mlxtran", level=1)
#' setSettings(project="theophylline.mlxtran", level=9)
#' setSettings(project="theophylline.mlxtran", new.project="theophyllineNew.mlxtran", level=7)
#' 
#' setSettings(level=7);  suppose that a project has been already loaded using the MlxConnectors
#' }
#' @export

setSettings  <- function(project=NULL, new.project=NULL, level=5) {
  
  if (!is.null(project)) {
    if(!file.exists(project)){
      message(paste0("ERROR: project '", project, "' does not exists"))
      return(invisible(FALSE))}
    
#    initializeMlxConnectors(software = "monolix")
    loadProject(project) 
    
    if (is.null(new.project)) 
      new.project <- project
  }
  
  n <- 9
  
  level <- round(level)
  if (level<1 | level>n)
    stop(paste0("level should be an integer beween 1 and ",n), call. = FALSE)
  
  
  setGeneralSettings(autochains=TRUE)
  setPopulationParameterEstimationSettings(exploratoryautostop=TRUE,
                                           smoothingautostop=TRUE)
  
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
  
  g1 <- getGeneralSettings()
  g2 <- getPopulationParameterEstimationSettings()
  g3 <- getLogLikelihoodEstimationSettings()
  g4 <- getConditionalDistributionSamplingSettings()
  g5 <- getStandardErrorEstimationSettings()
  g6 <- getConditionalModeEstimationSettings()
  
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
  setGeneralSettings(g1)
  setPopulationParameterEstimationSettings(g2)
  setLogLikelihoodEstimationSettings(g3)
  setConditionalDistributionSamplingSettings(g4)
  setStandardErrorEstimationSettings(g5)
  setConditionalModeEstimationSettings(g6)
  
  if (!is.null(new.project)) 
    saveProject(projectFile = new.project)
  else
    saveProject()
  
}