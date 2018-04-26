#> settings("warfpkpd1.mlxtran",1) #very fast but not  precise
#> settings("warfpkpd1.mlxtran",5) #indermediate
#> settings("warfpkpd1.mlxtran",9) #slow but very accurate

setSettings  <- function(project=NULL, final.project=NULL, level=5, new=FALSE) {
  n <- 9
  
  level <- round(level)
  if (level<1 | level>n)
    stop(paste0("level should be an integer beween 1 and ",n), call. = FALSE)
  
  if (!is.null(project)) 
    loadProject(project)
  
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
}
