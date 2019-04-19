.onLoad <- function(libname,pkgname){
  # packinfo <- utils::installed.packages()
  # if (!is.element("lixoftConnectors", packinfo[,1])){
  #   .error("You need to install the lixoftConnectors package in order to use Rsmlx")
  #   return(invisible(FALSE))
  # }
#  library(lixoftConnectors)
#  initializeLixoftConnectors(software="monolix", force=TRUE)
}
# 
.onAttach <- function(libname,pkgname){
  startMessage()
}

startMessage <- function() {
  packageStartupMessage("R speaks Monolix!")
}


