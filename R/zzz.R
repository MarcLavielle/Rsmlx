.onLoad <- function(libname, pkgname){
  # if (length(find.package("lixoftConnectors", quiet = TRUE))) {
  #   op <- options()
  #   op$lixoft_notificationOptions$info <- 1
  #   options(op)
  #   lixoftConnectors::initializeLixoftConnectors(software = "monolix", force = TRUE)
  # }
}

.onAttach <- function(libname,pkgname){
  startMessage()
}

startMessage <- function() {
  # if (!length(find.package("lixoftConnectors", quiet = TRUE))) {
  #   packageStartupMessage(paste0(
  #     "You need to install the lixoftConnectors package in order to use Rsmlx \n
  #     See http://monolix.lixoft.com/monolix-api/lixoftconnectors_installation/ for more information"
  #   ))
  # }
  packageStartupMessage("R speaks Monolix!")
}


