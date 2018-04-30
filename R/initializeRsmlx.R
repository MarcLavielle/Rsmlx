.onLoad <- function(libname,pkgname){
  #packageStartupMessage("This is Rsmlx package, enjoy it!")
  #initializeMlxConnectors(software = "monolix")
  # assign("MLXCONNECTORS_PKG_NAME", pkgname, envir = MlxEnvironment)
  # assign("MLXCONNECTORS_PKG_DIR", libname, envir = MlxEnvironment)
}
# 
.onAttach <- function(libname,pkgname){
  packageStartupMessage("You need to initialize the R connectors for Monolix with the command:")
  packageStartupMessage('initializeMlxConnectors(software = "monolix")')
  # #  packageStartupMessage("This Rsmlx package, enjoy it!")
# #  assign("MLXCONNECTORS_PKG_NAME", pkgname, envir = MlxEnvironment)
# #  assign("MLXCONNECTORS_PKG_DIR", libname, envir = MlxEnvironment)
}

