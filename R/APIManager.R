#' Initialize MlxConnectors API
#' 
#' Initialize MlxConnectors API for a given software
#' @param software (\emph{character}) Name of the software to be loaded : "monolix"\cr
#' @param mlxDirectory (\emph{character}) [optional] Path to installation directory of the Lixoft suite.
#' If no path is given, the one written in the lixoft.ini file is used.
#' @return A boolean equaling TRUE if the initialization has been successful and FALSE if not.
#' @examples
#' \dontrun{
#' initializeMlxConnectors(software = "monolix", mlxDirectory = "/path/to/mlxRuntime/")
#' }
#' @export
initializeMlxConnectors <- function (software, mlxDirectory = "") {
  if (is.character(software) == FALSE){
    .error("Unexpected type encountered for the \"software\" field. Please give a string corresponding to the name of the Lixoft software to be used.")
    return(invisible(FALSE))
  }  
  if (is.character(mlxDirectory) == FALSE){
    .error("Unexpected type encountered for the \"mlxDirectory\" field. Please give a string corresponding to the path to the Lixoft Suite installation directory.")
    return(invisible(FALSE))
  }
  
  status = .setMlxDirectory(software,mlxDirectory)
  if (status == FALSE){return(invisible(FALSE))}
  
  mlxLibraryInfo = .checkMlxLibraryInformations()
  if (mlxLibraryInfo$status == FALSE){return(invisible(FALSE))}
  
  
  # check if the mlxConnectors library is already loaded :
  loadedDLLs <- getLoadedDLLs();
  currentLoadedLibPath = loadedDLLs[[get("MLXCONNECTORS_LIB_NAME", envir = MlxEnvironment)]];
  if (!is.null(currentLoadedLibPath)){
    if (currentLoadedLibPath[["path"]] == mlxLibraryInfo$path) {
      .info(paste0("The library ", get("MLXCONNECTORS_LIB_NAME", envir = MlxEnvironment), " (\"", get("MLXCONNECTORS_LIB_PATH", envir = MlxEnvironment), "\") is already loaded -> nothing to be done."));
      return(invisible(TRUE));
    }
    else {
      .info(paste0("A library ", get("MLXCONNECTORS_LIB_NAME", envir = MlxEnvironment), " from an other software version (\"", get("MLXCONNECTORS_LIB_PATH", envir = MlxEnvironment), "\") is already loaded -> it will be unloaded before loading the new one."))
      .unloadMlxConnectorsLibrary()
    }
  }
  
  # set environment :
  Sys.setenv( LIXOFT_HOME = get("MLX_DIRECTORY", envir = MlxEnvironment))
  
  OS = .getOS()
  if (OS == "Windows") {
    runtimeDir = gsub("/", "\\\\", get("MLX_DIRECTORY", envir = MlxEnvironment))
    Sys.setenv(PATH = paste0( runtimeDir,"\\tools\\MinGW\\bin;",
                              runtimeDir,"\\...\\tools\\MinGW\\bin;",
                              get("SYSTEM_PATH", envir = MlxEnvironment) ) )
  }
  
  # load library :
  assign("MLXCONNECTORS_LIB_NAME", mlxLibraryInfo$name, envir = MlxEnvironment)
  assign("MLXCONNECTORS_LIB_PATH", mlxLibraryInfo$path, envir = MlxEnvironment)
  
  PATH <- Sys.getenv("PATH")
  
  if (OS == "Unix") {
    Sys.setenv( LD_LIBRARY_PATH = paste0( dirname(get("MLXCONNECTORS_LIB_PATH", envir = MlxEnvironment)),":",PATH), collapse = "" )
  } else if (OS == "Apple") {
    # path is set during installation
  } else if (OS == "Windows") {
    Sys.setenv(PATH = paste0( gsub("/", "\\\\", dirname(get("MLXCONNECTORS_LIB_PATH", envir = MlxEnvironment))),";",PATH) )
  }
  
  mlxDll = dyn.load(get("MLXCONNECTORS_LIB_PATH", envir = MlxEnvironment))
  .dynLibs(c(.dynLibs(), list(mlxDll)))
  
  Sys.setenv(PATH = PATH)
  
  
  # initialize session :
  arguments = list(software = software, libpath = get("MLXCONNECTORS_LIB_PATH", envir = MlxEnvironment))
  output = .processRequest("session", "initialize", arguments, "synchronous", type = "STATUS")
  
  if (!output) # unload library if initialization failed
    .unloadMlxConnectorsLibrary()
  
  return(invisible(output))
}

.setMlxDirectory <- function(software, mlxDirectory = ""){
  if (mlxDirectory == ""){
    .info(paste0("The path to ",software," installation directory has not been given."))
    OS = .getOS()
    if (OS == "Unix") {
      lixoftIniPath <- paste0(Sys.getenv("HOME"),"/lixoft/lixoft.ini")
    } else if (OS == "Apple") {
      lixoftIniPath <- paste0(Sys.getenv("HOME"),"/lixoft/lixoft.ini")
    } else if (OS == "Windows") {
      lixoftIniPath <- paste0(Sys.getenv("HOMEDRIVE"),Sys.getenv("HOMEPATH"),"\\lixoft\\lixoft.ini")
    } else {
      .error("Unknown OS.")
      return(invisible(FALSE))
    }
    
    if (.exists(lixoftIniPath) == FALSE){
      .error(paste0("Impossible to find an initialization file for Lixoft Suite.\n-> The Lixoft API will not been initialized."))
      return(invisible(FALSE))
    }
    lines = readLines(lixoftIniPath)
    for (i in 1:length(lines)){
      res = regexpr("monolixSuite=",lines[[i]])
      if (res != -1){
        mlxDirectory = substring(lines[[i]],res + attr(res,"match.length"))
      }
    }
    
    if (mlxDirectory == ""){
      .error(paste0("No setup informations for the \"",software,"\" software in the initialization file \"",lixoftIniPath,"\".\n-> The Lixoft API will not been initialized."))
      return(invisible(FALSE))
    } else {
      .info(paste0("The directory specified in the initialization file of the Lixoft Suite (located at \"",lixoftIniPath,"\") will be used by default.\n-> \"",mlxDirectory,"\""))
    }
  }
  
  if (.exists(mlxDirectory) == FALSE) {
    .error(paste0("The directory \"",mlxDirectory,"\" does not exist."))
    return(invisible(FALSE))
  } else {
    assign("MLX_DIRECTORY", mlxDirectory, envir = MlxEnvironment)
    return(invisible(TRUE))
  }
}

.checkMlxLibraryInformations <- function(){
  OS = .getOS()
  if (OS == "Unix") {
    mlxConnectorsLibName <- "libmlxConnectors"
    libSuffix <- ".so"
  }
  else if (OS == "Apple") {
    mlxConnectorsLibName <- "libmlxConnectors"
    libSuffix <- ".so" # ".dylib" won't be loaded by .C()
  }
  else if (OS == "Windows") {
    mlxConnectorsLibName <- "mlxConnectors"
    libSuffix <- ".dll"
  }
  else {
    .error("Unknown OS.")  
    return(invisible(FALSE))
  }
  
  mlxConnectorsLibPath <- paste0(get("MLX_DIRECTORY", envir = MlxEnvironment), "/lib/", mlxConnectorsLibName, libSuffix, collapse = "")
  if (.exists(mlxConnectorsLibPath) == FALSE){
    .error(paste0("Impossible to find mlxConnectors library at \"",mlxConnectorsLibPath,"\"."))
    return(invisible(list(status = FALSE)))
  }
  return(invisible( list(status = TRUE, name = mlxConnectorsLibName, path = mlxConnectorsLibPath) ))
}

.unloadMlxConnectorsLibrary <- function(){
  output = .processRequest("session", "unload", list(void = ""), "synchronous", type = "STATUS")
  
  if (output){
    dyn.unload(get("MLXCONNECTORS_LIB_PATH", envir = MlxEnvironment))
    return(TRUE)
  }
  else{
    return(FALSE)
  }
}

.onUnload <- function(libpath){
  #  unload mlxConnectors library :
  loadedDLLs <- getLoadedDLLs();
  if (exists("MLXCONNECTORS_LIB_NAME", envir = MlxEnvironment) && !is.null(loadedDLLs[[get("MLXCONNECTORS_LIB_NAME", envir = MlxEnvironment)]])) {
    .unloadMlxConnectorsLibrary()
  }
  
  # reset environment :
  Sys.setenv( LIXOFT_HOME = "")
  
  OS = .getOS()
  if (OS == "Unix") {
    Sys.setenv( LD_LIBRARY_PATH = get("SYSTEM_PATH", envir = MlxEnvironment) )
  } else if (OS == "Apple") {
    # path is set during installation
  } else if (OS == "Windows") {
    Sys.setenv( PATH = get("SYSTEM_PATH", envir = MlxEnvironment) )
  }
}

.onLoad <- function(libname,pkgname){
  # assign("MLXCONNECTORS_PKG_NAME", pkgname, envir = MlxEnvironment)
  # assign("MLXCONNECTORS_PKG_DIR", libname, envir = MlxEnvironment)
}

.onAttach <- function(libname,pkgname){}

.onDetach <- function(libpath){}

.Last.lib <- function(libpath){}

#' Get information about MlxEnvironment object
#' @export
getMlxEnvInfo <- function(){
  r <- list(ls.str(envir = MlxEnvironment), parent.env(MlxEnvironment))
  return(r)
}
