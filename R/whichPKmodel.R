#' Find a Monolix PK model
#'
#' Return the path of the Monolix PK model defined by a list of parameter names
#' See https://monolix.lixoft.com/rsmlx/whichPKmodel/ for more details.
#' @param parameter a vector of PK parameter names
#' @param mlxPath path to Monolix install
#' @param pkPath path to the Monolix PK library
#' @param lib boolean to define if the absolute path is returned
#' @examples
#' \dontrun{
#' whichPKmodel(parameter=c("Tlag", "Tk0", "V", "Cl"))
#' }
#' @export

#session <- "C:/ProgramData/Lixoft/MonolixSuite2018R2"


whichPKmodel <- function(parameter, mlxPath = NULL, pkPath = NULL, lib = FALSE) {
  
  initRsmlx()
  
  parameter[parameter=="k"] <- "ke"
  parameter[parameter=="Q"] <- "Q2"
  parameter[parameter=="V"] <- "V1"
  parameter <- sort(parameter)

  d <- gsub("lib:", "", mlx.getLibraryModelName("pk"))
  lib <- TRUE

  id0 <- c(grep("alpha",d), grep("bolus_",d))
  d <- d0 <- d[-id0]
  d <- gsub("kk","kek",d)
  d <- gsub("k.txt","ke.txt",d)
  d <- gsub("QV","Q2V",d)
  d <- gsub("Vk","V1k",d)
  d <- gsub("VCl","V1Cl",d)
  d <- gsub("VV","V1V",d)
  
  plist <- c("Tlag", "Mtt", "Ktr", "ka", "Tk0", "V1", "V2", "V3", "k12", "k21", "k13", "k31", "Q2", "Q3", "Vm", "Km", "Cl", "ke")
  for (k in (1:length(d))) {
    lk <- plist[unlist(lapply(plist, function(x) grepl(x, d[k])))]
    if (identical(sort(lk), parameter)) 
      if (lib)
        return(paste0("lib:",d0[k]))
      else
        return(file.path(pkPath,d0[k]))
  }
  
  stop("The model was not found in the PK library", call.=FALSE)
}



mlx.path <- function(){
  myOS <- Sys.info()['sysname']; 
  lixoft.path <- {
    if (myOS == "Windows"){ 
      file.path(Sys.getenv("USERPROFILE"),"lixoft")
    } else {
      file.path(Sys.getenv("HOME"),"lixoft")
    }
  } 
  
  lixoft.ini  <- file.path(lixoft.path,"lixoft.ini")
  if (!file.exists(lixoft.ini))
    stop(paste0("\nThe file ",lixoft.ini," does not exists."), call.=FALSE)
  
  get_lixoft_path <- function(name, lines){
    rx <- sprintf( "%s=", name)
    line <- grep( rx, lines, fixed=TRUE, value=TRUE)
    if( length(line) ){
      normalizePath( gsub( rx, "", line[1L] ) )
    }
  }
  
  lines <- readLines(lixoft.ini)
  monolix.path <- get_lixoft_path("monolixSuite", lines)
  monolix.path <- gsub("\\\\","/",monolix.path)
  return(monolix.path)
}
