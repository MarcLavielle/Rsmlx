
#' Initialize Rsmlx library
#' 
#' Initialize Rsmlx library
#' @return A boolean equaling TRUE if the initialization has been successful and FALSE if not.
#' @examples
#' \dontrun{
#' initRsmlx()
#' }
#' @export
initRsmlx <- function(){
  packinfo <- utils::installed.packages()
  if (!is.element("lixoftConnectors", packinfo[,1]))
    stop("You need to install the lixoftConnectors package in order to use Rsmlx", call. = FALSE)
    
  
  lixoftConnectorsState <- getLixoftConnectorsState(quietly = TRUE)
  
  if (!is.null(lixoftConnectorsState)){

    if (lixoftConnectorsState$software == "monolix") {
      status=TRUE
    } else {
      status = initializeLixoftConnectors(software = "monolix")
    }
    
  } else {
    status = initializeLixoftConnectors(software = "monolix")
  }
  return(invisible(status))
  
}

