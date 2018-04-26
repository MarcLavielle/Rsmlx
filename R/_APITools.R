.getOS = function(){
  if (.Platform$OS.type == "windows") { 
    return("Windows")
  } else if (Sys.info()["sysname"] == "Darwin") {
    return("Apple") 
  } else if (.Platform$OS.type == "unix") { 
    return("Unix")
  } else {
    stop("Unknown OS")
  }
}

.isInteger = function(value){
  if (is.numeric(value) == FALSE){
    return(invisible(FALSE))
  }
  else if ( (value > .Machine$integer.max) || (value <= -.Machine$integer.max) ){
    .warning(paste0("Integer values should be between ",-(.Machine$integer.max+1)," and ",.Machine$integer.max,"."))
    return(invisible(FALSE))                               
  }
  else if ( value - as.integer(value) == 0 ){
    return(invisible(TRUE))
  }
  else {
    return(invisible(FALSE))
  }
}

.evalNumerics = function(x){
  tryCatch(expr = {
    if (!is.logical(x) && !is.na(as.numeric(x)))
      x <- as.numeric(x)
  }, warning = function(e){})
  return(invisible(x))
}

.makeRequest = function(applicationName, functionName, arguments, requestType, wait = TRUE){
  loadedDLLs <- getLoadedDLLs();
  tryCatch(expr = {
    if (requestType != "synchronous"){ # if the request is asynchronous, the user value for "wait" will be ignored and set to TRUE by default
      wait = TRUE
    }
    
    outputStruct <- invisible ( .C( "applyJson",
                                    .encodeString(applicationName),
                                    .encodeString(functionName),
                                    RJSONIO::toJSON( .encodeString(arguments), digits = 15, .level = -Inf),
                                    .encodeString(requestType),
                                    wait,
                                    output_R = "",
                                    PACKAGE = get("MLXCONNECTORS_LIB_NAME", envir = MlxEnvironment)
    ) )
    return(invisible(outputStruct$output_R))
  },
  
  error = function(err){
    if ( !exists("MLXCONNECTORS_LIB_PATH",envir = MlxEnvironment) ){
      .error("Unknown OS -> Impossible to identify MlxConnectors library.")
      return(invisible())      
    } else if ( is.null(loadedDLLs[[get("MLXCONNECTORS_LIB_PATH",envir = MlxEnvironment)]]) ){
      .warning("The API has not been initialized yet. Please call initializeMlxConnectors() to make MlxConnectors package functions available.")
      return(invisible())     
    }
  },
  
  warning = function(war){
    if ( !exists("MLXCONNECTORS_LIB_PATH", envir = MlxEnvironment) || get("MLXCONNECTORS_LIB_PATH", envir = MlxEnvironment) == "" ){
      .warning("Unknown OS -> Impossible to identify MlxConnectors library.")
      return(invisible())      
    } else if ( is.null(loadedDLLs[[get("MLXCONNECTORS_LIB_PATH", envir = MlxEnvironment)]]) ){
      .warning("The API has not been initialized yet. Please call initializeMlxConnectors() to make MlxConnectors package functions available.")
      return(invisible())     
    }
  }
  ) 
}

.processRequest = function(applicationName, functionName, arguments, requestType, wait = TRUE,
                           type = "", ...){
  jsonOutput = .makeRequest(applicationName, functionName, arguments, requestType, wait)
  decodedOutput = .decodeFromJSON(jsonOutput, type)
  return(invisible(decodedOutput))
}

.checkOutputFormat = function(jsonObj){
  outputKeys = names(jsonObj)
  if ( "output" %in% outputKeys && "errors" %in% outputKeys && "warnings" %in% outputKeys && "information" %in% outputKeys ){
    return(invisible(TRUE))
  } else {
    .error("Bad json format encountered.")
    return(invisible(FALSE))
  }
}

.encodeNumericData = function(data){
  for (i in 1:length(data)){
    if (is.nan(data[[i]])){
      data[[i]] = "NaN"
    }
    else if (data[[i]] == Inf){
      data[[i]] = "Infinity"
    }
    else if (data[[i]] == -Inf){
      data[[i]] = "-Infinity"
    }
  }
  return(invisible(data))
}

.encodeString = function(struct){
  if ( is.list(struct) == TRUE && length(struct)>0 ){
    for (i in 1:length(struct)){
      if (is.character(struct[[i]])){
        struct[[i]] = enc2utf8(struct[[i]])
      }
      else if ( is.list(struct[[i]]) || is.numeric(struct[[i]]) ){
        struct[[i]] = .encodeString(struct[[i]])
      }
    }
  }
  else if (is.character(struct) == TRUE){
    struct = enc2utf8(struct)
  }
  else if (is.numeric(struct) == TRUE){
    struct = .encodeNumericData(struct)
  }
  return(invisible(struct))
}

.decodeFromJSON = function(jsonObj, type){
  tryCatch({
    decodedOutput = (fromJSON(jsonObj,digits = 15))
    
    if (.checkOutputFormat(decodedOutput)){
      if (length(decodedOutput$information) != 0){
        for (i in 1:length(decodedOutput$information))
          if (length(decodedOutput$information[[i]]) > 0)
            .info(decodedOutput$information[[i]],FALSE)
      }
      
      hasWarnings = length(decodedOutput$warnings) > 0
      if (hasWarnings){
        for (i in 1:length(decodedOutput$warnings))
          if (length(decodedOutput$warnings[[i]]) > 0)
            .warning(decodedOutput$warnings[[i]],FALSE)
      }
      
      hasErrors = length(decodedOutput$errors) > 0
      if (hasErrors){
        for (i in 1:length(decodedOutput$errors))
          if (length(decodedOutput$errors[[i]]) > 0)
            .error(decodedOutput$errors[[i]],FALSE)
      }
      
      if (type == "STATUS")
        return(invisible(!hasErrors))
      
      else {
        if (hasErrors || (hasWarnings && is.null(decodedOutput$output)))
          return(invisible())
        else
          return(invisible(decodedOutput$output)) 
      }
    }
    else {
      return(invisible())
    }
  },
  error = function(err){
    return(invisible())
  })
  
}

.error = function(str, addSource = FALSE, source = NULL){
  if (addSource == FALSE){
    source = ""
  } else if (is.null(source)){
    source = paste0("\n(from : ",deparse(sys.call(-1)),")")
  } else {
    source = parse0("\n(from : ",source,")")
  }
  message(paste0("[ERROR] ",str,source))
}

.warning = function(str, addSource = FALSE, source = NULL){
  if (addSource == FALSE){
    source = ""
  } else if (is.null(source)){
    source = paste0("\n(from : ",deparse(sys.call(-1)),")")
  } else {
    source = parse0("\n(from : ",source,")")
  }
  message(paste0("[WARNING] ",str,source))
}

.info = function(str, addSource = FALSE, source = NULL){
  if (addSource == FALSE){
    source = ""
  } else if (is.null(source)){
    source = paste0("\n(from : ",deparse(sys.call(-1)),")")
  } else {
    source = parse0("\n(from : ",source,")")
  }
  writeLines(paste0("[INFO] ",str,source)) 
}

.exists = function(path){
  i = nchar(path)
  while ( (substr(path,i,i) == "/" || substr(path,i,i) == "\\") && i>0 )
    i = i-1
  
  path = substr(path,1,i)
  return(invisible(file.exists(path)))
}