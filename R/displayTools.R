###: Display Mlx API Structures
###: 
###: [\bold{Tools}][\emph{Display}] \cr 
###: Display the structures retrieved by the MlxConnectors library in a user-friendly way.
###: @param structure [\emph{miscellanous}] The data structure to be displayed.
###: @examples 
###: \dontrun{
###: settings = getGeneralSettings()
###: mlxDisplay(settings)
###: }
###: @export
mlxDisplay = function(structure){
  if (is.null(structure) == TRUE){
    return(invisible())
  }
  
  if (is.list(structure) == TRUE){
    if (length(structure) == 0){ return(invisible()) }
    
    else if (length(structure) == 1){ print(structure[[1]]) }
    
    else {
      table = c()
      for (i in 1:length(structure)){
        if (is.list(structure[[i]]) == FALSE){
          if (is.array(structure[[i]]) == TRUE){
            table = rbind(table,"[array]")
          }
          else if (length(structure[[i]]) > 1){
            fieldNames = names(structure[[i]])
            if (length(fieldNames) == 0){
              table = rbind(table,paste0("{ ", toString(structure[[i]])," }"))
            }
            else {
              tmp = paste0("{ \"",fieldNames[[1]],"\":",toString(structure[[i]][[1]]));
              for (j in 2:length(structure[[i]])){
                tmp = paste0(tmp,", \"",fieldNames[[j]],"\":",toString(structure[[i]][[j]]))
              }
              tmp = paste0(tmp," }")
              table = rbind(table,tmp)
            }
          }
          else {
            fieldName = names(structure[[i]])
            if (length(fieldName) == 0){
              table = rbind(table,toString(structure[[i]]))
            } else {
              table = rbind(table,paste0(fieldName,":",toString(structure[[i]])))
            }
          }
        } else {
          table = rbind(table,paste0("structure< ",toString(names(structure[[i]]), width = 100)," >"))
        }
      }
      colnames(table) <- c("")
      rownames(table) <- names(structure)
      print(table, quote = FALSE)    
    }
  }
  
  else if (is.array(structure) == TRUE){
    print(structure)
  }
  
  else{
    if (is.null(names(structure)) == FALSE){
      print(structure, quote = TRUE)
    }
    else {
      dim(structure) <- c(length(structure),1)
      rownames(structure) <- rep("",length(structure))
      colnames(structure) <- ""
      print(structure) 
    }
  }
}

.buildArray = function(structure){
  table = c()
  for (i in 1:length(structure)){
    table = cbind(table, structure[[i]])
  }
  
  tableSize = nrow(table)
  colnames(table) <- names(structure)
  rowNames = names(structure[[1]])
  if (length(rowNames) == 0){
    if (tableSize > 1){
      rowNames = c(1:tableSize)
    } else {
      rowNames = ""
    }
  } 
  rownames(table) <- rowNames
  
  return(invisible(table))
}

.buildMatrix = function(structure){
  fieldNames = names(structure)
  
  if ( !("data" %in% fieldNames && "rownumber" %in% fieldNames) ){ 
    return(invisible())
  }
  
  matrix = c()
  columnnumber = length(structure$data)/structure$rownumber
  for ( iColumn in 0:(columnnumber-1) ){ # ! data vector is column major ordered !
    matrix = cbind( matrix, structure$data[(1 + iColumn*structure$rownumber) : ((iColumn+1)*structure$rownumber)] )
  }
  
  if ("rownames" %in% fieldNames && length(structure$columnnames) == structure$rownumber){
    rownames(matrix) <- structure$columnnames
  } else {
    rownames(matrix) <- rep("",dim(matrix)[[1]])
  }
  
  if ("columnnames" %in% fieldNames && length(structure$columnnames) == columnnumber){
    colnames(matrix) <- structure$columnnames
  } else {
    colnames(matrix) <- rep("",dim(matrix)[[2]])
  }
  return(invisible(matrix))
}