# ============================== PROJECT MANAGEMENT ============================== #
# Load-Save ---------------------------------------------------------------------- #
###: Load project from file
###:
###: Load a project by parsing the mlxtran-formated file whose path has been given as an input.
###: WARNING: R is sensitive between '\' and '/', only '/' can be used
###: @param projectFile (\emph{character}) Path to the project file. Can be absolute or relative to the current working directory.
###: @examples
###: \dontrun{
###: loadProject("/path/to/project/file.mlxtran") for Linux platform
###: loadProject("C:/Users/path/to/project/file.mlxtran") for Windows platform
###: }
###: @seealso \code{\link{saveProject}}
###: @export
loadProject = function(projectFile){
  if (is.character(projectFile) == FALSE){
    .error("Unexpected type encountered. Please give a string corresponding to the path to the mlxtran-formated file to be loaded.")
    return(invisible(FALSE))
  }
  
  arguments = list( normalizePath(projectFile,mustWork = FALSE) )
  output = .processRequest("monolix", "loadproject", arguments, "synchronous", type = "STATUS")
  return(invisible(output))
}

###: Save current project
###:
###: Save the current project as an Mlxtran-formated file.
###: @param projectFile [optional](\emph{character}) Path where to save a copy of the current mlxtran model. Can be absolute or relative to the current working directory.
###: If no path is given, the file used to build the current configuration is updated.
###: @examples
###: \dontrun{
###: saveProject("/path/to/project/file.mlxtran") # save a copy of the model
###: saveProject() # update current model
###: }
###: @seealso \code{\link{newProject}} \code{\link{loadProject}}
###: @export
saveProject = function(projectFile = ""){
  if (is.character(projectFile) == FALSE){
    .error("Unexpected type encountered. Please give a string corresponding to the path to the mlxtran-formated file to be used to save the current project.")
    return(invisible(FALSE))
  }
  
  if (projectFile == "")
    output = .processRequest("monolix", "saveproject", "", "synchronous", type = "STATUS")
  else
    output = .processRequest("monolix", "saveasproject", projectFile, "synchronous", type = "STATUS")
  
  return(invisible(output))
}

###: Create new project
###:
###: Create a new empty project providing model and data specification. The data specification is:
###:    \itemize{
###:    \item dataFile (\emph{string}): path to the data file
###:    \item headerTypes (\emph{array<character>}): vector of headers
###:	  \item observationTypes (\emph{list}): a list giving the type of each observation present in the data file (if there is only one y-type, the corresponding observation name can be omitted)
###:	  \item nbSSDoses (\emph{int}): number of steady-state doses (if there is a SS column)
###:}
###: @param modelFile (\emph{character}) Path to the model file. Can be absolute or relative to the current working directory.
###: @param data (\emph{list}) Structure describing the data.
###: @examples
###: \dontrun{
###: newProject(data = list(dataFile = "/path/to/data/file.txt", 
###:                        headerTypes = c("IGNORE","OBSERVATION"), 
###:                        observationTypes = "continuous"),
###:            modelFile = "/path/to/model/file.txt")
###: }
###: @seealso \code{\link{newProject}} \code{\link{saveProject}}
###: @export
newProject = function(modelFile, data){
  if (is.character(modelFile) == FALSE){
    .error("Unexpected type encountered for the field \"model\". Please give a string corresponding to the path to a model file.")
    return(invisible(FALSE))
  }
  
  if ( !is.list(data) || is.null(data$dataFile) || is.null(data$headerTypes) || is.null(data$observationTypes) ){
    .error("Unexpected type encountered for the field \"data\". Please give a list containing at least the data file, the header types and the observation types.")
    return(invisible(FALSE))
  }
  
  if (is.character(data$dataFile) == FALSE){
    .error("Unexpected type encountered for the field \"dataFile\" from \"data\". Please give a string corresponding to the path to a data file.")
    return(invisible(FALSE))
  }
  
  modelFile <- normalizePath(modelFile, mustWork = FALSE)
  data$dataFile <- normalizePath(data$dataFile, mustWork = FALSE)
  
  arguments = list(modelFile, data$dataFile)
  output = .processRequest("monolix", "createproject", arguments, "synchronous", type = "STATUS")
  if (output == FALSE) return (invisible(FALSE))
  
  output = setData(data)
  return(invisible(output))
}
# -------------------------------------------------------------------------------- #

# Model -------------------------------------------------------------------------- #
###: Set structural model file
###:
###: Set the structural model.
###: @param modelFile (\emph{character}) Path to the model file. Can be absolute or relative to the current working directory.
###: @examples
###: \dontrun{
###: setStructuralModel("/path/to/model/file.txt")
###: }
###: @seealso \code{\link{getStructuralModel}}
###: @export
setStructuralModel = function(modelFile){
  if (is.character(modelFile) == FALSE){
    .error("Unexpected type encountered. Please give a string corresponding to the path to a model file.")
    return(invisible(FALSE))
  }
  
  modelFile <- normalizePath(modelFile, mustWork = FALSE)
  output = .processRequest("monolix", "setmodel", modelFile, "synchronous", type = "STATUS")
  return(invisible(output))
}

###: Get structural model file
###:
###: Get the model file for the structural model used in the current project.
###: @return A string corresponding to the path to the structural model file.
###: @examples
###: \dontrun{
###: getStructuralModel() => "/path/to/model/inclusion/modelFile.txt"
###: }
###: @seealso \code{\link{setStructuralModel}}
###: @export
getStructuralModel = function(){
  arguments = list()
  output = .processRequest("monolix", "getmodel", arguments, "asynchronous")
  return(output)
}
# -------------------------------------------------------------------------------- #

# Data --------------------------------------------------------------------------- #
###: Set project data
###:
###: Set project data giving a data file and specifying headers and observations types.
###: @param dataFile (\emph{character}): Path to the data file. Can be absolute or relative to the current working directory.
###: @param headerTypes (\emph{array<character>}): A collection of header types. 
###: The possible header types are: "id", "time", "observation", "amount", "contcov", "catcov", "occ", "evid", "mdv", "obsid", "cens", "limit", "regressor","admid", "rate", "tinf", "ss", "ii", "addl", "date"
###: Notice that these are not the types displayed in the interface, these one are shortcuts.
###: @param observationTypes (\emph{list}): A list giving the type of each observation present in the data file. If there is only one y-type, the corresponding observation name can be omitted.
###: The possible observation types are "continuous", "discrete", and "event"
###: @param nbSSDoses [optional](\emph{int}): Number of doses (if there is a SS column).
###: @examples
###: \dontrun{
###: setData(dataFile = "/path/to/data/file.txt", headerTypes = c("IGNORE","OBSERVATION"), observationTypes = "continuous")
###: setData(dataFile = "/path/to/data/file.txt", headerTypes = c("IGNORE","OBSERVATION","YTYPE"), observationTypes = list(Concentration = "continuous", Level = "discrete"))
###: }
###: @seealso \code{\link{getData}}
###: @export
setData = function(dataFile, headerTypes, observationTypes, nbSSDoses = NULL){
  
  if (nargs() == 1 && is.list(dataFile)){
    arguments = dataFile
    if (!is.null(arguments$observationNames) && !is.null(arguments$observationTypes) && is.null(names(arguments$observationTypes)) &&
        length(arguments$observationNames) == length(arguments$observationTypes)){
      names(arguments$observationTypes) <- arguments$observationNames
    }
  }
  else {
    arguments = list(dataFile = dataFile, headerTypes = headerTypes, observationTypes = observationTypes)
    if (!is.null(nbSSDoses))
      arguments$nbSSDoses = nbSSDoses
  }
  
  if (is.character(arguments$dataFile) == FALSE){
    .error("Unexpected type encountered for field \"dataFile\". Please give a string.")
    return(invisible(FALSE))
  }
  
  if (!is.vector(arguments$headerTypes) && !is.character(headerTypes)){
    .error("Unexpected type encountered for field \"headerTypes\". Please give a collection of strings.")
    return(invisible(FALSE))
  }
  
  if (is.vector(arguments$observationTypes) || is.list(arguments$observationTypes)){
    yTypes = names(arguments$observationTypes)
    
    if (is.null(yTypes) && length(yTypes) == 1){
      if (!is.character(arguments$observationTypes[[1]])){
        .error(paste0("Unexpected type encountered at position ",i,". Please give a string corresponding to the wanted observation type name of the preceding y-type."))
        return(invisible(FALSE))
      }
      arguments$observationTypes = arguments$observationTypes[[1]]
    }
    else if (length(arguments$observationTypes) > 1){
      for (i in 1:length(arguments$observationTypes)){
        if (yTypes[i] == ""){
          .error(paste0("No y-type name found at position ",i,"."))
          return(invisible(FALSE))
        }
        else if (is.character(arguments$observationTypes[[yTypes[i]]]) == FALSE){
          .error(paste0("Unexpected type encountered at position ",i,". Please give a string corresponding to the wanted observation type name of the preceding y-type."))
          return(invisible(FALSE))
        }
      }
    }
  }
  else if (!is.character(arguments$observationTypes) || length(arguments$observationTypes) != 1){
    .error("Unexpected type encountered for field \"observationTypes\". Please give a list of string indexed by observation model names.")
    return(invisible(FALSE))
  }
  
  if (!is.null(arguments$nbSSDoses)){
    if (.isInteger(arguments$nbSSDoses) == FALSE){
      .error("Unexpected type encountered for field \"nbSSDoses\". Please give an integer.")
      return(invisible(FALSE))
    }
  }
  
  arguments$dataFile <- normalizePath(arguments$dataFile, mustWork = FALSE)
  output = .processRequest("monolix", "setdata", arguments, "synchronous", type = "STATUS")
  return(invisible(output))
}

###: Get project data
###:
###: Get a description of the data used in the current project. Available informations are:
###:    \itemize{
###:    \item dataFile (\emph{string}): path to the data file
###:    \item header (\emph{array<character>}): vector of header names
###:    \item headerTypes (\emph{array<character>}): vector of header types
###:    \item observationNames (\emph{vector<string>}): vector of observation names
###:    \item observationTypes (\emph{vector<string>}): vector of observation types
###:\item nbSSDoses (\emph{int}) : number of doses (if there is a SS column)
###:}
###: @return A list describing project data.
###: @examples
###: \dontrun{
###: data = getData()
###: data
###: -> $dataFile
###:      "/path/to/data/file.txt"
###:    $header
###:      c("ID","TIME","CONC","SEX","OCC")
###:    $headerTypes
###:      c("ID","TIME","OBSERVATION","CATEGORICAL COVARIATE","IGNORE")
###:    $observationNames
###:      c("concentration")
###:    $observationTypes
###:      c(concentration = "continuous")
###: }
###: @seealso \code{\link{setData}}
###: @export
getData = function(){
  arguments = list()
  output = .processRequest("monolix", "getdata", arguments, "asynchronous")
  return(output)
}
# -------------------------------------------------------------------------------- #
# ================================================================================ #

