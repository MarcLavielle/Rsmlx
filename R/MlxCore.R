# ======================= MLXCORE ADVANCED FUNCTIONNALITIES ====================== #
# Prediction --------------------------------------------------------------------- #
###: Compute predictions from the structural model
###:
###: [\bold{MlxCore}][\emph{Prediction}] \cr 
###: Call the monolix prediction function to compute observation models values on observation times for each subject of a set of individuals.
###:
###: @param individualParameters Individual parameter values associated to each one of the individual parameters present in the project, for a set of subjects which must be 
###: coherent with the list of individuals ids passed in "individualIds" field (ie, this length of the subject set must be the sum of the subject number of all the 
###: individuals selected by the "individualIds" field).
###: This input field accepts a dataframe indexed by individual parameter names (columns) and subject indexes (rows).
###: @param individualIds [optional] {\emph{vector<int>}} Ids of the individuals for which observation models should be computed. By default, all the individuals present in the project are considered.
###: @return For each prediction names, a vector giving the computed prediction at observation times for each subject.
###: @examples
###: \dontrun{
###: ids = c(1,4)
###: individualValuesForAllIndiv = getEstimatedIndividualParameters()$saem
###: 
###: predictions = computePredictions( individualParameters = individualValuesForAllIndiv[ids,],
###:                                   individualIds = ids )
###: 
###: predictions
###:   -> $Cc
###:        [3.8,6.75,...,3.4,5.1,...]
###:        |    id=1    |    id=4   |
###: }
###: @seealso \code{\link{getIndividualParameterModel}} \code{\link{getEstimatedIndividualParameters}}
###: @export
computePredictions = function(individualParameters,individualIds=NULL) {
  
  if (!is.data.frame(individualParameters)){
    message("ERROR [MlxCore.computePredictions]: Unexpected type encountered for \"individualParameters\" field. 
            Please give a list of individual parameters values indexed by individual parameter names." )
    return()  
  }
  individualParameters$id <- NULL
  
  if (is.null(individualIds)){
    individualIds = 1:nrow(individualParameters)
  } else if (!is.vector(individualIds)){
    message("ERROR [MlxCore.computePredictions]: Unexpected type encountered for \"individualIds\" field. Please provide a vector of integers.")
    return()
  }
  else {
    for (i in 1:length(individualIds)){
      if (!.isInteger(individualIds[[i]])){
        message("ERROR [MlxCore.computePredictions]: Unexpected type encountered for \"individualIds\" field. Please provide a vector of integers.")
        return()
      }
    }
  }

  individualValues = as.list(individualParameters[,1:dim(individualParameters)[[2]]])
  arguments = list(individualValues,individualIds)
  output = .processRequest("monolix", "computeF", arguments, "asynchronous")
  
  tryCatch(expr = {
    for (iModel in 1:length(output))  # for each observation model
      output[[iModel]] = unlist(output[[iModel]])
  },
  
  error = function(err){
    output = list()
    message("ERROR [MlxCore.computePredictions]: Bad output encountered.")
  })
  
  return(output)
}
# -------------------------------------------------------------------------------- #
# ================================================================================ #