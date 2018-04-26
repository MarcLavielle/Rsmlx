# ================================ COVARIATE MODEL ============================== #
###: Get covariates information
###: 
###: Get the name, the type and the values of the covariates present in the project.
###: @return A list containing the following fields :
###: \itemize{
###: \item name : (\emph{vector<string>}) covariate names
###: \item type : (\emph{vector<string>}) covariate types. Existing types are "continuous", "continuoustransformed", "categorical", "categoricaltransformed" and "latent".
###: \item modalityNumber : (\emph{vector<int>}) number of modalities (for latent covariates only)
###: \item covariate : a data frame giving the values of continuous and categorical covariates for each subject.
###: Latent covariate values exist only if they have been estimated, ie if the covariate is used and if the population parameters have been estimated.
###: Call \code{\link{getEstimatedIndividualParameters}} to retrieve them.
###: }
###: @examples
###: \dontrun{
###: info = getCovariateInformation()
###: info
###:   -> $name
###:      c("sex","wt","lcat")
###:   -> $type
###:      c(sex = "categorical", wt = "continuous", lcat = "latent")
###:   -> $modalityNumber
###:      c(lcat = 2)
###:   -> $covariate
###:      id   sex    wt
###:       1    M   66.7
###:       .    .      .
###:       N    F   59.0
###: }
###: @export
getCovariateInformation = function(){
  arguments = list()
  output = .processRequest("monolix", "getcovariatesinformation", arguments, "asynchronous")
  
  if (!is.null(output) && length(output) > 0){
    if (!is.null(output$covariate) && length(output$covariate) > 0){
      for (i in 1:length(output$covariate)){
        if (length(output$covariate[[i]]) > 0)
          output$covariate[[i]] <- .evalNumerics(output$covariate[[i]])
      }
      output$covariate <- as.data.frame( append(output$ids,output$covariate) )
    }
    else
      output$covariate <- NULL
    
    output$ids <- NULL
  }
  else
    output <- NULL
  
  return(output)
}

###: Add continuous transformed covariate
###: 
###: Create a new continuous covariate by transforming an existing one. Transformed covariates cannot be use to produce new covariates.
###: Call \code{\link{getCovariateInformation}} to know which covariates can be transformed.
###: @param ... A list of comma-separated pairs \{transformedCovariateName = (\emph{string})"transformation"\}
###: @examples
###: \dontrun{
###: addContinuousTransformedCovariate( tWt2 = "3*exp(Wt)"  )
###: }
###: @seealso \code{\link{getCovariateInformation}} \code{\link{removeCovariate}}
###: @export
addContinuousTransformedCovariate = function(...){
  arguments = list(...)
  transformedCovNames = names(arguments)
  
  covariates = list()
  
  if (length(transformedCovNames) == 0){
    if (length(arguments) > 0){
      .error("Unexpected type encountered. Please give a list of strings indexed by transformed covariate names.")
      return(invisible(FALSE)) 
    } else {
      return(invisible(TRUE))
    }
  }
  else {
    for (i in 1:length(transformedCovNames)){
      if (transformedCovNames[i] == ""){
        .error(paste0("No covariate name found at position ",i,"."))
        return(invisible(FALSE))
      }
      else {
        transCov = arguments[[transformedCovNames[i]]]
        if (is.character(transCov) == FALSE){
          .error(paste0("Unexpected type encountered at position ",i,". Please give a string corresponding to the formula used to create the new covariate." ))
          return(invisible(FALSE)) 
        }
        else {
          covariates = append(covariates, list(list(name = transformedCovNames[i], formula = transCov)))
        }
      }
    }
    
    output = .processRequest("monolix", "addcontinuoustransformedcovariate", covariates, "synchronous", type = "STATUS")
    return(invisible(output))
  }
}

###: Add categorical transformed covariate
###: 
###: Create a new categorical covariate by transforming an existing one. Transformed covariates cannot be use to produce new covariates.
###: Call \code{\link{getCovariateInformation}} to know which covariates can be transformed.
###: @param ... A list of comma-separated pairs \{transformedCovariateName = \{ from = (\emph{array<(string)>})["basicCovariateNames"], transformed = (\emph{array<array<string>>})"transformation"\} \}
###: @examples
###: \dontrun{
###: addCategoricalTransformedCovariate( Country2 = list( reference = "A1", from = "Country",
###:                                                   transformed = list( A1 = c("A","B"), A2 = c("C") ) ) ) 
###: }
###: @seealso \code{\link{getCovariateInformation}} \code{\link{removeCovariate}}
###: @export
addCategoricalTransformedCovariate = function(...){
  arguments = list(...)
  transformedCovNames = names(arguments)
  res <- vector("list", length(transformedCovNames))
  errorType = ". Please give a list, indexed by transformed covariate names, of lists containing the name of the original covariate (\"from\" = {string}), the reference (\"reference\" = {string}) and the transformation to be applied on (\"transformed\" = {array<array<string>>})."
  
  if (length(transformedCovNames) == 0){
    if (length(arguments) > 0){
      .error("Unexpected type encountered. Please give a list of indexed by transformed covariate names.")
      return(invisible(FALSE)) 
    } else {
      return(invisible(TRUE))
    }
  }
  else {
    for (i in 1:length(transformedCovNames)){
      
      if (transformedCovNames[i] == ""){
        .error(paste0("No covariate name found at position ",i,"."))
        return(invisible(FALSE))
      }
      else {
        transCov = arguments[[transformedCovNames[i]]]
        if (is.list(transCov) == FALSE || "reference" %in% names(transCov)== FALSE || "from" %in% names(transCov) == FALSE || "transformed" %in% names(transCov) == FALSE ||
            is.character(transCov[["reference"]]) == FALSE || is.character(transCov[["from"]]) == FALSE || is.list(transCov[["transformed"]]) == FALSE ){
          .error(paste0("Unexpected type encountered at position ",i,errorType))
          return(invisible(FALSE)) 
        }
        else {
          transformation = transCov[["transformed"]]
          transformationName = names(transformation)
          
          res[[i]] = list("name" = transformedCovNames[i], "reference" = transCov[["reference"]], "origin" = transCov[["from"]], "groups" = vector("list", length(transformationName)))
          
          for (j in 1:length(transformationName)){
            modalityName <- transformationName[j]
            if (modalityName == ""){
              .error(paste0("No modality name found at position ",j,"."))
              return(invisible(FALSE))
            }
            if (is.vector(transformation[[j]]) == FALSE){
              .error(paste0("Unexpected type encountered at position ",i,errorType))
              return(invisible(FALSE)) 
            }
            else {
              res[[i]][["groups"]][[j]] = list("name"= modalityName, "modalities" = vector("list", length(transformation[[modalityName]])))
              
              for (k in 1:length(transformation[[modalityName]])){
                if (is.character((transformation[[modalityName]])[[k]]) == FALSE ){
                  .error(paste0("Unexpected type encountered at position ",i,errorType))
                  return(invisible(FALSE)) 
                }
                res[[i]][["groups"]][[j]][["modalities"]][[k]] <- (transformation[[modalityName]])[[k]]
              }
            }
          }
        }
      }
    }
    
    output = .processRequest("monolix", "adddiscretetransformedcovariate", res, "synchronous", type = "STATUS")
    return(invisible(output))
  }
}

###: Add mixture to the covariate model
###:  
###: Add a new latent covariate to the current model giving its name and its modality number.
###: @param ... A list of comma-separated pairs \{latentCovariateName = (\emph{int})modalityNumber\}
###: @examples
###: \dontrun{
###: addMixture(lcat = 2)
###: }
###: @seealso \code{\link{getCovariateInformation}} \code{\link{removeCovariate}}
###: @export
addMixture = function(...){
  arguments = list(...)
  if (is.null(names(arguments)) && length(arguments) == 1){
    arguments <-arguments[[1]]
  }
  
  mixtures = names(arguments)
  if (length(mixtures) == 0){
    if (length(arguments) > 0){
      .error("Unexpected type encountered. Please give a list of integers indexed by latent covariate names.")
      return(invisible(FALSE)) 
    } else {
      return(invisible(TRUE))
    }
  }
  for (j in 1:length(mixtures)){
    if (mixtures[[j]] == "" || (.isInteger(arguments[[mixtures[j]]]) == FALSE)){
      .error(paste0("Unexpected type encountered at position ",j,". Please give a list of comma-separated pairs {latentCovariateName = (int)modalityNumber}."))
      return(invisible(FALSE)) 
    }
  }

  output = .processRequest("monolix", "addmixtures", arguments, "synchronous", type = "STATUS")
  return(invisible(output))
}

###: Remove covariate
###: 
###: Remove some of the transformed covariates (discrete and continuous) and/or latent covariates.
###: Call \code{\link{getCovariateInformation}} to know which covariates can be removed.
###: @param ... A list of covariate names.
###: @examples
###: \dontrun{
###: removeCovariate("tWt","lcat1")
###: }
###: @seealso \code{\link{getCovariateInformation}} \code{\link{addContinuousTransformedCovariate}} \code{\link{addCategoricalTransformedCovariate}}
###: \code{\link{addMixture}}
###: @export
removeCovariate = function(...){
  arguments = list(...)
  
  if (length(arguments) == 0){
    return(invisible(TRUE))
  }
  else {
    for (i in 1:length(arguments)){
      if (is.character(arguments[[i]]) == FALSE){
        .error(paste0("Unexpected type encountered at position ",i,". Please give strings corresponding to a list of transformed covariate names."))
        return(invisible(FALSE))    
      }
    }

    output = .processRequest("monolix", "removetransformedcovariate", arguments, "synchronous", type = "STATUS")
    return(invisible(output))
  }
}
# ================================================================================ #