#' Loglikelihood profiling  function
#'
#' Allow to fix a parameter, run the scenario (Lin or IS) and find the threshold
#' @param project project name
#' @param parameters a vector of parameters to search. [Optional] By default, all the parameters are run.
#' @param settings a list of settings of the llp. It contains
#' \itemize{
#' \item method, method for the Log-Likelihood calculation, it can be "lin" or "is" (default value is linearisation)
#' \item dLLthreshold, the threshold of -2LL (default value at 3.841)
#' \item nbMaxIterations, the maximum number of iterations to find the dLLthreshold (default value at 10)
#' \item method, the method of calculation of the Log-Likelihood. It can be "lin" or "is" (default value is 'lin')
#' \item toldLL, the tolerance in terms of -2LL for the search of the threshold (default value is 1e-3)
#' \item tolParam, the relative tolerance for the parameter (1e-2)
#' }
#' @examples
#' \dontrun{
#' llp(project) => run the llp for the project for all parameters
#' llp(project, settings = list(method = 'is')) => run the llp for the project for all parameters with the LogLikelihood computed with the Importance Sampling method
#' llp(project, settings = list(maxNbIterations = 3)) => run the llp for the project for all parameters at most 3 iterations to search the target on both directions
#' llp(project, parameters = c("Cl_pop","V_pop")) => run the llp for Cl_pop and V_pop (if those are population parameters)
#' llp(project, parameters = c("Cl_pop","V_pop"), settings = list(dLLthreshold = 2)) => run the llp for Cl_pop and V_pop (if those are population parameters) with a dLLthreshold of 2
#' }
#' @export
llp <-function(project, parameters = NULL, settings=NULL){

  ###################################################################################
  # Initialization
  ###################################################################################
  # Check and initialize the settings
  if(!file.exists(project)){
    message(paste0("ERROR: project : ", project,' does not exists'));
    return(invisible(FALSE))}

  if(!is.null(settings)){
    if(!.checkLLPinput(inputName = "settings", inputValue = settings)){return(invisible(FALSE))}
  }
  if(is.null(settings$maxNbIterations)){ maxNbIterations <- 10 }else{ maxNbIterations <- settings$maxNbIterations }
  if(is.null(settings$dLLthreshold)){ dLLthreshold <- 3.841 }else{ dLLthreshold <- settings$dLLthreshold }
  if(is.null(settings$method)){ method <- 'lin' }else{ method <- settings$method }
  if(is.null(settings$toldLL)){ toldLL <- .1 }else{ toldLL <- settings$toldLL }
  if(is.null(settings$tolParam)){ tolParam <- .01 }else{ tolParam <- settings$tolParam}

  ##################################################################################
  # Project initialization
  ##################################################################################
  loadProject(project)
  # Check if linearization is possible
  if(method=="lin"){useLinearization <- T}else{useLinearization <- F}
  # Check if linearization is possible
  for(indexObservationModel in 1:length(getObservationInformation()$name)){useLinearization <- useLinearization & (getObservationInformation()$type[[indexObservationModel]]=="continuous")}

  # Define the scenario
  currentScenario <- getScenario()
  currentScenario$linearization <- useLinearization
  currentScenario$tasks <-as.list(currentScenario$tasks)
  currentScenario$tasks$logLikelihoodEstimation <- T
  currentScenario$tasks$standardErrorEstimation <- F
  currentScenario$tasks$plots <- F
  if(useLinearization){
    currentScenario$tasks$conditionalDistributionSampling <- F
    currentScenario$tasks$conditionalModeEstimation <-T
  }else{
    currentScenario$tasks$conditionalDistributionSampling <- T
    currentScenario$tasks$conditionalModeEstimation <-F
  }
  # Set the scenario
  setScenario(currentScenario)

  # Set the location of the saved informations
  expDirectory <- getProjectSettings()$directory
  expDirectoryLlp <- paste0(expDirectory,"/llp/")
  dir.create(expDirectoryLlp, showWarnings = F)
  setProjectSettings(directory = expDirectoryLlp)
  modelName <- substr(basename(project), 1, nchar(basename(project))-8)
  file_llp <- paste0(expDirectoryLlp,modelName,"_llp.csv")

  # Define the parameters to look at
  populationParameters <- getPopulationParameterInformation()$name[which(!(getPopulationParameterInformation()$method=="FIXED"))]
  if(!is.null(parameters)){
    if(!.checkLLPinput(inputName = "parameters", inputValue = parameters)){return(invisible(FALSE))}
    populationParameters <- intersect(parameters, populationParameters)
  }

  #################################################################################
  # Get the reference LL and the initial conditions
  #################################################################################
  validParameterIndex <- match(x = populationParameters, table=getPopulationParameterInformation()$name)
  referenceParameterValues<- getPopulationParameterInformation()$initialValue[validParameterIndex]
  LLref <- .getLL()
  LLtarget <- LLref + dLLthreshold
  setInitialEstimatesToLastEstimates()
  referenceParameterValues <- getPopulationParameterInformation()$initialValue[validParameterIndex]
  names(referenceParameterValues) <- getPopulationParameterInformation()$name[validParameterIndex]
  # See if the results were not already computed
  llpCompute <- TRUE
  # if(file.exists(file_llp)){
  #   llp_res <- read.table(file_llp,header = T, sep = ",")
  #   savedParamNames <- unique(llp_res$name)
  #   # Check if the number of parameter is similar
  #   if(compareParamName(projectParamNames=names(popParam),paramNames=savedParamNames)){
  #     llpCompute <- FALSE
  #     llp <- list()
  #     for(indexParam in 1:length(savedParamNames)){
  #       ind_Param <- which(llp_res$name==savedParamNames[indexParam])
  #       llp[[indexParam]] <- data.frame(param = llp_res$param[ind_Param],
  #                                       paramT = llp_res$paramT[ind_Param],
  #                                       LL = llp_res$LL[ind_Param],
  #                                       name = llp_res$name[ind_Param],
  #                                       paramInit = llp_res$paramInit[ind_Param],
  #                                       thresh = llp_res$thresh[ind_Param],
  #                                       tol = llp_res$tol[ind_Param],
  #                                       scenario = llp_res$scenario[ind_Param])
  #     }
  #     cat("LLP was previously computed \n")
  #   }
  # }

  ###############################################################################
  # Define  the output figures
  ###############################################################################
  nbFig <- length(populationParameters)
  x_NbFig <- ceiling(max(sqrt(nbFig),1)); y_NbFig <- ceiling(nbFig/x_NbFig)
  par(mfrow = c(x_NbFig, y_NbFig), oma = c(0, 3, 1, 1), mar = c(3, 1, 1, 3), mgp = c(2, 1, 0), xpd = NA)
  ###############################################################################
  # Compute if necessary
  ###############################################################################
  if(llpCompute){

    #########################################################################
    # Get the profile LL on all parameters
    #########################################################################
    llp <- list()
    for(indexParam in 1:length(populationParameters)){
      # Initialization
      paramName <- populationParameters[indexParam]
      paramValue <- as.double(referenceParameterValues[indexParam])
      cat("/**********************************************************************/ \n",paste0("LL search on ",paramName,'\n'))

      LLup <- LLdown <- vector(length = maxNbIterations+1)
      LLup[1] <- LLdown[1] <- LLref
      paramValueUp <- paramValueDown <- vector(length = maxNbIterations+1)
      paramValueUp[1] <- paramValueDown[1] <- paramValue

      # Search for increasing parameter values
      for(indexIter in 2:(maxNbIterations+1)){
        # Parameter estimation
        paramValueUp[indexIter] <- .paramVariationUpdate(paramName, indexIter, vect_x = paramValueUp[1:(indexIter-1)], vect_y = LLup[1:(indexIter-1)], y_target = LLtarget, bIncreasing = TRUE)
        cat(paste0("Upper bound search / Iteration n ",toString(indexIter)," / ",paramName, " = ",toString(floor(1000*paramValueUp[indexIter])/1000),'\n'))
        LLup[indexIter] <- .getLL(referenceParameterValues, paramName, paramValue = paramValueUp[indexIter])
        LL_err <- abs(LLup[indexIter]-LLtarget); dx_rerr <- abs((paramValueUp[indexIter]-paramValueUp[indexIter-1])/paramValue)
        if((LL_err<toldLL)||(dx_rerr<tolParam)){break}
      }
      indexIterUp <- indexIter

      # Search for decreasing parameter values
      for(indexIter in 2:(maxNbIterations+1)){
        # Parameter estimation
        paramValueDown[indexIter] <- .paramVariationUpdate(paramName, indexIter, vect_x = paramValueDown[1:(indexIter-1)], vect_y = LLdown[1:(indexIter-1)], y_target = LLtarget, bIncreasing = FALSE)
        cat(paste0("Lower bound search / Iteration n?",toString(indexIter)," / ",paramName, " = ",toString(floor(1000*paramValueDown[indexIter])/1000),'\n'))
        LLdown[indexIter] <- .getLL(referenceParameterValues, paramName, paramValue = paramValueDown[indexIter])
        LL_err <- abs(LLdown[indexIter]-LLtarget); dx_rerr <- abs((paramValueDown[indexIter]-paramValueDown[indexIter-1])/paramValue)
        if((LL_err<toldLL)||(dx_rerr<tolParam)){break}
      }
      indexIterDown <- indexIter

      # Concatenate informaton and sort
      paramEval <- c(paramValue, paramValueUp[2:indexIterUp], paramValueDown[2:indexIterDown])
      llEval <- c(LLref, LLup[2:indexIterUp], LLdown[2:indexIterDown])
      sortParam <-sort(paramEval,index.return=TRUE)


      llp[[indexParam]] <- data.frame(param = sortParam$x, paramInit = paramValue,
                                      LL = llEval[sortParam$ix], name = paramName,
                                      thresh = dLLthreshold, tolParam = tolParam, toldLL = toldLL, useLinearization = useLinearization)

      ###############################################################################
      # Save a summary
      ###############################################################################
      if(indexParam==1){
        write.table(llp[[indexParam]],file=file_llp,quote = FALSE,
                    append=FALSE,sep = ',',col.names = TRUE, row.names  = FALSE)
      }else{
        write.table(llp[[indexParam]],file=file_llp,quote = FALSE,
                    append=TRUE,sep = ',',col.names = FALSE,row.names  = FALSE)
      }
      ###############################################################################
      # Plot the result
      ###############################################################################
      .plotFigure(llp2plot = llp[[indexParam]])
      .displaySummary(llp2display = llp[[indexParam]])

    }
  }

  ###############################################################################
  # Replot the result
  ###############################################################################
  par(mfrow = c(x_NbFig, y_NbFig), oma = c(0, 3, 1, 1),  mar = c(3, 1, 1, 3), mgp = c(2, 1, 0),  xpd = NA)
  for(indexFigure in 1:nbFig){
    .plotFigure(llp2plot = llp[[indexFigure]])
  }

  ###############################################################################
  # Display a summary
  ###############################################################################
  for(indexParam in 1:nbFig){
    cat("/**********************************************************************/ \n")
    .displaySummary(llp2display = llp[[indexParam]])
  }
  return(llp)
}

#########################################################################
# Function of the evaluation of the LL
#########################################################################
.getLL <- function(referenceParameterValues=NULL, paramName=NULL, paramValue=NULL, seed=NULL){
  # Renilitialize the project with initial values
  if(!is.null(referenceParameterValues)){
    for(indexParam in 1:length(referenceParameterValues)){
      parameterName <- names(referenceParameterValues)[indexParam]
      # Unfix all the parameters and reinit them
      eval(parse(text=paste0('setPopulationParameterInformation(',parameterName,' = list(method = "MLE", initialValue = as.double(referenceParameterValues[indexParam])))')))
    }
  }
  # Fix te parameter to a value
  if(!is.null(paramName)&!is.null(paramValue)){
    eval(parse(text=paste0('setPopulationParameterInformation(',paramName,' = list(method = "FIXED", initialValue =',as.double(paramValue),'))')))
  }
  # Run the scenario and get the LL
  runScenario(TRUE)
  return(as.double(getEstimatedLogLikelihood()[[1]][1]))
}

#########################################################################
# Function of the evaluation of the parameter estimation
#########################################################################
.paramVariationUpdate <-function(paramName, indexIter, vect_x, vect_y, y_target, bIncreasing){
  parameterLaw <- tolower(.getTransform(paramName)[[1]])
  if(bIncreasing){gain = 1}else{gain = -1}
  if(indexIter == 2){
    # Step 1 / look around the proposed values
    out <- .getValueTransformedReferential(paramName, paramValue = vect_x[1], offset = gain*.2)
  }else if(indexIter>2){
    # Step 2 / Find the optima of the parabola
    v_x <- switch(parameterLaw,
                  "lognormal" = log(vect_x) - log(vect_x[1]),
                  "logitnormal" = log(vect_x/(1-vect_x)) - log(vect_x[1]/(1-vect_x[1])),
                  vect_x[1:length(vect_y)] - vect_x[1]
    )
    v_y <- vect_y - vect_y[1]
    v_y <- pmax(v_y, 0)
    polyVal <- lm(formula = v_y ~ I(v_x^2)-1)
    coeff <- max(1,polyVal$coefficients[[1]])
    offset <- sqrt(abs(vect_y[1]-y_target)/coeff)

    if(indexIter == 3){
      # Make a first initial guess based on the parabola
      out <- .getValueTransformedReferential(paramName, paramValue = vect_x[1], offset = gain*offset)
    }else{
      # Newton iterations
      f <- vect_y[length(vect_y)]-y_target
      df <- 2*coeff*gain*offset
      x_prop <- vect_x[length(vect_y)]-f*df/(1+df^2)

      # Sort the previous results and see if it is in the good bounds
      # Else wise, we get the value by dichotomy
      sortX <- sort(vect_x,index.return=TRUE)
      xVal <- sortX$x
      yVal <- vect_y[sortX$ix]

      if(bIncreasing){
        ind_min <- max(which(yVal<y_target))
        if(length(which(yVal>y_target))>0){
          ind_max <- min(which(yVal>y_target))
        }else{
          ind_max <- length(xVal)
        }
      }else{
        ind_max <- min(which(yVal<y_target))
        if(length(which(yVal>y_target))>0){
          ind_min <- max(which(yVal>y_target))
        }else{
          ind_min <- 1
        }
      }

      if(bIncreasing){
        if(yVal[ind_max]<y_target){# Last point not sufficient
          out <- .getValueTransformedReferential(paramName, paramValue = xVal[ind_max], offset = gain*.2)
        }else{
          if((x_prop>xVal[ind_min])&(x_prop<xVal[ind_max])){# Newton method
            out <- x_prop
          }else{# Dichotomy method
            param_a <- (yVal[ind_max]-yVal[ind_min])/(xVal[ind_max]-xVal[ind_min])
            param_b <- yVal[ind_max]-param_a*xVal[ind_max]
            out <- (y_target - param_b)/param_a
          }
        }
      }else{
        if(yVal[ind_min]<y_target){# Last point not sufficient
          out <- .getValueTransformedReferential(paramName, paramValue = xVal[ind_min], offset = gain*.2)
        }else{
          if((x_prop>xVal[ind_min])&(x_prop<xVal[ind_max])){# Newton method
            out <- x_prop
          }else{# Dichotomy method
            param_a <- (yVal[ind_max]-yVal[ind_min])/(xVal[ind_max]-xVal[ind_min])
            param_b <- yVal[ind_max]-param_a*xVal[ind_max]
            out <- (y_target - param_b)/param_a
          }
        }
      }

    }
  }
  return(out)
}

###################################################################################
# Function that display a summary
###################################################################################
.displaySummary <- function(llp2display){
  xVal <- llp2display$param; yVal <- llp2display$LL
  xlabel <- llp2display$name[1]
  LLtarget <- yVal[which(xVal == llp2display$paramInit[1])]+llp2display$thresh[1]
  toldLL <- llp2display$toldLL[1]
  val_ref <- llp2display$paramInit[1]
  indexMin <- which(yVal==min(yVal))
  diffTarget <- yVal-LLtarget

  if(diffTarget[1]< -toldLL){
    CI_min <- "]-Inf"
    se_min <- "]-Inf"
    rse_min <- "]-Inf"
  }else{
    val <- xVal[which(abs(diffTarget[1:indexMin])==min(abs(diffTarget[1:indexMin])))]
    CI_min <- paste0('[',toString(floor(val*1000)/1000))
    se_min <- paste0('[',toString(floor((val - val_ref)*1000)/1000))
    rse_min <- paste0('[',toString(floor(100*(val - val_ref)/val_ref*1000)/1000))
  }
  if(diffTarget[length(yVal)]< -toldLL){
    CI_max <- "+Inf["
    se_max <- "+Inf["
    rse_max <- "+Inf["
  }else{
    val <- xVal[which(abs(diffTarget[indexMin:length(diffTarget)])==min(abs(diffTarget[indexMin:length(diffTarget)])))+indexMin-1]
    CI_max <- paste0(toString(floor(val*1000)/1000),']')
    se_max <- paste0(toString(floor((val-val_ref)*1000)/1000),']')
    rse_max <- paste0(toString(floor(100*(val-val_ref)/val_ref*1000)/1000),']')
  }
  cat("Parameter ", toString(xlabel), "\nValue ",toString(floor(val_ref*1000)/1000),"\nCI = ",CI_min,",",CI_max," \n")
  cat("SE = ",se_min,",",se_max,"\n")
  cat("RSE = ",rse_min,",",rse_max,"\n")
}

###################################################################################
# Function that plot the results
###################################################################################
.plotFigure <- function(llp2plot){
  xVal <- llp2plot$param; yVal <- llp2plot$LL
  xlabel <- llp2plot$name[1]
  LLtarget <- yVal[which(xVal == llp2plot$paramInit[1])]+llp2plot$thresh[1]
  ymin <- min(min(yVal),LLtarget)-1
  ymax <- max(yVal)+1
  xmax <- max(xVal)+.01*mean(abs(xVal))
  xmin <- min(xVal)-.01*mean(abs(xVal))

  plot(xVal,yVal,type="b",col="blue",ylab = "-2LL",ylim = c(ymin,ymax),xlab = xlabel , xlim=c(xmin,xmax))
  lines(xVal, 0*xVal+LLtarget, col="black", type="l", lty=3)
  lines(x=(0*seq(0,1,length.out = 100)+llp2plot$paramInit[1]), y=seq(ymin, max(llp2plot$LL),length.out = 100), col="black", type="l", lty=3)
}

###################################################################################
# Check the inputs
###################################################################################
.checkLLPinput = function(inputName, inputValue){
  isValid = TRUE
  inputName = tolower(inputName)
  if(inputName == tolower("parameters")){
    if(is.vector(inputValue) == FALSE){
      message("ERROR: Unexpected type encountered. parameters must be a vector")
      isValid = FALSE
    }else if(length(intersect(x = getPopulationParameterInformation()$name, y=inputValue))==0){
      message("ERROR: parameters have no valid population parameters in its definition.")
      isValid = FALSE
    }
  }else if(inputName == tolower("settings")){
    if(is.list(inputValue) == FALSE){
      message("ERROR: Unexpected type encountered. settings must be a list")
      isValid = FALSE
    }else {
      for (i in 1:length(inputValue)){
        if(!.checkLLPsettings(settingName = names(inputValue)[i], settingValue = inputValue[[i]])){
          isValid = FALSE
        }
      }
    }
  }
  return(invisible(isValid))
}

.checkLLPsettings = function(settingName, settingValue){
  isValid = TRUE
  settingName = tolower(settingName)
  if(settingName == tolower("maxNbIterations")){
    if((is.double(settingValue) == FALSE)&&(is.integer(settingValue) == FALSE)){
      message("ERROR: Unexpected type encountered. The maximum number of interations must be an integer.")
      isValid = FALSE
    }else{
      if(!(as.integer(settingValue) == settingValue)){
        message("ERROR: Unexpected type encountered. The maximum number of interations must be an integer.")
        isValid = FALSE
      }else if(settingValue<1){
        message("ERROR:the maximum number of replicates must be a strictly positive integer.")
        isValid = FALSE
      }
    }
  }else if(settingName == tolower("dLLthreshold")){
    if(is.double(settingValue) == FALSE){
      message("ERROR: Unexpected type encountered. The dLLthreshold must be a double.")
      isValid = FALSE
    }else{
      if(settingValue<0){
        message("ERROR: dLLthreshold must be strictly positive.")
        isValid = FALSE
      }
    }
  }else if(settingName == tolower("method")){
    if(is.character(settingValue) == FALSE){
      message("ERROR: Unexpected type encountered. The method must be a character.")
      isValid = FALSE
    }else{
      if(length(intersect(tolower(settingValue), c('lin','is')))==0){
        message("ERROR: the type must be a character in the following list: 'lin', 'is'.")
        isValid = FALSE
      }
    }
  }else if(settingName == tolower("toldLL")){
    if(is.double(settingValue) == FALSE){
      message("ERROR: Unexpected type encountered. toldLL must be a double.")
      isValid = FALSE
    }else{
      if(settingValue<0){
        message("ERROR: toldLL must be strictly positive.")
        isValid = FALSE
      }
    }
  }else if(settingName == tolower("tolParam")){
    if(is.double(settingValue) == FALSE){
      message("ERROR: Unexpected type encountered. tolParam must be a double.")
      isValid = FALSE
    }else{
      if(settingValue<0){
        message("ERROR: tolParam must be strictly positive.")
        isValid = FALSE
      }
    }
  }else{
    warning("WARNING: ",settingName,' is not a valid setting')
  }
  return(isValid)
}

#################################################################################################"""
# Get the distribution of the parameter
#################################################################################################"""
.getTransform <- function(paramName){
  # Get the individual transformation
  paramTransform <- NULL
  if(length(grep(pattern="_pop", x = paramName))>0){# This is a population parameter
    indexParam <- which(sub(pattern="_pop",replacement = "",x = paramName)==getIndividualParameterModel()$name)
    paramTransform <- getIndividualParameterModel()$distribution[indexParam]
  }else if(length(grep(pattern="omega", x = paramName))>0){
    paramTransform <- "logNormal"
  }else if(length(grep(pattern="gamma", x = paramName))>0){
    paramTransform <- "logNormal"
  }else if(length(grep(pattern="beta", x = paramName))>0){
    paramTransform <- "normal"
  }else{
    paramTransform <- "logNormal"
  }

  return(paramTransform)
}

#################################################################################################"""
# Get the value in the transformed referential
#################################################################################################"""
.getValueTransformedReferential <- function(paramName, paramValue, offset=NULL){

  paramTransform <- .getTransform(paramName)
  if(is.null(offset)){
    offset <- 0
  }
  offset[which(is.na(offset))] <- 0

  if(paramTransform == "normal"){
    outValue <- paramValue + offset
  }else if(paramTransform == "logNormal"){
    outValue <- exp(log(paramValue) + offset)
  }else if(paramTransform == "logitNormal"){
    tValue <- log(paramValue/(1-paramValue)) + offset
    outValue <- exp(tValue)/(1+exp(tValue))
  }else{
    outValue <- paramValue + offset
  }
  return(outValue)
}
