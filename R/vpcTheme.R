#' Create a customized VPC theme
#' 
#' @param update list containing the plot elements to be updated. Run `create_vpc_theme()` with no arguments to show an overview of default plot elements.
#' 
#' @details
#' This function creates a theme that customizes how the VPC looks, i.e. legend, colors, fills, transparencies, linetypes an sizes, etc. The following arguments can be specified in the input list:
#' \itemize{
#'  \item{obs_color} color for observations points
#'  \item{obs_size} size for observation points
#'  \item{obs_shape} shape for observation points
#'  \item{obs_legend} legend name for observation points (if NA, no legend for obs)
#'  \item{cens_color} color for censored points
#'  \item{cens_shape} shape for censored points
#'  \item{cens_size} size for censored points
#'  \item{cens_legend} legend name for censored points (if NA, no legend for cens)
#'  \item{emp_perc_color} color for empirical percentiles
#'  \item{emp_perc_linetype} linetype for empirical percentiles
#'  \item{emp_perc_size} size for empirical percentiles
#'  \item{emp_perc_legend} legend name for empirical percentiles
#'  \item{theo_perc_color} color for theorical percentiles
#'  \item{theo_perc_linetype} linetype for theorical percentiles
#'  \item{theo_perc_size} size for theorical percentiles
#'  \item{theo_perc_legend} legend name for theorical percentiles
#'  \item{theo_pi_median_color} color for prediction interval median
#'  \item{theo_pi_median_fill} fill for prediction interval median
#'  \item{theo_pi_median_alpha} alpha for prediction interval median
#'  \item{theo_pi_median_legend} legend name for prediction interval median
#'  \item{theo_pi_perc_color} color for prediction interval percentiles
#'  \item{theo_pi_perc_fill} fill for prediction interval percentiles
#'  \item{theo_pi_perc_alpha} alpha for prediction interval percentiles
#'  \item{theo_pi_perc_legend} legend name for prediction interval percentiles
#'  \item{outlier_areas_color} color for outlier areas
#'  \item{outlier_areas_fill} fill for outlier areas
#'  \item{outlier_areas_perc_alpha} alpha for outlier areas
#'  \item{outlier_areas_perc_legend} legend name for outlier areas
#'  \item{outlier_dots_color} color for outliers dots
#'  \item{outlier_dots_size} size for outliers dots
#'  \item{outlier_dots_shape} shape for outliers dots
#'  \item{outlier_dots_legend} legend name for outliers dots 
#'  \item{bins_color} color for bin limits
#'  \item{bins_linetype} linetype for bin limits
#'  \item{bins_size} size for bin limits
#'  \item{bins_legend} legend name for bin limits
#' }
#' @return A list with vpc theme specifiers
#' @export
#' @examples 
#' theme <- new_vpc_theme(update = list(
#'   obs_color = "red",
#'   cens_color = rgb(70, 130, 180, maxColorValue = 255),
#'   theo_pi_median_alpha = 0.3
#' ))
createVpcTheme <- function (update = NULL) {
  default <- structure(list(  
    obs_color = rgb(70, 130, 180, maxColorValue = 255),
    obs_size = 1,
    obs_shape = 16,
    obs_legend = "Observed data",
    
    cens_color = rgb(212, 66, 66, maxColorValue = 255),
    cens_size = 1,
    cens_shape = 16,
    cens_legend = "Censored data",
    
    emp_perc_color = rgb(70, 130, 180, maxColorValue = 255),
    emp_perc_linetype = "solid",
    emp_perc_size = 1,
    emp_perc_point = TRUE,
    emp_perc_legend = "Empirical percentiles",
    
    theo_perc_color = rgb(0, 0, 0, maxColorValue = 255),
    theo_perc_linetype = "dashed",
    theo_perc_size = 1,
    theo_perc_point = TRUE,
    theo_perc_legend = "Theorical percentiles",
    
    theo_pi_median_color = NA,
    theo_pi_median_fill = rgb(255, 0, 0, maxColorValue = 255),
    theo_pi_median_alpha = 0.2,
    theo_pi_median_legend = "Prediction interval median",
    
    theo_pi_perc_color = NA,
    theo_pi_perc_fill = rgb(0, 128, 255, maxColorValue = 255),
    theo_pi_perc_alpha = 0.3,
    theo_pi_perc_legend = "Prediction interval percentiles",
    
    outlier_areas_color = NA,
    outlier_areas_fill = rgb(255, 0, 0, maxColorValue = 255),
    outlier_areas_alpha = 1,
    outlier_areas_legend = "Outlier Areas",
    
    outlier_dots_color = rgb(255, 0, 0, maxColorValue = 255),
    outlier_dots_legend = "Outlier dots",
    outlier_dots_size = 3,
    outlier_dots_shape = 1,
    
    bins_color = rgb(255, 0, 0, maxColorValue = 255),
    bins_linetype = "solid",
    bins_size = 0.2,
    bins_legend = NA
    
  ), class = "vpc_theme")
  theme <- default
  n <- names(theme)
  if(is.null(update)) {
    #    stop(paste0("Please specify a list with plot elements to update. Available elements: \n  - ", paste(n, collapse="\n  - ")))
    return(theme)
  }
  if(!is.null(update) & length(names(update)) > 0) {
    for(i in seq(names(update))) {
      if(names(update)[i] %in% n) {
        theme[[names(update)[i]]] <- update[[names(update)[i]]]
      } else {
        warning(paste0("`", names(update)[i],"` is not recognized as a plot element, ignoring."))
      }
    }
  }
  return(theme)
}

#' `continuous_theme()` returns a default VPC theme for continuous_theme data
#' @export
#' @rdname createVpcTheme
continuous_theme <- function() {
  return(createVpcTheme())
}

#' `discrete_theme()` returns a default VPC theme for discrete data
#' @export
#' @rdname createVpcTheme
discrete_theme <- function() {
  update <- list(
    emp_perc_legend = "Empirical probability",
    emp_perc_point = FALSE, 
    theo_perc_legend = "Theoretical probability",
    theo_perc_point = FALSE,
    theo_pi_median_legend = "Prediction interval",
    theo_pi_median_fill = rgb(0, 128, 255, maxColorValue = 255),
    theo_pi_median_alpha = 0.3
  )
  return(createVpcTheme(update))
}

#' `event_theme()` returns a default VPC theme for event data
#' @export
#' @rdname createVpcTheme
event_theme <- function() {
  update <- list(
    emp_perc_legend = "Empirical",
    emp_perc_point = FALSE, 
    theo_perc_legend = "Predicted median",
    theo_perc_point = FALSE,
    theo_pi_median_legend = "Prediction interval",
    theo_pi_median_fill = rgb(0, 128, 255, maxColorValue = 255),
    theo_pi_median_alpha = 0.3
  )
  return(createVpcTheme(update))
}