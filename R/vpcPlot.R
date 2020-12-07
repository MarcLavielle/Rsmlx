#' Plot VPC
#'
#' @param vpcData (dataframe) dataframe with bins, vpc empirical and theoretical data.
#' @param obsData (dataframe) dataframe with observed data.
#' @param obsName (str) Name of observation in dataset header.
#' @param timeName (str) Name of time in dataset header.
#' @param curvesDisplay (\emph{vector(string)}) (\emph{optional}) (\emph{discrete data}) List of data curves to display in VPC plots.
#' (by default curvesDisplay = c("empiricalData", "predictionInterval"))
#' You must choose among the following curves:
#' \itemize{
#' \item \code{observedData} (\emph{continuous data}) Observed data. Specific to continuous data.
#' \item \code{censoredData} (\emph{continuous data}) Censored data. Specific to continuous data.
#' \item \code{empiricalData} Empirical data: 
#' empirical percentiles for continuous data; empirical probability for discrete data;
#' empirical curve for event data
#' \item \code{theoreticalData} TheoreticalData data. 
#' predicted percentiles for continuous data; theoretical probability for discrete data;
#' empiricalCurve for event data
#' \item \code{predictionInterval} Prediction interval.
#' \item \code{outliersDots} Add red dots indicating empirical percentiles that are outside prediction intervals.
#' \item \code{outlierAreas} Add red areas indicating empirical percentiles that are outside prediction intervals.
#' }
#' @param legend (\emph{bool}) (\emph{optional}) Add/remove legend (default FALSE).
#' @param grid (\emph{bool}) (\emph{optional}) Add/remove grid (default FALSE).
#' @param xlogScale (\emph{bool}) (\emph{optional}) Add/remove log scale for x axis (default FALSE).
#' @param ylogScale (\emph{bool}) (\emph{optional}) Add/remove log scale for x axis (default FALSE).
#' @param linearInterpolation (\emph{bool}) (\emph{optional}) If TRUE set piece wise display for prediction intervals,
#' else show bins as rectangular (default TRUE).
#' @param xlab (\emph{string}) (\emph{optional}) x-axis label. (`timeName` by default).
#' @param ylab (\emph{string}) (\emph{optional}) y-axis label. (`obsName` by default).
#' @param binLimits (\emph{bool}) (\emph{optional}) Add/remove vertical lines on the scatter plots to indicate the bins (default FALSE).
#' @param theme (\emph{vpc_theme}) (\emph{optional}) theme to be used in VPC. Expects list of class vpc_theme created with function createVpcTheme()
#' @return a ggplot2 object
#' @importFrom ggplot2 ggplot element_rect element_line geom_ribbon geom_point
#'             geom_rect geom_line theme aes xlab ylab facet_wrap facet_grid
#'             alpha scale_color_manual scale_shape_manual scale_fill_manual
#'             guide_legend geom_vline scale_x_log10 scale_y_log10 element_blank
#'             unit
#' @export
#' @examples
#' \dontrun{
#'   stats <- vpcStats(project="RsmlxDemo1.mlxtran")
#'   p <- plotVpc(stats$vpcPercentiles, stats$observations,
#'                stats$obsName, stats$timeName,
#'                curvesDisplay = c("empiricalData", "predictionInterval", "outlierAreas"),
#'                grid = TRUE, ylab = "concentration", xlab = "time (in hour)")
#' }
#' 
#' @seealso \link{vpcStats} \link{createVpcTheme} \link{vpc}
plotVpc <- function(vpcData, obsData, obsName, timeName,
                    curvesDisplay = c("empiricalData", "predictionInterval"),
                    legend = FALSE, grid = FALSE, xlogScale = FALSE, ylogScale = FALSE,
                    binLimits = FALSE, linearInterpolation = TRUE,
                    xlab = NULL, ylab = NULL, theme = NULL) {

  # Check arguments ------------------------------------------------------------
  allowedCurves <- c("observedData", "censoredData", "empiricalData",
                     "theoreticalData", "predictionInterval",
                     "outliersDots", "outlierAreas")
  curvesDisplay <- check_display(curvesDisplay, "curvesDisplay", allowedCurves)
  displayList <- get_display_list(curvesDisplay, allowedCurves)
  
  check_bool(legend, "legend")
  check_bool(grid, "grid")
  check_bool(xlogScale, "xlogScale")
  check_bool(ylogScale, "ylogScale")
  check_bool(binLimits, "binLimits")
  check_bool(linearInterpolation, "linearInterpolation")
  
  if (is.null(xlab)) xlab <- timeName
  if (is.null(ylab)) ylab <- obsName

  if(is.null(theme)) theme <- createVpcTheme()
  check_object_class(theme, "theme", "vpc_theme")
  
  # transform dataset (normalize continuous, discrete and event names)
  vpcData <- .prepareVpcData(vpcData)
  
  ## PLot ----------------------------------------------------------------------
  censored <- bins_middles <- bins_start <- bins_stop <- NULL
  empirical_median <- empirical_lower <- empirical_upper <- NULL
  theoretical_median_median <- theoretical_lower_median <- theoretical_upper_median <- NULL
  theoretical_median_piLower <- theoretical_upper_piLower <- theoretical_lower_piLower <- NULL
  theoretical_median_piUpper <- theoretical_upper_piUpper <- theoretical_lower_piUpper <- NULL
  p <- ggplot(vpcData) + xlab(xlab) + ylab(ylab) +
    theme(
      panel.background = element_rect(fill = "white", colour = "white",
                                      size = 2, linetype = "solid"),
      axis.line.x.bottom = element_line(color = "black", size = 0.3),
      axis.line.y.left = element_line(color = "black", size = 0.3),
      strip.background = element_rect(colour="white", fill="white", 
                                      size=1.5, linetype="solid"))
  wraps <- c()
  pos <- "top"
  if (length(unique(vpcData$category)) > 1) {
    wraps <- c(wraps, "category")
    vpcData["category"] <- factor(vpcData$category, levels = unique(vpcData$category))
    pos <- "left"
  }
  if (length(unique(vpcData$split)) > 1) {
    wraps <- c(wraps, "split")
    vpcData["split"] <- factor(vpcData$split, levels = unique(vpcData$split))
  }

  if (length(wraps) > 0) {
    p <- p + facet_wrap(
      as.formula(paste(".~", paste(wraps, collapse = " + "))),
      nrow=2, scales = "free", strip.position = pos
    )
    if ("category" %in% wraps)
      p <- p + ylab(NULL) + theme(strip.placement = "outside")
  }
  legendProp <- c("id", "label", "color", "linetype", "shape", "fill", "alpha")
  legendData <- as.data.frame(matrix(ncol = length(legendProp), nrow = 0, dimnames = list(NULL, legendProp)))
  
  if (grid) {
    p <- p + theme(
      panel.grid.major = element_line(size = 0.1, linetype = 'solid',
                                      colour = "grey"),
      panel.grid.minor = element_line(size = 0.1, linetype = 'solid',
                                      colour = "grey")
    )
  }
  
  if (displayList$observedData) {
    p <- p + geom_point(
      data = subset(obsData, censored == 0),
      aes(x = get(timeName), y = get(obsName), color = "c1", shape = "c1"),
      size = theme$obs_size
    )
    legendData = rbind(
      legendData,
      list(id = "c1", label = theme$obs_legend, color = theme$obs_color,
           linetype = NA, shape = theme$obs_shape, fill = NA, alpha = NA
    ))
  }
  
  if (displayList$censoredData) {
    if (nrow(subset(obsData, censored != 0)) > 0) {
      p <- p + geom_point(
        data = subset(obsData, censored != 0),
        aes(x = get(timeName), y = get(obsName), color = "c2", shape = "c2"),
        size = theme$cens_size
      )
      legendData = rbind(
        legendData,
        list(id = "c2", label = theme$cens_legend, color = theme$cens_color,
             linetype = NA, shape = theme$cens_shape, fill = NA, alpha = NA
        ))
    }
  }
  
  if (displayList$empiricalData) {
    if (is.element("empirical_median", names(vpcData))) {
      p <- p +
        geom_line(
          aes(x = bins_middles, y = empirical_median, color = "c3"),
          linetype = theme$emp_perc_linetype
        )
      if (theme$emp_perc_point)
        p <- p +
          geom_point(aes(x = bins_middles, y = empirical_median), color = theme$emp_perc_color)
    }
    if (all(is.element(c("empirical_lower", "empirical_upper"), names(vpcData)))) {
      p <- p +
        geom_line(
          aes(x = bins_middles, y = empirical_lower, color = "c3"),
          linetype = theme$emp_perc_linetype
        ) +
        geom_line(
          aes(x = bins_middles, y = empirical_upper, color = "c3"),
          linetype = theme$emp_perc_linetype
        )
      if (theme$emp_perc_point)
        p <- p +
        geom_point(aes(x = bins_middles, y = empirical_upper), color = theme$emp_perc_color) +
        geom_point(aes(x = bins_middles, y = empirical_lower), color = theme$emp_perc_color)
    }
    if (any(is.element(c("empirical_lower", "empirical_upper", "empirical_median"), names(vpcData)))) {
      legendData = rbind(
        legendData,
        list(id = "c3", label = theme$emp_perc_legend, 
             color = theme$emp_perc_color, linetype = theme$emp_perc_linetype,
             shape = NA, fill = NA, alpha = NA)
      )
    }
  }
  if (displayList$theoreticalData) {
    if (is.element("theoretical_median_median", names(vpcData))) {
      p <- p +
        geom_line(
          aes(x = bins_middles, y = theoretical_median_median, color = "c4"),
          linetype = theme$theo_perc_linetype
        )
      if (theme$theo_perc_point)
        p <- p +
          geom_point(aes(x = bins_middles, y = theoretical_median_median), color = theme$theo_perc_color)
    }
    if (all(is.element(c("theoretical_lower_median", "theoretical_upper_median"), names(vpcData)))) {
      p <- p +
        geom_line(
          aes(x = bins_middles, y = theoretical_lower_median, color = "c4"),
          linetype = theme$theo_perc_linetype
        ) +
        geom_line(
          aes(x = bins_middles, y = theoretical_upper_median, color = "c4"),
          linetype = theme$theo_perc_linetype
        )
      if (theme$theo_perc_point)
        p <- p +
          geom_point(aes(x = bins_middles, y = theoretical_upper_median), color = theme$theo_perc_color) +
          geom_point(aes(x = bins_middles, y = theoretical_lower_median), color = theme$theo_perc_color)
    }
    legendData = rbind(
      legendData,
      list(id = "c4", label = theme$theo_perc_legend,
           color = theme$theo_perc_color, linetype = theme$theo_perc_linetype,
           shape = NA, fill = NA, alpha = NA
    ))
  }
  if (displayList$predictionInterval) {
    extremesDf <- do.call(rbind, by(
      vpcData, vpcData$split,
      function(df) {
        df <- df[order(df$bins_start),]
        rowf <- df[1,]
        rowl <- utils::tail(df, n=1)
        rowf$bins_middles <- rowf$bins_start
        rowl$bins_middles <- rowl$bins_stop
        df <- rbind(rowf, df, rowl)
        return(df)
      }
    ))
    if (is.element("theoretical_median_median", names(vpcData))) {
      if (!linearInterpolation) {
        p <- p +
          geom_rect(
            aes(xmin = bins_start, xmax = bins_stop,
                ymin = theoretical_median_piLower, ymax = theoretical_median_piUpper,
                fill="c5"
            ),
            color = alpha(theme$theo_pi_median_color, theme$theo_pi_median_alpha)
          )
      } else {
        p <- p +
          geom_ribbon(
            data = extremesDf,
            aes(x = bins_middles, ymin = theoretical_median_piLower,
                ymax = theoretical_median_piUpper, fill="c5"
            ),
            color = alpha(theme$theo_pi_median_color, theme$theo_pi_median_alpha)
          )
      }
      legendData = rbind(
        legendData,
        list(id = "c5", label = theme$theo_pi_median_legend,
             fill = alpha(theme$theo_pi_median_fill, theme$theo_pi_median_alpha),
             linetype = NA, shape = NA, alpha = theme$theo_pi_median_alpha,
             color = NA
        ))

    }
    if (all(is.element(c("theoretical_lower_median", "theoretical_upper_median"), names(vpcData)))) {
      if (!linearInterpolation) {
        p <- p +
          geom_rect(
            aes(xmin = bins_start, xmax = bins_stop,
                ymin = theoretical_lower_piLower, ymax = theoretical_lower_piUpper,
                fill="c6"
            ),
            color = alpha(theme$theo_pi_perc_color, theme$theo_pi_perc_alpha),
          ) +
            geom_rect(
              aes(xmin = bins_start, xmax = bins_stop,
                  ymin = theoretical_upper_piLower, ymax = theoretical_upper_piUpper,
                  fill="c6"
              ),
              color = alpha(theme$theo_pi_perc_color, theme$theo_pi_perc_alpha),
            )
      } else {
        p <- p +
          geom_ribbon(
            data = extremesDf,
            aes(x = bins_middles, ymin = theoretical_lower_piLower,
                ymax = theoretical_lower_piUpper, fill="c6"
            ),
            color = alpha(theme$theo_pi_perc_color, theme$theo_pi_perc_alpha),
          ) +
          geom_ribbon(
            data = extremesDf,
            aes(x = bins_middles, ymin = theoretical_upper_piLower,
                ymax = theoretical_upper_piUpper, fill="c6"
            ),
            color = alpha(theme$theo_pi_perc_color, theme$theo_pi_perc_alpha),
          )
      }
      legendData = rbind(
        legendData,
        list(id = "c6", label = theme$theo_pi_perc_legend,
             fill = alpha(theme$theo_pi_perc_fill, theme$theo_pi_perc_alpha),
             linetype = NA, shape = NA, alpha = theme$theo_pi_perc_alpha,
             color = NA
        ))
      
    }
  }

  if (displayList$outlierAreas & linearInterpolation) {
    # interpolation
    vpcDataInterp <- do.call(rbind, by(
      vpcData,
      subset(vpcData, select = intersect(c("split", "category"), names(vpcData))),
      function(df) {
        yNames <- names(df)[!(names(df) %in% c("category", "split") | grepl("bins", names(df)))]
        r <- subset(df, select = c("bins_middles", "split", yNames))
        if (nrow(df) > 1)
          r <- .interpolate(df, "bins_middles", yNames)
        r$split <- unique(df$split)
        if (!is.null(df$category)) r$category <- unique(df$category)
        return(r)
      }
    ))
    if (is.element("theoretical_median_median", names(vpcData))) {
      p <- p +
        geom_ribbon(
          data = vpcDataInterp,
          aes(x = bins_middles,
              ymin = pmin(theoretical_median_piUpper, empirical_median),
              ymax = empirical_median, fill = "c7"
          ),
          color = alpha(theme$outlier_areas_color, theme$outlier_areas_alpha)
        ) +
        geom_ribbon(
          data = vpcDataInterp,
          aes(x = bins_middles,
              ymin = pmin(empirical_median, theoretical_median_piLower),
              ymax = theoretical_median_piLower, fill = "c7"
          ),
          color = alpha(theme$outlier_areas_color, theme$outlier_areas_alpha)
        )
    }
    if (all(is.element(c("theoretical_lower_median", "theoretical_upper_median"), names(vpcData)))) {
      p <- p +
        geom_ribbon(
          data = vpcDataInterp,
          aes(x = bins_middles,
              ymin = pmin(theoretical_lower_piUpper, empirical_lower),
              ymax = empirical_lower, fill = "c7"
          ),
          color = alpha(theme$outlier_areas_color, theme$outlier_areas_alpha)
        ) +
        geom_ribbon(
          data = vpcDataInterp,
          aes(x = bins_middles,
              ymin = pmin(empirical_lower, theoretical_lower_piLower),
              ymax = theoretical_lower_piLower, fill = "c7"
          ),
          color = alpha(theme$outlier_areas_color, theme$outlier_areas_alpha)
        ) +
        geom_ribbon(
          data = vpcDataInterp,
          aes(x = bins_middles,
              ymin = pmin(theoretical_upper_piUpper, empirical_upper),
              ymax = empirical_upper, fill = "c7"
          ),
          color = alpha(theme$outlier_areas_color, theme$outlier_areas_alpha)
        ) +
        geom_ribbon(
          data = vpcDataInterp,
          aes(x = bins_middles,
              ymin = pmin(empirical_upper, theoretical_upper_piLower),
              ymax = theoretical_upper_piLower, fill = "c7"
          ),
          color = alpha(theme$outlier_areas_color, theme$outlier_areas_alpha)
        )
    }
    legendData = rbind(
      legendData,
      list(id = "c7", label = theme$outlier_areas_legend,
           fill = alpha(theme$outlier_areas_fill, theme$outlier_areas_alpha),
           linetype = NA, shape = NA, alpha = theme$outlier_areas_alpha,
           color = NA
      ))
  }
  
  if (displayList$outliersDots) {
    dataOutliersLower <- subset(
      vpcData,
      theoretical_lower_piUpper < empirical_lower | theoretical_lower_piLower > empirical_lower
    )
    dataOutliersUpper <- subset(
      vpcData,
      theoretical_upper_piUpper < empirical_upper | theoretical_upper_piLower > empirical_upper
    )
    dataOutliersMedian <- subset(
      vpcData,
      theoretical_median_piUpper < empirical_median | theoretical_median_piLower > empirical_median
    )
    if (nrow(dataOutliersLower) > 0)
      p <- p +
      geom_point(
        data = dataOutliersLower,
        aes(x = bins_middles, y = empirical_lower, color = "c8", shape = "c8"),
        size = theme$outlier_dots_size
      ) 
    
    if (nrow(dataOutliersMedian) > 0)
      p <- p +
      geom_point(
        data = dataOutliersMedian,
        aes(x = bins_middles, y = empirical_median, color = "c8", shape = "c8"),
        size = theme$outlier_dots_size
      )
    
    if (nrow(dataOutliersUpper) > 0)
      p <- p +
      geom_point(
        data = dataOutliersUpper,
        aes(x = bins_middles, y = empirical_upper, color = "c8", shape = "c8"),
        size = theme$outlier_dots_size
      ) 
    
    if (any(c(nrow(dataOutliersLower), nrow(dataOutliersMedian), nrow(dataOutliersMedian)) > 0))
      legendData = rbind(
        legendData,
        list(id = "c8", label = theme$outlier_dots_legend, color = theme$outlier_dots_color,
             linetype = NA, shape = theme$outlier_dots_shape, fill = NA, alpha = NA
        ))
  }
  
  if (binLimits) {
    datavlines = do.call(rbind, by(
      subset(vpcData, select = c("bins_middles", "split")),
      vpcData$split,
      function(df) data.frame(bins_middles = unique(df$bins_middles), split = unique(df$split))
    ))
    p <- p + geom_vline(
      data = datavlines,
      aes(xintercept = bins_middles, color = "c9"),
      linetype = theme$bins_linetype,
      size = theme$bins_size
    )
    legendData = rbind(
      legendData,
      list(id = "c9", label = theme$bins_legend, color = theme$bins_color,
           linetype = theme$bins_linetype, shape = theme$bins_shape, fill = NA,
           alpha = NA
      ))
  }
  
  if (ylogScale) p <- p + scale_y_log10()
  if (xlogScale) p <- p + scale_x_log10()
  
  # legend
  colors <- stats::setNames(legendData[is.na(legendData$fill),]$color, legendData[is.na(legendData$fill),]$id)
  shapes <- stats::setNames(as.double(legendData[is.na(legendData$fill),]$shape), legendData[is.na(legendData$fill),]$id)
  labels <- stats::setNames(legendData[is.na(legendData$fill),]$label, legendData[is.na(legendData$fill),]$id)
  linetypes <- stats::setNames(legendData[is.na(legendData$fill),]$linetype, legendData[is.na(legendData$fill),]$id)
  linetypes[is.na(linetypes)] <- "blank"
  alphas <- as.double(legendData[!is.na(legendData$fill),]$alpha)
  p <- p +
    scale_color_manual(
      name = "", values = colors,
      labels = labels,
      guide=guide_legend(override.aes = list(shape = shapes, linetype = linetypes, fill = NA))
    ) +
    scale_shape_manual(values = shapes, guide=F)

  fills <- stats::setNames(legendData[!is.na(legendData$fill),]$fill, legendData[!is.na(legendData$fill),]$id)
  p <- p + scale_fill_manual(
    name = "", values = fills, labels = legendData[!is.na(legendData$fill),]$label,
    guide = guide_legend(title = "legend", override.aes = list(alpha = alphas))
  )
  if (!legend) {
    p <- p + theme(legend.position = "none")
  } else {
    p <- p + theme(
      legend.key = element_rect(fill = "white",colour = "white"),
      legend.title = element_blank(), legend.key.size = unit(0.5, "cm"),
      legend.position=c(.8,.75)
    )
  }
  return(p)
}

.interpolate <- function(data, xName, yNames, xInterp = NULL) {
  if (is.null(xInterp)) {
    x_range <- range(data[[xName]])
    x <- seq(x_range[1], x_range[2], length.out = 10000)
  } else {
    x <- xInterp
  }
  dfOut <- NULL
  for (i in seq(1, length(yNames))) {
    df <- as.data.frame(stats::approx(data[[xName]], data[[yNames[i]]], x))
    names(df) <- c(xName, yNames[i])
    if (is.null(dfOut)) dfOut <- df else dfOut <- merge(dfOut, df, by = c(xName))
  }
  return(dfOut)
}

.prepareVpcData <- function(vpcData) {
  vpcData <- .renameColumns(vpcData, "propCategory_empirical", "empirical_median")
  vpcData <- .renameColumns(vpcData, "propCategory_median", "theoretical_median_median")
  vpcData <- .renameColumns(vpcData, "propCategory_piLower", "theoretical_median_piLower")
  vpcData <- .renameColumns(vpcData, "propCategory_piUpper", "theoretical_median_piUpper")
  vpcData <- .renameColumns(vpcData, "survivalFunction", "empirical_median")
  vpcData <- .renameColumns(vpcData, "survivalFunction_median", "theoretical_median_median")
  vpcData <- .renameColumns(vpcData, "survivalFunction_piLower", "theoretical_median_piLower")
  vpcData <- .renameColumns(vpcData, "survivalFunction_piUpper", "theoretical_median_piUpper")
  vpcData <- .renameColumns(vpcData, "averageEventNumber", "empirical_median")
  vpcData <- .renameColumns(vpcData, "averageEventNumber_median", "theoretical_median_median")
  vpcData <- .renameColumns(vpcData, "averageEventNumber_piLower", "theoretical_median_piLower")
  vpcData <- .renameColumns(vpcData, "averageEventNumber_piUpper", "theoretical_median_piUpper")
  vpcData <- .renameColumns(vpcData, "time", "bins_middles")
  if (!is.element("bins_start", names(vpcData))) {
    vpcData$bins_start <- vpcData$bins_middles
  }
  if (!is.element("bins_stop", names(vpcData))) {
    vpcData$bins_stop <- vpcData$bins_middles
  }
  return(vpcData)
}

get_display_list <- function(display, allowedCurves) {
  displayList <- list()
  for (curve in allowedCurves) {
    if (curve %in% display) {
      displayList[curve] <- TRUE
    } else {
      displayList[curve] <- FALSE
    }
  }
  return(displayList)
}

check_display <- function(display, argname, allowedCurves) {
  if (!is.vector(display))
    stop("`", argname, "` must be a vector.", call. = FALSE)
  disaply <- check_in_vector(display, argname, allowedCurves, type = "warning")
  display <- display[display %in% allowedCurves]
  return(display)
}
