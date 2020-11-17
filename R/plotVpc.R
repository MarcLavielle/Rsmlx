plotVpc <- function(vpcData, obsData, obsName, timeName, settings, theme = NULL, minLog = 0.0) {
  # transform dataset (normalize continuous, discrete and event names)
  vpcData <- .prepareVpcData(vpcData)
  if(is.null(theme) || (class(theme) != "vpc_theme")) {
    theme <- createVpcTheme()
  }
  p <- ggplot(vpcData) + xlab(settings$xlab) + ylab(settings$ylab) +
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
  
  if (settings$grid) {
    p <- p + theme(
      panel.grid.major = element_line(size = 0.1, linetype = 'solid',
                                      colour = "grey"),
      panel.grid.minor = element_line(size = 0.1, linetype = 'solid',
                                      colour = "grey")
    )
  }
  
  if (settings$observedData) {
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
  
  if (settings$censoredData) {
    color <- rgb(212, 66, 66, maxColorValue = 255)
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
  
  if (settings$empiricalStats) {
    if (is.element("empirical_median", names(vpcData))) {
      p <- p +
        geom_line(
          aes(x = bins_middle, y = empirical_median, color = "c3"),
          linetype = theme$emp_perc_linetype
        )
      if (theme$emp_perc_point)
        p <- p +
          geom_point(aes(x = bins_middle, y = empirical_median), color = theme$emp_perc_color)
    }
    if (all(is.element(c("empirical_lower", "empirical_upper"), names(vpcData)))) {
      p <- p +
        geom_line(
          aes(x = bins_middle, y = empirical_lower, color = "c3"),
          linetype = theme$emp_perc_linetype
        ) +
        geom_line(
          aes(x = bins_middle, y = empirical_upper, color = "c3"),
          linetype = theme$emp_perc_linetype
        )
      if (theme$emp_perc_point)
        p <- p +
        geom_point(aes(x = bins_middle, y = empirical_upper), color = theme$emp_perc_color) +
        geom_point(aes(x = bins_middle, y = empirical_lower), color = theme$emp_perc_color)
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
  if (settings$theoricalMedian) {
    if (is.element("theoretical_median_median", names(vpcData))) {
      p <- p +
        geom_line(
          aes(x = bins_middle, y = theoretical_median_median, color = "c4"),
          linetype = theme$theo_perc_linetype
        )
      if (theme$theo_perc_point)
        p <- p +
          geom_point(aes(x = bins_middle, y = theoretical_median_median), color = theme$theo_perc_color)
    }
    if (all(is.element(c("theoretical_lower_median", "theoretical_upper_median"), names(vpcData)))) {
      p <- p +
        geom_line(
          aes(x = bins_middle, y = theoretical_lower_median, color = "c4"),
          linetype = theme$theo_perc_linetype
        ) +
        geom_line(
          aes(x = bins_middle, y = theoretical_upper_median, color = "c4"),
          linetype = theme$theo_perc_linetype
        )
      if (theme$theo_perc_point)
        p <- p +
          geom_point(aes(x = bins_middle, y = theoretical_upper_median), color = theme$theo_perc_color) +
          geom_point(aes(x = bins_middle, y = theoretical_lower_median), color = theme$theo_perc_color)
    }
    legendData = rbind(
      legendData,
      list(id = "c4", label = theme$theo_perc_legend,
           color = theme$theo_perc_color, linetype = theme$theo_perc_linetype,
           shape = NA, fill = NA, alpha = NA
    ))
  }
  if (settings$predictionInterval) {
    extremesDf <- do.call(rbind, by(
      vpcData, vpcData$split,
      function(df) {
        df <- df[order(df$bins_start),]
        rowf <- df[1,]
        rowl <- tail(df, n=1)
        rowf$bins_middle <- rowf$bins_start
        rowl$bins_middle <- rowl$bins_stop
        df <- rbind(rowf, df, rowl)
        return(df)
      }
    ))
    if (is.element("theoretical_median_median", names(vpcData))) {
      if (!settings$linearInterpolation) {
        p <- p +
          geom_rect(
            aes(xmin = bins_start, xmax = bins_stop,
                ymin = theoretical_median_lower, ymax = theoretical_median_upper,
                fill="c5"
            ),
            color = alpha(theme$theo_pi_median_color, theme$theo_pi_median_alpha)
          )
      } else {
        p <- p +
          geom_ribbon(
            data = extremesDf,
            aes(x = bins_middle, ymin = theoretical_median_lower,
                ymax = theoretical_median_upper, fill="c5"
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
      if (!settings$linearInterpolation) {
        p <- p +
          geom_rect(
            aes(xmin = bins_start, xmax = bins_stop,
                ymin = theoretical_lower_lower, ymax = theoretical_lower_upper,
                fill="c6"
            ),
            color = alpha(theme$theo_pi_perc_color, theme$theo_pi_perc_alpha),
          ) +
            geom_rect(
              aes(xmin = bins_start, xmax = bins_stop,
                  ymin = theoretical_upper_lower, ymax = theoretical_upper_upper,
                  fill="c6"
              ),
              color = alpha(theme$theo_pi_perc_color, theme$theo_pi_perc_alpha),
            )
      } else {
        p <- p +
          geom_ribbon(
            data = extremesDf,
            aes(x = bins_middle, ymin = theoretical_lower_lower,
                ymax = theoretical_lower_upper, fill="c6"
            ),
            color = alpha(theme$theo_pi_perc_color, theme$theo_pi_perc_alpha),
          ) +
          geom_ribbon(
            data = extremesDf,
            aes(x = bins_middle, ymin = theoretical_upper_lower,
                ymax = theoretical_upper_upper, fill="c6"
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

  if (settings$outlierAreas & settings$linearInterpolation) {
    # interpolation
    vpcDataInterp <- do.call(rbind, by(
      vpcData,
      subset(vpcData, select = intersect(c("split", "category"), names(vpcData))),
      function(df) {
        yNames <- names(df)[!(names(df) %in% c("category", "split") | grepl("bins", names(df)))]
        r <- subset(df, select = c("bins_middle", "split", yNames))
        if (nrow(df) > 1)
          r <- .interpolate(df, "bins_middle", yNames)
        r$split <- unique(df$split)
        if (!is.null(df$category)) r$category <- unique(df$category)
        return(r)
      }
    ))
    if (is.element("theoretical_median_median", names(vpcData))) {
      p <- p +
        geom_ribbon(
          data = vpcDataInterp,
          aes(x = bins_middle,
              ymin = pmin(theoretical_median_upper, empirical_median),
              ymax = empirical_median, fill = "c7"
          ),
          color = alpha(theme$outlier_areas_color, theme$outlier_areas_alpha)
        ) +
        geom_ribbon(
          data = vpcDataInterp,
          aes(x = bins_middle,
              ymin = pmin(empirical_median, theoretical_median_lower),
              ymax = theoretical_median_lower, fill = "c7"
          ),
          color = alpha(theme$outlier_areas_color, theme$outlier_areas_alpha)
        )
    }
    if (all(is.element(c("theoretical_lower_median", "theoretical_upper_median"), names(vpcData)))) {
      p <- p +
        geom_ribbon(
          data = vpcDataInterp,
          aes(x = bins_middle,
              ymin = pmin(theoretical_lower_upper, empirical_lower),
              ymax = empirical_lower, fill = "c7"
          ),
          color = alpha(theme$outlier_areas_color, theme$outlier_areas_alpha)
        ) +
        geom_ribbon(
          data = vpcDataInterp,
          aes(x = bins_middle,
              ymin = pmin(empirical_lower, theoretical_lower_lower),
              ymax = theoretical_lower_lower, fill = "c7"
          ),
          color = alpha(theme$outlier_areas_color, theme$outlier_areas_alpha)
        ) +
        geom_ribbon(
          data = vpcDataInterp,
          aes(x = bins_middle,
              ymin = pmin(theoretical_upper_upper, empirical_upper),
              ymax = empirical_upper, fill = "c7"
          ),
          color = alpha(theme$outlier_areas_color, theme$outlier_areas_alpha)
        ) +
        geom_ribbon(
          data = vpcDataInterp,
          aes(x = bins_middle,
              ymin = pmin(empirical_upper, theoretical_upper_lower),
              ymax = theoretical_upper_lower, fill = "c7"
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
  
  if (settings$outliersDots) {
    dataOutliersLower <- subset(
      vpcData,
      theoretical_lower_upper < empirical_lower | theoretical_lower_lower > empirical_lower
    )
    dataOutliersUpper <- subset(
      vpcData,
      theoretical_upper_upper < empirical_upper | theoretical_upper_lower > empirical_upper
    )
    dataOutliersMedian <- subset(
      vpcData,
      theoretical_median_upper < empirical_median | theoretical_median_lower > empirical_median
    )
    if (nrow(dataOutliersLower) > 0)
      p <- p +
      geom_point(
        data = dataOutliersLower,
        aes(x = bins_middle, y = empirical_lower, color = "c8", shape = "c8"),
        size = theme$outlier_dots_size
      ) 
    
    if (nrow(dataOutliersMedian) > 0)
      p <- p +
      geom_point(
        data = dataOutliersMedian,
        aes(x = bins_middle, y = empirical_median, color = "c8", shape = "c8"),
        size = theme$outlier_dots_size
      )
    
    if (nrow(dataOutliersUpper) > 0)
      p <- p +
      geom_point(
        data = dataOutliersUpper,
        aes(x = bins_middle, y = empirical_upper, color = "c8", shape = "c8"),
        size = theme$outlier_dots_size
      ) 
    
    if (any(c(nrow(dataOutliersLower), nrow(dataOutliersMedian), nrow(dataOutliersMedian)) > 0))
      legendData = rbind(
        legendData,
        list(id = "c8", label = theme$outlier_dots_legend, color = theme$outlier_dots_color,
             linetype = NA, shape = theme$outlier_dots_shape, fill = NA, alpha = NA
        ))
  }
  
  if (settings$binLimits) {
    color <- rgb(255, 0, 0, maxColorValue = 255)
    datavlines = do.call(rbind, by(
      subset(stats, select = c("bins_middle", "split")),
      stats$split,
      function(df) data.frame(bins_middle = unique(df$bins_middle), split = unique(df$split))
    ))
    p <- p + geom_vline(
      data = datavlines,
      aes(xintercept = bins_middle, color = "c9"),
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
  
  if (settings$ylogScale) p <- p + scale_y_log10()
  if (settings$xlogScale) p <- p + scale_x_log10()
  
  # legend
  colors <- setNames(legendData[is.na(legendData$fill),]$color, legendData[is.na(legendData$fill),]$id)
  shapes <- setNames(as.double(legendData[is.na(legendData$fill),]$shape), legendData[is.na(legendData$fill),]$id)
  labels <- setNames(legendData[is.na(legendData$fill),]$label, legendData[is.na(legendData$fill),]$id)
  linetypes <- setNames(legendData[is.na(legendData$fill),]$linetype, legendData[is.na(legendData$fill),]$id)
  linetypes[is.na(linetypes)] <- "blank"
  alphas <- as.double(legendData[!is.na(legendData$fill),]$alpha)
  p <- p +
    scale_color_manual(
      name = "", values = colors,
      labels = labels,
      guide=guide_legend(override.aes = list(shape = shapes, linetype = linetypes, fill = NA))
    ) +
    scale_shape_manual(values = shapes, guide=F)

  fills <- setNames(legendData[!is.na(legendData$fill),]$fill, legendData[!is.na(legendData$fill),]$id)
  p <- p + scale_fill_manual(
    name = "", values = fills, labels = legendData[!is.na(legendData$fill),]$label,
    guide = guide_legend(title = "legend", override.aes = list(alpha = alphas))
  )
  if (!settings$legend) {
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

.prepareVpcData <- function(vpcData) {
  vpcData <- .renameColumns(vpcData, "propCategory_empirical", "empirical_median")
  vpcData <- .renameColumns(vpcData, "propCategory_median", "theoretical_median_median")
  vpcData <- .renameColumns(vpcData, "propCategory_lower", "theoretical_median_lower")
  vpcData <- .renameColumns(vpcData, "propCategory_upper", "theoretical_median_upper")
  vpcData <- .renameColumns(vpcData, "survivalFunction", "empirical_median")
  vpcData <- .renameColumns(vpcData, "survivalFunction_median", "theoretical_median_median")
  vpcData <- .renameColumns(vpcData, "survivalFunction_lower", "theoretical_median_lower")
  vpcData <- .renameColumns(vpcData, "survivalFunction_upper", "theoretical_median_upper")
  vpcData <- .renameColumns(vpcData, "averageEventNumber", "empirical_median")
  vpcData <- .renameColumns(vpcData, "averageEventNumber_median", "theoretical_median_median")
  vpcData <- .renameColumns(vpcData, "averageEventNumber_lower", "theoretical_median_lower")
  vpcData <- .renameColumns(vpcData, "averageEventNumber_upper", "theoretical_median_upper")
  vpcData <- .renameColumns(vpcData, "time", "bins_middle")
  if (!is.element("bins_start", names(vpcData))) {
    vpcData$bins_start <- vpcData$bins_middle
  }
  if (!is.element("bins_stop", names(vpcData))) {
    vpcData$bins_stop <- vpcData$bins_middle
  }
  return(vpcData)
}