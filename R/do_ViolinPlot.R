#' Wrapper for \link[Seurat]{VlnPlot}.
#'
#' @inheritParams doc_function
#' @param plot_boxplot \strong{\code{\link[base]{logical}}} | Whether to plot a Box plot inside the violin or not.
#' @param pt.size \strong{\code{\link[base]{numeric}}} | Size of points in the Violin plot.
#' @param y_cut  \strong{\code{\link[base]{numeric}}} | Vector with the values in which the Violins should be cut. Only works for one feature.
#' @param line_width \strong{\code{\link[base]{numeric}}} | Width of the lines drawn in the plot. Defaults to 1.
#' @param boxplot_width \strong{\code{\link[base]{numeric}}} | Width of the boxplots. Defaults to 0.2.
#' @param share.y.lims \strong{\code{\link[base]{logical}}} | When querying multiple features, force the Y axis of all of them to be on the same range of values (this being the max and min of all features combined).
#' @param legend.ncol \strong{\code{\link[base]{numeric}}} | Number of columns in the legend.


#' @return A ggplot2 object containing a Violin Plot.
#' @export
#'
#' @example man/examples/examples_do_ViolinPlot.R
do_ViolinPlot <- function(sample,
                          features,
                          assay = NULL,
                          slot = NULL,
                          group.by = NULL,
                          split.by = NULL,
                          colors.use = NULL,
                          pt.size = 0,
                          line_width = 1,
                          y_cut = rep(NA, length(features)),
                          plot_boxplot = TRUE,
                          boxplot_width = 0.2,
                          legend.position = "none",
                          plot.title = NULL,
                          plot.subtitle = NULL,
                          plot.caption = NULL,
                          xlab = rep(NA, length(features)),
                          ylab = rep(NA, length(features)),
                          font.size = 14,
                          font.type = "sans",
                          rotate_x_axis_labels = 45,
                          plot.grid = TRUE,
                          grid.color = "grey75",
                          grid.type = "dashed",
                          flip = FALSE,
                          ncol = NULL,
                          share.y.lims = FALSE,
                          legend.title = NULL,
                          legend.ncol = NULL){
  check_suggests(function_name = "do_ViolinPlot")
  # Check if the sample provided is a Seurat object.
  check_Seurat(sample = sample)

  # Check the assay.
  out <- check_and_set_assay(sample = sample, assay = assay)
  sample <- out[["sample"]]
  assay <- out[["assay"]]
  # Check slot.
  slot <- check_and_set_slot(slot = slot)
  # Check logical parameters.
  logical_list <- list("plot_boxplot" = plot_boxplot,
                       "plot.grid" = plot.grid,
                       "flip" = flip,
                       "share.y.lims" = share.y.lims)
  check_type(parameters = logical_list, required_type = "logical", test_function = is.logical)
  # Check numeric parameters.
  numeric_list <- list("pt.size" = pt.size,
                       "y_cut" = y_cut,
                       "font.size" = font.size,
                       "line_width" = line_width,
                       "boxplot_width" = boxplot_width,
                       "rotate_x_axis_labels" = rotate_x_axis_labels,
                       "ncol" = ncol,
                       "legend.ncol" = legend.ncol)
  check_type(parameters = numeric_list, required_type = "numeric", test_function = is.numeric)
  # Check character parameters.
  character_list <- list("legend.position" = legend.position,
                         "features" = features,
                         "group.by" = group.by,
                         "colors.use" = colors.use,
                         "plot.title" = plot.title,
                         "plot.subtitle" = plot.subtitle,
                         "plot.caption" = plot.caption,
                         "xlab" = xlab,
                         "ylab" = ylab,
                         "font.type" = font.type,
                         "grid.color" = grid.color,
                         "grid.type" = grid.color,
                         "split.by" = split.by,
                         "legend.title" = legend.title)
  check_type(parameters = character_list, required_type = "character", test_function = is.character)

  `%>%` <- magrittr::`%>%`

  # Check X and Y labels.
  if (sum(is.na(xlab)) == length(features)){
    xlab <- rep("Groups", length(features))
  } else {
    assertthat::assert_that(length(xlab) == length(features),
                            msg = "Please provide as many values for xlab as features provided. Use NAs to skip.")
  }

  if (sum(is.na(ylab)) == length(features)){
    ylab <- features
  } else {
    assertthat::assert_that(length(ylab) == length(features),
                            msg = "Please provide as many values for ylab as features provided. Use NAs to skip.")
  }

  if (sum(is.na(y_cut)) != length(features)){
    assertthat::assert_that(length(y_cut) == length(features),
                            msg = "Please provide as many values for y_cut as features provided. Use NAs to skip.")
  }

  # Check the feature.
  features <- check_feature(sample = sample, features = features, permissive = TRUE)

  if (is.null(group.by)){
    if (is.null(colors.use)){
      colors.use <- generate_color_scale(levels(sample))
    } else {
      colors.use <- check_consistency_colors_and_names(sample = sample, colors = colors.use)
    }
  } else if (!(is.null(group.by))){
    if (is.null(colors.use)){
      if (is.factor(sample@meta.data[, group.by])){
        names.use <- levels(sample@meta.data[, group.by])
      } else {
        names.use <- sort(unique(sample@meta.data[, group.by]))
      }
      colors.use <- generate_color_scale(names.use)
    } else {
      colors.use <- check_consistency_colors_and_names(sample = sample, colors = colors.use, grouping_variable = group.by)
    }
  }
  check_colors(grid.color, parameter_name = "grid.color")

  check_parameters(parameter = font.type, parameter_name = "font.type")
  check_parameters(parameter = legend.position, parameter_name = "legend.position")
  check_parameters(parameter = grid.type, parameter_name = "grid.type")
  check_parameters(parameter = rotate_x_axis_labels, parameter_name = "rotate_x_axis_labels")

  list.plots <- list()
  counter <- 0

  # Get the feature limits.
  max_values <- c()
  min_values <- c()
  for(feature in features){
    max_values <- append(max_values, max(get_data_column(sample = sample, feature = feature, assay = assay, slot = slot)[, "feature"], na.rm = TRUE))
    min_values <- append(min_values, min(get_data_column(sample = sample, feature = feature, assay = assay, slot = slot)[, "feature"], na.rm = TRUE))
  }
  limits <- c(min(min_values), max(max_values))

  for (feature in features){
    counter <- counter + 1
    p <- get_data_column_in_context(sample = sample,
                                    feature = feature,
                                    assay = assay,
                                    slot = slot,
                                    group.by = group.by,
                                    split.by = split.by) %>%
         ggplot2::ggplot(mapping = ggplot2::aes(x = .data$group.by,
                                                y = .data$feature,
                                                fill = if (!is.null(split.by)){.data$split.by} else {.data$group.by})) +
         ggplot2::geom_violin(color = "black",
                              lwd = line_width)
    if (isTRUE(plot_boxplot)){
      assertthat::assert_that(is.null(split.by),
                              msg = "Boxplots are not implemented when split.by is set.")
      p <- p +
           ggplot2::geom_boxplot(fill = "white",
                                 color = "black",
                                 lwd = line_width,
                                 width = boxplot_width,
                                 outlier.shape = NA,
                                 fatten = 1)
    }
    p <- p +
         ggplot2::xlab(if (!is.na(xlab[counter])) {xlab[counter]} else {"Groups"}) +
         ggplot2::ylab(if (!is.na(ylab[counter])) {ylab[counter]} else {feature}) +
         ggplot2::labs(title = plot.title,
                       subtitle = plot.subtitle,
                       caption = plot.caption) +
         ggplot2::guides(fill = ggplot2::guide_legend(title = legend.title,
                                                      ncol = legend.ncol,
                                                      title.position = "top")) +
         ggplot2::theme_minimal(base_size = font.size) +
         ggplot2::theme(axis.text.x = ggplot2::element_text(color = "black",
                                                            face = "bold",
                                                            angle = get_axis_parameters(angle = rotate_x_axis_labels, flip = flip)[["angle"]],
                                                            hjust = get_axis_parameters(angle = rotate_x_axis_labels, flip = flip)[["hjust"]],
                                                            vjust = get_axis_parameters(angle = rotate_x_axis_labels, flip = flip)[["vjust"]]),
                        axis.text.y = ggplot2::element_text(face = "bold", color = "black"),
                        axis.title.y = ggplot2::element_text(face = "bold"),
                        axis.title.x = ggplot2::element_text(face = "bold"),
                        axis.line.x = if (isFALSE(flip)) {ggplot2::element_line(color = "black")} else if (isTRUE(flip)) {ggplot2::element_blank()},
                        axis.line.y = if (isTRUE(flip)) {ggplot2::element_line(color = "black")} else if (isFALSE(flip)) {ggplot2::element_blank()},
                        plot.title = ggplot2::element_text(face = "bold", hjust = 0),
                        plot.subtitle = ggplot2::element_text(hjust = 0),
                        plot.caption = ggplot2::element_text(hjust = 1),
                        plot.title.position = "plot",
                        panel.grid.major.x = ggplot2::element_blank(),
                        panel.grid.minor = ggplot2::element_blank(),
                        panel.grid.major.y = if (isTRUE(plot.grid)){ggplot2::element_line(color = grid.color, linetype = grid.type)} else {ggplot2::element_blank()},
                        text = ggplot2::element_text(family = font.type),
                        plot.caption.position = "plot",
                        legend.text = ggplot2::element_text(face = "bold"),
                        legend.position = legend.position,
                        legend.title = ggplot2::element_text(face = "bold"),
                        legend.justification = "center",
                        plot.margin = ggplot2::margin(t = 10, r = 10, b = 10, l = 10),
                        axis.ticks = ggplot2::element_line(color = "black"),
                        axis.line = ggplot2::element_line(color = "black"),
                        plot.background = ggplot2::element_rect(fill = "white", color = "white"),
                        panel.background = ggplot2::element_rect(fill = "white", color = "white"),
                        legend.background = ggplot2::element_rect(fill = "white", color = "white"))

    if (isTRUE(share.y.lims)){
      p <- p +
           ggplot2::ylim(limits)
    }

    if (!is.na(y_cut[counter])){
      p <- p +
           ggplot2::geom_hline(yintercept = y_cut[counter],
                               linetype = "longdash",
                               colour = "black",
                               size = 1,
                               na.rm = TRUE)
    }

    if (isTRUE(flip)){
      p <- p +
           ggplot2::coord_flip()
    }
    list.plots[[feature]] <- p
  }

  if (length(features) > 1){
    p <- patchwork::wrap_plots(list.plots, ncol = ncol)
  }

  return(p)
}
