#' Wrapper for \link[Seurat]{VlnPlot}.
#'
#' @inheritParams doc_function
#' @param plot_boxplot Logical. Whether to plot a Box plot inside the violin or not.
#' @param pt.size  Size of points in the Violin plot.
#' @param y_cut  \strong{\code{\link[base]{numeric}}} | Vector with the values in which the Violins should be cut. Only works for one feature.
#' @param line_width Integer. Width of the lines drawn in the plot. Defaults to 1.
#' @param boxplot_width Integer. Width of the boxplots. Defaults to 0.2.


#' @return A ggplot2 object containing a Violin Plot.
#' @export
#'
#' @example man/examples/examples_do_ViolinPlot.R
do_ViolinPlot <- function(sample,
                          feature,
                          assay = NULL,
                          slot = NULL,
                          group.by = NULL,
                          colors.use = NULL,
                          pt.size = 0,
                          line_width = 1,
                          y_cut = NULL,
                          plot_boxplot = TRUE,
                          boxplot_width = 0.2,
                          legend.position = "none",
                          plot.title = NULL,
                          plot.subtitle = NULL,
                          plot.caption = NULL,
                          xlab = "Groups",
                          ylab = feature,
                          font.size = 14,
                          font.type = "sans",
                          rotate_x_axis_labels = TRUE,
                          plot.grid = TRUE,
                          grid.color = "grey75",
                          grid.type = "dashed"){
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
                       "rotate_x_axis_labels" = rotate_x_axis_labels,
                       "plot.grid" = plot.grid)
  check_type(parameters = logical_list, required_type = "logical", test_function = is.logical)
  # Check numeric parameters.
  numeric_list <- list("pt.size" = pt.size,
                       "y_cut" = y_cut,
                       "font.size" = font.size,
                       "line_width" = line_width,
                       "boxplot_width" = boxplot_width)
  check_type(parameters = numeric_list, required_type = "numeric", test_function = is.numeric)
  # Check character parameters.
  character_list <- list("legend.position" = legend.position,
                         "feature" = feature,
                         "group.by" = group.by,
                         "colors.use" = colors.use,
                         "plot.title" = plot.title,
                         "plot.subtitle" = plot.subtitle,
                         "plot.caption" = plot.caption,
                         "xlab" = xlab,
                         "ylab" = ylab,
                         "font.type" = font.type,
                         "grid.color" = grid.color,
                         "grid.type" = grid.color)
  check_type(parameters = character_list, required_type = "character", test_function = is.character)

  # Check the feature.
  feature <- check_feature(sample = sample, features = feature, permissive = TRUE)

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

    p <- Seurat::VlnPlot(sample,
                         features = feature,
                         cols = colors.use,
                         group.by = group.by,
                         pt.size = pt.size) +
         ggplot2::xlab(xlab) +
         ggplot2::ylab(ylab) +
         ggplot2::labs(title = plot.title,
                       subtitle = plot.subtitle,
                       caption = plot.caption) +
         ggplot2::theme_minimal(base_size = font.size) +
         ggplot2::theme(axis.text.x = ggplot2::element_text(color = "black",
                                                            face = "bold",
                                                            angle = ifelse(isTRUE(rotate_x_axis_labels), 45, 0),
                                                            hjust = ifelse(isTRUE(rotate_x_axis_labels), 1, 0.5),
                                                            vjust = ifelse(isTRUE(rotate_x_axis_labels), 1, 1)),
                        axis.text.y = ggplot2::element_text(face = "bold", color = "black"),
                        axis.title.y = ggplot2::element_text(face = "bold"),
                        axis.title.x = ggplot2::element_text(face = "bold"),
                        axis.line.y = ggplot2::element_blank(),
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

    # Modify line width of violin plots.
    p$layers[[1]]$aes_params$size <- line_width
    # Modify color of the line.
    p$layers[[1]]$aes_params$colour <- "black"

    if (plot_boxplot == TRUE){
      p <- p &
           ggplot2::geom_boxplot(lwd = line_width, width = boxplot_width, fill = "white", outlier.colour = NA, color = "black", fatten = 1)
    }

    if (!(is.null(y_cut))){
      p <- p +
           ggplot2::geom_hline(yintercept = y_cut,
                               linetype = "longdash",
                               colour = "black",
                               size = 1)
    }
  return(p)
}
