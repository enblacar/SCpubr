#' Create ridge plots.
#'
#' This function computes ridge plots based on the ggridges packages.
#'
#' @param sample Your Seurat object.
#' @param feature Character. Feature to plot.
#' @param group.by Character. Metadata variable to group the values by.
#' @param split.by Character. Metadata variable to split the values by.
#' @param assay Character. Assay to retrieve data from. Defaults to SCT.
#' @param slot Character. Slot from the assay to retrieve data from. Defaults to data.
#' @param continuous Logical. Whether we want coloring based on the continuous scale (feature) or categorical scale (groups).
#' @param legend.title Character. Title for the legend.
#' @param legend.type Character. Type of legend to display. One of: normal, colorbar, colorsteps.
#' @param legend.position Character. Position of the legend in the plot. Will only work if legend is set to TRUE.
#' @param legend.framewidth,legend.tickwidth Numeric. Width of the lines of the box in the legend.
#' @param legend.framecolor,legend.tickcolor Character. Color of the lines of the box in the legend.
#' @param legend.length,legend.width Numeric. Length and width of the legend. Will adjust automatically depending on legend side.
#' @param font.size Numeric. Overall font size of the plot.
#' @param font.type Character. Font family for the plot: sans, mono, serif.
#' @param rotate_x_axis_labels Logical. Whether to rotate X axis labels.
#' @param font.size Numeric. Overall font size of the plot.
#' @param font.type Character. Font family for the plot: sans, mono, serif.
#' @param rotate_x_axis_labels Logical. Whether to rotate X axis labels.
#' @param plot_legend. Logical. Whether to plot the legend or not.
#' @param colors.use Character. Named vector of colors to use. Has to match the unique values of group.by or color.by (if used) when scale_type is set to categorical.
#' @param na.value Character. Color for NAs.
#' @param plot.title,plot.subtitle,plot.caption Character. Title, Subtitle and caption to use in the plot.
#' @param xlab,ylab Character. Titles for the X and Y axis.
#' @param compute_quantiles Logical. Whether to compute quantiles of the distribution and color the ridge plots by them.
#' @param compute_custom_quantiles Logical. Whether to compute custom quantiles.
#' @param quantiles Numeric vector of quantiles.
#' @param compute_distribution_tails Logical. Whether to compute distribution tails and color them.
#' @param prob_tails Numeric. The accumulated probability that the tails should contain.
#' @param color_by_probabilities Logical. Whether to color the ridges depending on the probability.
#' @param viridis_direction Numeric. Either 1 or -1. Controls how the gradient of viridis scale is formed.
#' @param alpha Numeric. How transparent ridges are.
#' @return A ggplot2 object.
#' @export
#'
#' @example /man/examples/examples_do_EnrichmentHeatmap.R
do_RidgePlot <- function(sample,
                         feature,
                         group.by = NULL,
                         split.by = NULL,
                         assay = "SCT",
                         slot = "data",
                         continuous_scale = FALSE,
                         legend.title = NULL,
                         legend.position = NULL,
                         legend.width = 1,
                         legend.length = 20,
                         legend.framewidth = 1.5,
                         legend.tickwidth = 1.5,
                         legend.framecolor = "grey50",
                         legend.tickcolor = "white",
                         legend.type = "colorbar",
                         colors.use = NULL,
                         font.size = 14,
                         font.type = "sans",
                         rotate_x_axis_labels = FALSE,
                         plot.legend = TRUE,
                         plot.title = NULL,
                         plot.subtitle = NULL,
                         plot.caption = NULL,
                         xlab = NULL,
                         ylab = NULL,
                         compute_quantiles = FALSE,
                         compute_custom_quantiles = FALSE,
                         quantiles = c(0.25, 0.50, 0.75),
                         compute_distribution_tails = FALSE,
                         prob_tails = 0.025,
                         color_by_probabilities = TRUE,
                         viridis_direction = -1,
                         alpha = 1){
  `%>%` <- purrr::`%>%`

  # Checks for packages.
  check_suggests(function_name = "do_RidgePlot")
  # Check if the sample provided is a Seurat object.
  check_Seurat(sample = sample)

  # Check logical parameters.
  logical_list <- list("continuous_scale" = continuous_scale,
                       "rotate_x_axis_labels" = rotate_x_axis_labels,
                       "plot.legend" = plot.legend,
                       "compute_quantiles" = compute_quantiles,
                       "compute_custom_quantiles" = compute_custom_quantiles,
                       "compute_distribution_tails" = compute_distribution_tails,
                       "color_by_probabilities" = color_by_probabilities)
  check_type(parameters = logical_list, required_type = "logical", test_function = is.logical)
  # Check numeric parameters.
  numeric_list <- list("legend.width" = legend.width,
                       "legend.length" = legend.length,
                       "legend.framewidth" = legend.framewidth,
                       "legend.tickwidth" = legend.tickwidth,
                       "font.size" = font.size,
                       "quantiles" = quantiles,
                       "prob_tails" = prob_tails,
                       "viridis_direction" = viridis_direction,
                       "alpha" = alpha)
  check_type(parameters = numeric_list, required_type = "numeric", test_function = is.numeric)
  # Check character parameters.
  character_list <- list("feature" = feature,
                         "group.by" = group.by,
                         "split.by" = split.by,
                         "assay" = assay,
                         "slot" = slot,
                         "legend.title" = legend.title,
                         "legend.position" = legend.position,
                         "legend.framecolor" = legend.framecolor,
                         "legend.tickcolor" = legend.tickcolor,
                         "legend.type" = legend.type,
                         "colors.use" = colors.use,
                         "font.type" = font.type,
                         "plot.title" = plot.title,
                         "plot.subtitle" = plot.subtitle,
                         "plot.caption" = plot.caption,
                         "xlab" = xlab,
                         "ylab" = ylab)

  check_colors(legend.tickcolor, parameter_name = "legend.tickcolor")
  check_colors(legend.framecolor, parameter_name = "legend.framecolor")
  if (!is.null(colors.use)){check_colors(colors.use, parameter_name = "colors.use")}

  if (is.null(legend.position)){
    legend.position <- ifelse(isTRUE(continuous_scale), "bottom", "none")
  }

  data <- get_data_column_in_context(sample = sample,
                                     feature = feature,
                                     assay = "SCT",
                                     slot = "data",
                                     group.by = group.by,
                                     split.by = split.by)
  if (isTRUE(continuous_scale)){
    if (isFALSE(compute_quantiles)){
      p <- data %>%
           ggplot2::ggplot(mapping = ggplot2::aes(x = .data$feature,
                                                  y = .data$group.by,
                                                  fill = ..x..)) +
           ggridges::geom_density_ridges_gradient(color = "black",
                                                  size = 1.25,
                                                  alpha = alpha) +
           ggplot2::scale_fill_viridis_c(option = "G",
                                         direction = viridis_direction,
                                         name = feature)

      p <- modify_continuous_legend(p = p,
                                    legend.aes = "fill",
                                    legend.type = legend.type,
                                    legend.position = legend.position,
                                    legend.length = legend.length,
                                    legend.width = legend.width,
                                    legend.framecolor = legend.framecolor,
                                    legend.tickcolor = legend.tickcolor,
                                    legend.framewidth = legend.framewidth,
                                    legend.tickwidth = legend.tickwidth)
    } else if (isTRUE(compute_quantiles)){
      if (isTRUE(compute_custom_quantiles)){
        labels <- c()
        for (i in seq_along(quantiles)){
          if (i == 1){
            labels <- c(labels, paste0("[0 , ", quantiles[i], "["))
          } else if (i == length(quantiles)){
            labels <- c(labels, paste0("]", quantiles[i], ", 1]"))
          } else {
            labels <- c(labels, paste0("]", quantiles[i - 1], ", ", quantiles[i], "]"))
            labels <- c(labels, paste0("]", quantiles[i], ", ", quantiles[i + 1], "]"))
          }
        }
        p <- data %>%
             ggplot2::ggplot(mapping = ggplot2::aes(x = .data$feature,
                                                    y = .data$group.by,
                                                    fill = ..quantile..)) +
             ggridges::stat_density_ridges(color = "black",
                                           size = 1.25,
                                           quantile_lines = TRUE,
                                           calc_ecdf = TRUE,
                                           geom = "density_ridges_gradient",
                                           quantiles = quantiles,
                                           alpha = alpha) +
             ggplot2::scale_fill_manual(values = viridis::viridis(n = length(quantiles) + 1, option = "G", direction = viridis_direction),
                                        name = ifelse(is.null(legend.title), "Probability", legend.title),
                                        labels = unique(labels))
      } else if (isTRUE(compute_distribution_tails)){
        p <- data %>%
             ggplot2::ggplot(mapping = ggplot2::aes(x = .data$feature,
                                                    y = .data$group.by,
                                                    fill = ..quantile..)) +
             ggridges::stat_density_ridges(color = "black",
                                           size = 1.25,
                                           quantile_lines = TRUE,
                                           calc_ecdf = TRUE,
                                           geom = "density_ridges_gradient",
                                           quantiles = c(0 + prob_tails, 1 - prob_tails),
                                           alpha = alpha) +
             ggplot2::scale_fill_manual(values = c("#134074", "grey75", "#721313"),
                                        labels = c(paste0("]0 , ", 0 + prob_tails, "]"),
                                                   paste0("]", 0 + prob_tails , ", ",  1 - prob_tails, "]"),
                                                   paste0("]", 1 - prob_tails , ", 1]")),
                                        name = ifelse(is.null(legend.title), "Probability", legend.title))
      } else if (isTRUE(color_by_probabilities)){
        p <- data %>%
             ggplot2::ggplot(mapping = ggplot2::aes(x = .data$feature,
                                                    y = .data$group.by,
                                                    fill = 0.5 - abs(0.5 - ..ecdf..))) +
             ggridges::stat_density_ridges(color = "black",
                                           size = 1.25,
                                           calc_ecdf = TRUE,
                                           geom = "density_ridges_gradient",
                                           alpha = alpha) +
             ggplot2::scale_fill_viridis_c(option = "G",
                                           name = "Tail probability",
                                           direction = viridis_direction)
        p <- modify_continuous_legend(p = p,
                                      legend.title = legend.title,
                                      legend.aes = "fill",
                                      legend.type = legend.type,
                                      legend.position = legend.position,
                                      legend.length = legend.length,
                                      legend.width = legend.width,
                                      legend.framecolor = legend.framecolor,
                                      legend.tickcolor = legend.tickcolor,
                                      legend.framewidth = legend.framewidth,
                                      legend.tickwidth = legend.tickwidth)
      }
    }

  } else if (isFALSE(continuous_scale)){
    p <- data %>%
         ggplot2::ggplot(mapping = ggplot2::aes(x = .data$feature,
                                                y = .data$group.by,
                                                fill = .data$group.by)) +
         ggridges::geom_density_ridges(color = "black",
                                       size = 1.25,
                                       alpha = alpha) +
         ggplot2::scale_fill_manual(values = if (is.null(colors.use)) {generate_color_scale(if (is.null(group.by)){levels(sample)} else {if(is.factor(sample@meta.data[, group.by])){levels(sample@meta.data[, group.by])} else {unique(sample@meta.data[, group.by])}})} else {colors.use},
                                    name = legend.title)
  }

  if (!is.null(split.by)){
    # Facet.
    p <- p +
         ggplot2::facet_grid( ~ .data$split.by)
  }


  p <- p +
       ggplot2::labs(title = plot.title,
                     subtitle = plot.subtitle,
                     caption = plot.caption) +
       ggplot2::xlab(if (is.null(xlab)) {feature} else (xlab)) +
       ggplot2::ylab(if (is.null(ylab)) {"Groups"} else (ylab)) +
       ggplot2::theme_minimal(base_size = font.size) +
       ggplot2::theme(axis.title = ggplot2::element_text(color = "black",
                                                         face = "bold"),
                      axis.line.x = ggplot2::element_line(color = "black"),
                      axis.text.x = ggplot2::element_text(color = "black",
                                                          face = "bold",
                                                          angle = ifelse(isTRUE(rotate_x_axis_labels), 90, 0),
                                                          hjust = ifelse(isTRUE(rotate_x_axis_labels), 1, 0.5),
                                                          vjust = ifelse(isTRUE(rotate_x_axis_labels), 0.5, 1)),
                      axis.text.y = ggplot2::element_text(color = "black", face = "bold"),
                      axis.ticks = ggplot2::element_line(color = "black"),
                      panel.grid.major = ggplot2::element_blank(),
                      plot.title.position = "plot",
                      plot.title = ggplot2::element_text(face = "bold", hjust = 0),
                      plot.subtitle = ggplot2::element_text(hjust = 0),
                      plot.caption = ggplot2::element_text(hjust = 1),
                      panel.grid = ggplot2::element_blank(),
                      text = ggplot2::element_text(family = font.type),
                      plot.caption.position = "plot",
                      legend.text = ggplot2::element_text(face = "bold"),
                      legend.position = ifelse(isTRUE(plot.legend), legend.position, "none"),
                      legend.title = ggplot2::element_text(face = "bold"),
                      legend.justification = "center",
                      plot.margin = ggplot2::margin(t = 10, r = 10, b = 10, l = 10),
                      plot.background = ggplot2::element_rect(fill = "white", color = "white"),
                      panel.background = ggplot2::element_rect(fill = "white", color = "white"),
                      legend.background = ggplot2::element_rect(fill = "white", color = "white"),
                      strip.text = ggplot2::element_text(color = "black", face = "bold"))


  return(p)
}
