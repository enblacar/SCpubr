#' Create ridge plots.
#'
#' This function computes ridge plots based on the \pkg{ggridges} package.
#'
#' @inheritParams doc_function
#' @param colors.use \strong{\code{\link[base]{character}}} | Named vector of colors to use. Has to match the unique values of group.by or color.by (if used) when scale_type is set to categorical.
#' @param compute_quantiles \strong{\code{\link[base]{logical}}} | Whether to compute quantiles of the distribution and color the ridge plots by them.
#' @param compute_custom_quantiles \strong{\code{\link[base]{logical}}} | Whether to compute custom quantiles.
#' @param quantiles \strong{\code{\link[base]{numeric}}} | Numeric vector of quantiles.
#' @param compute_distribution_tails \strong{\code{\link[base]{logical}}} | Whether to compute distribution tails and color them.
#' @param prob_tails \strong{\code{\link[base]{numeric}}} | The accumulated probability that the tails should contain.
#' @param color_by_probabilities \strong{\code{\link[base]{logical}}} | Whether to color the ridges depending on the probability.
#' @param continuous_scale \strong{\code{\link[base]{logical}}} | Whether to color the ridges depending on a categorical or continuous scale.
#' @return A ggplot2 object.
#' @export
#'
#' @example /man/examples/examples_do_RidgePlot.R
do_RidgePlot <- function(sample,
                         feature,
                         group.by = NULL,
                         split.by = NULL,
                         assay = "SCT",
                         slot = "data",
                         continuous_scale = FALSE,
                         legend.title = NULL,
                         legend.ncol = NULL,
                         legend.nrow = NULL,
                         legend.byrow = FALSE,
                         legend.position = NULL,
                         legend.width = 1,
                         legend.length = 20,
                         legend.framewidth = 0.5,
                         legend.tickwidth = 0.5,
                         legend.framecolor = "grey50",
                         legend.tickcolor = "white",
                         legend.type = "colorbar",
                         colors.use = NULL,
                         font.size = 14,
                         font.type = "sans",
                         rotate_x_axis_labels = 45,
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
                         color_by_probabilities = FALSE,
                         viridis_color_map = "G",
                         viridis_direction = 1,
                         plot.grid = TRUE,
                         grid.color = "grey75",
                         grid.type = "dashed",
                         flip = FALSE){
  check_suggests(function_name = "do_RidgePlot")
  `%>%` <- magrittr::`%>%`

  # Check if the sample provided is a Seurat object.
  check_Seurat(sample = sample)

  # Check logical parameters.
  logical_list <- list("continuous_scale" = continuous_scale,
                       "compute_quantiles" = compute_quantiles,
                       "compute_custom_quantiles" = compute_custom_quantiles,
                       "compute_distribution_tails" = compute_distribution_tails,
                       "color_by_probabilities" = color_by_probabilities,
                       "plot.grid" = plot.grid,
                       "flip" = flip,
                       "legend.nrow" = legend.nrow)
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
                       "rotate_x_axis_labels" = rotate_x_axis_labels,
                       "legend.ncol" = legend.ncol,
                       "legend.nrow" = legend.nrow)
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
                         "ylab" = ylab,
                         "viridis_color_map" = viridis_color_map,
                         "grid.color" = grid.color,
                         "grid.type" = grid.type)

  check_colors(legend.tickcolor, parameter_name = "legend.tickcolor")
  check_colors(legend.framecolor, parameter_name = "legend.framecolor")

  check_parameters(parameter = font.type, parameter_name = "font.type")
  check_parameters(parameter = legend.type, parameter_name = "legend.type")
  if (!is.null(legend.position)){check_parameters(parameter = legend.position, parameter_name = "legend.position")}
  check_parameters(parameter = viridis_direction, parameter_name = "viridis_direction")
  check_parameters(parameter = viridis_color_map, parameter_name = "viridis_color_map")
  check_parameters(parameter = grid.type, parameter_name = "grid.type")
  check_parameters(parameter = rotate_x_axis_labels, parameter_name = "rotate_x_axis_labels")

  if (!is.null(colors.use)){check_colors(colors.use, parameter_name = "colors.use")}

  if (is.null(legend.position)){
    legend.position <- ifelse(isTRUE(continuous_scale), "bottom", "none")
  }

  data <- get_data_column_in_context(sample = sample,
                                     feature = feature,
                                     assay = assay,
                                     slot = slot,
                                     group.by = group.by,
                                     split.by = split.by)
  if (isTRUE(continuous_scale)){
    if (isFALSE(compute_quantiles)){
      p <- data %>%
           ggplot2::ggplot(mapping = ggplot2::aes(x = .data$feature,
                                                  y = .data$group.by,
                                                  fill = ggplot2::after_stat(x))) +
           ggridges::geom_density_ridges_gradient(color = "black",
                                                  size = 1.25) +
           ggplot2::scale_fill_viridis_c(option = viridis_color_map,
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
                                                    fill = ggplot2::after_stat(quantile))) +
             ggridges::stat_density_ridges(color = "black",
                                           size = 1.25,
                                           quantile_lines = TRUE,
                                           calc_ecdf = TRUE,
                                           geom = "density_ridges_gradient",
                                           quantiles = quantiles) +
             ggplot2::scale_fill_manual(values = viridis::viridis(n = length(quantiles) + 1, option = viridis_color_map, direction = viridis_direction),
                                        name = ifelse(is.null(legend.title), "Probability", legend.title),
                                        labels = unique(labels)) +
             ggplot2::guides(fill = ggplot2::guide_legend(title = ifelse(is.null(legend.title), "Probability", legend.title),
                                                          title.position = "top",
                                                          title.hjust = 0.5,
                                                          ncol = legend.ncol,
                                                          nrow = legend.nrow,
                                                          byrow = legend.byrow))
      } else if (isTRUE(compute_distribution_tails)){
        p <- data %>%
             ggplot2::ggplot(mapping = ggplot2::aes(x = .data$feature,
                                                    y = .data$group.by,
                                                    fill = ggplot2::after_stat(quantile))) +
             ggridges::stat_density_ridges(color = "black",
                                           size = 1.25,
                                           quantile_lines = TRUE,
                                           calc_ecdf = TRUE,
                                           geom = "density_ridges_gradient",
                                           quantiles = c(0 + prob_tails, 1 - prob_tails)) +
             ggplot2::scale_fill_manual(values = c("#134074", "grey75", "#721313"),
                                        labels = c(paste0("]0 , ", 0 + prob_tails, "]"),
                                                   paste0("]", 0 + prob_tails, ", ",  1 - prob_tails, "]"),
                                                   paste0("]", 1 - prob_tails, ", 1]")),
                                        name = ifelse(is.null(legend.title), "Probability", legend.title))  +
             ggplot2::guides(fill = ggplot2::guide_legend(title = ifelse(is.null(legend.title), "Probability", legend.title),
                                                          title.position = "top",
                                                          title.hjust = 0.5,
                                                          ncol = legend.ncol))
      } else if (isTRUE(color_by_probabilities)){
        p <- data %>%
             ggplot2::ggplot(mapping = ggplot2::aes(x = .data$feature,
                                                    y = .data$group.by,
                                                    fill = 0.5 - abs(0.5 - ggplot2::after_stat(ecdf)))) +
             ggridges::stat_density_ridges(color = "black",
                                           size = 1.25,
                                           calc_ecdf = TRUE,
                                           geom = "density_ridges_gradient") +
             ggplot2::scale_fill_viridis_c(option = viridis_color_map,
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
                                       size = 1.25) +
         ggplot2::scale_fill_manual(values = if (is.null(colors.use)) {generate_color_scale(if (is.null(group.by)){levels(sample)} else {if(is.factor(sample@meta.data[, group.by])){levels(sample@meta.data[, group.by])} else {unique(sample@meta.data[, group.by])}})} else {colors.use},
                                    name = legend.title) +
         ggplot2::guides(fill = ggplot2::guide_legend(title = legend.title,
                                                      title.position = "top",
                                                      title.hjust = 0.5,
                                                      ncol = legend.ncol))
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
                      axis.line.y = if (isFALSE(flip)) {ggplot2::element_line(color = "black")} else if (isTRUE(flip)) {ggplot2::element_blank()},
                      axis.line.x = if (isTRUE(flip)) {ggplot2::element_line(color = "black")} else if (isFALSE(flip)) {ggplot2::element_blank()},
                      axis.text.x = ggplot2::element_text(color = "black",
                                                          face = "bold",
                                                          angle = get_axis_parameters(angle = rotate_x_axis_labels, flip = flip)[["angle"]],
                                                          hjust = get_axis_parameters(angle = rotate_x_axis_labels, flip = flip)[["hjust"]],
                                                          vjust = get_axis_parameters(angle = rotate_x_axis_labels, flip = flip)[["vjust"]]),
                      axis.text.y = ggplot2::element_text(color = "black", face = "bold"),
                      axis.ticks = ggplot2::element_line(color = "black"),
                      panel.grid.major = ggplot2::element_blank(),
                      panel.grid.minor = ggplot2::element_blank(),
                      plot.title.position = "plot",
                      plot.title = ggplot2::element_text(face = "bold", hjust = 0),
                      plot.subtitle = ggplot2::element_text(hjust = 0),
                      plot.caption = ggplot2::element_text(hjust = 1),
                      panel.grid.major.y = ggplot2::element_blank(),
                      panel.grid.major.x = if (isTRUE(plot.grid)){ggplot2::element_line(color = grid.color, linetype = grid.type)} else {ggplot2::element_blank()},
                      text = ggplot2::element_text(family = font.type),
                      plot.caption.position = "plot",
                      legend.text = ggplot2::element_text(face = "bold"),
                      legend.position = legend.position,
                      legend.title = ggplot2::element_text(face = "bold"),
                      legend.justification = "center",
                      plot.margin = ggplot2::margin(t = 10, r = 10, b = 10, l = 10),
                      plot.background = ggplot2::element_rect(fill = "white", color = "white"),
                      panel.background = ggplot2::element_rect(fill = "white", color = "white"),
                      legend.background = ggplot2::element_rect(fill = "white", color = "white"),
                      strip.text = ggplot2::element_text(color = "black", face = "bold"))

  if (isTRUE(flip)){
    p <- p +
         ggplot2::coord_flip()
  }


  return(p)
}
