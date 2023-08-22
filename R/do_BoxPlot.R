
#' Generate Box Plots.
#'
#' @inheritParams doc_function
#' @inheritParams ggsignif::geom_signif
#'
#' @param boxplot.line.color \strong{\code{\link[base]{character}}} | Color of the borders of the boxplots if use_silhouette is FALSE.
#' @param outlier.color \strong{\code{\link[base]{character}}} | Color of the outlier dots.
#' @param outlier.alpha \strong{\code{\link[base]{numeric}}} | Alpha applied to the outliers.
#' @param boxplot.linewidth \strong{\code{\link[base]{numeric}}} | Width of the lines in the boxplots. Also controls the lines of the tests applied if use_test is set to true.
#' @param boxplot.width \strong{\code{\link[base]{numeric}}} | Width of the boxplots.
#' @param order \strong{\code{\link[base]{logical}}} | Whether to order the boxplots by average values. Can not be used alongside split.by.
#' @param use_silhouette \strong{\code{\link[base]{logical}}} | Whether to color the borders of the boxplots instead of the inside area.
#' @param use_test \strong{\code{\link[base]{logical}}} | Whether to apply a statistical test to a given pair of elements. Can not be used alongside split.by.
#'
#' @return A ggplot2 object.
#' @export
#'
#' @example /man/examples/examples_do_BoxPlot.R

do_BoxPlot <- function(sample,
                       feature,
                       group.by = NULL,
                       split.by = NULL,
                       assay = NULL,
                       slot = "data",
                       font.size = 14,
                       font.type = "sans",
                       axis.text.x.angle = 45,
                       colors.use = NULL,
                       na.value = "grey75",
                       plot.title = NULL,
                       plot.subtitle = NULL,
                       plot.caption = NULL,
                       xlab = NULL,
                       ylab = NULL,
                       legend.title = NULL,
                       legend.title.position = "top",
                       legend.position = "bottom",
                       boxplot.line.color = "black",
                       outlier.color = "black",
                       outlier.alpha = 0.5,
                       boxplot.linewidth = 0.5,
                       boxplot.width = NULL,
                       plot.grid = TRUE,
                       grid.color = "grey75",
                       grid.type = "dashed",
                       flip = FALSE,
                       order = FALSE,
                       use_silhouette = FALSE,
                       use_test = FALSE,
                       comparisons = NULL,
                       test = "wilcox.test",
                       map_signif_level = TRUE,
                       plot.title.face = "bold",
                       plot.subtitle.face = "plain",
                       plot.caption.face = "italic",
                       axis.title.face = "bold",
                       axis.text.face = "plain",
                       legend.title.face = "bold",
                       legend.text.face = "plain"){
  # Add lengthy error messages.
  withr::local_options(.new = list("warning.length" = 8170))

  check_suggests(function_name = "do_BoxPlot")
  # Check if the sample provided is a Seurat object.
  check_Seurat(sample = sample)

  # Check the assay.
  out <- check_and_set_assay(sample = sample, assay = assay)
  sample <- out[["sample"]]
  assay <- out[["assay"]]
  # Check slot.
  slot <- check_and_set_slot(slot = slot)
  # Check logical parameters.
  logical_list <- list("flip" = flip,
                       "plot.grid" = plot.grid,
                       "order" = order,
                       "use_silhouette" = use_silhouette,
                       "map_signif_level" = map_signif_level)
  check_type(parameters = logical_list, required_type = "logical", test_function = is.logical)
  # Check numeric parameters.
  numeric_list <- list("font.size" = font.size,
                       "outlier.alpha" = outlier.alpha,
                       "boxplot.linewidth" = boxplot.linewidth,
                       "boxplot.width" = boxplot.width,
                       "axis.text.x.angle" = axis.text.x.angle)
  check_type(parameters = numeric_list, required_type = "numeric", test_function = is.numeric)
  # Check character parameters.
  character_list <- list("feature" = feature,
                         "group.by" = group.by,
                         "split.by" = split.by,
                         "assay" = assay,
                         "slot" = slot,
                         "font.type" = font.type,
                         "colors.use" = colors.use,
                         "plot.title" = plot.title,
                         "plot.subtitle" = plot.subtitle,
                         "plot.caption" = plot.caption,
                         "xlab" = xlab,
                         "ylab" = ylab,
                         "legend.title" = legend.title,
                         "legend.title.position" = legend.title.position,
                         "legend.position" = legend.position,
                         "boxplot.line.color" = boxplot.line.color,
                         "outlier.color" = outlier.color,
                         "grid.color" = grid.color,
                         "grid.type" = grid.type,
                         "comparisons" = comparisons,
                         "test" = test,
                         "plot.title.face" = plot.title.face,
                         "plot.subtitle.face" = plot.subtitle.face,
                         "plot.caption.face" = plot.caption.face,
                         "axis.title.face" = axis.title.face,
                         "axis.text.face" = axis.text.face,
                         "legend.title.face" = legend.title.face,
                         "legend.text.face" = legend.text.face)
  check_type(parameters = character_list, required_type = "character", test_function = is.character)

  # Check the feature.
  feature <- check_feature(sample = sample, features = feature, permissive = TRUE)

  `%>%` <- magrittr::`%>%`
  
  
  check_colors(na.value, parameter_name = "na.value")
  check_colors(boxplot.line.color, parameter_name = "boxplot.line.color")
  check_colors(outlier.color, parameter_name = "outlier.color")
  check_colors(grid.color, parameter_name = "grid.color")

  check_parameters(parameter = font.type, parameter_name = "font.type")
  check_parameters(parameter = legend.position, parameter_name = "legend.position")
  check_parameters(parameter = grid.type, parameter_name = "grid.type")
  check_parameters(parameter = axis.text.x.angle, parameter_name = "axis.text.x.angle")
  check_parameters(plot.title.face, parameter_name = "plot.title.face")
  check_parameters(plot.subtitle.face, parameter_name = "plot.subtitle.face")
  check_parameters(plot.caption.face, parameter_name = "plot.caption.face")
  check_parameters(axis.title.face, parameter_name = "axis.title.face")
  check_parameters(axis.text.face, parameter_name = "axis.text.face")
  check_parameters(legend.title.face, parameter_name = "legend.title.face")
  check_parameters(legend.text.face, parameter_name = "legend.text.face")
  
  if (is.null(legend.title)){
    if (is.null(split.by)){
      if (is.null(group.by)) {
        legend.title <- "Groups"
      } else {
        legend.title <- group.by
      }
    } else {
      legend.title <- split.by
    }
  }
  
  # Check group.by.
  out <- check_group_by(sample = sample,
                        group.by = group.by,
                        is.heatmap = FALSE)
  sample <- out[["sample"]]
  group.by <- out[["group.by"]]
  
  if (is.null(colors.use)){
    if (is.null(split.by)){
      colors.use <- generate_color_scale(names_use = if (is.factor(sample@meta.data[, group.by])) {
        levels(sample@meta.data[, group.by])
      } else {
        sort(unique(sample@meta.data[, group.by]))
      })
    } else {
      colors.use <- generate_color_scale(names_use = if (is.factor(sample@meta.data[, split.by])) {levels(sample@meta.data[, split.by])} else {sort(unique(sample@meta.data[, split.by]))})
    }
  } else {
    check_colors(colors.use, parameter_name = "colors.use")
    check_consistency_colors_and_names(sample = sample, colors = colors.use, grouping_variable = ifelse(!is.null(split.by), split.by, group.by))
  }

  data <- get_data_column_in_context(sample,
                                     feature = feature,
                                     assay = assay,
                                     slot = slot,
                                     group.by = group.by,
                                     split.by = split.by)
  if (isTRUE(order) & is.null(split.by)){
    data <- data %>%
      dplyr::mutate("group.by" = factor(as.character(.data[["group.by"]]),
                                        levels = {data %>%
                                                  tibble::as_tibble() %>%
                                                  dplyr::group_by(.data[["group.by"]]) %>%
                                                  dplyr::summarise("median" = stats::median(.data[["feature"]], na.rm = TRUE)) %>%
                                                  dplyr::arrange(if(base::isFALSE(flip)){dplyr::desc(.data[["median"]])} else {.data[["median"]]}) %>%
                                                  dplyr::pull(.data[["group.by"]]) %>%
                                                  as.character()}))
  }
  if (isTRUE(order)){
    assertthat::assert_that(is.null(split.by),
                            msg = paste0(add_cross(), crayon_body("Parameter "),
                                         crayon_key("split.by"),
                                         crayon_body(" cannot be used alonside "),
                                         crayon_key("order"),
                                         crayon_body(".")))
  }

  if (!is.null(split.by)){
    assertthat::assert_that(base::isFALSE(order),
                            msg = paste0(add_cross(), crayon_body("Parameter "),
                                         crayon_key("split.by"),
                                         crayon_body(" cannot be used alonside "),
                                         crayon_key("order"),
                                         crayon_body(".")))
  }

  if (isTRUE(use_silhouette) & is.null(split.by)){
    p <- data %>%
         ggplot2::ggplot(mapping = ggplot2::aes(x = .data[["group.by"]],
                                                y = .data[["feature"]],
                                                color = .data[["group.by"]])) +
         ggplot2::scale_color_manual(values = colors.use, na.value = na.value) +
         ggplot2::geom_boxplot(outlier.color = outlier.color,
                               outlier.alpha = outlier.alpha,
                               width = boxplot.width,
                               lwd = boxplot.linewidth,
                               fatten = 1,
                               na.rm = TRUE)   +
         ggplot2::guides(color = ggplot2::guide_legend(title = legend.title,
                                                       title.position = legend.title.position,
                                                       title.hjust = 0.5))
  } else if (isTRUE(use_silhouette) & !is.null(split.by)){
    stop(paste0(add_cross(), crayon_body("Parameter "), crayon_key("use_silhouette"),  crayon_body("can not be used alongside "), crayon_key("split.by"), crayon_body(".")), call. = FALSE)
  } else if (base::isFALSE(use_silhouette)){
    if (is.null(split.by)){
      p <- data %>%
           ggplot2::ggplot(mapping = ggplot2::aes(x = .data[["group.by"]],
                                                  y = .data[["feature"]],
                                                  fill = .data[["group.by"]]))
    } else {
      p <- data %>%
           ggplot2::ggplot(mapping = ggplot2::aes(x = .data[["group.by"]],
                                                  y = .data[["feature"]],
                                                  fill = .data[["split.by"]]))
    }
    p <- p +
         ggplot2::scale_fill_manual(values = colors.use, na.value = na.value) +
         ggplot2::geom_boxplot(color = boxplot.line.color,
                               outlier.color = outlier.color,
                               outlier.alpha = outlier.alpha,
                               width = boxplot.width,
                               lwd = boxplot.linewidth,
                               fatten = 1,
                               key_glyph = "rect",
                               na.rm = TRUE)   +
         ggplot2::guides(fill = ggplot2::guide_legend(title = legend.title,
                                                      title.position = legend.title.position,
                                                      title.hjust = 0.5))
  }

   p <- p +
        ggplot2::labs(title = plot.title,
                      subtitle = plot.subtitle,
                      caption = plot.caption) +
        ggplot2::xlab(if (is.null(xlab)) {"Groups"} else (xlab)) +
        ggplot2::ylab(if (is.null(ylab)) {feature} else (ylab)) +
        ggplot2::theme_minimal(base_size = font.size) +
        ggplot2::theme(axis.title = ggplot2::element_text(color = "black",
                                                          face = axis.title.face),
                       axis.line.x = if (base::isFALSE(flip)) {ggplot2::element_line(color = "black")} else if (isTRUE(flip)) {ggplot2::element_blank()},
                       axis.line.y = if (isTRUE(flip)) {ggplot2::element_line(color = "black")} else if (base::isFALSE(flip)) {ggplot2::element_blank()},
                       axis.text.x = ggplot2::element_text(color = "black",
                                                           face = axis.text.face,
                                                           angle = get_axis_parameters(angle = axis.text.x.angle, flip = flip)[["angle"]],
                                                           hjust = get_axis_parameters(angle = axis.text.x.angle, flip = flip)[["hjust"]],
                                                           vjust = get_axis_parameters(angle = axis.text.x.angle, flip = flip)[["vjust"]]),
                       axis.text.y = ggplot2::element_text(color = "black", face = axis.text.face),
                       axis.ticks = ggplot2::element_line(color = "black"),
                       panel.grid.major = ggplot2::element_blank(),
                       panel.grid.major.y = if (base::isFALSE(flip)) {if (isTRUE(plot.grid)){ggplot2::element_line(color = grid.color, linetype = grid.type)}} else if (isTRUE(flip)) {ggplot2::element_blank()},
                       panel.grid.major.x = if (isTRUE(flip)) {if (isTRUE(plot.grid)){ggplot2::element_line(color = grid.color, linetype = grid.type)}} else if (base::isFALSE(flip)) {ggplot2::element_blank()},
                       plot.title.position = "plot",
                       plot.title = ggplot2::element_text(face = plot.title.face, hjust = 0),
                       plot.subtitle = ggplot2::element_text(face = plot.subtitle.face, hjust = 0),
                       plot.caption = ggplot2::element_text(face = plot.caption.face, hjust = 1),
                       legend.text = ggplot2::element_text(face = legend.text.face),
                       legend.title = ggplot2::element_text(face = legend.title.face),
                       panel.grid = ggplot2::element_blank(),
                       text = ggplot2::element_text(family = font.type),
                       plot.caption.position = "plot",
                       legend.position = legend.position,
                       legend.justification = "center",
                       plot.margin = ggplot2::margin(t = 10, r = 10, b = 10, l = 10),
                       plot.background = ggplot2::element_rect(fill = "white", color = "white"),
                       panel.background = ggplot2::element_rect(fill = "white", color = "white"),
                       legend.background = ggplot2::element_rect(fill = "white", color = "white"),
                       strip.text =ggplot2::element_text(color = "black", face = "bold"))

  if (isTRUE(flip)){
    p <- p + ggplot2::coord_flip()
  }

  if (isTRUE(use_test) & is.null(split.by)){
    if (!(is.null(comparisons))){
      p <- p +
           ggsignif::geom_signif(comparisons = comparisons,
                                 map_signif_level = map_signif_level,
                                 test = test,
                                 color = "black",
                                 size = boxplot.linewidth,
                                 textsize = font.size - 8,
                                 family = font.type,
                                 fontface = "bold")
    } else {
      stop(paste0(add_cross(), crayon_body("Please provide the pair of groups to test.")), call. = FALSE)
    }
  } else if (isTRUE(use_test) & !is.null(split.by)){
    stop(paste0(add_cross(), crayon_body("Tests can not be made if "), crayon_key("split.by"),  crayon_body(" is set.")), call. = FALSE)
  }
   

  return(p)
}
