
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
                       assay = "SCT",
                       slot = "data",
                       font.size = 14,
                       font.type = "sans",
                       rotate_x_axis_labels = TRUE,
                       colors.use = NULL,
                       na.value = "grey75",
                       plot.title = NULL,
                       plot.subtitle = NULL,
                       plot.caption = NULL,
                       xlab = NULL,
                       ylab = NULL,
                       legend.title = if (is.null(split.by)){if (is.null(group.by)) {"Groups"} else {group.by}} else {split.by},
                       legend.title.position = "top",
                       legend.position = if (is.null(split.by)) {"none"} else {"bottom"},
                       boxplot.line.color = "black",
                       outlier.color = "black",
                       outlier.alpha = 0.5,
                       boxplot.linewidth = 1,
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
                       map_signif_level = TRUE){

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
                       "rotate_x_axis_labels" = rotate_x_axis_labels,
                       "plot.grid" = plot.grid,
                       "order" = order,
                       "use_silhouette" = use_silhouette,
                       "map_signif_level" = map_signif_level)
  check_type(parameters = logical_list, required_type = "logical", test_function = is.logical)
  # Check numeric parameters.
  numeric_list <- list("font.size" = font.size,
                       "outlier.alpha" = outlier.alpha,
                       "boxplot.linewidth" = boxplot.linewidth,
                       "boxplot.width" = boxplot.width)
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
                         "test" = test)
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


  if (is.null(group.by)){
    sample[["group.by"]] <- Seurat::Idents(sample)
    group.by <- "group.by"
  } else {
    sample[["group.by"]] <- sample@meta.data[, group.by]
    group.by <- "group.by"
  }

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
    check_consistency_colors_and_names(sample = sample, colors = colors.use, grouping_variable = group.by)
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
                                                  dplyr::summarise("mean" = mean(.data[["feature"]])) %>%
                                                  dplyr::arrange(if(isFALSE(flip)){dplyr::desc(.data[["mean"]])} else {.data[["mean"]]}) %>%
                                                  dplyr::pull(.data[["group.by"]]) %>%
                                                  as.character()}))
  }
  if (isTRUE(order)){
    assertthat::assert_that(is.null(split.by),
                            msg = "Parameter order can not be used alongside split.by.")
  }

  if (!is.null(split.by)){
    assertthat::assert_that(isFALSE(order),
                            msg = "Parameter order can not be used alongside split.by.")
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
                               key_glyph = "rect")
  } else if (isTRUE(use_silhouette) & !is.null(split.by)){
    stop("Parameter use_silhouetter can not be used alongside split.by.", call. = FALSE)
  } else if (isFALSE(use_silhouette)){
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
                               key_glyph = "rect")
  }

   p <- p +
        ggplot2::labs(title = plot.title,
                      subtitle = plot.subtitle,
                      caption = plot.caption) +
        ggplot2::xlab(if (is.null(xlab)) {"Groups"} else (xlab)) +
        ggplot2::ylab(if (is.null(ylab)) {feature} else (ylab))  +
        ggplot2::guides(fill = ggplot2::guide_legend(title = legend.title,
                                                     title.position = legend.title.position,
                                                     title.hjust = 0.5)) +
        ggplot2::theme_minimal(base_size = font.size) +
        ggplot2::theme(axis.title = ggplot2::element_text(color = "black",
                                                          face = "bold"),
                       axis.line.x = if (isFALSE(flip)) {ggplot2::element_line(color = "black")} else if (isTRUE(flip)) {ggplot2::element_blank()},
                       axis.line.y = if (isTRUE(flip)) {ggplot2::element_line(color = "black")} else if (isFALSE(flip)) {ggplot2::element_blank()},
                       axis.text.x = ggplot2::element_text(color = "black",
                                                           face = "bold",
                                                           angle = ifelse(isTRUE(rotate_x_axis_labels), 45, 0),
                                                           hjust = ifelse(isTRUE(rotate_x_axis_labels), 1, 0.5),
                                                           vjust = ifelse(isTRUE(rotate_x_axis_labels), 1, 1)),
                       axis.text.y = ggplot2::element_text(color = "black", face = "bold"),
                       axis.ticks = ggplot2::element_line(color = "black"),
                       panel.grid.major = ggplot2::element_blank(),
                       panel.grid.major.y = if (isFALSE(flip)) {if (isTRUE(plot.grid)){ggplot2::element_line(color = grid.color, linetype = grid.type)}} else if (isTRUE(flip)) {ggplot2::element_blank()},
                       panel.grid.major.x = if (isTRUE(flip)) {if (isTRUE(plot.grid)){ggplot2::element_line(color = grid.color, linetype = grid.type)}} else if (isFALSE(flip)) {ggplot2::element_blank()},
                       plot.title.position = "plot",
                       plot.title = ggplot2::element_text(face = "bold", hjust = 0),
                       plot.subtitle = ggplot2::element_text(hjust = 0),
                       plot.caption = ggplot2::element_text(hjust = 1),
                       panel.grid = ggplot2::element_blank(),
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
      stop("Please provide the pair of groups to test.", call. = FALSE)
    }
  } else if (isTRUE(use_test) & !is.null(split.by)){
    stop("Tests can not be made if split.by is set.", call. = FALSE)
  }
  return(p)
}
