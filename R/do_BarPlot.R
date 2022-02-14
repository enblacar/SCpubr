#' Wrapper for computing publication ready bar plots.
#'
#' @param sample  Seurat object.
#' @param features  Main variable in the bar plot - This will define the height of the bars. Must be a metadata variable stored in `object@meta.data`.
#' @param group.by  Secondary variable to group the bar plot for - This will define the total number of bars. Must be a metadata variable stored in `object@meta.data`.
#' @param order.by  One of the unique values in group.by that will be used to reorder the bars according to its descending proportion.
#' @param labels.order  Vector of labels to explicitly state the order or the bars. Better fed after knowing the order. This is meant to be used to tweak the end plot if needed.
#' @param position  Either "fill" or "stack." Position "fill" will generate a bar plot with one column and the proportions of values for each group inside, while "stack" plots the bars together.
#' @param xlab  Title for the X axis.
#' @param ylab  Title for the Y axis.
#' @param colors.use  Palette of colors to use. It must match the group.by variable in terms of length and names.
#' @param legend Whether to plot the legend.
#' @param legend.title  Logical stating whether the legend title is shown or not.
#' @param legend.position  Position of the legend in the plot.
#' @param legend.ncol  Number of columns in the legend.
#' @param legend.text.size  Font size of the legend labels.
#' @param legend.title.size  Fantasize of the legend title.
#' @param legend.icon.size  Size of the icons in legend.
#' @param legend.position  Position of the legend in the plot. Will only work if legend is set to TRUE.
#' @param axis.text.size  Font size for axis text.
#' @param axis.title.size  Font size for axis title.
#' @param plot.title.size  Font size for the plot title.
#' @param legend.byrow  Logical stating whether the legend is filled by row or not.
#' @param plot.title  Title to use in the plot.
#' @param horizontal Whether to plot the Bar plot horizontally.
#' @param verbose Use warnings.
#'
#' @return A ggplot2 object containing a Bar plot.
#' @export
#'
#' @examples
#' \dontrun{
#' TBD
#' }
do_BarPlot <- function(sample,
                       features,
                       group.by = NULL,
                       labels.order = NULL,
                       order.by = NULL,
                       position = "stack",
                       xlab = NULL,
                       ylab = NULL,
                       plot.title = NULL,
                       legend = TRUE,
                       legend.position = "right",
                       legend.title = FALSE,
                       legend.ncol = 1,
                       legend.text.size = 12,
                       legend.title.size = 12,
                       axis.text.size = 16,
                       axis.title.size = 16,
                       plot.title.size = 18,
                       legend.icon.size = 4,
                       legend.byrow = FALSE,
                       colors.use = NULL,
                       horizontal = TRUE,
                       verbose = TRUE){
    # Checks for packages.
    check_suggests(function_name = "do_BarPlot")

    # Define pipe operator internally.
    `%>%` <- purrr::`%>%`

    # Check logical parameters.
    logical_list <- list("legend" = legend,
                         "legend.title" = legend.title,
                         "legend.byrow" = legend.byrow,
                         "horizontal" = horizontal)
    check_type(parameters = logical_list, required_type = "logical", test_function = is.logical)
    # Check numeric parameters.
    numeric_list <- list("legend.icon.size" = legend.icon.size,
                         "axis.text.size" = axis.text.size,
                         "axis.title.size" = axis.title.size,
                         "plot.title.size" = plot.title.size,
                         "legend.text.size" = legend.text.size,
                         "legend.title.size" = legend.title.size,
                         "legend.ncol" = legend.ncol)
    check_type(parameters = numeric_list, required_type = "numeric", test_function = is.numeric)
    # Check character parameters.
    character_list <- list("legend.position" = legend.position,
                           "features" = features,
                           "group.by" = group.by,
                           "colors.use" = colors.use,
                           "plot.title" = plot.title,
                           "xlab" = xlab,
                           "ylab" = ylab,
                           "order.by" = order.by,
                           "position" = position)
    check_type(parameters = character_list, required_type = "character", test_function = is.character)

    # If no color scale is provided, generate a custom one.
    if (is.null(colors.use)){
      if (is.null(group.by)){
        # Generate a color palette equal to the number of identities in the seurat object.
        colors.use <- generate_color_scale(names_use = unique(sample[[]][, feature]))
      } else if (!is.null(group.by)) {
        # Generate a color palette equal to the number of unique values in group.by variable.
        colors.use <- generate_color_scale(names_use = unique(sample[[]][, group.by]))
      }
    }

    for (feature in features){
      # Enforce the features to be part of the metadata.
      check_feature(sample = sample, features = feature, enforce_check = "metadata", enforce_parameter = "features")

      if (is.null(group.by)){
        factor_levels <- compute_factor_levels(sample = sample, feature = feature)
        if (isTRUE(horizontal)){factor_levels <- rev(factor_levels)}
        if (verbose){
          if (isTRUE(legend)){warning("Recommended settings without using group.by is to set legend to FALSE.")}
          if (position != "stack"){warning("Recommended settings without using group.by is to set position to 'stack'.")}
        }
        p <- sample@meta.data %>%
          dplyr::select(!!rlang::sym(feature)) %>%
          dplyr::group_by(!!rlang::sym(feature)) %>%
          dplyr::summarise(n = dplyr::n()) %>%
          dplyr::arrange(dplyr::desc(.data$n)) %>%
          dplyr::mutate(x_values = as.factor(!!(rlang::sym(feature)))) %>%
          dplyr::mutate(x_values = factor(.data$x_values, levels = factor_levels)) %>%
          ggplot2::ggplot(mapping = ggplot2::aes(x = .data$x_values, y = .data$n, fill = .data$x_values)) +
          ggplot2::geom_bar(position = position, stat="identity", width = 1,
                            colour="black",
                            size = 1) +
          ggpubr::theme_pubr(legend = legend.position) +
          ggplot2::scale_fill_manual(values = colors.use, na.value = "grey75") +
          ggplot2::theme(axis.title.x = ggplot2::element_text(size = axis.title.size, face = "bold"),
                         axis.title.y = ggplot2::element_text(size = axis.title.size, face = "bold"),
                         axis.text = ggplot2::element_text(size = axis.text.size, face = "bold"),
                         legend.text = ggplot2::element_text(size = legend.text.size, face = "bold"),
                         legend.title = ggplot2::element_text(size = legend.title.size, face = "bold"),
                         plot.title = ggplot2::element_text(size = plot.title.size, face = "bold", hjust = 0.5)) +
          ggplot2::guides(color = ggplot2::guide_legend(ncol = legend.ncol,
                                                        byrow = legend.byrow,
                                                        override.aes = list(size = legend.icon.size)))
      } else {
        check_feature(sample = sample, features = group.by, enforce_check = "metadata", enforce_parameter = "group.by")
        # Check the order of labels.
        if (!(is.null(labels.order))){
          factor_levels <- labels.order
        } else {
          if (!is.null(order.by)){
            if (!(order.by %in% unique(sample@meta.data[, group.by]))){
              stop("Parameter order.by (", order.by, ") not present in the unique values of parameter group.by (", group.by, ").")
            }
            factor_levels <- compute_factor_levels(sample = sample, feature = feature, group.by = group.by, order.by = order.by)
          } else {
            factor_levels <- compute_factor_levels(sample = sample, feature = feature, group.by = group.by)
          }
        }
        if (isTRUE(horizontal)){factor_levels <- rev(factor_levels)}
        if (verbose){
          if (isFALSE(legend)){warning("Recommended settings when using group.by is to set legend to TRUE.")}
          if (position != "fill"){warning("Recommended settings when using group.by is to set position to 'fill'.")}
        }

        p <- sample@meta.data %>%
          dplyr::select(!!rlang::sym(feature), !!rlang::sym(group.by)) %>%
          dplyr::group_by(!!rlang::sym(group.by), !!rlang::sym(feature)) %>%
          dplyr::summarise(n = dplyr::n()) %>%
          dplyr::arrange(dplyr::desc(.data$n)) %>%
          dplyr::mutate(x_values = as.factor(!!(rlang::sym(feature)))) %>%
          dplyr::mutate(x_values = factor(.data$x_values, levels = factor_levels)) %>%
          ggplot2::ggplot(mapping = ggplot2::aes(x = .data$x_values, y = .data$n, fill = !!rlang::sym(group.by))) +
          ggplot2::geom_bar(position = position, stat="identity", width = 1,
                            colour="black",
                            size = 1) +
          ggpubr::theme_pubr(legend = legend.position) +
          ggplot2::scale_fill_manual(values = colors.use, na.value = "grey75") +
          ggplot2::theme(axis.title.x = ggplot2::element_text(size = axis.title.size, face = "bold"),
                         axis.title.y = ggplot2::element_text(size = axis.title.size, face = "bold"),
                         axis.text = ggplot2::element_text(size = axis.text.size, face = "bold"),
                         legend.text = ggplot2::element_text(size = legend.text.size, face = "bold"),
                         legend.title = ggplot2::element_text(size = legend.title.size, face = "bold"),
                         plot.title = ggplot2::element_text(size = plot.title.size, face = "bold", hjust = 0.5)) +
          ggplot2::guides(color = ggplot2::guide_legend(ncol = legend.ncol,
                                                        byrow = legend.byrow,
                                                        override.aes = list(size = legend.icon.size)))
      }

      # Add X axis label.
      if (!is.null(xlab)){
        p <- p & ggplot2::xlab(xlab)
      } else {
        p <- p & ggplot2::xlab("")
      }
      # Add Y axis label.
      if (!is.null(ylab)){
        p <- p & ggplot2::ylab(ylab)
      } else {
        p <- p & ggplot2::ylab("")
      }
      # Add plot title.
      if (!is.null(plot.title)){
        p <- p + ggplot2::ggtitle(plot.title)
      }
      # Whether to flip the axis or not.
      if (isTRUE(horizontal)){
        p <- p + ggplot2::coord_flip()
      }
      # Whether to plot the legend.
      if (isFALSE(legend)){
        p <- p + Seurat::NoLegend()
      }
      # Whether to remove legend title.
      if (isFALSE(legend.title)){
        p <- p + ggpubr::rremove("legend.title")
      }
    }
    # Return the plot.
    return(p)

}
