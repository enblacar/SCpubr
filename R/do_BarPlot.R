#' Wrapper for computing publication ready bar plots.
#'
#' @param sample  Seurat object.
#' @param var.to.plot  Main variable in the bar plot. Example: seurat_clusters
#' @param group.by  Secondary variable to group the bar plot for. Example: orig.ident
#' @param order.by  Value in var.to.plot to order the items in group.by for.
#' @param labels.order  Vector of labels to explicitly state the order or the bars.
#' @param position  Either "fill" or "stack." Position "fill" will generate a bar plot with one column and the proportions of values for each group inside, while "stack" plots the bars together.
#' @param xlab  Title for the X axis.
#' @param ylab  Title for the Y axis.
#' @param colors.use  Palette of colors to use. It must match the group.by variable in terms of length and names.
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
#'
#' @return A ggplot2 object containing a Bar plot.
#' @export
#'
#' @examples
#' \dontrun{
#' TBD
#' }
do_BarPlot <- function(sample,
                       var.to.plot,
                       group.by = NULL,
                       labels.order = NULL,
                       order.by = NULL,
                       position = "stack",
                       xlab = "",
                       ylab = "",
                       plot.title = "",
                       legend.position = "bottom",
                       legend.title = FALSE,
                       legend.ncol = 1,
                       legend.text.size = 22,
                       legend.title.size = 22,
                       axis.text.size = 22,
                       axis.title.size = 22,
                       plot.title.size = 24,
                       legend.icon.size = 4,
                       legend.byrow = FALSE,
                       colors.use = NULL,
                       horizontal = TRUE){
    # Checks for packages.
    check_suggests(function_name = "do_BarPlot")

    # Define pipe operator internally.
    `%>%` <- purrr::`%>%`

    # If no color scale is provided, generate a custom one.
    if (is.null(colors.use)){
        if (is.null(group.by)){
            # Generate a color palette equal to the number of identities in the seurat object.
            colors.use <- colortools::setColors("#2874A6", length(levels(sample)))
            names(colors.use) <- levels(sample)
        } else {
            # Generate a color palette equal to the number of unique values in group.by variable.
            colors.use <-  colortools::setColors("#2874A6", length(unique(sample[[]][, group.by])))
            names(colors.use) <- unique(sample[[]][, group.by])
        }
    }

    # Bar plot without grouping variables.
    if (is.null(group.by)){
        if (!(is.null(labels.order))){
            factor_levels <- labels.order
        }
        p <- sample@meta.data %>%
            dplyr::select(!!rlang::sym(var.to.plot)) %>%
            dplyr::group_by(!!rlang::sym(var.to.plot)) %>%
            dplyr:: summarise(n = dplyr::n()) %>%
            dplyr::mutate(x_values = as.factor(!!rlang::sym(var.to.plot))) %>%
            ggplot2::ggplot(mapping = ggplot2::aes(x = .data$x_values, y = .data$n, fill = .data$x_values)) +
            ggplot2::geom_bar(position = position, stat="identity", width = 1,
                              colour="black",
                              size = 1) +
            ggpubr::theme_pubr(legend = legend.position) +
            ggplot2::scale_fill_manual(values = colors.use, na.value = "grey75") +
            ggplot2::xlab(xlab) +
            ggplot2::ylab(ylab) +
            ggplot2::ggtitle(plot.title) +
            ggplot2::theme(axis.title.x = ggplot2::element_text(size = axis.title.size, face = "bold"),
                           axis.title.y = ggplot2::element_text(size = axis.title.size, face = "bold"),
                           axis.text = ggplot2::element_text(size = axis.text.size, face = "bold"),
                           legend.text = ggplot2::element_text(size = legend.text.size, face = "bold"),
                           legend.title = ggplot2::element_text(size = legend.title.size, face = "bold"),
                           plot.title = ggplot2::element_text(size = plot.title.size, face = "bold", hjust = 0.5)) +
            ggplot2::guides(color = ggplot2::guide_legend(ncol = legend.ncol,
                                                          byrow = legend.byrow,
                                                          override.aes = list(size = legend.icon.size)))

    # Barplot grouping by the split.by variable.
    } else {

        if (!(is.null(order.by))){
            factor_levels <- sample@meta.data %>%
                dplyr::select(!!rlang::sym(var.to.plot), !!rlang::sym(group.by)) %>%
                dplyr::group_by(!!rlang::sym(group.by), !!rlang::sym(var.to.plot)) %>%
                dplyr::summarise(n = dplyr::n()) %>%
                dplyr::mutate(x_value = !!rlang::sym(group.by)) %>%
                dplyr::filter(.data$x_value == order.by) %>%
                dplyr::mutate(num_cells = {sample@meta.data %>% dplyr::select(!!rlang::sym(var.to.plot)) %>% dplyr::group_by(!!rlang::sym(var.to.plot)) %>% dplyr::summarise(n = dplyr::n()) %>%
                        dplyr::filter(!!rlang::sym(var.to.plot) %in% unique(sample@meta.data[, c(group.by, var.to.plot)][sample@meta.data[, c(group.by, var.to.plot)][, group.by] == order.by, ][, var.to.plot])) %>%
                        dplyr::pull(.data$n)}) %>%
                dplyr::mutate(frac = .data$n/.data$num_cells) %>%
                dplyr::arrange(dplyr::desc(.data$frac)) %>%
                dplyr::pull(!!rlang::sym(group.by))

            total_levels <- unique(sample[[group.by]])

            if (length(factor_levels) != length(total_levels)){
                factor_levels <- c(factor_levels, total_levels[!(total_levels %in% factor_levels)])
            }
            factor_levels <- rev(factor_levels)

        } else {
            data <- sample@meta.data %>%
                dplyr::select(!!rlang::sym(var.to.plot), !!rlang::sym(group.by)) %>%
                dplyr::group_by(!!rlang::sym(group.by), !!rlang::sym(var.to.plot)) %>%
                dplyr::summarise(n = dplyr::n()) %>%
                dplyr::arrange(dplyr::desc(.data$n)) %>%
                dplyr::mutate(x_values = as.factor(!!(rlang::sym(var.to.plot))))

            factor_levels <- rev(sort(unique(data$x_values)))
        }

        if (!is.null(labels.order)){
            factor_levels <- labels.order
        }

        p <- sample@meta.data %>%
            dplyr::select(!!rlang::sym(var.to.plot), !!rlang::sym(group.by)) %>%
            dplyr::group_by(!!rlang::sym(group.by), !!rlang::sym(var.to.plot)) %>%
            dplyr::summarise(n = dplyr::n()) %>%
            dplyr::arrange(dplyr::desc(.data$n)) %>%
            dplyr::mutate(x_values = as.factor(!!(rlang::sym(var.to.plot)))) %>%
            dplyr::mutate(x_values = factor(.data$x_values, levels = factor_levels)) %>%
            ggplot2::ggplot(mapping = ggplot2::aes(x = .data$x_values, y = .data$n, fill = !!rlang::sym(group.by))) +
            ggplot2::geom_bar(position = position, stat="identity", width = 1,
                              colour="black",
                              size = 1) +
            ggpubr::theme_pubr(legend = legend.position) +
            ggplot2::scale_fill_manual(values = colors.use, na.value = "grey75") +
            ggplot2::xlab(xlab) +
            ggplot2::ylab(ylab) +
            ggplot2::ggtitle(plot.title) +
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
    # Whether to flip the axis or not.
    if (horizontal == TRUE){
        p <- p + ggplot2::coord_flip()
    }

    # Whether to remove legend title.
    if (legend.title == FALSE){
        p <- p + ggpubr::rremove("legend.title")
    }
    # Return the plot.
    return(p)

}
