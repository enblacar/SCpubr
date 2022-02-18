#' Wrapper for \link[Seurat]{VlnPlot}.
#'
#'
#' @param sample  Seurat object.
#' @param features Features to represent.
#' @param assay Assay to use. Defauls to the current assay.
#' @param slot Data slot to use. Character. Only one of: counts, data, scale.data. Defaults to "data".
#' @param group.by  Variable you want the cells to be colored for.
#' @param split.by  CURRENTLY NOT WORKING. Split into as many plots as unique values in the variable provided.
#' @param colors.use  Named vector of colors matching the unique values in either the current identities in the seurat object of the unique values in group.by or split.by.
#' @param plot_boxplot Logical. Whether to plot a boxplot inside the violin or not.
#' @param legend  Whether to plot the legend or not.
#' @param legend.position Position of the legend in the plot. Will only work if legend is set to TRUE.
#' @param plot.title  Title to use in the plot.
#' @param individual.titles Titles for each feature if needed. Either NULL or a vector of equal length of features.
#' @param pt.size  Size of points in the VlnPlot.
#' @param xlab  Title for the X axis.
#' @param ylab  Title for the Y axis.
#' @param fontsize Base fontsize for the figure.
#' @param y_cut  Vector with the values in which the Violins should be cut. Only works for one feature.
#' @param legend.ncol  Number of columns in the legend.
#' @param ncol Numeric. Number of columns to arrange multiple plots into.
#' @param rotate_x_labels Logical. Whether to rotate X axis labels to horizontal or not. If multiple features, a vector of logicals of the same length.


#' @return A ggplot2 object containing a Violin Plot.
#' @export
#'
#' @examples
#' \dontrun{
#' TBD
#' }
do_VlnPlot <- function(sample,
                       features,
                       assay = NULL,
                       slot = NULL,
                       group.by = NULL,
                       split.by = NULL,
                       legend = FALSE,
                       colors.use = NULL,
                       pt.size = 0,
                       y_cut = NULL,
                       plot_boxplot = TRUE,
                       legend.position = "bottom",
                       plot.title = NULL,
                       individual.titles = NULL,
                       xlab = NULL,
                       ylab = NULL,
                       fontsize = 16,
                       ncol = NULL,
                       legend.ncol = 3,
                       rotate_x_labels = NULL){
    # Checks for packages.
    check_suggests(function_name = "do_VlnPlot")
    # Check the assay.
    out <- check_and_set_assay(sample = sample, assay = assay)
    sample <- out[["sample"]]
    assay <- out[["assay"]]
    # Check slot.
    slot <- check_and_set_slot(slot = slot)
    # Check logical parameters.
    logical_list <- list("legend" = legend,
                         "plot_boxplot" = plot_boxplot,
                         "rotate_x_labels" = rotate_x_labels)
    check_type(parameters = logical_list, required_type = "logical", test_function = is.logical)
    # Check numeric parameters.
    numeric_list <- list("pt.size" = pt.size,
                         "y_cut" = y_cut,
                         "fontsize" = fontsize,
                         "legend.ncol" = legend.ncol,
                         "ncol" = ncol)
    check_type(parameters = numeric_list, required_type = "numeric", test_function = is.numeric)
    # Check character parameters.
    character_list <- list("legend.position" = legend.position,
                           "features" = features,
                           "group.by" = group.by,
                           "split.by" = split.by,
                           "colors.use" = colors.use,
                           "plot.title" = plot.title,
                           "xlab" = xlab,
                           "ylab" = ylab)
    check_type(parameters = character_list, required_type = "character", test_function = is.character)

    # Check the feature.
    check_feature(sample = sample, features = features)

    # Define fontsize parameters.
    plot.title.fontsize <- fontsize + 2
    axis.text.fontsize <- fontsize
    axis.title.fontsize <- fontsize + 1
    legend.text.fontsize <- fontsize - 5
    legend.title.fontsize <- fontsize - 4

    # Disable split.by.
    if (!(is.null(split.by))){stop("This option is currently not available.", call. = F)}
    # Check for y_cut and only having 1 feature.
    if (!(is.null(y_cut))){
      if(length(features) != length(y_cut)){
      stop('Total number of y_cut values does not match the number of features provided.', call. = F)
      }
    }

    # Check for y_cut and only having 1 feature.
    if (!(is.null(rotate_x_labels))){
      if(length(features) != length(rotate_x_labels)){
        stop('Total number of rotate_x_labels values does not match the number of features provided.', call. = F)
      }
    }

    # If group.by and split.by are NULL.
    if (is.null(group.by) & is.null(split.by)){
      if (is.null(colors.use)){
        colors.use <- generate_color_scale(levels(sample))
      } else {
        check_consistency_colors_and_names(sample = sample, colors = colors.use)
      }
    # If group.by is not NULL but split.by is NULL.
    } else if (!(is.null(group.by)) & is.null(split.by)){
      if (is.null(colors.use)){
        names.use <- sort(unique(sample@meta.data[, group.by]))
        if (is.factor(names.use)){names.use <- levels(names.use)}
        colors.use <- generate_color_scale(names.use)
      } else {
        check_consistency_colors_and_names(sample = sample, colors = colors.use, grouping_variable = group.by)
      }
    # If group by is NULL but split.by is not NULL.
    } else if (!(is.null(split.by)) & is.null(group.by)){
      if (is.null(colors.use)){
        names.use <- sort(unique(sample@meta.data[, split.by]))
        if (is.factor(names.use)){names.use <- levels(names.use)}
        colors.use <- generate_color_scale(names.use)
      } else {
        check_consistency_colors_and_names(sample = sample, colors = colors.use, grouping_variable = split.by)
      }
    } else if (!(is.null(split.by)) & !(is.null(group.by))){stop("Either group.by or split.by has to be NULL.", call. = F)}

    if (!(is.null(individual.titles)) & length(individual.titles) != length(features)){
      stop("The length of individual titles has to be equal to the number of features.", call. = F)
    }
    counter <- 0
    list.plots <- list()
    for (feature in features) {
      counter <- counter + 1
      if (!is.null(y_cut)){y_cut_select <- y_cut[counter]}
      if (!is.null(rotate_x_labels)){x_label_select <- rotate_x_labels[counter]}
      # Check the value for y_cut_select is on range.
      if (!(is.null(y_cut))){
        if (!(is.na(y_cut_select))){
          check_limits(sample = sample, feature = feature, value_name = "y_cut", value = y_cut_select)
        }
      }
      plot <- Seurat::VlnPlot(sample,
                              features = feature,
                              cols = colors.use,
                              group.by = group.by,
                              split.by = split.by,
                              pt.size = pt.size) &
        ggpubr::theme_pubr(legend = legend.position) &
        ggpubr::rremove("x.title") &
        ggplot2::theme(axis.text.x = ggplot2::element_text(size = axis.text.fontsize, angle = 90, vjust = 0.5, hjust = 1, face = "bold"),
                       axis.text.y = ggplot2::element_text(size = axis.text.fontsize, face = "bold"),
                       axis.title = ggplot2::element_text(face = "bold", size = axis.title.fontsize),
                       legend.text = ggplot2::element_text(size = legend.text.fontsize, hjust = 0, face = "bold"),
                       legend.title = ggplot2::element_text(size = legend.title.fontsize, face = "bold"),
                       plot.title = ggplot2::element_text(size = plot.title.fontsize, face = "bold", hjust = 0.5)) &
        ggplot2::guides(fill = ggplot2::guide_legend(ncol = legend.ncol))

      if (!is.null(xlab)){
        plot <- plot & ggplot2::xlab(xlab)
      }
      if (!is.null(ylab)){
        plot <- plot & ggplot2::ylab(ylab)
      }
      if (plot_boxplot == TRUE){
        plot <- plot &
          ggplot2::geom_boxplot(width = 0.2, fill = "white", outlier.colour = NA)
      }

      if (legend == FALSE){
        plot <- plot & Seurat::NoLegend()
      }

      if (!(is.null(rotate_x_labels))){
        if (isTRUE(x_label_select)){
          plot <- plot & ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 0))
        }
      }

      if (!(is.null(y_cut))){
          if (!(is.na(y_cut_select))){
            plot <- plot &
              ggplot2::geom_hline(yintercept = y_cut_select, linetype = "dashed", colour = "black", size = 1, alpha = 0.5)
          }
      }

      if (!(is.null(individual.titles))){
        if (!(is.na(individual.titles[counter]))){
          plot <- plot + ggplot2::ggtitle(individual.titles[counter])
        }
      }
      list.plots[[feature]] <- plot
    }
    plot <- patchwork::wrap_plots(list.plots, ncol = ncol)

    if (!(is.null(plot.title))){
      plot <- plot + patchwork::plot_annotation(title = plot.title,
                                                theme = ggplot2::theme(plot.title = ggplot2::element_text(size = plot.title.fontsize + 2,
                                                                                                          face = "bold",
                                                                                                          hjust = 0.5)))
    }
    return(plot)
}
