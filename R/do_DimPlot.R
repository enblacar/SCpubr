#' Wrapper for \link[Seurat]{DimPlot}.
#'
#'
#' @param sample Seurat object.
#' @param reduction Reduction to use. Can be the canonical ones such as "umap", "pca", or any custom ones, such as "diffusion". If you are unsure about which reductions you have, use `Seurat::Reductions(sample)`. Defaults to "umap" if present or to the last computed reduction if the argument is not provided.
#' @param group.by Variable you want the cells to be colored for.
#' @param split.by Split into as many plots as unique values in the variable provided.
#' @param colors.use Vector of named HEX values to color the cells. It has to match the number of unique values in either `Seurat::Idents(sample)` or the group.by or split.by variable. For split.by, a single color can be provided and each panel will be colored by it.
#' @param label Whether to plot the cluster labels in the UMAP. The cluster labels will have the same color as the cluster colors.
#' @param cells.highlight Vector of cells for which the DimPlot should focus into. The rest of the cells will be grayed out.
#' @param shuffle Whether to shuffle the cells or not, so that they are not plotted cluster-wise. Recommended.
#' @param pt.size Point size of the cells.
#' @param sizes.highlight Point size of highlighted cells using cells.highlight parameter.
#' @param legend Whether to plot the legend or not.
#' @param legend.title Logical stating whether the legend title is shown or not.
#' @param legend.ncol Number of columns in the legend.
#' @param fontsize Base fontsize of the figure.
#' @param legend.icon.size Size of the icons in legend.
#' @param legend.position Position of the legend in the plot. Will only work if legend is set to TRUE.
#' @param legend.byrow Logical stating whether the legend is filled by row or not.
#' @param plot.title Title to use in the plot.
#' @param ncol Number of columns used in the arrangement of the output plot using "split.by" parameter.
#' @param dims Vector of 2 numerics indicating the dimensions to plot out of the selected reduction. Defaults to c(1, 2) if not specified.
#' @param repel Whether to repel the labels if label is set to TRUE.
#' @param raster Whether to raster the resulting plot. This is recommendable if plotting a lot of cells.
#' @param label.color HEX code for the color of the text in the labels if label is set to TRUE.
#'
#' @return  A ggplot2 object containing a DimPlot.
#' @export
#'
#' @example man/examples/examples_do_DimPlot.R
do_DimPlot <- function(sample,
                       reduction = NULL,
                       group.by = NULL,
                       split.by = NULL,
                       colors.use = NULL,
                       shuffle = TRUE,
                       pt.size = 0.5,
                       label = FALSE,
                       label.color = "black",
                       repel = TRUE,
                       cells.highlight = NULL,
                       sizes.highlight = 0.5,
                       legend = TRUE,
                       ncol = NULL,
                       plot.title = NULL,
                       legend.title = FALSE,
                       legend.position = "right",
                       legend.ncol = 1,
                       legend.icon.size = 4,
                       legend.byrow = FALSE,
                       raster = FALSE,
                       dims = c(1, 2),
                       fontsize = 14){
    # Checks for packages.
    check_suggests(function_name = "do_DimPlot")
    # Check the reduction.
    reduction <- check_and_set_reduction(sample = sample, reduction = reduction)
    # Check the dimensions.
    dimensions <- check_and_set_dimensions(sample = sample, reduction = reduction, dims = dims)
    # Check logical parameters.
    logical_list <- list("label" = label,
                         "repel" = repel,
                         "shuffle" = shuffle,
                         "legend" = legend,
                         "legend.title" = legend.title,
                         "legend.byrow" = legend.byrow,
                         "raster" = raster)
    check_type(parameters = logical_list, required_type = "logical", test_function = is.logical)
    # Check numeric parameters.
    numeric_list <- list("pt.size" = pt.size,
                         "sizes.highlight" = sizes.highlight,
                         "legend.ncol" = legend.ncol,
                         "fontsize" = fontsize,
                         "legend.icon.size" = legend.icon.size,
                         "ncol" = ncol)
    check_type(parameters = numeric_list, required_type = "numeric", test_function = is.numeric)
    # Check character parameters.
    character_list <- list("legend.position" = legend.position,
                           "plot.title" = plot.title,
                           "cells.highlight" = cells.highlight)
    check_type(parameters = character_list, required_type = "character", test_function = is.character)

    # Checks to ensure proper function.
    if (!(is.null(split.by)) & !(is.null(group.by))){stop("Either group.by or split.by has to be NULL.")}
    if (!(is.null(cells.highlight)) & !(is.null(group.by))){stop("Either group.by or cells.highlight has to be NULL.")}
    if (!(is.null(cells.highlight)) & !(is.null(split.by))){stop("Either split.by or cells.highlight has to be NULL.")}

    # Check for label.color.
    check_colors(label.color, parameter_name = "label.color")

    # Automatically generate colors.
    if (is.null(colors.use)){
      colors.use <- {
        # When everything is NULL.
        if (is.null(group.by) & is.null(split.by) & is.null(cells.highlight)){
          generate_color_scale(levels(sample))
        # When everything is NULL but group.by.
        } else if (!(is.null(group.by)) & is.null(split.by) & is.null(cells.highlight)){
          data.use <- sample[[]][, group.by, drop = F]
          names.use <- if (is.factor(data.use[, 1])){levels(data.use[, 1])} else {sort(unique(data.use[, 1]))}
          generate_color_scale(names.use)
        # When everything is NULL but split.by.
        } else if (is.null(group.by) & !(is.null(split.by)) & is.null(cells.highlight)){
          data.use <- sample[[]][, split.by, drop = F]
          names.use <- if (is.factor(data.use[, 1])){levels(data.use[, 1])} else {sort(unique(data.use[, 1]))}
          generate_color_scale(names.use)
        } else if (is.null(group.by) & is.null(split.by) & !(is.null(cells.highlight))){
          colors.use <- "#0A305F"
        }
      }
    # But, if the user has provided some.
    } else{
      # Check that the provided values are valid color representations.
      check_colors(colors.use, parameter_name = "colors.use")
      # When everything is NULL.
      if (is.null(group.by) & is.null(split.by) & is.null(cells.highlight)){
        check_consistency_colors_and_names(sample = sample, colors = colors.use)
        # When everything is NULL but group.by.
      } else if (!(is.null(group.by)) & is.null(split.by) & is.null(cells.highlight)){
        check_consistency_colors_and_names(sample = sample, colors = colors.use, grouping_variable = group.by)
        # When everything is NULL but split.by.
      } else if (is.null(group.by) & !(is.null(split.by)) & is.null(cells.highlight)){
        if (length(colors.use) != 1){
          check_consistency_colors_and_names(sample = sample, colors = colors.use, grouping_variable = split.by)
        }
      }
    }

    # Define fontsize parameters.
    plot.title.fontsize <- fontsize + 2
    axis.text.fontsize <- fontsize
    axis.title.fontsize <- fontsize + 1
    legend.text.fontsize <- fontsize - 2
    legend.title.fontsize <- fontsize - 2

    # If the UMAP does not need to be split in multiple panes (default case).
    if (is.null(cells.highlight) & is.null(split.by)){
        p <- Seurat::DimPlot(sample,
                                  reduction = reduction,
                                  label = label,
                                  dims = dims,
                                  repel = ifelse(is.null(label) == TRUE, NULL, TRUE),
                                  label.box = ifelse(is.null(label) == TRUE, NULL, TRUE),
                                  label.color = ifelse(is.null(label) == TRUE, NULL, "black"),
                                  shuffle = TRUE,
                                  pt.size = pt.size,
                                  group.by = group.by,
                                  cols = colors.use,
                                  raster = raster,
                                  ncol = ncol) +
            ggpubr::theme_pubr(legend = legend.position) +
            ggplot2::theme(plot.title = ggplot2::element_text(size = plot.title.fontsize, face = "bold", hjust = 0.5),
                           legend.text = ggplot2::element_text(size = legend.text.fontsize, face = "bold"),
                           legend.title = ggplot2::element_text(size = legend.title.fontsize, face = "bold")) +
            ggplot2::guides(color = ggplot2::guide_legend(ncol = legend.ncol,
                                                          byrow = legend.byrow,
                                                          override.aes = list(size = legend.icon.size)))

    }
    # If the UMAP has to be split in multiple panes.
    else if (is.null(cells.highlight) & !(is.null(split.by))){
        # If the user provided multiple highlighting colors.
        multiple_colors <- ifelse(length(colors.use) > 1, TRUE, FALSE)
        # List to store each individual plots.
        list.plots <- list()
        # Recover all metadata.
        data <- sample[[]]
        # Retrieve the metadata column belonging to the split.by parameter.
        data.use <- data[, split.by, drop = F]
        # Retrieve the plotting order, keep factor levels if the column is a factor.
        plot_order <- if (is.factor(data.use[, 1])){levels(data.use[, 1])} else {sort(unique(data.use[, 1]))}
        # Iterate over each unique value in split.by parameter.
        for (iteration in plot_order){
            # Retrieve the cells that do belong to the iteration's split.by value.
            cells.highlight <- rownames(data.use)[which(data.use == iteration)]
            p <- Seurat::DimPlot(sample,
                                      reduction = reduction,
                                      dims = dims,
                                      cells.highlight = cells.highlight,
                                      sizes.highlight = sizes.highlight,
                                      pt.size = pt.size,
                                      raster = raster,
                                      ncol = ncol) +
                ggplot2::ggtitle(iteration) +
                ggpubr::theme_pubr(legend = legend.position) +
                ggplot2::scale_color_manual(labels = c("Unselected", "Selected"),
                                            values = c("grey75", ifelse(multiple_colors == TRUE, colors.use[[iteration]], colors.use)))  +
                ggplot2::theme(plot.title = ggplot2::element_text(size = plot.title.fontsize, face = "bold", hjust = 0.5),
                               legend.text = ggplot2::element_text(size = legend.text.fontsize, face = "bold"),
                               legend.title = ggplot2::element_text(size = legend.title.fontsize, face = "bold")) +
                ggplot2::guides(color = ggplot2::guide_legend(ncol = legend.ncol,
                                                              byrow = legend.byrow,
                                                              override.aes = list(size = legend.icon.size)))
            list.plots[[iteration]] <- p
        }
        # Assemble individual plots as a patch.
        p <- patchwork::wrap_plots(list.plots, ncol = ncol)
    }


    # If the user wants to highlight some of the cells.
    else if (!(is.null(cells.highlight))) {
        p <- Seurat::DimPlot(sample,
                                  reduction = reduction,
                                  cells.highlight = cells.highlight,
                                  sizes.highlight = sizes.highlight,
                                  dims = dims,
                                  pt.size = pt.size,
                                  raster = raster,
                                  ncol = ncol) +
            ggpubr::theme_pubr(legend = legend.position) +
            ggplot2::scale_color_manual(labels = c("Unselected", "Selected"),
                                        values = c("grey", colors.use))  +
            ggplot2::theme(plot.title = ggplot2::element_text(size = plot.title.fontsize, face = "bold", hjust = 0.5),
                           legend.text = ggplot2::element_text(size = legend.text.fontsize, face = "bold"),
                           legend.title = ggplot2::element_text(size = legend.title.fontsize, face = "bold")) +
            ggplot2::guides(color = ggplot2::guide_legend(ncol = legend.ncol,
                                                          byrow = legend.byrow,
                                                          override.aes = list(size = legend.icon.size)))
    }


    # General additions to all kind of plots.
    if (!is.null(plot.title)){
      if (!(is.null(split.by))){
        p <- p + patchwork::plot_annotation(title = plot.title,
                                            theme = ggplot2::theme(plot.title = ggplot2::element_text(size = plot.title.fontsize + 1,
                                                                                                      face = "bold",
                                                                                                      hjust = 0.5)))
      } else {
        p <- p + ggplot2::ggtitle(plot.title)
      }
    }
    # Legend treatment.
    if (legend == FALSE){
        p <- p & ggpubr::rremove("legend.title") & ggpubr::rremove("legend")
    } else if (legend == TRUE) {
        p <- p &  ggpubr::rremove("legend.title")
    }

    if (legend.title == FALSE) {
        p <- p & ggpubr::rremove("legend.title")
    }
    # For embeddings that are umap of tsne, we remove all axes..
    if (reduction %in% c("umap", "tsne")){
        # if dims is first and then second.
        if (sum(dims == c(1, 2)) == 2){
          p <- p & Seurat::NoAxes()
        } else {
          labels <- colnames(sample@reductions[[reduction]][[]])[dims]
          p <- p &
            ggpubr::rremove("axis") &
            ggpubr::rremove("axis.text") &
            ggpubr::rremove("ticks") &
            ggplot2::theme(axis.title.x = ggplot2::element_text(size = axis.title.fontsize, face = "bold"),
                           axis.title.y = ggplot2::element_text(size = axis.title.fontsize, face = "bold", angle = 90)) &
            ggplot2::xlab(labels[1]) & ggplot2::ylab(labels[2])
        }
    # For diffusion maps, we do want to keep at least the axis titles so that we know which DC are we plotting.
    } else {
        labels <- colnames(sample@reductions[[reduction]][[]])[dims]
        p <- p &
            ggpubr::rremove("axis") &
            ggpubr::rremove("axis.text") &
            ggpubr::rremove("ticks") &
            ggplot2::theme(axis.title.x = ggplot2::element_text(size = axis.title.fontsize, face = "bold"),
                           axis.title.y = ggplot2::element_text(size = axis.title.fontsize, face = "bold", angle = 90)) &
            ggplot2::xlab(labels[1]) & ggplot2::ylab(labels[2])
    }
    # Label treatment.
    if (label == TRUE && is.null(cells.highlight)){
        p$layers[[2]]$aes_params$fontface <- "bold"
    }
    # Return the final plot.
    return(p)
}

