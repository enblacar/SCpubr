#' Wrapper for \link[Seurat]{DimPlot}.
#'
#'
#' @param sample Seurat object.
#' @param reduction Reduction to use. Can be the canonical ones such as "umap", "pca", or any custom ones, such as "diffusion". If you are unsure about which reductions you have, use `Seurat::Reductions(sample)`. Defaults to "umap" if present or to the last computed reduction if the argument is not provided.
#' @param group.by Variable you want the cells to be colored for.
#' @param split.by Split into as many plots as unique values in the variable provided.
#' @param colors.use Vector of named HEX values to color the cells. It has to match the number of unique values in either `Seurat::Idents(sample)` or the group.by variable.
#' @param colors.split Vector of named HEX values to color the cells in a split DimPlot. It has to match the unique values in split.by parameter.
#' @param label Whether to plot the cluster labels in the UMAP. The cluster labels will have the same color as the cluster colors.
#' @param cells.highlight Vector of cells for which the DimPlot should focus into. The rest of the cells will be grayed out.
#' @param colors.highlight HEX color code to use with the highlighted cells.
#' @param shuffle Whether to shuffle the cells or not, so that they are not plotted cluster-wise. Recommended.
#' @param pt.size Point size of the cells.
#' @param sizes.highlight Point size of highlighted cells using cells.highlight parameter.
#' @param legend Whether to plot the legend or not.
#' @param legend.title Logical stating whether the legend title is shown or not.
#' @param legend.ncol Number of columns in the legend.
#' @param legend.text.size Font size of the legend labels.
#' @param legend.title.size Font size of the legend title.
#' @param legend.icon.size Size of the icons in legend.
#' @param legend.position Position of the legend in the plot. Will only work if legend is set to TRUE.
#' @param legend.byrow Logical stating whether the legend is filled by row or not.
#' @param plot.title Title to use in the plot.
#' @param ncol Number of columns used in the arrangement of the output plot using "split.by" parameter.
#' @param dims Vector of 2 numerics indicating the dimensions to plot out of the selected reduction. Defaults to c(1, 2) if not specified.
#' @param repel Whether to repel the labels if label is set to TRUE.
#' @param raster Whether to raster the resulting plot. This is recommendable if plotting a lot of cells.
#' @param label.color HEX code for the color of the text in the labels if label is set to TRUE.
#' @param ... Other parameters for \link[Seurat]{DimPlot}.
#'
#' @return  A ggplot2 object containing a DimPlot.
#' @export
#'
#' @examples
#' \dontrun{
#' TBD
#' }
do_DimPlot <- function(sample,
                       reduction = NULL,
                       label = FALSE,
                       label.color = "black",
                       repel = TRUE,
                       shuffle = TRUE,
                       pt.size = 0.5,
                       sizes.highlight = 0.5,
                       group.by = NULL,
                       split.by = NULL,
                       colors.split = "#0A305F",
                       cells.highlight = NULL,
                       colors.highlight = "#0A305F",
                       legend = TRUE,
                       legend.title = FALSE,
                       legend.position = "bottom",
                       legend.ncol = 4,
                       legend.text.size = 14,
                       legend.title.size = 14,
                       legend.icon.size = 4,
                       legend.byrow = FALSE,
                       colors.use = NULL,
                       plot.title = "",
                       ncol = NULL,
                       raster = FALSE,
                       dims = c(1, 2),
                       ...){
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
                         "legend.text.size" = legend.text.size,
                         "legend.title.size" = legend.title.size,
                         "legend.icon.size" = legend.icon.size,
                         "ncol" = ncol)
    check_type(parameters = numeric_list, required_type = "numeric", test_function = is.numeric)
    # Check character parameters.
    character_list <- list("legend.position" = legend.position,
                           "plot.title" = plot.title,
                           "cells.highlight" = cells.highlight)
    check_type(parameters = character_list, required_type = "character", test_function = is.character)

    # Checks to ensure proper function.
    # Check whether the names of colors.use match the unique values in group.by or whether the number of colors is lower to the number of unique values in group.by.
    if (!(is.null(group.by)) & !(is.null(colors.use))){check_consistency_colors_and_names(sample = sample, colors = colors.use, groping_variable = group.by)}
    # Check whether the names of colors.split match the unique values in split.by or whether the number of colors is lower to the number of unique values in split.by.
    if (!(is.null(split.by)) & colors.split != "#0A305F" & length(colors.split) > 1){check_consistency_colors_and_names(sample = sample, colors = colors.split, groping_variable = split.by)}
    # Check for colors.highlight.
    if (colors.highlight != "#0A305F" & !(is.null(colors.highlight))){check_colors(colors.highlight, parameter_name = "colors.highlight")}
    # Check for colors.use.
    if (!is.null(colors.use)){check_colors(colors.use, parameter_name = "colors.use")}
    # Check for colors.split.
    if (!(is.null(colors.split)) & colors.split != "#0A305F" & colors.split != TRUE){check_colors(colors.split, parameter_name = "colors.split")}
    # Check for label.color.
    check_colors(label.color, parameter_name = "label.color")

    # Automatically generate color palettes when the user has not defined one.
    # If the user does not want to highlight any cells (Regular case.).
    if (is.null(cells.highlight) & is.null(colors.use)){
        # If no special grouping is set up, DimPlot defaults back to Seurat::Idents(sample) or levels(sample).
        if (is.null(group.by)){colors.use <- generate_color_scale(levels(sample))} else {colors.use <- generate_color_scale(unique(sample[[]][, group.by]))}
    }

    # If the user wants different coloring but has not provided a vector of colors, then resort to the default coloring.
    if (colors.split == TRUE){
      # Generate a vector of colors equal to the number of identities in the sample.
      data.use <- sample[[]][, split.by, drop = F]
      names.use <- if (is.factor(data.use[, 1])){levels(data.use[, 1])} else {sort(unique(data.use[, 1]))}
      colors.split <- generate_color_scale(names.use)
    }


    # If the UMAP does not need to be split in multiple panes (default case).
    if (is.null(split.by)){
        p.umap <- Seurat::DimPlot(sample,
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
                                  ncol = ncol,
                                  ...) +
            ggpubr::theme_pubr(legend = legend.position) +
            ggplot2::ggtitle(plot.title) +
            ggplot2::theme(plot.title = ggplot2::element_text(face = "bold", hjust = 0.5),
                           legend.text = ggplot2::element_text(size = legend.text.size, face = "bold"),
                           legend.title = ggplot2::element_text(size = legend.title.size, face = "bold")) +
            ggplot2::guides(color = ggplot2::guide_legend(ncol = legend.ncol,
                                                          byrow = legend.byrow,
                                                          override.aes = list(size = legend.icon.size)))
        # If a custom color scale is provided, this line's aim is to reorder the legend labels into alphabetical order.
        if (!(is.null(colors.use))){
            p.umap <- p.umap + ggplot2::scale_color_manual(values = colors.use, breaks = sort(names(colors.use)))
        }

    }
    # If the UMAP has to be split in multiple panes.
    else if (!(is.null(split.by))){
        # If the user provided multiple highlighting colors.
        multiple_colors <- ifelse(length(colors.split) > 1, TRUE, FALSE)
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
            p.umap <- Seurat::DimPlot(sample,
                                      reduction = reduction,
                                      dims = dims,
                                      cells.highlight = cells.highlight,
                                      sizes.highlight = sizes.highlight,
                                      pt.size = pt.size,
                                      raster = raster,
                                      ncol = ncol,
                                      ...) +
                ggplot2::ggtitle(iteration) +
                ggpubr::theme_pubr(legend = legend.position) +
                ggplot2::scale_color_manual(labels = c("Unselected", "Selected"),
                                            values = c("grey75", ifelse(multiple_colors == TRUE, colors.split[[iteration]], colors.split)))  +
                ggplot2::theme(plot.title = ggplot2::element_text(face = "bold", hjust = 0.5),
                               legend.text = ggplot2::element_text(size = legend.text.size, face = "bold"),
                               legend.title = ggplot2::element_text(size = legend.title.size, face = "bold")) +
                ggplot2::guides(color = ggplot2::guide_legend(ncol = legend.ncol,
                                                              byrow = legend.byrow,
                                                              override.aes = list(size = legend.icon.size)))
            list.plots[[iteration]] <- p.umap
        }
        # Assemble individual plots as a patch.
        p.umap <- patchwork::wrap_plots(list.plots, ncol = ncol)
    }


    # If the user wants to highlight some of the cells.
    else if (!(is.null(cells.highlight))) {
        p.umap <- Seurat::DimPlot(sample,
                                  reduction = reduction,
                                  cells.highlight = cells.highlight,
                                  sizes.highlight = sizes.highlight,
                                  dims = dims,
                                  pt.size = pt.size,
                                  raster = raster,
                                  ncol = ncol,
                                  ...) +
            ggplot2::ggtitle(plot.title) +
            ggpubr::theme_pubr(legend = legend.position) +
            ggplot2::scale_color_manual(labels = c("Unselected", "Selected"),
                                        values = c("grey", colors.highlight))  +
            ggplot2::theme(plot.title = ggplot2::element_text(face = "bold", hjust = 0.5),
                           legend.text = ggplot2::element_text(size = legend.text.size, face = "bold"),
                           legend.title = ggplot2::element_text(size = legend.title.size, face = "bold")) +
            ggplot2::guides(color = ggplot2::guide_legend(ncol = legend.ncol,
                                                          byrow = legend.byrow,
                                                          override.aes = list(size = legend.icon.size)))
    }


    # General additions to all kind of plots.
    # Legend treatment.
    if (legend == FALSE){
        p.umap <- p.umap & ggpubr::rremove("legend.title") & ggpubr::rremove("legend")
    } else if (legend == TRUE) {
        p.umap <- p.umap &  ggpubr::rremove("legend.title")
    }

    if (legend.title == FALSE) {
        p.umap <- p.umap & ggpubr::rremove("legend.title")
    }
    # For embeddings that are not diffusion maps, we remove all axes..
    if (reduction != "diffusion"){
        p.umap <- p.umap & Seurat::NoAxes()
    # For diffusion maps, we do want to keep at least the axis titles so that we know which DC are we plotting.
    } else {
        p.umap <- p.umap &
            ggpubr::rremove("axis") &
            ggpubr::rremove("axis.text") &
            ggpubr::rremove("ticks") &
            ggplot2::theme(axis.title.x = ggplot2::element_text(size = 14, face = "bold"),
                           axis.title.y = ggplot2::element_text(size = 14, face = "bold")) &
            ggplot2::xlab(paste0("DC_", dims[1])) & ggplot2::ylab(paste0("DC_", dims[2]))
    }
    # Label treatment.
    if (label == TRUE && is.null(cells.highlight)){
        p.umap$layers[[2]]$aes_params$fontface <- "bold"
    }
    # Return the final plot.
    return(p.umap)
}

