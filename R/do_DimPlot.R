#' Wrapper for \link[Seurat]{DimPlot}.
#'
#'
#' @param sample Seurat object.
#' @param reduction Reduction to use. Can be the canonical ones such as "umap", "pca", or any custom ones, such as "diffusion". If you are unsure about which reductions you have, use `Seurat::Reductions(sample)`. Defaults to "umap" if present or to the last computed reduction if the argument is not provided.
#' @param group.by Variable you want the cells to be colored for.
#' @param split.by Split into as many plots as unique values in the variable provided.
#' @param colors.use Vector of named HEX values to color the cells. It has to match the number of unique values in either `Seurat::Idents(sample)` or the group.by or split.by variable. For split.by, a single color can be provided and each panel will be colored by it.
#' @param label Whether to plot the cluster labels in the UMAP. The cluster labels will have the same color as the cluster colors.
#' @param cells.highlight,idents.highlight Vector of cells/identities for which the DimPlot should focus into. The rest of the cells will be grayed out. Both parameters can be used at the same time.
#' @param idents.keep Vector of identities to keep. This will effectively set the rest of the cells that do not match the identities provided to NA, therefore coloring them according to na.value parameter.
#' @param shuffle Whether to shuffle the cells or not, so that they are not plotted cluster-wise. Recommended.
#' @param order Vector of identities to be plotted. Either one with all identities or just some, which will be plotted last.
#' @param pt.size Point size of the cells.
#' @param sizes.highlight Point size of highlighted cells using cells.highlight parameter.
#' @param legend Whether to plot the legend or not.
#' @param legend.title Logical stating whether the legend title is shown or not.
#' @param legend.title.position Character stating where to place the title of the legend.
#' @param legend.ncol,legend.nrow Number of columns/rows in the legend.
#' @param fontsize Base fontsize of the figure.
#' @param legend.icon.size Size of the icons in legend.
#' @param legend.position Position of the legend in the plot. Will only work if legend is set to TRUE.
#' @param legend.byrow Logical stating whether the legend is filled by row or not.
#' @param plot.title Title to use in the plot.
#' @param ncol Number of columns used in the arrangement of the output plot using "split.by" parameter.
#' @param dims Vector of 2 numerics indicating the dimensions to plot out of the selected reduction. Defaults to c(1, 2) if not specified.
#' @param repel Whether to repel the labels if label is set to TRUE.
#' @param raster Whether to raster the resulting plot. This is recommendable if plotting a lot of cells.
#' @param raster.dpi Numeric. Pixel resolution for rasterized plots. Defaults to 512, as per default `Seurat::DimPlot()` behavior.
#' @param label.color HEX code for the color of the text in the labels if label is set to TRUE.
#' @param na.value Color value for NA.
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
                       order = NULL,
                       pt.size = 0.5,
                       label = FALSE,
                       label.color = "black",
                       repel = TRUE,
                       cells.highlight = NULL,
                       idents.highlight = NULL,
                       idents.keep = NULL,
                       sizes.highlight = 0.5,
                       legend = TRUE,
                       ncol = NULL,
                       plot.title = NULL,
                       legend.title = FALSE,
                       legend.position = "right",
                       legend.title.position = "top",
                       legend.ncol = NULL,
                       legend.nrow = NULL,
                       legend.icon.size = 4,
                       legend.byrow = FALSE,
                       raster = FALSE,
                       raster.dpi = 2048,
                       dims = c(1, 2),
                       fontsize = 14,
                       na.value = "grey75"){
    # Checks for packages.
    check_suggests(function_name = "do_DimPlot")
    # Check if the sample provided is a Seurat object.
    check_Seurat(sample = sample)
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
                         "legend.nrow" = legend.nrow,
                         "fontsize" = fontsize,
                         "legend.icon.size" = legend.icon.size,
                         "ncol" = ncol,
                         "raster.dpi" = raster.dpi)
    check_type(parameters = numeric_list, required_type = "numeric", test_function = is.numeric)
    # Check character parameters.
    character_list <- list("legend.position" = legend.position,
                           "plot.title" = plot.title,
                           "cells.highlight" = cells.highlight,
                           "idents.keep" = idents.keep,
                           "order" = order,
                           "na.value" = na.value,
                           "idents.highlight" = idents.highlight,
                           "legend.title.position" = legend.title.position)
    check_type(parameters = character_list, required_type = "character", test_function = is.character)

    # Checks to ensure proper function.
    group_by_and_split_by_used <- !(is.null(split.by)) & !(is.null(group.by))
    group_by_and_highlighting_cells <- (!(is.null(cells.highlight)) | !(is.null(idents.highlight))) & !(is.null(group.by))
    split_by_and_highlighting_cells <- (!(is.null(cells.highlight)) | !(is.null(idents.highlight))) & !(is.null(split.by))
    order_and_shuffle_used <- !(is.null(order)) & isTRUE(shuffle)
    if (group_by_and_split_by_used){stop("Either group.by or split.by has to be NULL.", call. = F)}
    if (group_by_and_highlighting_cells){stop("Either group.by or cells.highlight has to be NULL.", call. = F)}
    if (split_by_and_highlighting_cells){stop("Either split.by or cells.highlight has to be NULL.", call. = F)}
    if (order_and_shuffle_used){warning("Setting up a custom order while 'shuffle = TRUE' might result in unexpected behaviours.\nPlease consider using it alongside 'shuffle = FALSE'.", call. = FALSE)}
    # Check for label.color.
    check_colors(label.color, parameter_name = "label.color")
    check_colors(na.value, parameter_name = "na.value")
    if (length(na.value) != 1){stop("Please provide only one color to na.value.", call. = FALSE)}

    # Check for raster and pt.size.
    if (isTRUE(raster) & pt.size < 1){
      warning("Setting raster = TRUE and pt.size < 1 will result in the cells being ploted as a cross. This behaviour can not be modified, but setting pt.size to 1 or higher solves it. For DimPlots, optimized values would be pt.size = 3 and raster.dpi = 2048.", call. = F)
    }
    # Automatically generate colors.
    # If the user has provided some colors.
    if (is.null(colors.use)){
      colors.use <- {
        # Default parameters.
        default_parameters <- is.null(group.by) & is.null(split.by) & is.null(cells.highlight) & is.null(idents.highlight)
        # Group.by was used.
        group_by_is_used <- !(is.null(group.by)) & is.null(split.by) & is.null(cells.highlight) & is.null(idents.highlight)
        # Split.by was used.
        split_by_is_used <- is.null(group.by) & !(is.null(split.by)) & is.null(cells.highlight) & is.null(idents.highlight)
        # Cells.highlight or idents.highlight was used.
        highlighting_cells <- is.null(group.by) & is.null(split.by) & (!(is.null(cells.highlight)) | !(is.null(idents.highlight)))
        if (default_parameters){
          generate_color_scale(levels(sample))
        } else if (group_by_is_used){
          data.use <- sample[[]][, group.by, drop = F]
          names.use <- if (is.factor(data.use[, 1])){levels(data.use[, 1])} else {sort(unique(data.use[, 1]))}
          generate_color_scale(names.use)
        } else if (split_by_is_used){
          data.use <- sample[[]][, split.by, drop = F]
          names.use <- if (is.factor(data.use[, 1])){levels(data.use[, 1])} else {sort(unique(data.use[, 1]))}
          generate_color_scale(names.use)
        } else if (highlighting_cells){
          colors.use <- "#0A305F"
        }
      }
    # But, if the user has provided some.
    } else {
      # Check that the provided values are valid color representations.
      check_colors(colors.use, parameter_name = "colors.use")
      # All set to null.
      default_parameters <- is.null(group.by) & is.null(split.by) & is.null(cells.highlight) & is.null(idents.highlight)
      # Group.by was used.
      group_by_is_used <- !(is.null(group.by)) & is.null(split.by) & is.null(cells.highlight) & is.null(idents.highlight)
      # Split.by was used.
      split_by_is_used <- is.null(group.by) & !(is.null(split.by)) & is.null(cells.highlight) & is.null(idents.highlight)
      # When either cells.highlight or idents.highlight was used.
      highlighting_cells <- is.null(group.by) & is.null(split.by) & (!(is.null(cells.highlight)) | !(is.null(idents.highlight)))
      # When everything is NULL.
      if (default_parameters){
        colors.use <- check_consistency_colors_and_names(sample = sample, colors = colors.use)
        # When everything is NULL but group.by.
      } else if (group_by_is_used){
        colors.use <- check_consistency_colors_and_names(sample = sample, colors = colors.use, grouping_variable = group.by)
        # When everything is NULL but split.by.
      } else if (split_by_is_used){
        if (length(colors.use) != 1){
          colors.use <- check_consistency_colors_and_names(sample = sample, colors = colors.use, grouping_variable = split.by)
        }
      } else if (highlighting_cells){
        if (length(colors.use) > 1){
          stop("Provide only one color if cells.highlight or idents.highlight is used.", call. = F)
        }
      }
    }

    # Set cells to NA according to idents.keep.
    # If the user does not want to highlight cells or split by identities but wants to remove some identities.
    idents_keep_used <- is.null(cells.highlight) & is.null(idents.highlight) & !(is.null(idents.keep))
    if (idents_keep_used){
      # If both group.by and split.by are not used.
      group_by_and_split_by_are_null <- is.null(group.by) & is.null(split.by)
      # If group.by is used.
      group_by_is_used <- !(is.null(group.by)) & is.null(split.by)
      # If
      split_by_is_used <- is.null(group.by) & !(is.null(split.by))
      if (group_by_and_split_by_are_null){
        # Check that idents.keep matches the values.
        if (isFALSE(length(idents.keep) == sum(idents.keep %in% levels(sample)))){
          stop("All the values in idents.keep must be in levels(sample).", call. = F)
        }
        Seurat::Idents(sample)[!(Seurat::Idents(sample) %in% idents.keep)] <- NA
        # Generate the new color scale
        colors.use <- check_consistency_colors_and_names(sample = sample, colors = colors.use)
      # If not, check with the values in group.by.
      } else if (group_by_is_used) {
        # Check that idents.keep matches the values.
        if (isFALSE(length(idents.keep) == sum(idents.keep %in% unique(sample@meta.data[, group.by])))){
          stop("All the values in idents.keep must be in the group.by variable provided.", call. = F)
        }
        sample@meta.data[, group.by][!(sample@meta.data[, group.by] %in% idents.keep)] <- NA
        # Generate the new color scale
        colors.use <- check_consistency_colors_and_names(sample = sample, colors = colors.use, grouping_variable = group.by)
      # If split.by is used intead.
      } else if (split_by_is_used){
        if (isFALSE(length(idents.keep) == sum(idents.keep %in% unique(sample@meta.data[, split.by])))){
          stop("All the values in idents.keep must be in the split.by variable provided.", call. = F)
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
    not_highlighting_and_not_split_by <- is.null(cells.highlight) & is.null(idents.highlight) & is.null(split.by)
    split_by_used <- is.null(cells.highlight) & is.null(idents.highlight) & !(is.null(split.by))
    highlighting_cells <- !(is.null(cells.highlight)) | !(is.null(idents.highlight))
    if (not_highlighting_and_not_split_by){
        p <- Seurat::DimPlot(sample,
                                  reduction = reduction,
                                  label = label,
                                  dims = dims,
                                  repel = ifelse(is.null(label) == TRUE, NULL, TRUE),
                                  label.box = ifelse(is.null(label) == TRUE, NULL, TRUE),
                                  label.color = ifelse(is.null(label) == TRUE, NULL, "black"),
                                  na.value = na.value,
                                  shuffle = shuffle,
                                  order = order,
                                  pt.size = pt.size,
                                  group.by = group.by,
                                  cols = colors.use,
                                  raster = raster,
                                  raster.dpi = c(raster.dpi, raster.dpi),
                                  ncol = ncol) &
            ggpubr::theme_pubr(legend = legend.position) &
            ggplot2::theme(plot.title = ggplot2::element_text(size = plot.title.fontsize, face = "bold", hjust = 0.5),
                           legend.text = ggplot2::element_text(size = legend.text.fontsize, face = "bold"),
                           legend.title = ggplot2::element_text(size = legend.title.fontsize, face = "bold")) &
            ggplot2::guides(color = ggplot2::guide_legend(ncol = legend.ncol,
                                                          nrow = legend.nrow,
                                                          byrow = legend.byrow,
                                                          override.aes = list(size = legend.icon.size),
                                                          title.position = legend.title.position))
        if (!(is.null(group.by))){
          # Remove automatic title.
          p <- p & ggplot2::ggtitle("")
        }
    }
    # If the UMAP has to be split in multiple panes.
    else if (split_by_used){
        # If the user provided multiple highlighting colors.
        multiple_colors <- ifelse(length(colors.use) > 1, TRUE, FALSE)
        # List to store each individual plots.
        list.plots <- list()
        # Recover all metadata.
        data.use <- sample@meta.data[, split.by, drop = F]
        # Retrieve the plotting order, keep factor levels if the column is a factor.
        plot_order <- if (is.factor(data.use[, 1])){levels(data.use[, 1])} else {sort(unique(data.use[, 1]))}
        # If idents.keep is used, subset to only these values.
        if (!(is.null(idents.keep))){
          plot_order <- if (is.factor(data.use[, 1])){levels(data.use[, 1])[levels(data.use[, 1]) %in% idents.keep]} else {sort(unique(data.use[, 1])[unique(data.use[, 1]) %in% idents.keep])}
          # If the user wants more than one color.
          if (isTRUE(multiple_colors)){
            colors.use <- colors.use[names(colors.use) %in% idents.keep]
          }
        }
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
                                 raster.dpi = c(raster.dpi, raster.dpi),
                                 ncol = ncol) &
                  ggplot2::ggtitle(iteration) &
                  ggpubr::theme_pubr(legend = legend.position)
            p <- add_scale(p = p,
                           function_use = ggplot2::scale_color_manual(labels = c("Not selected", iteration),
                                                                      values = c(na.value, ifelse(multiple_colors == TRUE, colors.use[[iteration]], colors.use))),
                           scale = "color") &
                 ggplot2::theme(plot.title = ggplot2::element_text(size = plot.title.fontsize, face = "bold", hjust = 0.5),
                                legend.text = ggplot2::element_text(size = legend.text.fontsize, face = "bold"),
                                legend.title = ggplot2::element_text(size = legend.title.fontsize, face = "bold")) &
                 ggplot2::guides(color = ggplot2::guide_legend(ncol = legend.ncol,
                                                               nrwo = legend.nrow,
                                                               byrow = legend.byrow,
                                                               override.aes = list(size = legend.icon.size),
                                                               title.position = legend.title.position))
            list.plots[[iteration]] <- p
        }
        # Assemble individual plots as a patch.
        p <- patchwork::wrap_plots(list.plots, ncol = ncol)
    }


    # If the user wants to highlight some of the cells.
    else if (highlighting_cells){
      # Compute the cells to highlight.
      if (is.null(idents.highlight) & !(is.null(cells.highlight))){
        # Only if cells.highlight parameters is used.
        cells.use <- cells.highlight
      } else if (!(is.null(idents.highlight)) & is.null(cells.highlight)){
        # Only if idents.highlight parameter is used.
        # Check if the provided identities are part of the active identities in the object.
        check_identity(sample = sample, identities = idents.highlight)
        cells.use <- names(Seurat::Idents(sample)[Seurat::Idents(sample) %in% idents.highlight])
      } else if (!(is.null(idents.highlight)) & !(is.null(cells.highlight))){
        # Check if the provided identities are part of the active identities in the object.
        check_identity(sample = sample, identities = idents.highlight)
        # Both idents.highlight and cells.highlight are used.
        cells.1 <- cells.highlight
        cells.2 <- names(Seurat::Idents(sample)[Seurat::Idents(sample) %in% idents.highlight])
        cells.use <- unique(c(cells.1, cells.2))
      }
        p <- Seurat::DimPlot(sample,
                             reduction = reduction,
                             cells.highlight = cells.use,
                             sizes.highlight = sizes.highlight,
                             dims = dims,
                             pt.size = pt.size,
                             raster = raster,
                             raster.dpi = c(raster.dpi, raster.dpi),
                             ncol = ncol) &
             ggpubr::theme_pubr(legend = legend.position)
        p <- add_scale(p = p,
                       function_use = ggplot2::scale_color_manual(labels = c("Not selected", "Selected cells"),
                                                                  values = c(na.value, colors.use)),
                       scale = "color") &
             ggplot2::theme(plot.title = ggplot2::element_text(size = plot.title.fontsize, face = "bold", hjust = 0.5),
                            legend.text = ggplot2::element_text(size = legend.text.fontsize, face = "bold"),
                            legend.title = ggplot2::element_text(size = legend.title.fontsize, face = "bold")) &
             ggplot2::guides(color = ggplot2::guide_legend(ncol = legend.ncol,
                                                           nrow = legend.nrow,
                                                           byrow = legend.byrow,
                                                           override.aes = list(size = legend.icon.size),
                                                           title.position = legend.title.position))
    }


    # General additions to all kind of plots.
    if (!is.null(plot.title)){
      if (!(is.null(split.by))){
        p <- p & patchwork::plot_annotation(title = plot.title,
                                            theme = ggplot2::theme(plot.title = ggplot2::element_text(size = plot.title.fontsize + 1,
                                                                                                      face = "bold",
                                                                                                      hjust = 0.5)))
      } else {
        p <- p & ggplot2::ggtitle(plot.title)
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

