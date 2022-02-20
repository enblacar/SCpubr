#' Wrapper for \link[Seurat]{FeaturePlot}.
#'
#'
#' @param sample Seurat object.
#' @param assay Assay to use. Defauls to the current assay.
#' @param reduction Reduction to use. Can be the canonical ones such as "umap", "pca", or any custom ones, such as "diffusion". If you are unsure about which reductions you have, use `Seurat::Reductions(sample)`. Defaults to "umap" if present or to the last computed reduction if the argument is not provided.
#' @param slot Data slot to use. Character. Only one of: counts, data, scale.data. Defaults to "data".
#' @param features Features to plot. It can be a single one or a vector of multiple features. Similar behavior as with \link[Seurat]{FeaturePlot}.
#' @param pt.size Point size.
#' @param split.by Split the plot in as many unique values stored in the provided metadata column.
#' @param split.by.idents Vector of identities to plot. The gradient scale will also be subset to only the values of such identities.
#' @param legend Whether to plot the legend or not.
#' @param legend.position Position of the legend in the plot. Will only work if legend is set to TRUE.
#' @param plot.title Title for the plot.
#' @param fontsize Base fontsize of the figure.
#' @param ncol Number of columns to use in the arrangement of the output if more than one feature is queried to the function.
#' @param cells.highlight Vector of cells for which the FeaturePlot should focus into. The rest of the cells will be grayed out.
#' @param idents.highlight Vector of identities that the FeaturePlot should focus into. Has to match the current Seurat identities in `Seurat::Idents(sample)`.
#' @param dims Vector of 2 numerics indicating the dimensions to plot out of the selected reduction. Defaults to c(1, 2) if not specified.
#' @param viridis_color_map Character. A capital letter from A to H or the scale name as in \link[viridis]{scale_fill_viridis}.
#' @param verbose Whether to show warnings.
#' @param individual.titles Titles for each feature if needed. Either NULL or a vector of equal length of features.
#' @return  A ggplot2 object containing a Feature Plot.
#' @export
#'
#' @examples
#' \dontrun{
#' TBD
#' }
do_FeaturePlot <- function(sample,
                           assay = NULL,
                           features,
                           reduction = NULL,
                           pt.size = 0.5,
                           legend = TRUE,
                           split.by = NULL,
                           split.by.idents = NULL,
                           slot = NULL,
                           legend.position = "right",
                           plot.title = NULL,
                           fontsize = 14,
                           ncol = NULL,
                           cells.highlight = NULL,
                           idents.highlight = NULL,
                           dims = c(1, 2),
                           viridis_color_map = "D",
                           verbose = TRUE,
                           individual.titles = NULL){
    # Checks for packages.
    check_suggests(function_name = "do_FeaturePlot")
    # Check the assay.
    out <- check_and_set_assay(sample = sample, assay = assay)
    sample <- out[["sample"]]
    assay <- out[["assay"]]
    # Check the reduction.
    reduction <- check_and_set_reduction(sample = sample, reduction = reduction)
    # Check the dimensions.
    dimensions <- check_and_set_dimensions(sample = sample, reduction = reduction, dims = dims)
    # Check logical parameters.
    logical_list <- list("legend" = legend,
                         "verbose" = verbose)
    check_type(parameters = logical_list, required_type = "logical", test_function = is.logical)
    # Check numeric parameters.
    numeric_list <- list("pt.size" = pt.size,
                         "ncol" = ncol,
                         "fontsize" = fontsize)
    check_type(parameters = numeric_list, required_type = "numeric", test_function = is.numeric)
    # Check character parameters.
    character_list <- list("legend.position" = legend.position,
                           "features" = features,
                           "cells.highlight" = cells.highlight,
                           "idents.highlight" = idents.highlight,
                           "slot" = slot,
                           "split.by" = split.by,
                           "plot.title" = plot.title,
                           "split.by.idents" = split.by.idents,
                           "viridis_color_map" = viridis_color_map,
                           "individual.titles" = individual.titles)
    check_type(parameters = character_list, required_type = "character", test_function = is.character)

    # Check slot.
    slot <- check_and_set_slot(slot = slot)

    # Check split.by is on metadata.
    if (!(is.null(split.by))){check_feature(sample = sample, features = split.by, enforce_check = "metadata", enforce_parameter = "split.by")}

    # Check individual titles.
    if (!(is.null(individual.titles))){
      if(length(features) != length(individual.titles)){
        stop('Total number of individual titles does not match the number of features provided.', call. = F)
      }
    }

    # Check viridis_color_map.
    check_viridis_color_map(viridis_color_map = viridis_color_map, verbose = verbose)

    # Define fontsize parameters.
    plot.title.fontsize <- fontsize + 2
    axis.text.fontsize <- fontsize
    axis.title.fontsize <- fontsize + 1
    legend.text.fontsize <- fontsize - 4
    legend.title.fontsize <- fontsize - 4

    # Regular FeaturePlot.
    check <- is.null(split.by) & is.null(cells.highlight) & is.null(idents.highlight)
    if (check){
        # Check if the feature is actually in the object.
        check_feature(sample = sample, features = features)
        p <- Seurat::FeaturePlot(sample,
                                 features,
                                 reduction = reduction,
                                 order = T,
                                 dims = dims,
                                 pt.size = pt.size,
                                 ncol = ncol) &
            Seurat::NoAxes() &
            viridis::scale_color_viridis(na.value = "grey75", option = viridis_color_map) &
            ggplot2::theme(plot.title = ggplot2::element_text(size = plot.title.fontsize, face = "bold", hjust = 0.5),
                           legend.text = ggplot2::element_text(size = legend.text.fontsize, face = "bold", hjust = 1),
                           legend.position = legend.position)
        # Special patches for diffusion maps: Adding "DC" labels to the axis.
        if (reduction == "diffusion"){
            p <- p & ggplot2::xlab(paste0("DC_", dims[1])) & ggplot2::ylab(paste0("DC_", dims[2]))
        }

    # Modified FeaturePlot including only a subset of cells.
    } else {
      # Check if the feature is actually in the object.
      dim_colnames <- check_feature(sample = sample, features = features, dump_reduction_names = TRUE)
        # Get the subset of wanted cells according to the combination of idents.highlight and cells.highlight parameters.
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
        } else if (!(is.null(split.by))){
            if (is.null(split.by.idents)){
              cells.use <- colnames(sample)
            } else {
              if (sum(duplicated(split.by.idents)) != 0){
                message("Found and removed duplicated values in split.by.idents.")
                split.by.idents <- split.by.idents[!duplicated(split.by.idents)]
              }
              cells.use <- names(Seurat::Idents(sample)[Seurat::Idents(sample) %in% split.by.idents])
            }
        }
        # Plots are generated independently if more than one feature is provided.
        list.plots <- list()
        count_iteration <- 0
        for (feature in features){
            # A "dummy" metadta column is generated using the values of the selected feature.
            if (feature %in% colnames(sample@meta.data)) {
              sample$dummy <- sample@meta.data[, feature]
            } else if (feature %in% rownames(sample)){
              sample$dummy <- Seurat::GetAssayData(object = sample, slot = slot)[feature, ]
            } else if (feature %in% dim_colnames){
              for(red in Seurat::Reductions(object = sample)){
                if (feature %in% colnames(sample@reductions[[red]][[]])){
                  red.feature <- red
                  sample$dummy <- sample@reductions[[red.feature]][[]][, feature]
                }
              }
            }
            # Assign NAs to the values corresponding to the cells  not selected.
            sample$dummy[!(names(sample$dummy) %in% cells.use)] <- NA

            if (is.null(split.by)){
              p.loop <- Seurat::FeaturePlot(sample,
                                            "dummy",
                                            reduction = reduction,
                                            order = T,
                                            dims = dims,
                                            pt.size = pt.size) +
                # This is actually a "fake cell" with alpha 0 (invisible) sitting in the left-down corner, which helps adding a new legend label for the grayed out not selected (NS) cells.
                ggplot2::geom_point(mapping = ggplot2::aes(x = min(Seurat::Embeddings(sample, reduction)[, 1]),
                                                           y = min(Seurat::Embeddings(sample, reduction)[, 2]),
                                                           fill = "NS"),
                                    alpha = 0) +
                Seurat::NoAxes() +
                viridis::scale_color_viridis(na.value = "grey75", option = viridis_color_map) +
                ggplot2::ggtitle(feature) +
                ggplot2::theme(plot.title = ggplot2::element_text(size = plot.title.fontsize, face = "bold", hjust = 0.5),
                               legend.text = ggplot2::element_text(size = legend.text.fontsize, face = "bold"),
                               legend.position = legend.position,
                               legend.spacing = ggplot2::unit(0.01, "cm")) +
                # Override aesthetic of the legend, providing the desired gray color.
                ggplot2::guides(fill = ggplot2::guide_legend("", override.aes = list(color = "grey75",
                                                                                     alpha = 1)))
            } else {
              # Recover all metadata.
              data <- sample[[]]
              # Retrieve the metadata column belonging to the split.by parameter.
              data.use <- data[, split.by, drop = F]
              # Retrieve the plotting order, keep factor levels if the column is a factor.
              if (is.null(split.by.idents)){
                plot_order <- if (is.factor(data.use[, 1])){levels(data.use[, 1])} else {sort(unique(data.use[, 1]))}
              } else {
                plot_order <- sort(split.by.idents)
              }
              list.plots.split.by <- list()
              count_plot <- 0
              count_iteration <- count_iteration + 1
              if (legend.position != "right"){
                if (verbose){warning("Using split.by parameter will force the legend of each individual plot to the right.")}
              }
              limits <- c(min(sample$dummy), max(sample$dummy))
              for (iteration in plot_order){
                count_plot <- count_plot + 1
                # Create a second dummy variable.
                sample$dummy2 <- sample$dummy
                cells.iteration <- sample[[]][, split.by] == iteration
                sample$dummy2[!(cells.iteration)] <- NA
                p.loop <- Seurat::FeaturePlot(sample,
                                              "dummy2",
                                              reduction = reduction,
                                              order = T,
                                              dims = dims,
                                              pt.size = pt.size) +
                  # This is actually a "fake cell" with alpha 0 (invisible) sitting in the left-down corner, which helps adding a new legend label for the grayed out not selected (NS) cells.
                  ggplot2::geom_point(mapping = ggplot2::aes(x = min(Seurat::Embeddings(sample, reduction)[, 1]),
                                                             y = min(Seurat::Embeddings(sample, reduction)[, 2]),
                                                             fill = "NS"),
                                      alpha = 0) +
                  Seurat::NoAxes() +
                  ggplot2::scale_color_viridis_c(limits = limits, na.value = "grey75", option = viridis_color_map) +
                  ggplot2::ggtitle(feature) +
                  ggplot2::ggtitle(ifelse(count_iteration == 1, iteration, "")) +
                  ggplot2::ylab(ifelse(count_plot == 1, feature, "")) +
                  ggplot2::theme(plot.title = ggplot2::element_text(size = plot.title.fontsize, face = "bold", hjust = 0.5),
                                 axis.title.y = ggplot2::element_text(size = axis.title.fontsize, face = "bold", vjust = 0.5),
                                 legend.text = ggplot2::element_text(size = legend.text.fontsize, face = "bold"),
                                 legend.position = ifelse(legend.position != "right", "right", legend.position),
                                 legend.spacing = ggplot2::unit(0.01, "cm")) +
                  # Override aesthetic of the legend, providing the desired gray color.
                  ggplot2::guides(fill = ggplot2::guide_legend("", override.aes = list(color = "grey75",
                                                                                       alpha = 1)))
                # Only keep the legend of the last plot, as we are setting the same color scale.
                if (count_plot != length(plot_order)) {p.loop <- p.loop + Seurat::NoLegend()}

                list.plots.split.by[[iteration]] <- p.loop
              }
              # Assemble individual plots as a patch.
              p.loop <- patchwork::wrap_plots(list.plots.split.by, ncol = length(plot_order))
            }

            # Patch for diffusion maps.
            if (reduction == "diffusion"){
                # Add "DC" labels.
                p.loop <- p.loop & ggplot2::xlab(paste0("DC_", dims[1])) & ggplot2::ylab(paste0("DC_", dims[2]))
            }
            # Add the plot to the list.
            list.plots[[feature]] <- p.loop
        }
        # Generate the final plot with patchwork and use the "ncol" parameter value for the number of colums.
        if (is.null(split.by)){
          p <- patchwork::wrap_plots(list.plots, ncol = ncol)
        } else {
          if (!(is.null(ncol))){message("Ignoring 'ncol' parameter as split.by parameter is used. Number of columns equal to each unique value in split.by.")}
          p <- patchwork::wrap_plots(list.plots, ncol = 1)
        }


        # Patch for the case in which features only contains one element and a custom plot title is provided.
        # Basically, as this is a "patchwork" object, the way the title has to be set is different than using "ggplot2::ggtitle()".
        if (!is.null(plot.title) & length(features) == 1){
            p[[1]]$labels$title <- plot.title
            p <- p[[1]]
        }
    }

    # Add custom title.
    if (!is.null(plot.title)){
        if (length(features) > 1){
            p <- p + patchwork::plot_annotation(title = plot.title,
                                                theme = ggplot2::theme(plot.title = ggplot2::element_text(size = plot.title.fontsize + 1,
                                                                                                          face = "bold",
                                                                                                          hjust = 0.5)))
        } else {
            p <- p + ggplot2::ggtitle(plot.title)
        }
    }

    # Add individual titles.
    if (!is.null(individual.titles)){
      for (counter in seq(1,length(features))){
        if (!(is.na(individual.titles[counter]))){
          p[[counter]]$labels$title <- individual.titles[counter]
        }
      }
    }

    # Remove legend.
    if (legend == FALSE){
        p <- p + ggpubr::rremove("legend")
    }

    # Further patch for diffusion maps.
    if (reduction == "diffusion"){
        # Fix the axis scale so that the highest and lowest values are in the range of the DCs (previously was around +-1.5, while DCs might range to +-0.004 or so).
        p <- p &
            ggplot2::xlim(c(min(sample@reductions$diffusion[[]][, paste0("DC_", dims[1])]),
                          max(sample@reductions$diffusion[[]][, paste0("DC_", dims[1])]))) &
            ggplot2::ylim(c(min(sample@reductions$diffusion[[]][, paste0("DC_", dims[2])]),
                          max(sample@reductions$diffusion[[]][, paste0("DC_", dims[2])]))) &
            # Remove axis elements so that the axis title is the only thing left.
            ggpubr::rremove("axis") &
            ggpubr::rremove("axis.text") &
            ggpubr::rremove("ticks") &
            ggplot2::theme(axis.title.x = ggplot2::element_text(size = axis.title.fontsize, face = "bold"),
                                                                                                                                                                                                                                                                                                                                                                                                                                                            axis.title.y = ggplot2::element_text(size = 14, face = "bold")) &
            ggplot2::theme(axis.title.y = ggplot2::element_text(angle = 90))
    }
    # Return the plot.
    return(p)
}
