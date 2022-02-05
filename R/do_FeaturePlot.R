#' Wrapper for \link[Seurat]{FeaturePlot}.
#'
#'
#' @param sample Seurat object.
#' @param assay Assay to use.
#' @param reduction Reduction to use. Can be the canonical ones such as "umap", "pca", or any custom ones, such as "diffusion". If you are unsure about which reductions you have, use `Seurat::Reductions(sample)`.
#' @param features Features to plot. It can be a single one or a vector of multiple features. Similar behavior as with \link[Seurat]{FeaturePlot}.
#' @param pt.size Point size.
#' @param legend Whether to plot the legend or not.
#' @param legend.position Position of the legend in the plot. Will only work if legend is set to TRUE.
#' @param scale.begin Value to where you want the continuous color scale to start. It has to be in the range of the feature's values.
#' @param scale.end Value to where you want the continuous color scale to end. It has to be in the range of the feature's values.
#' @param plot.title Title for the plot.
#' @param ncol Number of columns to use in the arrangement of the output if more than one feature is queried to the function.
#' @param cells.highlight Vector of cells for which the FeaturePlot should focus into. The rest of the cells will be grayed out.
#' @param idents.highlight Vector of identities that the FeaturePlot should focus into. Has to match the current Seurat identities in `Seurat::Idents(sample)`.
#' @param dims Vector of 2 dimensions to use. Defaults to first and second dimensions.
#' @return  A ggplot2 object containing a Feature Plot.
#' @export
#'
#' @examples
#' \dontrun{
#' TBD
#' }
do_FeaturePlot <- function(sample,
                           assay = "SCT",
                           features,
                           reduction = "umap",
                           pt.size = 0.5,
                           legend = TRUE,
                           legend.position = "right",
                           plot.title = NULL,
                           ncol = NULL,
                           scale.begin = NULL,
                           scale.end = NULL,
                           cells.highlight = NULL,
                           idents.highlight = NULL,
                           dims = c(1, 2)){

    # Regular FeaturePlot.
    if (is.null(cells.highlight) & is.null(idents.highlight)){
        p <- Seurat::FeaturePlot(sample,
                                 features,
                                 reduction = reduction,
                                 order = T,
                                 dims = dims,
                                 pt.size = pt.size,
                                 ncol = ncol) &
            Seurat::NoAxes() &
            viridis::scale_color_viridis(na.value = "grey75") &
            ggplot2::theme(plot.title = ggplot2::element_text(face = "bold", hjust = 0.5),
                           legend.text = ggplot2::element_text(size = 10, face = "bold", hjust = 1),
                           legend.position = legend.position)
        # Special patches for diffusion maps: Adding "DC" labels to the axis.
        if (reduction == "diffusion"){
            p <- p & ggplot2::xlab(paste0("DC_", dims[1])) & ggplot2::ylab(paste0("DC_", dims[2]))
        }
        # Set the range of the color scale to scale.begin and scale.end parameters.
        if (!is.null(scale.begin) | !is.null(scale.end)){
            p <- p & ggplot2::scale_color_continuous(type = "viridis", limits = c(scale.begin, scale.end))
        }
    # Modified FeaturePlot including only a subset of cells.
    } else {
        # Get the subset of wanted cells according to the combination of idents.highlight and cells.highlight parameters.
        if (is.null(idents.highlight) & !(is.null(cells.highlight))){
            # Only if cells.highlight parameters is used.
            cells.use <- cells.highlight
        } else if (!(is.null(idents.highlight)) & is.null(cells.highlight)){
            # Only if idents.highlight parameter is used.
            cells.use <- names(Seurat::Idents(sample)[Seurat::Idents(sample) %in% idents.highlight])
        } else if (!(is.null(idents.highlight)) & !(is.null(cells.highlight))){
            # Both idents.highlight and cells.highlight are used.
            cells.1 <- cells.highlight
            cells.2 <- names(Seurat::Idents(sample)[Seurat::Idents(sample) %in% idents.highlight])
            cells.use <- unique(c(cells.1, cells.2))
        }
        # Plots are generated independently if more than one feature is provided.
        list.plots <- list()
        for (feature in features){
            # A "dummy" metadta column is generated using the values of the selected feature.
            # If the feature is in the metadata columns.
            if (feature %in% colnames(sample[[]])){
                sample$dummy <- sample[[]][, feature]
            # If the feature is a gene in the current active assay.
            } else {
                sample$dummy <- sample@assays[[assay]]@data[feature, ]
            }
            # Assign NAs to the values corresponding to the cells  not selected.
            sample$dummy[!(names(sample$dummy) %in% cells.use)] <- NA

            p.loop <- Seurat::FeaturePlot(sample,
                                          "dummy",
                                          reduction = reduction,
                                          order = T,
                                          dims = dims,
                                          pt.size = pt.size) +
                # This is actually a "fake cell" with alpha 0 (invisible), which helps adding a new legend label for the grayed out not selected (NS) cells.
                ggplot2::geom_point(mapping = ggplot2::aes(x = min(Seurat::Embeddings(sample, reduction)[, 1]),
                                                           y = min(Seurat::Embeddings(sample, reduction)[, 2]),
                                                           fill = "NS"),
                                    alpha = 0) +
                Seurat::NoAxes() +
                viridis::scale_color_viridis(na.value = "grey75") +
                ggplot2::ggtitle(feature) +
                ggplot2::theme(plot.title = ggplot2::element_text(face = "bold", hjust = 0.5),
                               legend.text = ggplot2::element_text(size = 10, face = "bold"),
                               legend.position = legend.position) +
                # Override aesthetic of the legend, providing the desired gray color.
                ggplot2::guides(fill = ggplot2::guide_legend("", override.aes = list(color = "grey75",
                                                                                     alpha = 1)))
            # Patch for diffusion maps.
            if (reduction == "diffusion"){
                # Add "DC" labels.
                p.loop <- p.loop + ggplot2::xlab(paste0("DC_", dims[1])) + ggplot2::ylab(paste0("DC_", dims[2]))
            }
            # Set the color scale to the values provided for scale.being and scale.end.
            if (!is.null(scale.begin) | !is.null(scale.end)){
                p.loop <- p.loop & ggplot2::scale_color_continuous(type = "viridis", limits = c(scale.begin, scale.end))
            }
            # Add the plot to the list.
            list.plots[[feature]] <- p.loop
        }
        # Generate the final plot with patchwork and use the "ncol" parameter value for the number of colums.
        p <- patchwork::wrap_plots(list.plots, ncol = ncol)

        # Patch for the case in which features only contains one element and a custom plot title is provided.
        # Basically, as this is a "patchwork" object, the way the title has to be set is different than using "ggplot2::ggtitle()".
        if (!is.null(plot.title) & length(features) == 1){
            p[[1]]$labels$title <- plot.title
            p <- p[[1]]
        }
        # Remove the dummy variable.
        sample$dummy <- NULL

    }

    # Add custom title.
    if (!is.null(plot.title)){
        if (length(features) > 1){
            p <- p + patchwork::plot_annotation(title = plot.title,
                                                theme = ggplot2::theme(plot.title = ggplot2::element_text(size = 18,
                                                                                                          face = "bold",
                                                                                                          hjust = 0.5)))
        } else {
            p <- p + ggplot2::ggtitle(plot.title)
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
            ggplot2::theme(axis.title.x = ggplot2::element_text(size = 14, face = "bold"),
                                                                                                                                                                                                                                                                                                                                                                                                                                                            axis.title.y = ggplot2::element_text(size = 14, face = "bold")) &
            ggplot2::theme(axis.title.y = ggplot2::element_text(angle = 90))
    }
    # Return the plot.
    return(p)
}
