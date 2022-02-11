#' Wrapper for \link[Seurat]{FeaturePlot}.
#'
#'
#' @param sample Seurat object.
#' @param assay Assay to use. Defauls to the current assay.
#' @param reduction Reduction to use. Can be the canonical ones such as "umap", "pca", or any custom ones, such as "diffusion". If you are unsure about which reductions you have, use `Seurat::Reductions(sample)`. Defaults to "umap" if present or to the last computed reduction if the argument is not provided.
#' @param slot Data slot to use. Character. Only one of: counts, data, scale.data. Defaults to "data".
#' @param features Features to plot. It can be a single one or a vector of multiple features. Similar behavior as with \link[Seurat]{FeaturePlot}.
#' @param pt.size Point size.
#' @param legend Whether to plot the legend or not.
#' @param legend.position Position of the legend in the plot. Will only work if legend is set to TRUE.
#' @param plot.title Title for the plot.
#' @param ncol Number of columns to use in the arrangement of the output if more than one feature is queried to the function.
#' @param cells.highlight Vector of cells for which the FeaturePlot should focus into. The rest of the cells will be grayed out.
#' @param idents.highlight Vector of identities that the FeaturePlot should focus into. Has to match the current Seurat identities in `Seurat::Idents(sample)`.
#' @param dims Vector of 2 numerics indicating the dimensions to plot out of the selected reduction. Defaults to c(1, 2) if not specified.
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
                           slot = NULL,
                           legend.position = "right",
                           plot.title = NULL,
                           ncol = NULL,
                           cells.highlight = NULL,
                           idents.highlight = NULL,
                           dims = c(1, 2)){
    # Checks for packages.
    check_suggests(function_name = "do_FeaturePlot")
    # Check the assay.
    out <- check_and_set_assay(sample, assay = assay)
    sample <- out[["sample"]]
    assay <- out[["assay"]]
    # Check the reduction.
    reduction <- check_and_set_reduction(sample = sample, reduction = reduction)
    # Check the dimensions.
    dimensions <- check_and_set_dimensions(sample = sample, reduction = reduction, dims = dims)
    # Check logical parameters.
    logical_list <- list("legend" = legend)
    check_type(parameters = logical_list, required_type = "logical", test_function = is.logical)
    # Check numeric parameters.
    numeric_list <- list("pt.size" = pt.size)
    check_type(parameters = numeric_list, required_type = "numeric", test_function = is.numeric)
    if(!(is.null(ncol))){check_type(parameters = list("ncol" = ncol), required_type = "numeric", test_function = is.numeric)}
    # Check character parameters.
    character_list <- list("legend.position" = legend.position,
                           "plot.title" = plot.title,
                           "features" = features)
    check_type(parameters = character_list, required_type = "character", test_function = is.character)
    if(!(is.null(cells.highlight))){check_type(parameters = list("cells.highlight" = cells.highlight), required_type = "chracter", test_function = is.character)}
    if(!(is.null(idents.highlight))){check_type(parameters = list("idents.highlight" = idents.highlight), required_type = "chracter", test_function = is.character)}
    if(!(is.null(slot))){check_type(parameters = list("slot" = slot), required_type = "chracter", test_function = is.character)}
    # Check slot.
    slot <- check_and_set_slot(slot = slot)

    # Regular FeaturePlot.
    if (is.null(cells.highlight) & is.null(idents.highlight)){
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
            viridis::scale_color_viridis(na.value = "grey75") &
            ggplot2::theme(plot.title = ggplot2::element_text(face = "bold", hjust = 0.5),
                           legend.text = ggplot2::element_text(size = 10, face = "bold", hjust = 1),
                           legend.position = legend.position)
        # Special patches for diffusion maps: Adding "DC" labels to the axis.
        if (reduction == "diffusion"){
            p <- p & ggplot2::xlab(paste0("DC_", dims[1])) & ggplot2::ylab(paste0("DC_", dims[2]))
        }

    # Modified FeaturePlot including only a subset of cells.
    } else {
        # Check if the feature is actually in the object.
        check_feature(sample = sample, features = features)
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
