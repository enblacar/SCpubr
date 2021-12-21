#' Wrapper for \link[Seurat]{DimPlot}.
#'
#'
#' @param sample # Seurat object.
#' @param reduction # Reduction to use. Can be the canonical ones such as "umap", "pca", or any custom ones, such as "diffusion". If you are unsure about which reductions you have, use `Seurat::Reductions(sample)`.
#' @param group.by # Variable you want the cells to be colored for.
#' @param split.by # Split into as many plots as unique values in the variable provided.
#' @param colors.use # Vector of named HEX values to color the cells. It has to match the number of unique values in either `Seurat::Idents(sample)` or the group.by variable.
#' @param cols.split # Vector of named HEX values to color the cells in a split DimPlot. It has to match the unique values in split.by parameter.
#' @param label # Whether to plot the cluster labels in the UMAP. The cluster labels will have the same color as the cluster colors.
#' @param cells.highlight # Vector of cells for which the DimPlot should focus into. The rest of the cells will be grayed out.
#' @param cols.highlight # HEX color code to use with the highlighted cells.
#' @param shuffle # Whether to shuffle the cells or not, so that they are not plotted cluster-wise. Recommended.
#' @param pt.size # Point size of the cells.
#' @param legend # Whether to plot the legend or not.
#' @param legend.title # Logical stating whether the legend title is shown or not.
#' @param legend.ncol # Number of columns in the legend.
#' @param legend.text.size # Fontsize of the legend labels.
#' @param legend.title.size # Fontisize of the legend title.
#' @param legend.icon.size # Size of the icons in legend.
#' @param legend.position # Position of the legend in the plot. Will only work if legend is set to TRUE.
#' @param legend.byrow # Logical stating whether the legend is filled by row or not.
#' @param plot.title # Title to use in the plot.
#' @param ncol # Number of columns used in the arrangement of the output plot using "split.by" parameter.
#'
#' @return
#' @export
#'
#' @examples
do_DimPlot <- function(sample,
                       reduction = "umap",
                       label = FALSE,
                       label.color = "black",
                       repel = TRUE,
                       shuffle = TRUE,
                       pt.size = 0.5,
                       group.by = NULL,
                       split.by = NULL,
                       cols.split = "#0A305F",
                       cells.highlight = NULL,
                       cols.highlight = "#0A305F",
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
                       dims = c(1, 2)){

    # Checks to ensure proper function.
    # Check whether the names of colors.use match the unique values in group.by or whether the number of colors is lower to the number of unique values in group.by.
    if (!(is.null(group.by)) & !(is.null(colors.use))){
        if (sum(names(colors.use) %in% unique(sample[[]][, group.by])) != length(unique(sample[[]][, group.by]))){
            stop('The names of the color vector provided to "colors.use" do not entirely match the unique values in "group.by" parameter.')
        }
        if (length(colors.use) != length(unique(sample[[]][, group.by]))){
            stop('The number of values provided to "colors.use" is lower than the unique values in "group.by" parameter.')
        }
    }
    # Check whether the names of cols.split match the unique values in split.by or whether the number of colors is lower to the number of unique values in split.by.
    if (!(is.null(split.by)) & cols.split != "#0A305F"){
        if (length(cols.split) > 1){
            if (sum(names(cols.split) %in% unique(sample[[]][, split.by])) != length(unique(sample[[]][, split.by]))){
                stop('The names of the color vector provided to "cols.highlight" do not entirely match the unique values in "split.by" parameter.')
            }
            if (length(cols.split) != length(unique(sample[[]][, split.by]))){
                stop('The number of values provided to "cols.split" is lower than the unique values in "split.by" parameter.')
            }
        }
    }

    # Check whether the user has provided only one color to cols.highlight.
    if (cols.highlight != "#0A305F"){
        # Check if the input is a character object or not.
        if (!(is.character(cols.highlight))){
            stop("The value for cols.highlight must be a character containing a HEX code.")
        }
    }
    # From: https://stackoverflow.com/a/13290832
    # Check that the input colors are valid color representations.
    areColors <- function(x) {
        sapply(x, function(X) {
            tryCatch(is.matrix(grDevices::col2rgb(X)),
                     error = function(e) FALSE)
        })
    }
    # Check for colors.use.
    if (!is.null(colors.use)){
        check <- areColors(colors.use)
        if (sum(check) != length(colors.use)){
            stop("Not all provided colors for colors.use are valid color representations.")
        }
    }
    # Check for cols.split.
    if (cols.split != "#0A305F"){
        check <- areColors(cols.split)
        if (sum(check) != length(cols.split)){
            stop("Not all provided colors for cols.split are valid color representations.")
        }
    }
    # Check for cols.highlight.
    if (sum(areColors(cols.highlight)) != length(cols.highlight)){
        stop("The value for cols.highlight is not a valid color representation.")
    }

    # If the user does not want to highlight any cells (Regular case.).
    if (is.null(cells.highlight)){

        # If no special color palette is provided by the user, generate a default one.
        if (is.null(colors.use)){
            # If no special grouping is set up, DimPlot defaults back to Seurat::Idents(sample) or levels(sample).
            if (is.null(group.by)){
                # Generate a vector of colors equal to the number of identities in the sample.
                colors.use <- colortools::wheel("#457b9d", length(levels(sample)))
                names(colors.use) <- levels(sample)
            } else {
                # Generate a vector of colors equal to the number of unique identities in the grouping variable.
                colors.use <-  colortools::wheel("#457b9d", length(unique(sample[[]][, group.by])))
                names(colors.use) <- unique(sample[[]][, group.by])
            }
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
                                      ncol = ncol) +
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
        # If the UMAP has to be split in multiple panes.
        } else {
            # If the user provided multiple highlighting colors.
            if (length(cols.split) > 1) {
                # List to store each individual plots.
                list.plots <- list()
                # Recover all metadata.
                data <- sample[[]]
                # Retrieve the metadata column belonging to the split.by parameter.
                data.use <- data[, split.by, drop = F]

                # Iterate over each unique value in split.by parameter.
                for (iteration in unique(data.use[, 1])){
                    # Retrieve the cells that do belong to the iteration's split.by value.
                    cells.highlight <- rownames(data.use[data.use[[split.by]] == iteration, , drop = F])
                    p.umap <- Seurat::DimPlot(sample,
                                              reduction = reduction,
                                              dims = dims,
                                              cells.highlight = cells.highlight,
                                              sizes.highlight = pt.size,
                                              pt.size = pt.size,
                                              raster = raster,
                                              ncol = ncol) +
                        ggplot2::ggtitle(iteration) +
                        ggpubr::theme_pubr(legend = legend.position) +
                        ggplot2::scale_color_manual(labels = c("Unselected", "Selected"),
                                                    values = c("grey75", cols.split[[iteration]]))  +
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
            # If the user did not provide a vector of colors, therefore using the default value in this function.
            } else {
                # List to store each individual plot.
                list.plots <- list()
                # Retrieve metadta from sample.
                data <- sample[[]]
                # Retrieve only the metadata belonging to split.by parameter.
                data.use <- data[, split.by, drop = F]
                # Iterate over each unique value in split.by parameter.
                for (iteration in unique(data.use[, 1])){
                    # Recover the cells for which the iteration value of split.by is true.
                    cells.highlight <- rownames(data.use[data.use[[split.by]] == iteration, , drop = F])
                    p.umap <- Seurat::DimPlot(sample,
                                              reduction = reduction,
                                              dims = dims,
                                              cells.highlight = cells.highlight,
                                              sizes.highlight = pt.size,
                                              pt.size = pt.size,
                                              raster = raster,
                                              ncol = ncol) +
                        ggplot2::ggtitle(iteration) +
                        ggpubr::theme_pubr(legend = legend.position) +
                        ggplot2::scale_color_manual(labels = c("Unselected", "Selected"),
                                                    values = c("grey75", cols.split))  +
                        ggplot2::theme(plot.title = ggplot2::element_text(face = "bold", hjust = 0.5),
                                       legend.text = ggplot2::element_text(size = legend.text.size, face = "bold"),
                                       legend.title = ggplot2::element_text(size = legend.title.size, face = "bold")) +
                        ggplot2::guides(color = ggplot2::guide_legend(ncol = legend.ncol,
                                                                      byrow = legend.byrow,
                                                                      override.aes = list(size = legend.icon.size)))
                    list.plots[[iteration]] <- p.umap
                }
                # Assemble individual panes together.
                p.umap <- patchwork::wrap_plots(list.plots, ncol = ncol)
            }
        }
    # If the user wants to highlight some of the cells.
    } else {
        p.umap <- Seurat::DimPlot(sample,
                                  reduction = reduction,
                                  cells.highlight = cells.highlight,
                                  sizes.highlight = pt.size,
                                  dims = dims,
                                  pt.size = pt.size,
                                  raster = raster,
                                  ncol = ncol) +
            ggplot2::ggtitle(plot.title) +
            ggpubr::theme_pubr(legend = legend.position) +
            ggplot2::scale_color_manual(labels = c("Unselected", "Selected"),
                                        values = c("grey", cols.highlight))  +
            ggplot2::theme(plot.title = ggplot2::element_text(face = "bold", hjust = 0.5),
                           legend.text = ggplot2::element_text(size = legend.text.size, face = "bold"),
                           legend.title = ggplot2::element_text(size = legend.title.size, face = "bold")) +
            ggplot2::guides(color = ggplot2::guide_legend(ncol = legend.ncol,
                                                          byrow = legend.byrow,
                                                          override.aes = list(size = legend.icon.size)))
    }

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

