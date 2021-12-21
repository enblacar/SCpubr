#' Wrapper for Nebulosa::plot_density in Seurat.
#'
#'
#' @param sample # Seurat object.
#' @param reduction # Reduction to use. Can be the canonical ones such as "umap", "pca", or any custom ones, such as "diffusion". If you are unsure about which reductions you have, use `Seurat::Reductions(sample)`.
#' @param features # Features to plot density for. It can be a single one or a vector of multiple features.
#' @param joint # Whether to plot different features as joint density.
#' @param legend.position # Position of the legend in the plot.
#' @param plot.title # Title to use in the plot.
#'
#' @return
#' @export
#'
#' @examples
do_Nebulosa_plot <- function(sample,
                             features,
                             joint = FALSE,
                             plot.title = NULL,
                             legend.position = "right",
                             reduction = "umap"){
    # Plot a density plot using Nebulosa package.
    plot <- Nebulosa::plot_density(object = sample,
                                   features = features,
                                   joint = joint,
                                   reduction = reduction) &
            Seurat::NoAxes() &
            ggplot2::theme(plot.title = ggplot2::element_text(face = "bold", hjust = 0.5),
                           legend.text = ggplot2::element_text(size = 10, face = "bold"),
                           legend.title = ggplot2::element_text(size = 12, face = "bold"),
                           legend.position = legend.position)
    # Add a title.
    if (!(is.null(plot.title))){
        if (length(features) == 1){
            plot <- plot + ggplot2::ggtitle(plot.title)
        } else {
            plot <- plot + patchwork::plot_annotation(title = plot.title,
                                                      theme = ggplot2::theme(plot.title = element_text(size = 18,
                                                                                                       face = "bold",
                                                                                                       hjust = 0.5)))
        }
    }
    return(plot)
}
