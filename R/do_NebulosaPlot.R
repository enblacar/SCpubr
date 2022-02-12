#' Wrapper for Nebulosa::plot_density in Seurat.
#'
#'
#' @param sample Seurat object.
#' @param reduction Reduction to use. Can be the canonical ones such as "umap", "pca", or any custom ones, such as "diffusion". If you are unsure about which reductions you have, use `Seurat::Reductions(sample)`.
#' @param features Features to plot density for. It can be a single one or a vector of multiple features.
#' @param slot Slot to retrieve the data from. Defaults to "data".
#' @param size Size of the dots.
#' @param combine Whether to create a single plot out of multiple features.
#' @param method Kernel density estimation method. Either "ks" or "wkde" or both. See \link[Nebulosa]{plot_density} for more details.
#' @param dims Vector of 2 dims to plot the data. By default, first and second from the specified reduction.
#' @param joint Whether to plot different features as joint density.
#' @param pal Viridis palette to use. Use the names.
#' @param shape Shape of the geometry (ggplot number).
#' @param legend Whether to plot the legend or not. Logical.
#' @param legend.position Position of the legend in the plot.
#' @param plot.title Title to use in the plot.
#'
#' @return  A ggplot2 object containing a Nebulosa plot.
#' @export
#'
#' @examples
#' \dontrun{
#' TBD
#' }
do_NebulosaPlot <- function(sample,
                             features,
                             slot = NULL,
                             dims = c(1, 2),
                             size = 1,
                             reduction = NULL,
                             combine = TRUE,
                             method = c("ks", "wkde"),
                             joint = FALSE,
                             plot.title = NULL,
                             pal = "viridis",
                             shape = 16,
                             legend = TRUE,
                             legend.position = "right"){
  # Checks for packages.
  check_suggests(function_name = "do_Nebulosa_plot")
  # Check the reduction.
  reduction <- check_and_set_reduction(sample = sample, reduction = reduction)
  # Check the dimensions.
  dimensions <- check_and_set_dimensions(sample = sample, reduction = reduction, dims = dims)
  # Check logical parameters.
  logical_list <- list("legend" = legend,
                       "combine" = combine,
                       "joint" = joint)
  check_type(parameters = logical_list, required_type = "logical", test_function = is.logical)
  # Check numeric parameters.
  numeric_list <- list("size" = size,
                       "shape" = shape)
  check_type(parameters = numeric_list, required_type = "numeric", test_function = is.numeric)
  # Check character parameters.
  character_list <- list("legend.position" = legend.position,
                         "features" = features,
                         "pal" = pal,
                         "method" = method,
                         "plot.title" = plot.title,
                         "slot" = slot)
  check_type(parameters = character_list, required_type = "character", test_function = is.character)
  # Check slot.
  slot <- check_and_set_slot(slot = slot)

  # Check if the feature is actually in the object.
  check_feature(sample = sample, features = features)

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

    # Remove legend.
    if (legend == FALSE){
      plot <- plot + Seurat::NoLegend()
    }
    # Add a title.
    if (!(is.null(plot.title))){
        if (length(features) == 1){
            plot <- plot + ggplot2::ggtitle(plot.title)
        } else {
            plot <- plot + patchwork::plot_annotation(title = plot.title,
                                                      theme = ggplot2::theme(plot.title = ggplot2::element_text(size = 18,
                                                                                                                face = "bold",
                                                                                                                hjust = 0.5)))
        }
    }
    return(plot)
}
