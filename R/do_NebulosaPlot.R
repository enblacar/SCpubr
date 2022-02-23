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
#' @param return_only_joint Whether to only return the joint density panel. Logical.
#' @param shape Shape of the geometry (ggplot number).
#' @param legend Whether to plot the legend or not. Logical.
#' @param legend.position Position of the legend in the plot.
#' @param fontsize Base fontsize of the plot.
#' @param plot.title Title to use in the plot.
#' @param individual.titles Titles for each feature if needed. Either NULL or a vector of equal length of features.
#' @param viridis_color_map Character. A capital letter from A to H or the scale name as in \link[viridis]{scale_fill_viridis}.
#' @param verbose Whether to show warnings.
#'
#' @return  A ggplot2 object containing a Nebulosa plot.
#' @export
#'
#' @example /man/examples/examples_do_NebulosaPlot.R
do_NebulosaPlot <- function(sample,
                             features,
                             slot = NULL,
                             dims = c(1, 2),
                             size = 1,
                             reduction = NULL,
                             combine = TRUE,
                             method = c("ks", "wkde"),
                             joint = FALSE,
                             return_only_joint = NULL,
                             plot.title = NULL,
                             individual.titles = NULL,
                             shape = 16,
                             legend = TRUE,
                             fontsize = 14,
                             legend.position = "right",
                             viridis_color_map = "D",
                             verbose = TRUE){
  # Checks for packages.
  check_suggests(function_name = "do_Nebulosa_plot")
  # Check the reduction.
  reduction <- check_and_set_reduction(sample = sample, reduction = reduction)
  # Check the dimensions.
  dimensions <- check_and_set_dimensions(sample = sample, reduction = reduction, dims = dims)
  # Check logical parameters.
  logical_list <- list("legend" = legend,
                       "combine" = combine,
                       "joint" = joint,
                       "return_only_joint" = return_only_joint)
  check_type(parameters = logical_list, required_type = "logical", test_function = is.logical)
  # Check numeric parameters.
  numeric_list <- list("size" = size,
                       "shape" = shape)
  check_type(parameters = numeric_list, required_type = "numeric", test_function = is.numeric)
  # Check character parameters.
  character_list <- list("legend.position" = legend.position,
                         "features" = features,
                         "method" = method,
                         "plot.title" = plot.title,
                         "slot" = slot,
                         "individual.titles" = individual.titles)
  check_type(parameters = character_list, required_type = "character", test_function = is.character)
  # Check slot.
  slot <- check_and_set_slot(slot = slot)

  # Check if the feature is actually in the object.
  check_feature(sample = sample, features = features)

  # Check individual titles.
  if (!(is.null(individual.titles))){
    # If joint is set up.
    if (isTRUE(joint)){
      if (!(is.null(return_only_joint)) & isTRUE(return_only_joint)){
        stop("If return_only_joint is set to TRUE, then only one title is needed. Use plot.title instead.", call. = F)
      } else {
        if(length(features) + 1 != length(individual.titles)){
          stop('Total number of individual titles does not match the number of features provided + 1 (for the joint density).', call. = F)
        }
      }
    } else {
      if(length(features) != length(individual.titles)){
        stop('Total number of individual titles does not match the number of features provided.', call. = F)
      }
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

  # Plot a density plot using Nebulosa package.
    p <- Nebulosa::plot_density(object = sample,
                                features = features,
                                joint = joint,
                                reduction = reduction) &
         Seurat::NoAxes() &
         ggplot2::theme(plot.title = ggplot2::element_text(size = plot.title.fontsize, face = "bold", hjust = 0.5),
                        legend.text = ggplot2::element_text(size = legend.text.fontsize, face = "bold"),
                        legend.title = ggplot2::element_text(size = legend.title.fontsize, face = "bold"),
                        legend.position = legend.position)
    # Compute the total number of plots according to whether joint is set to TRUE or not.
    if (isTRUE(joint)){
      num_plots <- length(features) + 1
    } else {
      num_plots <- length(features)
    }
    p <- add_scale(p = p,
                   num_plots = num_plots,
                   scale = "color",
                   function_use = viridis::scale_color_viridis(na.value = "grey75",
                                                               option = viridis_color_map))
    # Remove legend.
    if (legend == FALSE){
      p <- p + Seurat::NoLegend()
    }
    # Add a title.
    if (!(is.null(plot.title))){
        if (length(features) == 1 | (!(is.null(return_only_joint)) & isTRUE(return_only_joint))){
            if (isTRUE(return_only_joint)){
              p <- p[[length(features) + 1]]
              p <- p & ggplot2::ggtitle(plot.title)
            } else {
              p <- p & ggplot2::ggtitle(plot.title)
            }
        } else {
            p <- p & patchwork::plot_annotation(title = plot.title,
                                                theme = ggplot2::theme(plot.title = ggplot2::element_text(size = plot.title.fontsize + 1,
                                                                                                          face = "bold",
                                                                                                          hjust = 0.5)))
        }
    }

    # Add individual titles.
    if (!is.null(individual.titles)){
      times <- length(features)
      if (isTRUE(joint)){times <- times + 1}
      for (counter in seq(1, times)){
        if (!(is.na(individual.titles[counter]))){
          p[[counter]]$labels$title <- individual.titles[counter]
        }
      }
    }

    return(p)
}
