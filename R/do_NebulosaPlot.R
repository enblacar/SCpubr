#' Wrapper for Nebulosa::plot_density in Seurat.
#'
#'
#' @param sample Seurat object.
#' @param reduction Reduction to use. Can be the canonical ones such as "umap", "pca", or any custom ones, such as "diffusion". If you are unsure about which reductions you have, use `Seurat::Reductions(sample)`.
#' @param features Features to plot density for. It can be a single one or a vector of multiple features.
#' @param slot Slot to retrieve the data from. Defaults to "data".
#' @param pt.size Size of the dots.
#' @param combine Whether to create a single plot out of multiple features.
#' @param method Kernel density estimation method. Either "ks" or "wkde" or both. See \link[Nebulosa]{plot_density} for more details.
#' @param dims Vector of 2 dims to plot the data. By default, first and second from the specified reduction.
#' @param joint Whether to plot different features as joint density.
#' @param return_only_joint Whether to only return the joint density panel. Logical.
#' @param shape Shape of the geometry (ggplot number).
#' @param legend Whether to plot the legend or not. Logical.
#' @param legend.type Character. Type of legend to display. One of: normal, colorbar, colorsteps.
#' @param legend.position Position of the legend in the plot.
#' @param legend.framewidth,legend.tickwidth Width of the lines of the box in the legend.
#' @param legend.framecolor,legend.tickcolor Color of the lines of the box in the legend.
#' @param legend.length,legend.width Length and width of the legend. Will adjust automatically depending on legend side.
#' @param font.size Base font.size of the plot.
#' @param font.type Character. Base font for the plot. One of mono, serif or sans.
#' @param plot.title,plot.subtitle,plot.caption Title to use in the plot.
#' @param individual.titles,individual.subtitles,individual.captions Titles, subtitles and captions for each feature if needed. Either NULL or a vector of equal length of features.
#' @param viridis_color_map Character. A capital letter from A to H or the scale name as in \link[viridis]{scale_fill_viridis}.
#' @param plot_cell_borders Logical. Whether to plot border around cells.
#' @param border.size Numeric. Width of the border of the cells.
#' @param border.color Character. Color to use for the border of the cells.
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
                             pt.size = 1,
                             reduction = NULL,
                             combine = TRUE,
                             method = c("ks", "wkde"),
                             joint = FALSE,
                             return_only_joint = NULL,
                             plot.title = NULL,
                             plot.subtitle = NULL,
                             plot.caption = NULL,
                             individual.titles = NULL,
                             individual.subtitles = NULL,
                             individual.captions = NULL,
                             shape = 16,
                             legend = TRUE,
                             legend.type = "colorbar",
                             legend.framewidth = 1.5,
                             legend.tickwidth = 1.5,
                             legend.length = 20,
                             legend.width = 1,
                             legend.framecolor = "grey50",
                             legend.tickcolor = "white",
                             font.size = 14,
                             font.type = "sans",
                             legend.position = "bottom",
                             plot_cell_borders = TRUE,
                             border.size = 2,
                             border.color = "black",
                             viridis_color_map = "D",
                             verbose = TRUE){
  # Check if the sample provided is a Seurat object.
  check_Seurat(sample = sample)

  # Checks for packages.
  check_suggests(function_name = "do_NebulosaPlot")
  # Check the reduction.
  reduction <- check_and_set_reduction(sample = sample, reduction = reduction)
  # Check the dimensions.
  dimensions <- check_and_set_dimensions(sample = sample, reduction = reduction, dims = dims)
  # Check logical parameters.
  logical_list <- list("legend" = legend,
                       "combine" = combine,
                       "joint" = joint,
                       "return_only_joint" = return_only_joint,
                       "plot_cell_borders" = plot_cell_borders)
  check_type(parameters = logical_list, required_type = "logical", test_function = is.logical)
  # Check numeric parameters.
  numeric_list <- list("pt.size" = pt.size,
                       "shape" = shape,
                       "legend.framewidth" = legend.framewidth,
                       "legend.tickwidth" = legend.tickwidth,
                       "legend.length" = legend.length,
                       "legend.width" = legend.width,
                       "dims" = dims,
                       "border.size" = border.size)
  check_type(parameters = numeric_list, required_type = "numeric", test_function = is.numeric)
  # Check character parameters.
  if (is.list(features)){
    warning("Features provided as a list. Unlisting the list. Please use a character vector next time.", call. = F)
    features <- unique(unlist(features))
  }
  character_list <- list("legend.position" = legend.position,
                         "features" = features,
                         "method" = method,
                         "plot.title" = plot.title,
                         "plot.subtitle" = plot.subtitle,
                         "plot.caption" = plot.caption,
                         "slot" = slot,
                         "individual.titles" = individual.titles,
                         "legend.framecolor" = legend.framecolor,
                         "legend.tickcolor" = legend.tickcolor,
                         "legend.type" = legend.type,
                         "font.type" = font.type,
                         "border.color" = border.color)
  check_type(parameters = character_list, required_type = "character", test_function = is.character)
  # Check slot.
  slot <- check_and_set_slot(slot = slot)

  # Check if the feature is actually in the object.
  features <- check_feature(sample = sample, features = features, permissive = TRUE)
  features <- remove_duplicated_features(features = features)

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

  # Check the colors provided to legend.framecolor and legend.tickcolor and border.color.
  check_colors(legend.framecolor, parameter_name = "legend.framecolor")
  check_colors(legend.tickcolor, parameter_name = "legend.tickcolor")
  check_colors(border.color, parameter_name = "border.color")

  # Check font.type.
  if (!(font.type %in% c("sans", "serif", "mono"))){
    stop("Please select one of the following for font.type: sans, serif, mono.", call. = F)
  }

  # Check the legend.type.
  if (!(legend.type %in% c("normal", "colorbar", "colorsteps"))){
    stop("Please select one of the following for legend.type: normal, colorbar, colorsteps.", call. = FALSE)
  }

  # Check the legend.position.
  if (!(legend.position %in% c("top", "bottom", "left", "right"))){
    stop("Please select one of the following for legend.position: top, bottom, left, right.", call. = FALSE)
  }

  # Define legend parameters.
  if (legend.position %in% c("top", "bottom")){
    legend.barwidth <- legend.length
    legend.barheight <- legend.width
  } else if (legend.position %in% c("left", "right")){
    legend.barwidth <- legend.width
    legend.barheight <- legend.length
  }

  # Plot a density plot using Nebulosa package.
    p <- Nebulosa::plot_density(object = sample,
                                features = features,
                                joint = joint,
                                reduction = reduction,
                                dims = dims) &
         ggplot2::theme_minimal(base_size = font.size) &
         ggplot2::theme(axis.title = ggplot2::element_blank(),
                        axis.text = ggplot2::element_blank(),
                        plot.title = ggtext::element_markdown(face = "bold", hjust = 0),
                        plot.subtitle = ggtext::element_markdown(hjust = 0),
                        plot.caption = ggtext::element_markdown(hjust = 1),
                        plot.title.position = "plot",
                        panel.grid = ggplot2::element_blank(),
                        text = ggplot2::element_text(family = font.type),
                        plot.caption.position = "plot",
                        legend.text = ggplot2::element_text(face = "bold"),
                        legend.position = legend.position,
                        legend.title = ggplot2::element_text(face = "bold"),
                        legend.justification = "center",
                        plot.margin = ggplot2::margin(t = 10, r = 10, b = 10, l = 10),
                        panel.grid.major = ggplot2::element_blank(),
                        plot.background = ggplot2::element_rect(fill = "white", color = "white"),
                        panel.background = ggplot2::element_rect(fill = "white", color = "white"),
                        legend.background = ggplot2::element_rect(fill = "white", color = "white"))

    # Compute the total number of plots according to whether joint is set to TRUE or not.
    if (isTRUE(joint)){
      num_plots <- length(features) + 1
    } else {
      num_plots <- length(features)
    }
    p <- add_scale(p = p,
                   num_plots = num_plots,
                   scale = "color",
                   function_use = ggplot2::scale_color_viridis_c(na.value = "grey75",
                                                                 option = viridis_color_map))

    for (plot_num in seq(1:num_plots)){
      # Set size of dots.
      p[[plot_num]]$layers[[1]]$aes_params$size <- pt.size

      if (num_plots == 1){
        p <- modify_continuous_legend(p = p,
                                      legend.aes = "color",
                                      legend.type = legend.type,
                                      legend.position = legend.position,
                                      legend.length = legend.length,
                                      legend.width = legend.width,
                                      legend.framecolor = legend.framecolor,
                                      legend.tickcolor = legend.tickcolor,
                                      legend.framewidth = legend.framewidth,
                                      legend.tickwidth = legend.tickwidth)
      } else {
        p[[plot_num]] <- modify_continuous_legend(p = p[[plot_num]],
                                                  legend.aes = "color",
                                                  legend.type = legend.type,
                                                  legend.position = legend.position,
                                                  legend.length = legend.length,
                                                  legend.width = legend.width,
                                                  legend.framecolor = legend.framecolor,
                                                  legend.tickcolor = legend.tickcolor,
                                                  legend.framewidth = legend.framewidth,
                                                  legend.tickwidth = legend.tickwidth)
      }

    }

    # For embeddings that are umap of tsne, we remove all axes..
    if (reduction %in% c("umap", "tsne")){
      # if dims is first and then second.
      if (sum(dims == c(1, 2)) == 2){
        p <- p &
          ggplot2::theme(axis.title.x = ggplot2::element_blank(),
                         axis.title.y = ggplot2::element_blank(),
                         axis.text = ggplot2::element_blank(),
                         axis.ticks = ggplot2::element_blank(),
                         axis.line = ggplot2::element_blank())
      } else {
        labels <- colnames(sample@reductions[[reduction]][[]])[dims]
        p <- p &
          ggplot2::theme(axis.text = ggplot2::element_blank(),
                         axis.ticks = ggplot2::element_blank(),
                         axis.line = ggplot2::element_blank(),
                         axis.title = ggplot2::element_text(face = "bold", hjust = 0.5)) &
          ggplot2::xlab(labels[1]) &
          ggplot2::ylab(labels[2])
      }
      # For diffusion maps, we do want to keep at least the axis titles so that we know which DC are we plotting.
    } else {
      labels <- colnames(sample@reductions[[reduction]][[]])[dims]
      p <- p &
        ggplot2::theme(axis.text = ggplot2::element_blank(),
                       axis.ticks = ggplot2::element_blank(),
                       axis.line = ggplot2::element_blank(),
                       axis.title = ggplot2::element_text(face = "bold", hjust = 0.5)) &
        ggplot2::xlab(labels[1]) &
        ggplot2::ylab(labels[2])
    }

    # Further patch for diffusion maps.
    if (reduction == "diffusion"){
      # Fix the axis scale so that the highest and lowest values are in the range of the DCs (previously was around +-1.5, while DCs might range to +-0.004 or so).
      p <-  suppressMessages({
        p &
          ggplot2::xlim(c(min(sample@reductions$diffusion[[]][, paste0("DC_", dims[1])]),
                          max(sample@reductions$diffusion[[]][, paste0("DC_", dims[1])]))) &
          ggplot2::ylim(c(min(sample@reductions$diffusion[[]][, paste0("DC_", dims[2])]),
                          max(sample@reductions$diffusion[[]][, paste0("DC_", dims[2])]))) &
          # Remove axis elements so that the axis title is the only thing left.
          ggplot2::theme(axis.text = ggplot2::element_blank(),
                         axis.ticks = ggplot2::element_blank(),
                         axis.line = ggplot2::element_blank(),
                         axis.title = ggplot2::element_text(face = "bold", hjust = 0.5))
      })

    }

    # Add cell borders.

    if (isTRUE(plot_cell_borders)){
      # Generate base layer.
      labels <- colnames(sample@reductions[[reduction]][[]])[dims]
      df <- data.frame(x = Seurat::Embeddings(sample, reduction = reduction)[, labels[1]],
                       y = Seurat::Embeddings(sample, reduction = reduction)[, labels[2]])

      base_layer <- ggplot2::geom_point(data = df, mapping = ggplot2::aes(x = .data$x,
                                                                          y = .data$y),
                                        colour = border.color,
                                        size = pt.size * border.size,
                                        show.legend = FALSE)
      if (num_plots > 1){
        for (plot_num in seq(1, num_plots)){
          # Add cell borders.
          if (isTRUE(plot_cell_borders)){
            p[[plot_num]]$layers <- append(base_layer, p[[plot_num]]$layers)
          }
        }
      } else {
        p$layers <- append(base_layer, p$layers)
      }

    }

    # Remove legend.
    if (legend == FALSE){
      p <- p +
           ggplot2::theme(legend.position = "none")
    }

    if (isTRUE(return_only_joint)){
      p <- p[[length(features) + 1]]
    }

    # Add a title.
    if (!(is.null(plot.title))){
        if (length(features) == 1 | (isTRUE(return_only_joint))){
          p <- p +
               ggplot2::labs(title = plot.title)
        } else {
            p <- p & patchwork::plot_annotation(title = plot.title)
        }
    }

    # Add a subtitle.
    if (!(is.null(plot.subtitle))){
      if (length(features) == 1 | (isTRUE(return_only_joint))){
        p <- p &
             ggplot2::labs(subtitle = plot.subtitle)
      } else {
        p <- p & patchwork::plot_annotation(subtitle = plot.subtitle)
      }
    }

    # Add a caption
    if (!(is.null(plot.caption))){
      if (length(features) == 1 | (isTRUE(return_only_joint))){
        p <- p &
             ggplot2::labs(caption = plot.caption)
      } else {
        p <- p & patchwork::plot_annotation(caption = plot.caption)
      }
    }


    # Add individual titles.
    if (!is.null(individual.titles) & !(isTRUE(return_only_joint))){
      times <- length(features)
      if (isTRUE(joint)){times <- times + 1}
      for (counter in seq(1, times)){
        if (!(is.na(individual.titles[counter]))){
          p[[counter]]$labels$title <- individual.titles[counter]
        }
      }
    }

    # Add individual subtitles.
    if (!is.null(individual.subtitles) & !(isTRUE(return_only_joint))){
      times <- length(features)
      if (isTRUE(joint)){times <- times + 1}
      for (counter in seq(1, times)){
        if (!(is.na(individual.subtitles[counter]))){
          p[[counter]]$labels$subtitle <- individual.subtitles[counter]
        }
      }
    }

    # Add individual titles.
    if (!is.null(individual.captions) & !(isTRUE(return_only_joint))){
      times <- length(features)
      if (isTRUE(joint)){times <- times + 1}
      for (counter in seq(1, times)){
        if (!(is.na(individual.captions[counter]))){
          p[[counter]]$labels$caption <- individual.captions[counter]
        }
      }
    }

    return(p)
}
