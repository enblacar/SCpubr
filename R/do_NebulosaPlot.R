#' Wrapper for Nebulosa::plot_density in Seurat.
#'
#' @inheritParams doc_function
#' @inheritParams Nebulosa::plot_density
#' @param combine \strong{\code{\link[base]{logical}}} | Whether to create a single plot out of multiple features.
#' @param joint \strong{\code{\link[base]{logical}}} | Whether to plot different features as joint density.
#' @param return_only_joint \strong{\code{\link[base]{logical}}} | Whether to only return the joint density panel.
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
                            return_only_joint = FALSE,
                            plot.title = NULL,
                            plot.subtitle = NULL,
                            plot.caption = NULL,
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
                            viridis_color_map = "G",
                            viridis_direction = 1,
                            verbose = TRUE,
                            na.value = "grey75",
                            plot.axes = FALSE){
  # Check if the sample provided is a Seurat object.
  check_Seurat(sample = sample)

  # Checks for packages.
  check_suggests(function_name = "do_NebulosaPlot")
  # Check the reduction.
  reduction <- check_and_set_reduction(sample = sample, reduction = reduction)
  # Check the dimensions.
  dimensions <- check_and_set_dimensions(sample = sample, reduction = reduction, dims = dims)
  # Check logical parameters.
  logical_list <- list("combine" = combine,
                       "joint" = joint,
                       "return_only_joint" = return_only_joint,
                       "plot_cell_borders" = plot_cell_borders,
                       "plot.axes" = plot.axes)
  check_type(parameters = logical_list, required_type = "logical", test_function = is.logical)
  # Check numeric parameters.
  numeric_list <- list("pt.size" = pt.size,
                       "legend.framewidth" = legend.framewidth,
                       "legend.tickwidth" = legend.tickwidth,
                       "legend.length" = legend.length,
                       "legend.width" = legend.width,
                       "dims" = dims,
                       "border.size" = border.size,
                       "viridis_direction" = viridis_direction)
  check_type(parameters = numeric_list, required_type = "numeric", test_function = is.numeric)
  # Check character parameters.
  if (is.list(features)){
    warning("Features provided as a list. Unlisting the list. Please use a character vector next time.", call. = FALSE)
    features <- unique(unlist(features))
  }
  character_list <- list("legend.position" = legend.position,
                         "features" = features,
                         "method" = method,
                         "plot.title" = plot.title,
                         "plot.subtitle" = plot.subtitle,
                         "plot.caption" = plot.caption,
                         "slot" = slot,
                         "legend.framecolor" = legend.framecolor,
                         "legend.tickcolor" = legend.tickcolor,
                         "legend.type" = legend.type,
                         "font.type" = font.type,
                         "border.color" = border.color,
                         "na.value" = na.value)
  check_type(parameters = character_list, required_type = "character", test_function = is.character)
  # Check slot.
  slot <- check_and_set_slot(slot = slot)

  # Check if the feature is actually in the object.
  features <- check_feature(sample = sample, features = features, permissive = TRUE)
  features <- remove_duplicated_features(features = features)

  # Check the colors provided to legend.framecolor and legend.tickcolor and border.color.
  check_colors(legend.framecolor, parameter_name = "legend.framecolor")
  check_colors(legend.tickcolor, parameter_name = "legend.tickcolor")
  check_colors(border.color, parameter_name = "border.color")
  check_colors(na.value, parameter_name = "na.value")

  check_parameters(parameter = font.type, parameter_name = "font.type")
  check_parameters(parameter = legend.type, parameter_name = "legend.type")
  check_parameters(parameter = legend.position, parameter_name = "legend.position")
  check_parameters(parameter = viridis_direction, parameter_name = "viridis_direction")
  check_parameters(parameter = viridis_color_map, parameter_name = "viridis_color_map")

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
         ggplot2::theme(plot.title = ggplot2::element_text(face = "bold", hjust = 0),
                        plot.subtitle = ggplot2::element_text(hjust = 0),
                        plot.caption = ggplot2::element_text(hjust = 1),
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
                   function_use = ggplot2::scale_color_viridis_c(na.value = na.value,
                                                                 option = viridis_color_map,
                                                                 direction = viridis_direction))

    for (plot_num in seq(1:num_plots)){
      # Set size of dots.
      p[[plot_num]]$layers[[1]]$aes_params$size <- pt.size

      if (legend.position != "none"){
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
    }

    # For embeddings that are umap of tsne, we remove all axes..
    if (reduction %in% c("umap", "tsne")){
      # if dims is first and then second.
      if (sum(dims == c(1, 2)) == 2){
        p <- p &
          ggplot2::theme(axis.title = if (isFALSE(plot.axes)){ggplot2::element_blank()} else {ggplot2::element_text(color = "black", face = "bold", hjust = 0.5)},
                         axis.text = if (isFALSE(plot.axes)){ggplot2::element_blank()} else {ggplot2::element_text(color = "black")},
                         axis.ticks = if (isFALSE(plot.axes)){ggplot2::element_blank()} else {ggplot2::element_line(color = "black")},
                         axis.line = if (isFALSE(plot.axes)){ggplot2::element_blank()} else {ggplot2::element_line(color = "black")})
      } else {
        labels <- colnames(sample@reductions[[reduction]][[]])[dims]
        p <- p &
          ggplot2::theme(axis.text = if (isFALSE(plot.axes)){ggplot2::element_blank()} else {ggplot2::element_text(color = "black")},
                         axis.ticks = if (isFALSE(plot.axes)){ggplot2::element_blank()} else {ggplot2::element_line(color = "black")},
                         axis.line = if (isFALSE(plot.axes)){ggplot2::element_blank()} else {ggplot2::element_line(color = "black")},
                         axis.title = ggplot2::element_text(face = "bold", hjust = 0.5, color = "black")) &
          ggplot2::xlab(labels[1]) &
          ggplot2::ylab(labels[2])
      }
      # For diffusion maps, we do want to keep at least the axis titles so that we know which DC are we plotting.
    } else {
      labels <- colnames(sample@reductions[[reduction]][[]])[dims]
      p <- p &
        ggplot2::theme(axis.text = if (isFALSE(plot.axes)){ggplot2::element_blank()} else {ggplot2::element_text(color = "black")},
                       axis.ticks = if (isFALSE(plot.axes)){ggplot2::element_blank()} else {ggplot2::element_line(color = "black")},
                       axis.line = if (isFALSE(plot.axes)){ggplot2::element_blank()} else {ggplot2::element_line(color = "black")},
                       axis.title = ggplot2::element_text(face = "bold", hjust = 0.5, color = "black")) &
        ggplot2::xlab(labels[1]) &
        ggplot2::ylab(labels[2])
    }

    # Further patch for diffusion maps.
    if (reduction == "diffusion"){
      labels <- colnames(sample@reductions[["diffusion"]][[]])[dims]
      # Fix the axis scale so that the highest and lowest values are in the range of the DCs (previously was around +-1.5, while DCs might range to +-0.004 or so).
      p <-  suppressMessages({
        p &
          ggplot2::xlim(c(min(sample@reductions$diffusion[[]][, labels[1]]),
                          max(sample@reductions$diffusion[[]][, labels[1]]))) &
          ggplot2::ylim(c(min(sample@reductions$diffusion[[]][, labels[2]]),
                          max(sample@reductions$diffusion[[]][, labels[2]]))) &
          # Remove axis elements so that the axis title is the only thing left.
          ggplot2::theme(axis.text = if (isFALSE(plot.axes)){ggplot2::element_blank()} else {ggplot2::element_text(color = "black")},
                         axis.ticks = if (isFALSE(plot.axes)){ggplot2::element_blank()} else {ggplot2::element_line(color = "black")},
                         axis.line = if (isFALSE(plot.axes)){ggplot2::element_blank()} else {ggplot2::element_line(color = "black")},
                         axis.title = ggplot2::element_text(face = "bold", hjust = 0.5, color = "black"))
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

    return(p)
}
