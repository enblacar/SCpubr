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
                            legend.framewidth = 0.5,
                            legend.tickwidth = 0.5,
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
                            viridis.palette = "G",
                            viridis.direction = 1,
                            verbose = TRUE,
                            na.value = "grey75",
                            plot.axes = FALSE,
                            number.breaks = 5,
                            use_viridis = FALSE,
                            sequential.palette = "YlGnBu",
                            sequential.direction = 1,
                            plot.title.face = "bold",
                            plot.subtitle.face = "plain",
                            plot.caption.face = "italic",
                            axis.title.face = "bold",
                            axis.text.face = "plain",
                            legend.title.face = "bold",
                            legend.text.face = "plain"){
  # Add lengthy error messages.
  withr::local_options(.new = list("warning.length" = 8170))
  
  `%>%` <- magrittr::`%>%`
  check_suggests(function_name = "do_NebulosaPlot")
  # Check if the sample provided is a Seurat object.
  check_Seurat(sample = sample)

  # Check the reduction.
  reduction <- check_and_set_reduction(sample = sample, reduction = reduction)
  # Check the dimensions.
  dims <- check_and_set_dimensions(sample = sample, reduction = reduction, dims = dims)
  # Check logical parameters.
  logical_list <- list("combine" = combine,
                       "joint" = joint,
                       "return_only_joint" = return_only_joint,
                       "plot_cell_borders" = plot_cell_borders,
                       "plot.axes" = plot.axes,
                       "use_viridis" = use_viridis)
  check_type(parameters = logical_list, required_type = "logical", test_function = is.logical)
  # Check numeric parameters.
  numeric_list <- list("pt.size" = pt.size,
                       "legend.framewidth" = legend.framewidth,
                       "legend.tickwidth" = legend.tickwidth,
                       "legend.length" = legend.length,
                       "legend.width" = legend.width,
                       "dims" = dims,
                       "border.size" = border.size,
                       "viridis.direction" = viridis.direction,
                       "number.breaks" = number.breaks,
                       "sequential.direction" = sequential.direction)
  check_type(parameters = numeric_list, required_type = "numeric", test_function = is.numeric)
  # Check character parameters.
  if (is.list(features)){
    warning(paste0(add_warning(), crayon_key("Features"),
                   crayon_body(" provided as a "),
                   crayon_key("list"),
                   crayon_body(". Unlisting it. Please use a "),
                   crayon_key("character vector")), call. = FALSE)
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
                         "na.value" = na.value,
                         "sequential.palette" = sequential.palette,
                         "plot.title.face" = plot.title.face,
                         "plot.subtitle.face" = plot.subtitle.face,
                         "plot.caption.face" = plot.caption.face,
                         "axis.title.face" = axis.title.face,
                         "axis.text.face" = axis.text.face,
                         "legend.title.face" = legend.title.face,
                         "legend.text.face" = legend.text.face)
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
  check_parameters(parameter = viridis.palette, parameter_name = "viridis.palette")
  check_parameters(parameter = number.breaks, parameter_name = "number.breaks")
  check_parameters(plot.title.face, parameter_name = "plot.title.face")
  check_parameters(plot.subtitle.face, parameter_name = "plot.subtitle.face")
  check_parameters(plot.caption.face, parameter_name = "plot.caption.face")
  check_parameters(axis.title.face, parameter_name = "axis.title.face")
  check_parameters(axis.text.face, parameter_name = "axis.text.face")
  check_parameters(legend.title.face, parameter_name = "legend.title.face")
  check_parameters(legend.text.face, parameter_name = "legend.text.face")
  check_parameters(viridis.direction, parameter_name = "viridis.direction")
  check_parameters(sequential.direction, parameter_name = "sequential.direction")
  
  colors.gradient <- compute_continuous_palette(name = ifelse(isTRUE(use_viridis), viridis.palette, sequential.palette),
                                                use_viridis = use_viridis,
                                                direction = ifelse(isTRUE(use_viridis), viridis.direction, sequential.direction),
                                                enforce_symmetry = FALSE)
  
  # Plot a density plot using Nebulosa package.
    p <- Nebulosa::plot_density(object = sample,
                                features = features,
                                joint = joint,
                                reduction = reduction,
                                dims = dims) &
         ggplot2::theme_minimal(base_size = font.size) &
         ggplot2::theme(plot.title = ggplot2::element_text(face = plot.title.face, hjust = 0),
                        plot.subtitle = ggplot2::element_text(face = plot.subtitle.face, hjust = 0),
                        plot.caption = ggplot2::element_text(face = plot.caption.face, hjust = 1),
                        legend.text = ggplot2::element_text(face = legend.text.face),
                        legend.title = ggplot2::element_text(face = legend.title.face),
                        plot.title.position = "plot",
                        panel.grid = ggplot2::element_blank(),
                        text = ggplot2::element_text(family = font.type),
                        plot.caption.position = "plot",
                        legend.position = legend.position,
                        legend.justification = "center",
                        plot.margin = ggplot2::margin(t = 0, r = 0, b = 0, l = 0),
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
    for (counter in seq_len(num_plots)){
      if (counter > length(features)){
        name.use <-  "Joint Density"
      } else {
        name.use <- "Density"
      }
      if (num_plots == 1){
        limits <- c(p$data[, "feature", drop = FALSE] %>% dplyr::arrange(.data$feature) %>% utils::head(1) %>% dplyr::pull(.data$feature),
                    p$data[, "feature", drop = FALSE] %>% dplyr::arrange(.data$feature) %>% utils::tail(1) %>% dplyr::pull(.data$feature))

        scale.setup <- compute_scales(sample = sample,
                                      feature = features,
                                      assay = NULL,
                                      reduction = NULL,
                                      slot = slot,
                                      number.breaks = number.breaks,
                                      min.cutoff = NA,
                                      max.cutoff = NA,
                                      flavor = "Seurat",
                                      enforce_symmetry = FALSE,
                                      from_data = TRUE,
                                      limits.use = limits)

        p <- add_scale(p = p,
                       function_use = ggplot2::scale_color_gradientn(colors = colors.gradient,
                                                                     na.value = na.value,
                                                                     name = name.use,
                                                                     breaks = scale.setup$breaks,
                                                                     labels = scale.setup$labels,
                                                                     limits = scale.setup$limits),
                       scale = "color")

      } else {
        limits <- c(p[[counter]]$data[, "feature", drop = FALSE] %>% dplyr::arrange(.data$feature) %>% utils::head(1) %>% dplyr::pull(.data$feature),
                    p[[counter]]$data[, "feature", drop = FALSE] %>% dplyr::arrange(.data$feature) %>% utils::tail(1) %>% dplyr::pull(.data$feature))

        scale.setup <- compute_scales(sample = sample,
                                      feature = features,
                                      assay = NULL,
                                      reduction = NULL,
                                      slot = slot,
                                      number.breaks = number.breaks,
                                      min.cutoff = NA,
                                      max.cutoff = NA,
                                      flavor = "Seurat",
                                      enforce_symmetry = FALSE,
                                      from_data = TRUE,
                                      limits.use = limits)
        p[[counter]] <- add_scale(p =  p[[counter]],
                                  function_use = ggplot2::scale_color_gradientn(colors = colors.gradient,
                                                                                na.value = na.value,
                                                                                name = name.use,
                                                                                breaks = scale.setup$breaks,
                                                                                labels = scale.setup$labels,
                                                                                limits = scale.setup$limits),
                                  scale = "color")
      }


      # Set size of dots.
      p[[counter]]$layers[[1]]$aes_params$size <- pt.size

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
          p[[counter]] <- modify_continuous_legend(p = p[[counter]],
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


    if (base::isFALSE(plot.axes)){
      p <- p &
           ggplot2::theme(axis.title = ggplot2::element_blank(),
                          axis.text = ggplot2::element_blank(),
                          axis.ticks = ggplot2::element_blank(),
                          axis.line = ggplot2::element_blank())
    } else {
      p <- p &
           ggplot2::theme(axis.title = ggplot2::element_text(face = axis.title.face),
                          axis.text = ggplot2::element_text(face = axis.text.face),
                          axis.ticks = ggplot2::element_line(color = "black"),
                          axis.line = ggplot2::element_line(color = "black"))
    }

    # Further patch for diffusion maps.
    if (stringr::str_starts(reduction, "diff|DIFF")){
      labels <- colnames(sample@reductions[["diffusion"]][[]])[dims]
      # Fix the axis scale so that the highest and lowest values are in the range of the DCs (previously was around +-1.5, while DCs might range to +-0.004 or so).
      p <-  suppressMessages({
        p &
          ggplot2::xlim(c(min(sample@reductions$diffusion[[]][, labels[1]]),
                          max(sample@reductions$diffusion[[]][, labels[1]]))) &
          ggplot2::ylim(c(min(sample@reductions$diffusion[[]][, labels[2]]),
                          max(sample@reductions$diffusion[[]][, labels[2]])))
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
