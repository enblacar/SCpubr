#' Wrapper for \link[Seurat]{FeaturePlot}.
#'
#' @inheritParams doc_function
#' @param split.by.idents \strong{\code{\link[base]{character}}} | Vector of identities to plot. The gradient scale will also be subset to only the values of such identities.
#' @param individual.titles,individual.subtitles,individual.captions \strong{\code{\link[base]{character}}} | Titles or subtitles. for each feature if needed. Either NULL or a vector of equal length of features.
#' @param order \strong{\code{\link[base]{logical}}} | Whether to order the cells based on expression.
#'
#' @return  A ggplot2 object containing a Feature Plot.
#' @export
#'
#' @example /man/examples/examples_do_FeaturePlot.R
do_FeaturePlot <- function(sample,
                           features,
                           assay = NULL,
                           reduction = NULL,
                           slot = NULL,
                           order = TRUE,
                           split.by = NULL,
                           split.by.idents = NULL,
                           cells.highlight = NULL,
                           idents.highlight = NULL,
                           dims = c(1, 2),
                           enforce_symmetry = FALSE,
                           pt.size = 1,
                           font.size = 14,
                           font.type = "sans",
                           legend.title = NULL,
                           legend.type = "colorbar",
                           legend.position = "bottom",
                           legend.framewidth = 1.5,
                           legend.tickwidth = 1.5,
                           legend.length = 20,
                           legend.width = 1,
                           legend.framecolor = "grey50",
                           legend.tickcolor = "white",
                           plot.title = NULL,
                           plot.subtitle = NULL,
                           plot.caption = NULL,
                           individual.titles = NULL,
                           individual.subtitles = NULL,
                           individual.captions = NULL,
                           ncol = NULL,
                           viridis_color_map = "G",
                           viridis_direction = 1,
                           raster = FALSE,
                           raster.dpi = 1024,
                           plot_cell_borders = TRUE,
                           border.size = 2,
                           border.color = "black",
                           na.value = "grey75",
                           verbose = TRUE,
                           plot.axes = FALSE){

  check_suggests(function_name = "do_FeaturePlot")
  # Check if the sample provided is a Seurat object.
  check_Seurat(sample = sample)
  # Check the assay.
  out <- check_and_set_assay(sample = sample, assay = assay)
  sample <- out[["sample"]]
  assay <- out[["assay"]]
  # Check the reduction.
  reduction <- check_and_set_reduction(sample = sample, reduction = reduction)
  # Check the dimensions.
  dimensions <- check_and_set_dimensions(sample = sample, reduction = reduction, dims = dims)
  # Check logical parameters.
  logical_list <- list("verbose" = verbose,
                       "raster" = raster,
                       "plot_cell_borders" = plot_cell_borders,
                       "order" = order,
                       "enforce_symmetry" = enforce_symmetry,
                       "plot.axes" = plot.axes)
  check_type(parameters = logical_list, required_type = "logical", test_function = is.logical)
  # Check numeric parameters.
  numeric_list <- list("pt.size" = pt.size,
                       "ncol" = ncol,
                       "font.size" = font.size,
                       "raster.dpi" = raster.dpi,
                       "legend.framewidth" = legend.framewidth,
                       "legend.tickwidth" = legend.tickwidth,
                       "legend.length" = legend.length,
                       "legend.width" = legend.width,
                       "border.size" = border.size,
                       "viridis_direction" = viridis_direction)
  check_type(parameters = numeric_list, required_type = "numeric", test_function = is.numeric)
  # Check character parameters.
  # Workaround for features.
  if (is.list(features)){
    warning("Features provided as a list. Unlisting the list. Please use a character vector next time.", call. = FALSE)
    features <- unique(unlist(features))
  }
  character_list <- list("legend.position" = legend.position,
                         "features" = features,
                         "cells.highlight" = cells.highlight,
                         "idents.highlight" = idents.highlight,
                         "slot" = slot,
                         "split.by" = split.by,
                         "plot.title" = plot.title,
                         "plot.subtitle" = plot.subtitle,
                         "plot.caption" = plot.caption,
                         "split.by.idents" = split.by.idents,
                         "viridis_color_map" = viridis_color_map,
                         "individual.titles" = individual.titles,
                         "individual.subtitles" = individual.subtitles,
                         "individual.captions" = individual.captions,
                         "legend.framecolor" = legend.framecolor,
                         "legend.tickcolor" = legend.tickcolor,
                         "legend.type" = legend.type,
                         "font.type" = font.type,
                         "border.color" = border.color,
                         "legend.title" = legend.title,
                         "na.value" = na.value)
  check_type(parameters = character_list, required_type = "character", test_function = is.character)

  # Check slot.
  slot <- check_and_set_slot(slot = slot)

  # Check split.by is on metadata.
  if (!(is.null(split.by))){check_feature(sample = sample, features = split.by, enforce_check = "metadata", enforce_parameter = "split.by")}

  # Check individual titles.
  if (length(features) > 1 & !is.null(individual.titles)){
    assertthat::assert_that(length(features) == length(individual.titles),
                            msg = 'Total number of individual titles does not match the number of features provided.')
  }

  if (length(features) > 1 & !is.null(individual.subtitles)){
    assertthat::assert_that(length(features) == length(individual.subtitles),
                            msg = 'Total number of individual subtitles does not match the number of features provided.')
  }

  if (length(features) > 1 & !is.null(individual.captions)){
    assertthat::assert_that(length(features) == length(individual.captions),
                            msg = 'Total number of individual captions does not match the number of features provided.')
  }

  check_colors(border.color, parameter_name = "border.color")
  check_colors(na.value, parameter_name = "na.value")
  check_colors(legend.framecolor, parameter_name = "legend.framecolor")
  check_colors(legend.tickcolor, parameter_name = "legend.tickcolor")

  check_parameters(parameter = font.type, parameter_name = "font.type")
  check_parameters(parameter = legend.type, parameter_name = "legend.type")
  check_parameters(parameter = legend.position, parameter_name = "legend.position")
  check_parameters(parameter = viridis_direction, parameter_name = "viridis_direction")
  check_parameters(parameter = viridis_color_map, parameter_name = "viridis_color_map")


  # Define legend parameters. Width and height values will change depending on the legend orientation.
  if (legend.position %in% c("top", "bottom")){
    legend.barwidth <- legend.length
    legend.barheight <- legend.width
  } else if (legend.position %in% c("left", "right")){
    legend.barwidth <- legend.width
    legend.barheight <- legend.length
  }

  # Check for raster and pt.size.
  if (isTRUE(raster) & pt.size < 1){
    warning("Setting raster = TRUE and pt.size < 1 will result in the cells being ploted as a cross. This behaviour can not be modified, but setting pt.size to 1 or higher solves it. For Feature plots, optimized values would be pt.size = 3 and raster.dpi = 2048.", call. = FALSE)
  }


  # Generate base layer.
  if (isTRUE(plot_cell_borders)){
    labels <- colnames(sample@reductions[[reduction]][[]])[dims]
    df <- data.frame(x = Seurat::Embeddings(sample, reduction = reduction)[, labels[1]],
                     y = Seurat::Embeddings(sample, reduction = reduction)[, labels[2]])

    if (isFALSE(raster)){
      base_layer <- ggplot2::geom_point(data = df, mapping = ggplot2::aes(x = .data$x,
                                                                          y = .data$y),
                                        colour = border.color,
                                        size = pt.size * border.size,
                                        show.legend = FALSE)
    } else if (isTRUE(raster)){
      base_layer <- scattermore::geom_scattermore(data = df,
                                                  mapping = ggplot2::aes(x = .data$x,
                                                                         y = .data$y),
                                                  color = border.color,
                                                  size = pt.size * border.size,
                                                  stroke = pt.size / 2,
                                                  show.legend = FALSE,
                                                  pointsize = pt.size * border.size,
                                                  pixels = c(raster.dpi, raster.dpi))
    }
  }

  # CONDITION: Regular FeaturePlot, under default parameters.
  default_parameters <- is.null(split.by) & is.null(cells.highlight) & is.null(idents.highlight)
  # If only default parameters are used.
  if (default_parameters){
    # Check if the feature is actually in the object.
    features <- check_feature(sample = sample,
                              features = features,
                              permissive = TRUE)
    # Remove duplicated features.
    features <- remove_duplicated_features(features = features)
    p <- Seurat::FeaturePlot(sample,
                             features,
                             slot = slot,
                             reduction = reduction,
                             order = order,
                             dims = dims,
                             pt.size = pt.size,
                             ncol = ncol,
                             raster = raster,
                             raster.dpi = c(raster.dpi, raster.dpi)) &
      # Remove Seurat::FeaturePlot() default plot title.
      ggplot2::ggtitle("")

    # Add color scales.
    num_plots <- length(features)
    for (counter in seq(1, num_plots)){
      if (num_plots == 1){
        if (isFALSE(enforce_symmetry)){
          p <- add_scale(p = p,
                         function_use = ggplot2::scale_color_viridis_c(na.value = na.value,
                                                                       option = viridis_color_map,
                                                                       direction = viridis_direction),
                         scale = "color")
        } else if (isTRUE(enforce_symmetry)){
          p.build <- ggplot2::ggplot_build(p)
          feature.select <- gsub("-", ".", features)
          limits <- c(min(p.build$plot$data[, feature.select]),
                      max(p.build$plot$data[, feature.select]))
          end_value <- max(abs(limits))
          p <- add_scale(p = p,
                         function_use = ggplot2::scale_color_gradientn(colors = c("#033270", "#4091C9", "#fdf0d5", "#c94040", "#65010C"),
                                                                       limits = c(-end_value, end_value),
                                                                       na.value = na.value),
                         scale = "color")
        }
      } else if (num_plots > 1){
        if (isFALSE(enforce_symmetry)){
          p[[counter]] <- add_scale(p = p[[counter]],
                                    function_use = ggplot2::scale_color_viridis_c(na.value = na.value,
                                                                                  option = viridis_color_map,
                                                                                  direction = viridis_direction),
                                    scale = "color")
        } else if (isTRUE(enforce_symmetry)){
          p.build <- ggplot2::ggplot_build(p[[counter]])
          feature.select <- gsub("-", ".",  features[counter])
          limits <- c(min(p.build$plot$data[, feature.select]),
                      max(p.build$plot$data[, feature.select]))
          end_value <- max(abs(limits))
          p[[counter]] <- add_scale(p = p[[counter]],
                                    function_use = ggplot2::scale_color_gradientn(colors = c("#033270", "#4091C9", "#fdf0d5", "#c94040", "#65010C"),
                                                                                  limits = c(-end_value, end_value),
                                                                                  na.value = na.value),
                                    scale = "color")
        }
      }
    }

    # Special patches for diffusion maps: Adding "DC" labels to the axis.
    if (reduction == "diffusion"){
      p <- p &
        ggplot2::xlab(paste0("DC_", dims[1])) &
        ggplot2::ylab(paste0("DC_", dims[2]))
    }

    # Add cell borders.
    if (isTRUE(plot_cell_borders)){
      counter <- 0
      for (feature in features){
        counter <- counter + 1
        p[[counter]]$layers <- append(base_layer, p[[counter]]$layers)
      }
    }

    # Modified FeaturePlot including only a subset of cells.
  } else {
    # Check if the feature is actually in the object.
    output_list <- check_feature(sample = sample,
                                 features = features,
                                 dump_reduction_names = TRUE,
                                 permissive = TRUE)
    features <- output_list[["features"]]
    dim_colnames <- output_list[["reduction_names"]]
    # Remove duplicated features.
    features <- remove_duplicated_features(features = features)
    # Get the subset of wanted cells according to the combination of idents.highlight and cells.highlight parameters.
    if (is.null(idents.highlight) & !(is.null(cells.highlight))){
      # Only if cells.highlight parameters is used.
      cells.use <- cells.highlight
    } else if (!(is.null(idents.highlight)) & is.null(cells.highlight)){
      # Only if idents.highlight parameter is used.
      # Check if the provided identities are part of the active identities in the object.
      check_identity(sample = sample, identities = idents.highlight)
      # Get the names of the cells to use.
      cells.use <- names(Seurat::Idents(sample)[Seurat::Idents(sample) %in% idents.highlight])
    } else if (!(is.null(idents.highlight)) & !(is.null(cells.highlight))){
      # Check if the provided identities are part of the active identities in the object.
      check_identity(sample = sample, identities = idents.highlight)
      # Both idents.highlight and cells.highlight are used.
      cells.1 <- cells.highlight
      cells.2 <- names(Seurat::Idents(sample)[Seurat::Idents(sample) %in% idents.highlight])
      # Get the names of the cells to use.
      cells.use <- unique(c(cells.1, cells.2))
      # If split.by is used.
    } else if (!(is.null(split.by))){
      # No identities selected by the user.
      if (is.null(split.by.idents)){
        cells.use <- colnames(sample)
        # Identitites selected by the user.
      } else {
        # Check if the identitites are duplicated. If so, remove them.
        if (sum(duplicated(split.by.idents)) != 0){
          message("Found and removed duplicated values in split.by.idents.")
          split.by.idents <- split.by.idents[!duplicated(split.by.idents)]
        }
        # Get the names of the cells to plot.
        cells.use <- names(Seurat::Idents(sample)[Seurat::Idents(sample) %in% split.by.idents])
      }
    }
    # Plots are generated independently if more than one feature is provided.
    list.plots <- list()
    # Counter depicting the feature used. Will increase each feature used.
    count_iteration <- 0
    # Iterate over the features.
    for (feature in features){
      # A "dummy" metadata column is generated using the values of the selected feature.
      ## Based on whether the feature is in the metadata of the object.
      if (feature %in% colnames(sample@meta.data)) {
        sample$dummy <- sample@meta.data[, feature]
        ## Or is a gene in the object.
      } else if (feature %in% rownames(sample)){
        sample$dummy <- Seurat::GetAssayData(object = sample, slot = slot)[feature, ]
        ## Or is a dimensional reduction component.
      } else if (feature %in% dim_colnames){
        # Iterate over each dimensional reduction in the object.
        for(red in Seurat::Reductions(object = sample)){
          # If the feature to plot is one of the dimensional reduction components.
          if (feature %in% colnames(sample@reductions[[red]][[]])){
            red.feature <- red
            sample$dummy <- sample@reductions[[red.feature]][[]][, feature]
          }
        }
      }
      # Assign NAs to the values corresponding to the cells  not selected.
      sample$dummy[!(names(sample$dummy) %in% cells.use)] <- NA
      # If split.by is not used.
      if (is.null(split.by)){
        feature.use <- "dummy"
        p.loop <- Seurat::FeaturePlot(sample,
                                      feature.use,
                                      reduction = reduction,
                                      slot = slot,
                                      order = order,
                                      dims = dims,
                                      pt.size = pt.size,
                                      raster = raster,
                                      raster.dpi = c(raster.dpi, raster.dpi))

        # Add scale.

        if (isFALSE(enforce_symmetry)){
          p.loop <- add_scale(p = p.loop,
                              function_use = ggplot2::scale_color_viridis_c(na.value = na.value,
                                                                            option = viridis_color_map,
                                                                            direction = viridis_direction),
                              scale = "color")
        } else if (isTRUE(enforce_symmetry)){
          p.build <- ggplot2::ggplot_build(p.loop)
          feature.select <- gsub("-", ".", feature.use)
          limits <- c(min(p.build$plot$data[, feature.select]),
                      max(p.build$plot$data[, feature.select]))
          end_value <- max(abs(limits))
          p.loop <- add_scale(p = p.loop,
                              function_use = ggplot2::scale_color_gradientn(colors = c("#033270", "#4091C9", "#fdf0d5", "#c94040", "#65010C"),
                                                                            limits = c(-end_value, end_value),
                                                                            na.value = na.value),
                              scale = "color")
        }
        p.loop <- p.loop +
          ggplot2::ggtitle("")

        # Add cell borders.
        if (isTRUE(plot_cell_borders)){
          p.loop$layers <- append(base_layer, p.loop$layers)
        }

      } else if (!(is.null(split.by))){
        # Recover all metadata.
        data <- sample[[]]
        # Retrieve the metadata column belonging to the split.by parameter.
        data.use <- data[, split.by, drop = FALSE]
        # Retrieve the plotting order, keep factor levels if the column is a factor.
        if (is.null(split.by.idents)){
          plot_order <- if (is.factor(data.use[, 1])){levels(data.use[, 1])} else {sort(unique(data.use[, 1]))}
        } else {
          plot_order <- sort(split.by.idents)
        }
        list.plots.split.by <- list()
        count_plot <- 0 # Will update for each unique value in split.by.
        count_iteration <- count_iteration + 1 # Will update for each feature.
        # Compute the limits of the variable.
        ## This should be replaced by a function in utils.R
        limits <- c(min(sample$dummy[!is.na(sample$dummy)]), max(sample$dummy[!is.na(sample$dummy)]))
        # For each value in split.by.
        for (iteration in plot_order){
          feature.use <- "dummy2"
          count_plot <- count_plot + 1
          # Create a second dummy variable storing the values of the feature but only for the selected subset.
          sample$dummy2 <- sample$dummy
          # Get the names of the cells used in this subset.
          cells.iteration <- sample@meta.data[, split.by] == iteration
          # Assign the cells that are not part of the iteration to NA.
          sample$dummy2[!(cells.iteration)] <- NA
          p.loop <- Seurat::FeaturePlot(sample,
                                        feature.use,
                                        slot = slot,
                                        reduction = reduction,
                                        order = order,
                                        dims = dims,
                                        pt.size = pt.size,
                                        raster = raster,
                                        raster.dpi = c(raster.dpi, raster.dpi))

          if (isFALSE(enforce_symmetry)){
            p.loop <- add_scale(p = p.loop,
                                function_use = ggplot2::scale_color_viridis_c(na.value = na.value,
                                                                              option = viridis_color_map,
                                                                              limits = limits,
                                                                              direction = viridis_direction),
                                scale = "color")
          } else if (isTRUE(enforce_symmetry)){
            p.build <- ggplot2::ggplot_build(p.loop)
            end_value <- max(abs(limits))
            p.loop <- add_scale(p = p.loop,
                                function_use = ggplot2::scale_color_gradientn(colors = c("#033270", "#4091C9", "#fdf0d5", "#c94040", "#65010C"),
                                                                              limits = c(-end_value, end_value),
                                                                              na.value = na.value),
                                scale = "color")
          }
          p.loop <- p.loop +
            ggplot2::ggtitle(iteration)

          if (legend.position != "none"){
            p.loop <- modify_continuous_legend(p = p.loop,
                                               legend.title = if (is.null(legend.title)){feature} else {legend.title},
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

          if (iteration != plot_order[length(plot_order)]){
            p.loop <- p.loop & Seurat::NoLegend()
          }

          # Add cell borders.
          if (isTRUE(plot_cell_borders)){
            p.loop$layers <- append(base_layer, p.loop$layers)
          }

          list.plots.split.by[[iteration]] <- p.loop
        }
        p.loop <- patchwork::wrap_plots(list.plots.split.by,
                                        ncol = ncol,
                                        guides = "collect") +
          patchwork::plot_annotation(title = ifelse(typeof(individual.titles) == "character", individual.titles[[count_iteration]], ""),
                                     subtitle = ifelse(typeof(individual.subtitles) == "character", individual.subtitles[[count_iteration]], ""),
                                     caption = ifelse(typeof(individual.captions) == "character", individual.captions[[count_iteration]], ""))

      }

      # Patch for diffusion maps.
      if (reduction == "diffusion"){
        # Add "DC" labels.
        p.loop <- p.loop &
          ggplot2::xlab(paste0("DC_", dims[1])) &
          ggplot2::ylab(paste0("DC_", dims[2]))
      } else if (reduction == "pca"){
        p.loop <- p.loop &
          ggplot2::xlab(paste0("PC_", dims[1])) &
          ggplot2::ylab(paste0("PC_", dims[2]))
      }
      # Add the plot to the list.
      list.plots[[feature]] <- p.loop
    }

    # Generate the final plot with patchwork and use the "ncol" parameter value for the number of columns.
    p <- patchwork::wrap_plots(list.plots, nrow = 1)

  }

  # Fix the extra space and add theme parameters.
  p <- p &
    ggplot2::theme_minimal(base_size = font.size) &
    ggplot2::theme(plot.margin = ggplot2::margin(t = 10, r = 10, b = 10, l = 10),
                   plot.title = ggplot2::element_text(face = "bold",
                                                      hjust = ifelse(!(is.null(split.by)), 0.5, 0)),
                   plot.subtitle = ggplot2::element_text(hjust = 0),
                   plot.caption = ggplot2::element_text(hjust = 1),
                   panel.grid = ggplot2::element_blank(),
                   plot.title.position = "plot",
                   plot.caption.position = "plot",
                   text = ggplot2::element_text(family = font.type),
                   legend.text = ggplot2::element_text(face = "bold"),
                   legend.position = legend.position,
                   legend.title = ggplot2::element_text(face = "bold"),
                   legend.justification = "center",
                   plot.background = ggplot2::element_rect(fill = "white", color = "white"),
                   panel.background = ggplot2::element_rect(fill = "white", color = "white"),
                   legend.background = ggplot2::element_rect(fill = "white", color = "white"))
  if (is.null(split.by) & legend.position != "none"){
    counter <- 0
    for (feature in features){
      counter <- counter + 1
      p[[counter]] <- modify_continuous_legend(p = p[[counter]],
                                               legend.title = legend.title,
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

  # Add custom title.
  if (!is.null(plot.title)){
    if (length(features) > 1 | !(is.null(split.by)) | !(is.null(cells.highlight)) | !(is.null(idents.highlight))){
      p <- p +
        patchwork::plot_annotation(title = plot.title)
    } else {
      p <- p +
        ggplot2::labs(title = plot.title)
    }
  }

  # Add custom subtitle.
  if (!is.null(plot.subtitle)){
    if (length(features) > 1 | !(is.null(split.by)) | !(is.null(cells.highlight)) | !(is.null(idents.highlight))){
      p <- p +
        patchwork::plot_annotation(subtitle = plot.subtitle)
    } else {
      p <- p +
        ggplot2::labs(subtitle = plot.subtitle)
    }
  }

  # Add custom caption
  if (!is.null(plot.caption)){
    if (length(features) > 1 | !(is.null(split.by)) | !(is.null(cells.highlight)) | !(is.null(idents.highlight))){
      p <- p +
        patchwork::plot_annotation(caption = plot.caption)
    } else {
      p <- p +
        ggplot2::labs(caption = plot.caption)
    }
  }

  # Add individual titles.
  if (!is.null(individual.titles)){
    for (counter in seq(1,length(features))){
      if (!(is.na(individual.titles[counter]))){
        if (is.null(split.by)){
          p[[counter]]$labels$title <- individual.titles[counter]
        }
      }
    }
  }

  # Add individual subtitles.
  if (!is.null(individual.subtitles)){
    for (counter in seq(1,length(features))){
      if (!(is.na(individual.subtitles[counter]))){
        if (is.null(split.by)){
          p[[counter]]$labels$subtitle <- individual.subtitles[counter]
        }
      }
    }
  }

  # Add individual captions
  if (!is.null(individual.captions)){
    for (counter in seq(1,length(features))){
      if (!(is.na(individual.captions[counter]))){
        if (is.null(split.by)){
          p[[counter]]$labels$caption <- individual.captions[counter]
        }
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
                       axis.line =if (isFALSE(plot.axes)){ggplot2::element_blank()} else {ggplot2::element_line(color = "black")})
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

  # Add legend titles.
  if (is.null(split.by) & is.null(legend.title)){
    counter <- 0
    for (feature in features){
      counter <- counter + 1
      p[[counter]]$guides$colour$title <- feature
    }
  }

  # Return the plot.
  return(p)
}
