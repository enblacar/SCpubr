#' Wrapper for \link[Seurat]{FeaturePlot}.
#'
#'
#' @param sample Seurat object.
#' @param features Features to plot. It can be a single one or a vector of multiple features. Similar behavior as with \link[Seurat]{FeaturePlot}.
#' @param assay Assay to use. Defaults to the current assay.
#' @param reduction Reduction to use. Can be the canonical ones such as "umap", "pca", or any custom ones, such as "diffusion". If you are unsure about which reductions you have, use `Seurat::Reductions(sample)`. Defaults to "umap" if present or to the last computed reduction if the argument is not provided.
#' @param slot Data slot to use. Character. Only one of: counts, data, scale.data. Defaults to "data".
#' @param split.by Split the plot in as many unique values stored in the provided metadata column.
#' @param split.by.idents Vector of identities to plot. The gradient scale will also be subset to only the values of such identities.
#' @param cells.highlight Vector of cells for which the FeaturePlot should focus into. The rest of the cells will be grayed out.
#' @param idents.highlight Vector of identities that the FeaturePlot should focus into. Has to match the current Seurat identities in `Seurat::Idents(sample)`.
#' @param dims Vector of 2 numerics indicating the dimensions to plot out of the selected reduction. Defaults to c(1, 2) if not specified.
#' @param pt.size Point size.
#' @param font.size Base font.size of the figure.
#' @param font.type Character. Base font for the plot. One of mono, serif or sans.
#' @param legend Whether to plot the legend or not.
#' @param legend.type Character. Type of legend to display. One of: normal, colorbar, colorsteps.
#' @param legend.position Position of the legend in the plot. Will only work if legend is set to TRUE.
#' @param legend.framewidth,legend.tickwidth Width of the lines of the box in the legend.
#' @param legend.framecolor,legend.tickcolor Color of the lines of the box in the legend.
#' @param legend.length,legend.width Length and width of the legend. Will adjust automatically depending on legend side.
#' @param plot.title,plot.subtitle,plot.caption Title for the plot or subtitle or caption.
#' @param individual.titles,individual.subtitles,individual.captions Titles or subtitles. for each feature if needed. Either NULL or a vector of equal length of features.
#' @param ncol Number of columns to use in the arrangement of the output if more than one feature is queried to the function.
#' @param viridis_color_map Character. A capital letter from A to H or the scale name as in \link[viridis]{scale_fill_viridis}.
#' @param raster Whether to raster the resulting plot. This is recommendable if plotting a lot of cells.
#' @param raster.dpi Numeric. Pixel resolution for rasterized plots. Defaults to 512, as per default `Seurat::FeaturePlot()` behavior.
#' @param verbose Whether to show warnings.

#' @return  A ggplot2 object containing a Feature Plot.
#' @export
#'
#' @example /man/examples/examples_do_FeaturePlot.R
do_FeaturePlot <- function(sample,
                           features,
                           assay = NULL,
                           reduction = NULL,
                           slot = NULL,
                           split.by = NULL,
                           split.by.idents = NULL,
                           cells.highlight = NULL,
                           idents.highlight = NULL,
                           dims = c(1, 2),
                           pt.size = 0.5,
                           font.size = 14,
                           font.type = "sans",
                           legend = TRUE,
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
                           viridis_color_map = "D",
                           raster = FALSE,
                           raster.dpi = 1024,
                           verbose = TRUE){


  # Checks for packages.
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
  logical_list <- list("legend" = legend,
                       "verbose" = verbose,
                       "raster" = raster)
  check_type(parameters = logical_list, required_type = "logical", test_function = is.logical)
  # Check numeric parameters.
  numeric_list <- list("pt.size" = pt.size,
                       "ncol" = ncol,
                       "font.size" = font.size,
                       "raster.dpi" = raster.dpi,
                       "legend.framewidth" = legend.framewidth,
                       "legend.tickwidth" = legend.tickwidth,
                       "legend.length" = legend.length,
                       "legend.width" = legend.width)
  check_type(parameters = numeric_list, required_type = "numeric", test_function = is.numeric)
  # Check character parameters.
  # Workaround for features.
  if (is.list(features)){
    warning("Features provided as a list. Unlisting the list. Please use a character vector next time.", call. = F)
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
                         "font.type" = font.type)
  check_type(parameters = character_list, required_type = "character", test_function = is.character)

  # Check slot.
  slot <- check_and_set_slot(slot = slot)

  # Check split.by is on metadata.
  if (!(is.null(split.by))){check_feature(sample = sample, features = split.by, enforce_check = "metadata", enforce_parameter = "split.by")}

  # Check individual titles.
  if (!(is.null(individual.titles))){
    if(length(features) != length(individual.titles)){
      stop('Total number of individual titles does not match the number of features provided.', call. = F)
    }
  }

  # Check individual subtitles.
  if (!(is.null(individual.subtitles))){
    if(length(features) != length(individual.subtitles)){
      stop('Total number of individual subtitles does not match the number of features provided.', call. = F)
    }
  }

  # Check individual captions
  if (!(is.null(individual.captions))){
    if(length(features) != length(individual.captions)){
      stop('Total number of individual captions does not match the number of features provided.', call. = F)
    }
  }

  # Check viridis_color_map.
  check_viridis_color_map(viridis_color_map = viridis_color_map, verbose = verbose)

  # Check the colors provided to legend.framecolor and legend.tickcolor.
  check_colors(legend.framecolor, parameter_name = "legend.framecolor")
  check_colors(legend.tickcolor, parameter_name = "legend.tickcolor")

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
    warning("Setting raster = TRUE and pt.size < 1 will result in the cells being ploted as a cross. This behaviour can not be modified, but setting pt.size to 1 or higher solves it. For Feature plots, optimized values would be pt.size = 3 and raster.dpi = 2048.", call. = F)
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
                             reduction = reduction,
                             order = T,
                             dims = dims,
                             pt.size = pt.size,
                             ncol = ncol,
                             raster = raster,
                             raster.dpi = c(raster.dpi, raster.dpi)) &
      # Remove Seurat::FeaturePlot() default plot title.
      ggplot2::ggtitle("")
    # Add viridis and supress warnings.
    p <- add_scale(p = p,
                   num_plots = length(features),
                   function_use = ggplot2::scale_color_viridis_c(na.value = "grey75",
                                                                 option = viridis_color_map),
                   scale = "color")
    # Special patches for diffusion maps: Adding "DC" labels to the axis.
    if (reduction == "diffusion"){
      p <- p &
        ggplot2::xlab(paste0("DC_", dims[1])) &
        ggplot2::ylab(paste0("DC_", dims[2]))
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
                                      order = T,
                                      dims = dims,
                                      pt.size = pt.size,
                                      raster = raster,
                                      raster.dpi = c(raster.dpi, raster.dpi))
        p.loop <- add_scale(p = p.loop,
                            num_plots = length(feature),
                            function_use = ggplot2::scale_color_viridis_c(na.value = "grey75",
                                                                          option = viridis_color_map),
                            scale = "color")
        p.loop <- p.loop +
          ggplot2::ggtitle("")

      } else if (!(is.null(split.by))){
        # Recover all metadata.
        data <- sample[[]]
        # Retrieve the metadata column belonging to the split.by parameter.
        data.use <- data[, split.by, drop = F]
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
                                        reduction = reduction,
                                        order = T,
                                        dims = dims,
                                        pt.size = pt.size,
                                        raster = raster,
                                        raster.dpi = c(raster.dpi, raster.dpi))
          p.loop <- add_scale(p = p.loop,
                              num_plots = length(feature),
                              function_use = ggplot2::scale_color_viridis_c(na.value = "grey75",
                                                                            option = viridis_color_map,
                                                                            limits = limits),
                              scale = "color")
          p.loop <- p.loop +
            ggplot2::ggtitle(iteration)

          if (legend.type == "normal"){
            p.loop <- p.loop +
                      ggplot2::guides(color = ggplot2::guide_colorbar(title = feature,
                                                                      title.position = "top",
                                                                      title.hjust = 0.5))
          } else if (legend.type == "colorbar"){
            p.loop <- p.loop +
                      ggplot2::guides(color = ggplot2::guide_colorbar(title = feature,
                                                                      title.position = "top",
                                                                      barwidth = legend.barwidth,
                                                                      barheight = legend.barheight,
                                                                      title.hjust = 0.5,
                                                                      ticks.linewidth = legend.tickwidth,
                                                                      frame.linewidth = legend.framewidth,
                                                                      frame.colour = legend.framecolor,
                                                                      ticks.colour = legend.tickcolor))
          } else if (legend.type == "colorsteps"){
            p.loop <- p.loop +
              ggplot2::guides(color = ggplot2::guide_colorsteps(title = feature,
                                                                title.position = "top",
                                                                barwidth = legend.barwidth,
                                                                barheight = legend.barheight,
                                                                title.hjust = 0.5,
                                                                ticks.linewidth = legend.tickwidth,
                                                                frame.linewidth = legend.framewidth,
                                                                frame.colour = legend.framecolor,
                                                                ticks.colour = legend.tickcolor))
          }

          if (iteration != plot_order[length(plot_order)]){
            p.loop <- p.loop & Seurat::NoLegend()
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
    ggplot2::coord_cartesian(expand = FALSE) &
    ggplot2::theme(plot.margin = ggplot2::margin(t = 10, r = 10, b = 10, l = 10),
                   plot.title = ggtext::element_markdown(face = "bold",
                                                         hjust = ifelse(!(is.null(split.by)), 0.5, 0)),
                   plot.subtitle = ggtext::element_markdown(hjust = 0),
                   plot.caption = ggtext::element_markdown(hjust = 1),
                   panel.grid = ggplot2::element_blank(),
                   plot.title.position = "plot",
                   plot.caption.position = "plot",
                   text = ggplot2::element_text(family = font.type),
                   legend.text = ggplot2::element_text(face = "bold"),
                   legend.position = legend.position,
                   legend.title = ggplot2::element_text(face = "bold"),
                   legend.justification = "center",
                   axis.title.x = ggplot2::element_text(face = "bold"),
                   axis.title.y = ggplot2::element_text(face = "bold",
                                                        angle = 90),
                   plot.background = ggplot2::element_rect(fill = "white", color = "white"),)
  if (is.null(split.by)){
    if (legend.type == "normal"){
      p <- p &
           ggplot2::guides(color = ggplot2::guide_colorbar(title.position = "top",
                                                           title.hjust = 0.5))
    } else if (legend.type == "colorbar"){
      p <- p &
           ggplot2::guides(color = ggplot2::guide_colorbar(title.position = "top",
                                                           barwidth = legend.barwidth,
                                                           barheight = legend.barheight,
                                                           title.hjust = 0.5,
                                                           ticks.linewidth = legend.tickwidth,
                                                           frame.linewidth = legend.framewidth,
                                                           frame.colour = legend.framecolor,
                                                           ticks.colour = legend.tickcolor))
    } else if (legend.type == "colorsteps"){
      p <- p &
           ggplot2::guides(color = ggplot2::guide_colorsteps(title.position = "top",
                                                             barwidth = legend.barwidth,
                                                             barheight = legend.barheight,
                                                             title.hjust = 0.5,
                                                             ticks.linewidth = legend.tickwidth,
                                                             frame.linewidth = legend.framewidth,
                                                             frame.colour = legend.framecolor,
                                                             ticks.colour = legend.tickcolor))
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



  # Remove legend.
  if (legend == FALSE){
    p <- p &
         ggplot2::theme(legend.position = "none")
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
                       axis.line = ggplot2::element_blank()) &
        ggplot2::xlab(labels[1]) &
        ggplot2::ylab(labels[2])
    }
    # For diffusion maps, we do want to keep at least the axis titles so that we know which DC are we plotting.
  } else {
    labels <- colnames(sample@reductions[[reduction]][[]])[dims]
    p <- p &
      ggplot2::theme(axis.text = ggplot2::element_blank(),
                     axis.ticks = ggplot2::element_blank(),
                     axis.line = ggplot2::element_blank()) &
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
                       axis.line = ggplot2::element_blank())
    })

  }

  # Add legend titles.
  if (is.null(split.by)){
    counter <- 0
    for (feature in features){
      counter <- counter + 1
      p[[counter]]$guides$colour$title <- feature
    }
  }

  # Return the plot.
  return(p)
}
