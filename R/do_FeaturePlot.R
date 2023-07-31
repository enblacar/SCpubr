#' Wrapper for \link[Seurat]{FeaturePlot}.
#'
#' @inheritParams doc_function
#' @param idents.keep \strong{\code{\link[base]{character}}} | Vector of identities to plot. The gradient scale will also be subset to only the values of such identities.
#' @param individual.titles,individual.subtitles,individual.captions \strong{\code{\link[base]{character}}} | Titles or subtitles. for each feature if needed. Either NULL or a vector of equal length of features.
#' @param order \strong{\code{\link[base]{logical}}} | Whether to order the cells based on expression.
#' @param group.by \strong{\code{\link[base]{character}}} | Metadata variable based on which cells are grouped. This will effectively introduce a big dot in the center of each cluster, colored using a categorical color scale or with the values provided by the user in \strong{\code{group.by.colors.use}}. It will also displays a legend.
#' @param group.by.legend \strong{\code{\link[base]{character}}} | Title for the legend when \strong{\code{group.by}} is used. Use \strong{\code{NA}} to disable it and \strong{\code{NULL}} to use the default column title provided in \strong{\code{group.by}}.
#' @param group.by.dot.size \strong{\code{\link[base]{numeric}}} | Size of the dots placed in the middle of the groups.
#' @param group.by.cell_borders \strong{\code{\link[base]{logical}}} | Plots another border around the cells displaying the same color code of the dots displayed with \strong{\code{group.by}}. Legend is shown always with alpha = 1 regardless of the alpha settings.
#' @param group.by.colors.use \strong{\code{\link[base]{character}}} | Colors to use for the group dots.
#' @param group.by.show.dots \strong{\code{\link[base]{logical}}} | Controls whether to place in the middle of the groups.
#' @param group.by.cell_borders.alpha \strong{\code{\link[base]{numeric}}} | Controls the transparency of the new borders drawn by \strong{\code{group.by.cell_borders}}.
#' @return  A ggplot2 object containing a Feature Plot.
#' @export
#'
#' @example /man/examples/examples_do_FeaturePlot.R
do_FeaturePlot <- function(sample,
                           features,
                           assay = NULL,
                           reduction = NULL,
                           slot = NULL,
                           order = FALSE,
                           group.by = NULL,
                           group.by.colors.use = NULL,
                           group.by.legend = NULL,
                           group.by.show.dots = TRUE,
                           group.by.dot.size = 8,
                           group.by.cell_borders = FALSE,
                           group.by.cell_borders.alpha = 0.1,
                           split.by = NULL,
                           idents.keep = NULL,
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
                           legend.framewidth = 0.5,
                           legend.tickwidth = 0.5,
                           legend.length = 20,
                           legend.width = 1,
                           legend.framecolor = "grey50",
                           legend.tickcolor = "white",
                           legend.ncol = NULL,
                           legend.nrow = NULL,
                           legend.byrow = FALSE,
                           plot.title = NULL,
                           plot.subtitle = NULL,
                           plot.caption = NULL,
                           individual.titles = NULL,
                           individual.subtitles = NULL,
                           individual.captions = NULL,
                           ncol = NULL,
                           use_viridis = FALSE,
                           viridis.palette = "G",
                           viridis.direction = 1,
                           raster = FALSE,
                           raster.dpi = 1024,
                           plot_cell_borders = TRUE,
                           border.size = 2,
                           border.color = "black",
                           border.density = 1,
                           na.value = "grey75",
                           verbose = TRUE,
                           plot.axes = FALSE,
                           min.cutoff = rep(NA, length(features)),
                           max.cutoff = rep(NA, length(features)),
                           plot_density_contour = FALSE,
                           contour.position = "bottom",
                           contour.color = "grey90",
                           contour.lineend = "butt",
                           contour.linejoin = "round",
                           contour_expand_axes = 0.25,
                           label = FALSE,
                           label.color = "black",
                           label.size = 4,
                           number.breaks = 5,
                           diverging.palette = "RdBu",
                           diverging.direction = -1,
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
  
  check_suggests(function_name = "do_FeaturePlot")
  # Check if the sample provided is a Seurat object.
  check_Seurat(sample = sample)
  # Check the assay.
  out <- check_and_set_assay(sample = sample, assay = assay)
  sample <- out[["sample"]]
  assay <- out[["assay"]]
  
  sample <- check_Assay5(sample, assay = assay)
  
  # Check the reduction.
  reduction <- check_and_set_reduction(sample = sample, reduction = reduction)
  # Check the dimensions.
  dims <- check_and_set_dimensions(sample = sample, reduction = reduction, dims = dims)
  # Check logical parameters.
  logical_list <- list("verbose" = verbose,
                       "raster" = raster,
                       "plot_cell_borders" = plot_cell_borders,
                       "order" = order,
                       "enforce_symmetry" = enforce_symmetry,
                       "plot.axes" = plot.axes,
                       "plot_density_contour" = plot_density_contour,
                       "label" = label,
                       "legend.byrow" = legend.byrow,
                       "group.by.cell_borders" = group.by.cell_borders,
                       "group.by.show.dots" = group.by.show.dots,
                       "use_viridis" = use_viridis)
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
                       "viridis.direction" = viridis.direction,
                       "min.cutoff" = min.cutoff,
                       "max.cutoff" = max.cutoff,
                       "contour_expand_axes" = contour_expand_axes,
                       "label.size" = label.size,
                       "number.breaks" = number.breaks,
                       "border.density" = border.density,
                       "legend.nrow" = legend.nrow,
                       "legend.ncol" = legend.ncol,
                       "group.by.dot.size" = group.by.dot.size,
                       "group.by.cell_borders.alpha" = group.by.cell_borders.alpha,
                       "sequential.direction" = sequential.direction,
                       "diverging.direction" = diverging.direction)
  check_type(parameters = numeric_list, required_type = "numeric", test_function = is.numeric)
  # Check character parameters.
  # Workaround for features.
  if (is.list(features)){
    warning(paste0(add_warning(), crayon_body("Features provided as a list. Unlisting the list. Please use a character vector next time."), call. = FALSE))
    features <- unique(unlist(features))
  }
  character_list <- list("legend.position" = legend.position,
                         "features" = features,
                         "cells.highlight" = cells.highlight,
                         "idents.highlight" = idents.highlight,
                         "slot" = slot,
                         "group.by" = group.by,
                         "split.by" = split.by,
                         "plot.title" = plot.title,
                         "plot.subtitle" = plot.subtitle,
                         "plot.caption" = plot.caption,
                         "idents.keep" = idents.keep,
                         "viridis.palette" = viridis.palette,
                         "individual.titles" = individual.titles,
                         "individual.subtitles" = individual.subtitles,
                         "individual.captions" = individual.captions,
                         "legend.framecolor" = legend.framecolor,
                         "legend.tickcolor" = legend.tickcolor,
                         "legend.type" = legend.type,
                         "font.type" = font.type,
                         "border.color" = border.color,
                         "legend.title" = legend.title,
                         "na.value" = na.value,
                         "contour.position" = contour.position,
                         "contour.color" = contour.color,
                         "contour.lineend" = contour.lineend,
                         "contour.linejoin" = contour.linejoin,
                         "label.color" = label.color,
                         "group.by.colors.use" = group.by.colors.use,
                         "diverging.palette" = diverging.palette,
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

  # Check split.by is on metadata.
  if (!(is.null(split.by))){check_feature(sample = sample, features = split.by, enforce_check = "metadata", enforce_parameter = "split.by")}

  # Check individual titles.
  if (length(features) > 1 & !is.null(individual.titles)){
    assertthat::assert_that(length(features) == length(individual.titles),
                            msg = paste0(add_cross(), crayon_body("The total number of "),
                                         crayon_key("individual titles"),
                                         crayon_body(" does not match the number of "),
                                         crayon_key("features"),
                                         crayon_body(" provided.")))
  }

  if (length(features) > 1 & !is.null(individual.subtitles)){
    assertthat::assert_that(length(features) == length(individual.subtitles),
                            msg = paste0(add_cross(), crayon_body("The total number of "),
                                         crayon_key("individual subtitles"),
                                         crayon_body(" does not match the number of "),
                                         crayon_key("features"),
                                         crayon_body(" provided.")))
  }

  if (length(features) > 1 & !is.null(individual.captions)){
    assertthat::assert_that(length(features) == length(individual.captions),
                            msg = paste0(add_cross(), crayon_body("The total number of "),
                                         crayon_key("individual captions"),
                                         crayon_body(" does not match the number of "),
                                         crayon_key("features"),
                                         crayon_body(" provided.")))
  }

  ## Check that the contour_expand_axes is between 0 and 1.
  assertthat::assert_that(contour_expand_axes <= 1,
                          msg = paste0(add_cross(), crayon_body("Please provide a value to "),
                                       crayon_key("countour_expand_axes"),
                                       crayon_body(" lower or equal than "),
                                       crayon_key("1"),
                                       crayon_body(".")))

  assertthat::assert_that(contour_expand_axes >= 0,
                          msg = paste0(add_cross(), crayon_body("Please provide a value to "),
                                       crayon_key("countour_expand_axes"),
                                       crayon_body(" higher or equal than "),
                                       crayon_key("0"),
                                       crayon_body(".")))


  check_colors(border.color, parameter_name = "border.color")
  check_colors(na.value, parameter_name = "na.value")
  check_colors(legend.framecolor, parameter_name = "legend.framecolor")
  check_colors(legend.tickcolor, parameter_name = "legend.tickcolor")
  check_colors(contour.color, parameter_name = "contour.color")
  check_colors(label.color, parameter_name = "label.color")
  if (!is.null(group.by)){
    if (!is.null(group.by.colors.use)){
      check_colors(group.by.colors.use, parameter_name = "group.by.colors.use")
      check_consistency_colors_and_names(sample = sample,
                                         colors = group.by.colors.use,
                                         grouping_variable = group.by)
    } else {
      data.use <- sample@meta.data[, group.by, drop = FALSE]
      # If the variable is a factor, use the levels as order. If not, order the values alphabetically.
      names.use <- if (is.factor(data.use[, 1])){levels(data.use[, 1])} else {sort(unique(data.use[, 1]))}
      # Generate the color scale to be used based on the unique values of group.by.
      group.by.colors.use <- generate_color_scale(names.use)
    }
  }



  check_parameters(parameter = font.type, parameter_name = "font.type")
  check_parameters(parameter = legend.type, parameter_name = "legend.type")
  check_parameters(parameter = legend.position, parameter_name = "legend.position")
  check_parameters(parameter = viridis.palette, parameter_name = "viridis.palette")
  check_parameters(parameter = contour.lineend, parameter_name = "contour.lineend")
  check_parameters(parameter = contour.linejoin, parameter_name = "contour.linejoin")
  check_parameters(parameter = contour.position, parameter_name = "contour.position")
  check_parameters(parameter = number.breaks, parameter_name = "number.breaks")
  check_parameters(parameter = border.density, parameter_name = "border.density")
  check_parameters(parameter = diverging.palette, parameter_name = "diverging.palette")
  check_parameters(parameter = sequential.palette, parameter_name = "sequential.palette")
  check_parameters(plot.title.face, parameter_name = "plot.title.face")
  check_parameters(plot.subtitle.face, parameter_name = "plot.subtitle.face")
  check_parameters(plot.caption.face, parameter_name = "plot.caption.face")
  check_parameters(axis.title.face, parameter_name = "axis.title.face")
  check_parameters(axis.text.face, parameter_name = "axis.text.face")
  check_parameters(legend.title.face, parameter_name = "legend.title.face")
  check_parameters(legend.text.face, parameter_name = "legend.text.face")
  check_parameters(viridis.direction, parameter_name = "viridis.direction")
  check_parameters(sequential.direction, parameter_name = "sequential.direction")
  check_parameters(diverging.direction, parameter_name = "diverging.direction")

  if (length(min.cutoff) != length(features)){
    warning(paste0(add_warning(), crayon_body("Please provide as many values to "),
                   crayon_key("min.cutoff"),
                   crayon_body(" as "),
                   crayon_key("features"),
                   crayon_body(" provided. The values will be used in order and, when outside of the range, no cutoffs will be applied.")), call. = FALSE)
  }
  
  if (length(max.cutoff) != length(features)){
    warning(paste0(add_warning(), crayon_body("Please provide as many values to "),
                   crayon_key("max.cutoff"),
                   crayon_body(" as "),
                   crayon_key("features"),
                   crayon_body(" provided. The values will be used in order and, when outside of the range, no cutoffs will be applied.")), call. = FALSE)
  }
  
  # Generate the continuous color palette.
  if (isTRUE(enforce_symmetry)){
    colors.gradient <- compute_continuous_palette(name = diverging.palette,
                                                  use_viridis = FALSE,
                                                  direction = diverging.direction,
                                                  enforce_symmetry = enforce_symmetry)
  } else {
    colors.gradient <- compute_continuous_palette(name = ifelse(isTRUE(use_viridis), viridis.palette, sequential.palette),
                                                  use_viridis = use_viridis,
                                                  direction = ifelse(isTRUE(use_viridis), viridis.direction, sequential.direction),
                                                  enforce_symmetry = enforce_symmetry)
  }
  
  
  # Generate base layer.
  if (isTRUE(plot_cell_borders)){
    out <- compute_umap_layer(sample = sample,
                              labels = colnames(sample@reductions[[reduction]][[]])[dims],
                              pt.size = pt.size,
                              dot.size = group.by.dot.size,
                              border.density = border.density,
                              border.size = border.size,
                              border.color = border.color,
                              raster = raster,
                              raster.dpi = raster.dpi,
                              reduction = reduction,
                              group.by = group.by,
                              split.by = split.by,
                              na.value = na.value,
                              alpha = group.by.cell_borders.alpha,
                              n = 100)
    base_layer <- out$base_layer
    na_layer <- out$na_layer
    if (!is.null(group.by)){
      center_layers <- out$center_layers
      center_layer_1 <- center_layers$center_layer_1
      center_layer_2 <- center_layers$center_layer_2
      color_layer <- out$color_layer
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

    if (utils::packageVersion("Seurat") >= "4.1.0"){
      p <- Seurat::FeaturePlot(sample,
                               features,
                               slot = slot,
                               reduction = reduction,
                               order = order,
                               dims = dims,
                               pt.size = pt.size,
                               ncol = ncol,
                               raster = raster,
                               raster.dpi = c(raster.dpi, raster.dpi),
                               min.cutoff = min.cutoff,
                               max.cutoff = max.cutoff,
                               label = label,
                               label.size = label.size,
                               label.color = label.color)
    } else { # nocov start
      p <- Seurat::FeaturePlot(sample,
                               features,
                               slot = slot,
                               reduction = reduction,
                               order = order,
                               dims = dims,
                               pt.size = pt.size,
                               ncol = ncol,
                               raster = raster,
                               min.cutoff = min.cutoff,
                               max.cutoff = max.cutoff,
                               label = label,
                               label.size = label.size,
                               label.color = label.color)
    } # nocov end
    p$layers[[length(p$layers)]]$aes_params$fontface <- "bold"
     p <- p &
          # Remove Seurat::FeaturePlot() default plot title.
          ggplot2::ggtitle("")

    # Add color scales.
    num_plots <- length(features)
    for (counter in seq(1, num_plots)){
      if (num_plots == 1){
        scale.setup <- compute_scales(sample = sample,
                                      feature = features,
                                      assay = assay,
                                      reduction = NULL,
                                      slot = slot,
                                      number.breaks = number.breaks,
                                      min.cutoff = min.cutoff,
                                      max.cutoff = max.cutoff,
                                      flavor = "Seurat",
                                      enforce_symmetry = enforce_symmetry)
        p <- add_scale(p = p,
                       function_use = ggplot2::scale_color_gradientn(colors = colors.gradient,
                                                                     na.value = na.value,
                                                                     name = legend.title,
                                                                     breaks = scale.setup$breaks,
                                                                     labels = scale.setup$labels,
                                                                     limits = scale.setup$limits),
                       scale = "color")
      } else if (num_plots > 1){
        
        feature.use <- features[counter]
        scale.setup <- compute_scales(sample = sample,
                                      feature = feature.use,
                                      assay = assay,
                                      reduction = NULL,
                                      slot = slot,
                                      number.breaks = number.breaks,
                                      min.cutoff = min.cutoff[counter],
                                      max.cutoff = max.cutoff[counter],
                                      flavor = "Seurat",
                                      enforce_symmetry = enforce_symmetry)

        p[[counter]] <- add_scale(p = p[[counter]],
                                  function_use = ggplot2::scale_color_gradientn(colors = colors.gradient,
                                                                                na.value = na.value,
                                                                                name = legend.title,
                                                                                breaks = scale.setup$breaks,
                                                                                labels = scale.setup$labels,
                                                                                limits = scale.setup$limits),
                                  scale = "color")
      }
    }

    # Special patches for diffusion maps: Adding "DC" labels to the axis.
    if (stringr::str_starts(reduction, "diff|DIFF")){
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
        if (!(is.null(group.by))){
          if (isTRUE(group.by.show.dots)){
            p[[counter]]$layers <- append(p[[counter]]$layers, c(center_layer_2, center_layer_1))
          }
          if (isTRUE(group.by.cell_borders)){
            p[[counter]]$layers <- append(color_layer, p[[counter]]$layers)
          }
          if (is.null(group.by.legend)){
            legend.name <- group.by
          } else {
            if (is.na(group.by.legend)){
              legend.name <- ""
            } else {
              legend.name <- group.by.legend
            }
          }
          p[[counter]] <- p[[counter]] +
                          ggplot2::scale_fill_manual(values = group.by.colors.use,
                                                     name = legend.name) +
                          ggplot2::guides(fill = ggplot2::guide_legend(title = legend.name,
                                                                       title.position = "top",
                                                                       title.hjust = 0.5,
                                                                       ncol = legend.ncol,
                                                                       nrow = legend.nrow,
                                                                       byrow = legend.byrow,
                                                                       override.aes = list(alpha = 1)))
        }
      }
    }
    if (isTRUE(plot_density_contour)){
      counter <- 0
      for (feature in features){
        counter <- counter + 1

        data <- ggplot2::ggplot_build(p[[counter]])

        density_layer <- ggplot2::stat_density_2d(data = data$data[[1]],
                                                  mapping = ggplot2::aes(x = .data$x,
                                                                         y = .data$y),
                                                  color = contour.color,
                                                  lineend = contour.lineend,
                                                  linejoin = contour.linejoin)
        if (contour.position == "bottom"){
          p[[counter]]$layers <- append(density_layer, p[[counter]]$layers)
        } else if (contour.position == "top"){
          p[[counter]]$layers <- append(p[[counter]]$layers, density_layer)
        }

        min_x <- min(data$data[[1]]$x) * (1 + contour_expand_axes)
        max_x <- max(data$data[[1]]$x) * (1 + contour_expand_axes)
        min_y <- min(data$data[[1]]$y) * (1 + contour_expand_axes)
        max_y <- max(data$data[[1]]$y) * (1 + contour_expand_axes)
        # Expand axes limits to allocate the new contours.
        suppressMessages({
          p[[counter]] <- p[[counter]] +
            ggplot2::xlim(c(min_x, max_x)) +
            ggplot2::ylim(c(min_y, max_y))
        })

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
      if (is.null(idents.keep)){
        cells.use <- colnames(sample)
        # Identitites selected by the user.
      } else {
        # Check if the identitites are duplicated. If so, remove them.
        if (sum(duplicated(idents.keep)) != 0){
          message(paste0(add_info(), crayon_body("Found and removed duplicated values in idents.keep.")))
          idents.keep <- idents.keep[!duplicated(idents.keep)]
        }
        # Get the names of the cells to plot.
        cells.use <- names(Seurat::Idents(sample)[Seurat::Idents(sample) %in% idents.keep])
      }
    }
    # Plots are generated independently if more than one feature is provided.
    list.plots <- list()
    # Counter depicting the feature used. Will increase each feature used.
    count_iteration <- 0
    cutoff.counter <- 0
    # Iterate over the features.
    for (feature in features){
      cutoff.counter <- cutoff.counter + 1
      min.cutoff.use <- min.cutoff[cutoff.counter]
      max.cutoff.use <- max.cutoff[cutoff.counter]

      # A "dummy" metadata column is generated using the values of the selected feature.
      ## Based on whether the feature is in the metadata of the object.
      if (feature %in% colnames(sample@meta.data)) {
        sample$dummy <- sample@meta.data[, feature]
        ## Or is a gene in the object.
      } else if (feature %in% rownames(sample)){
        sample$dummy <- .GetAssayData(sample = sample, slot = slot, assay = assay)[feature, ]
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

        sample.use <- sample[, cells.use]

        if (isTRUE(plot_cell_borders)){
          out <- compute_umap_layer(sample = sample.use,
                                    labels = colnames(sample.use@reductions[[reduction]][[]])[dims],
                                    pt.size = pt.size,
                                    border.density = border.density,
                                    border.size = border.size,
                                    border.color = border.color,
                                    raster = raster,
                                    raster.dpi = raster.dpi,
                                    reduction = reduction,
                                    group.by = group.by,
                                    split.by = split.by,
                                    n = 100)
          base_layer_subset <- out$base_layer
        }

        if (utils::packageVersion("Seurat") >= "4.1.0"){
          p.loop <- Seurat::FeaturePlot(sample.use,
                                        feature.use,
                                        reduction = reduction,
                                        slot = slot,
                                        order = order,
                                        dims = dims,
                                        pt.size = pt.size,
                                        raster = raster,
                                        raster.dpi = c(raster.dpi, raster.dpi),
                                        min.cutoff = min.cutoff.use,
                                        max.cutoff = max.cutoff.use,
                                        label = label,
                                        label.size = label.size,
                                        label.color = label.color)
        } else { # nocov start
          p.loop <- Seurat::FeaturePlot(sample.use,
                                        feature.use,
                                        reduction = reduction,
                                        slot = slot,
                                        order = order,
                                        dims = dims,
                                        pt.size = pt.size,
                                        raster = raster,
                                        min.cutoff = min.cutoff.use,
                                        max.cutoff = max.cutoff.use,
                                        label = label,
                                        label.size = label.size,
                                        label.color = label.color)
        } # nocov end
        p.loop$layers[[length(p.loop$layers)]]$aes_params$fontface <- "bold"

        # Add scale.

        scale.setup <- compute_scales(sample = sample,
                                      feature = feature,
                                      assay = assay,
                                      reduction = NULL,
                                      slot = slot,
                                      number.breaks = number.breaks,
                                      min.cutoff = min.cutoff.use,
                                      max.cutoff = max.cutoff.use,
                                      flavor = "Seurat",
                                      enforce_symmetry = enforce_symmetry)

        p.loop <- add_scale(p = p.loop,
                            function_use = ggplot2::scale_color_gradientn(colors = colors.gradient,
                                                                          na.value = na.value,
                                                                          name = legend.title,
                                                                          breaks = scale.setup$breaks,
                                                                          labels = scale.setup$labels,
                                                                          limits = scale.setup$limits),
                            scale = "color")
        
        p.loop <- p.loop +
                  ggplot2::ggtitle("")

        # Add cell borders.
        if (isTRUE(plot_cell_borders)){
          p.loop$layers <- append(base_layer_subset, p.loop$layers)
          p.loop$layers <- append(na_layer, p.loop$layers)
          p.loop$layers <- append(base_layer, p.loop$layers)
          
          suppressMessages({
            p.loop <- p.loop + 
              ggplot2::scale_x_continuous(limits = c(min(p.loop$layers[[1]]$data$x), 
                                                     max(p.loop$layers[[1]]$data$x))) + 
              ggplot2::scale_y_continuous(limits = c(min(p.loop$layers[[1]]$data$y), 
                                                     max(p.loop$layers[[1]]$data$y)))
          })
          
          if (!(is.null(group.by))){
            if (isTRUE(group.by.show.dots)){
              p.loop$layers <- append(p.loop$layers, c(center_layer_2, center_layer_1))
            }

            if (isTRUE(group.by.cell_borders)){
              p.loop$layers <- append(color_layer, p.loop$layers)
            }
            if (is.null(group.by.legend)){
              legend.name <- group.by
            } else {
              if (is.na(group.by.legend)){
                legend.name <- ""
              } else {
                legend.name <- group.by.legend
              }
            }
            p.loop <- p.loop +
              ggplot2::scale_fill_manual(values = group.by.colors.use,
                                         name = legend.name) +
              ggplot2::guides(fill = ggplot2::guide_legend(title = legend.name,
                                                           title.position = "top",
                                                           title.hjust = 0.5,
                                                           ncol = legend.ncol,
                                                           nrow = legend.nrow,
                                                           byrow = legend.byrow,
                                                           override.aes = list(alpha = 1)))
          }
        }

        if (isTRUE(plot_density_contour)){
          data <- ggplot2::ggplot_build(p.loop)

          density_layer <- ggplot2::stat_density_2d(data = data$data[[1]],
                                                    mapping = ggplot2::aes(x = .data$x,
                                                                           y = .data$y),
                                                    color = contour.color,
                                                    lineend = contour.lineend,
                                                    linejoin = contour.linejoin)
          if (contour.position == "bottom"){
            p.loop$layers <- append(density_layer, p.loop$layers)
          } else if (contour.position == "top"){
            p.loop$layers <- append(p.loop$layers, density_layer)
          }

          min_x <- min(data$data[[1]]$x) * (1 + contour_expand_axes)
          max_x <- max(data$data[[1]]$x) * (1 + contour_expand_axes)
          min_y <- min(data$data[[1]]$y) * (1 + contour_expand_axes)
          max_y <- max(data$data[[1]]$y) * (1 + contour_expand_axes)
          # Expand axes limits to allocate the new contours.
          suppressMessages({
            p.loop <- p.loop +
              ggplot2::xlim(c(min_x, max_x)) +
              ggplot2::ylim(c(min_y, max_y))
          })

        }

      } else if (!(is.null(split.by))){
        # Recover all metadata.
        data <- sample[[]]
        # Retrieve the metadata column belonging to the split.by parameter.
        data.use <- data[, split.by, drop = FALSE]
        # Retrieve the plotting order, keep factor levels if the column is a factor.
        if (is.null(idents.keep)){
          plot_order <- if (is.factor(data.use[, 1])){levels(data.use[, 1])} else {sort(unique(data.use[, 1]))}
        } else {
          plot_order <- sort(idents.keep)
        }
        list.plots.split.by <- list()
        count_plot <- 0 # Will update for each unique value in split.by.
        count_iteration <- count_iteration + 1 # Will update for each feature.
        # Compute the limits of the variable.
        ## This should be replaced by a function in utils.R
        limits <- c(min(sample$dummy[!is.na(sample$dummy)]), max(sample$dummy[!is.na(sample$dummy)]))
        if (!is.na(min.cutoff.use)){
          limits[1] <- min.cutoff.use
        }
        if (!is.na(max.cutoff.use)){
          limits[2] <- max.cutoff.use
        }
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

          sample.use <- sample[, cells.iteration]

          if (isTRUE(plot_cell_borders)){
            out <- compute_umap_layer(sample = sample.use,
                                      labels = colnames(sample.use@reductions[[reduction]][[]])[dims],
                                      pt.size = pt.size,
                                      border.density = border.density,
                                      border.size = border.size,
                                      border.color = border.color,
                                      raster = raster,
                                      raster.dpi = raster.dpi,
                                      reduction = reduction,
                                      group.by = group.by,
                                      split.by = split.by,
                                      n = 100)
            base_layer_subset <- out$base_layer
          }




          if (utils::packageVersion("Seurat") >= "4.1.0"){
            p.loop <- Seurat::FeaturePlot(sample.use,
                                          feature.use,
                                          slot = slot,
                                          reduction = reduction,
                                          order = order,
                                          dims = dims,
                                          pt.size = pt.size,
                                          raster = raster,
                                          raster.dpi = c(raster.dpi, raster.dpi),
                                          min.cutoff = min.cutoff.use,
                                          max.cutoff = max.cutoff.use,
                                          label = label,
                                          label.size = label.size,
                                          label.color = label.color)
          } else { # nocov start
            p.loop <- Seurat::FeaturePlot(sample.use,
                                          feature.use,
                                          slot = slot,
                                          reduction = reduction,
                                          order = order,
                                          dims = dims,
                                          pt.size = pt.size,
                                          raster = raster,
                                          min.cutoff = min.cutoff.use,
                                          max.cutoff = max.cutoff.use,
                                          label = label,
                                          label.size = label.size,
                                          label.color = label.color)
          } # nocov end
          p.loop$layers[[length(p.loop$layers)]]$aes_params$fontface <- "bold"


          scale.setup <- compute_scales(sample = sample,
                                        feature = feature,
                                        assay = assay,
                                        reduction = NULL,
                                        slot = slot,
                                        number.breaks = number.breaks,
                                        min.cutoff = min.cutoff.use,
                                        max.cutoff = max.cutoff.use,
                                        flavor = "Seurat",
                                        enforce_symmetry = enforce_symmetry)

          p.loop <- add_scale(p = p.loop,
                              function_use = ggplot2::scale_color_gradientn(colors = colors.gradient,
                                                                            na.value = na.value,
                                                                            name = legend.title,
                                                                            breaks = scale.setup$breaks,
                                                                            labels = scale.setup$labels,
                                                                            limits = scale.setup$limits),
                              scale = "color")

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
            p.loop$layers <- append(base_layer_subset, p.loop$layers)
            p.loop$layers <- append(na_layer, p.loop$layers)
            p.loop$layers <- append(base_layer, p.loop$layers)
    
            suppressMessages({
              p.loop <- p.loop + 
                ggplot2::scale_x_continuous(limits = c(min(p.loop$layers[[1]]$data$x), 
                                                       max(p.loop$layers[[1]]$data$x))) + 
                ggplot2::scale_y_continuous(limits = c(min(p.loop$layers[[1]]$data$y), 
                                                       max(p.loop$layers[[1]]$data$y)))
            })
            
            if (!(is.null(group.by))){
              if (isTRUE(group.by.show.dots)){
                p.loop$layers <- append(p.loop$layers, c(center_layer_2, center_layer_1))
              }
              if (isTRUE(group.by.cell_borders)){
                p.loop$layers <- append(color_layer, p.loop$layers)
              }
              if (is.null(group.by.legend)){
                legend.name <- group.by
              } else {
                if (is.na(group.by.legend)){
                  legend.name <- ""
                } else {
                  legend.name <- group.by.legend
                }
              }
              p.loop <- p.loop +
                ggplot2::scale_fill_manual(values = group.by.colors.use,
                                           name = legend.name) +
                ggplot2::guides(fill = ggplot2::guide_legend(title = legend.name,
                                                             title.position = "top",
                                                             title.hjust = 0.5,
                                                             ncol = legend.ncol,
                                                             nrow = legend.nrow,
                                                             byrow = legend.byrow,
                                                             override.aes = list(alpha = 1)))
            }
          }

          if (isTRUE(plot_density_contour)){
            data <- ggplot2::ggplot_build(p.loop)

            density_layer <- ggplot2::stat_density_2d(data = data$data[[1]],
                                                      mapping = ggplot2::aes(x = .data$x,
                                                                             y = .data$y),
                                                      color = contour.color,
                                                      lineend = contour.lineend,
                                                      linejoin = contour.linejoin)
            if (contour.position == "bottom"){
              p.loop$layers <- append(density_layer, p.loop$layers)
            } else if (contour.position == "top"){
              p.loop$layers <- append(p.loop$layers, density_layer)
            }

            min_x <- min(data$data[[1]]$x) * (1 + contour_expand_axes)
            max_x <- max(data$data[[1]]$x) * (1 + contour_expand_axes)
            min_y <- min(data$data[[1]]$y) * (1 + contour_expand_axes)
            max_y <- max(data$data[[1]]$y) * (1 + contour_expand_axes)
            # Expand axes limits to allocate the new contours.
            suppressMessages({
              p.loop <- p.loop +
                ggplot2::xlim(c(min_x, max_x)) +
                ggplot2::ylim(c(min_y, max_y))
            })

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
      if (stringr::str_starts(reduction, "diff|DIFF")){
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
       ggplot2::theme(plot.margin = ggplot2::margin(t = 0, r = 0, b = 0, l = 0),
                      plot.title = ggplot2::element_text(face = plot.title.face,
                                                         hjust = ifelse(!(is.null(split.by)), 0.5, 0)),
                      plot.subtitle = ggplot2::element_text(face = plot.subtitle.face, hjust = 0),
                      plot.caption = ggplot2::element_text(face = plot.caption.face, hjust = 1),
                      panel.grid = ggplot2::element_blank(),
                      plot.title.position = "plot",
                      plot.caption.position = "plot",
                      text = ggplot2::element_text(family = font.type),
                      legend.position = legend.position,
                      legend.text = ggplot2::element_text(face = legend.text.face),
                      legend.title = ggplot2::element_text(face = legend.title.face),
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
    for (counter in seq(1, length(features))){
      if (!(is.na(individual.titles[counter]))){
        if (is.null(split.by)){
          p[[counter]]$labels$title <- individual.titles[counter]
        }
      }
    }
  }

  # Add individual subtitles.
  if (!is.null(individual.subtitles)){
    for (counter in seq(1, length(features))){
      if (!(is.na(individual.subtitles[counter]))){
        if (is.null(split.by)){
          p[[counter]]$labels$subtitle <- individual.subtitles[counter]
        }
      }
    }
  }

  # Add individual captions
  if (!is.null(individual.captions)){
    for (counter in seq(1, length(features))){
      if (!(is.na(individual.captions[counter]))){
        if (is.null(split.by)){
          p[[counter]]$labels$caption <- individual.captions[counter]
        }
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
    labels <- colnames(sample@reductions[[reduction]][[]])[dims]
    # Fix the axis scale so that the highest and lowest values are in the range of the DCs (previously was around +-1.5, while DCs might range to +-0.004 or so).
    p <-  suppressMessages({
      p &
        ggplot2::xlim(c(min(sample@reductions[[reduction]][[]][, labels[1]], na.rm = TRUE),
                        max(sample@reductions[[reduction]][[]][, labels[1]], na.rm = TRUE))) &
        ggplot2::ylim(c(min(sample@reductions[[reduction]][[]][, labels[2]], na.rm = TRUE),
                        max(sample@reductions[[reduction]][[]][, labels[2]], na.rm = TRUE)))
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
