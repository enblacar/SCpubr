#' Wrapper for \link[Seurat]{DimPlot}.
#'
#' @inheritParams doc_function
#' @param idents.keep \strong{\code{\link[base]{character}}} | Vector of identities to keep. This will effectively set the rest of the cells that do not match the identities provided to NA, therefore coloring them according to na.value parameter.
#' @param shuffle \strong{\code{\link[base]{logical}}} | Whether to shuffle the cells or not, so that they are not plotted cluster-wise. Recommended.
#' @param order \strong{\code{\link[base]{character}}} | Vector of identities to be plotted. Either one with all identities or just some, which will be plotted last.
#' @param sizes.highlight \strong{\code{\link[base]{numeric}}} | Point size of highlighted cells using cells.highlight parameter.
#' @return  A ggplot2 object containing a DimPlot.
#' @export
#'
#' @example man/examples/examples_do_DimPlot.R
do_DimPlot <- function(sample,
                       reduction = NULL,
                       group.by = NULL,
                       split.by = NULL,
                       colors.use = NULL,
                       shuffle = TRUE,
                       order = NULL,
                       raster = FALSE,
                       pt.size = 1,
                       label = FALSE,
                       label.color = "black",
                       label.fill = "white",
                       label.size = 4,
                       label.box = TRUE,
                       repel = FALSE,
                       cells.highlight = NULL,
                       idents.highlight = NULL,
                       idents.keep = NULL,
                       sizes.highlight = 1,
                       ncol = NULL,
                       plot.title = NULL,
                       plot.subtitle = NULL,
                       plot.caption = NULL,
                       legend.title = NULL,
                       legend.position = "bottom",
                       legend.title.position = "top",
                       legend.ncol = NULL,
                       legend.nrow = NULL,
                       legend.icon.size = 4,
                       legend.byrow = FALSE,
                       raster.dpi = 2048,
                       dims = c(1, 2),
                       font.size = 14,
                       font.type = "sans",
                       na.value = "grey75",
                       plot_cell_borders = TRUE,
                       border.size = 2,
                       border.color = "black",
                       border.density = 1,
                       plot_marginal_distributions = FALSE,
                       marginal.type = "density",
                       marginal.size = 5,
                       marginal.group = TRUE,
                       plot.axes = FALSE,
                       plot_density_contour = FALSE,
                       contour.position = "bottom",
                       contour.color = "grey90",
                       contour.lineend = "butt",
                       contour.linejoin = "round",
                       contour_expand_axes = 0.25,
                       plot.title.face = "bold",
                       plot.subtitle.face = "plain",
                       plot.caption.face = "italic",
                       axis.title.face = "bold",
                       axis.text.face = "plain",
                       legend.title.face = "bold",
                       legend.text.face = "plain"){
  # Add lengthy error messages.
  withr::local_options(.new = list("warning.length" = 8170))
  
  check_suggests(function_name = "do_DimPlot")
  # Check if the sample provided is a Seurat object.
  check_Seurat(sample = sample)
  
  sample <- check_Assay5(sample)
  
  # Check the reduction.
  reduction <- check_and_set_reduction(sample = sample, reduction = reduction)
  # Check the dimensions.
  dims <- check_and_set_dimensions(sample = sample, reduction = reduction, dims = dims)
  # Check logical parameters.
  logical_list <- list("label" = label,
                       "repel" = repel,
                       "shuffle" = shuffle,
                       "legend.byrow" = legend.byrow,
                       "raster" = raster,
                       "plot_marginal_distributions" = plot_marginal_distributions,
                       "marginal.group" = marginal.group,
                       "plot_cell_borders" = plot_cell_borders,
                       "plot.axes" = plot.axes,
                       "plot_density_contour" = plot_density_contour,
                       "label.box" = label.box)
  check_type(parameters = logical_list, required_type = "logical", test_function = is.logical)
  # Check numeric parameters.
  numeric_list <- list("pt.size" = pt.size,
                       "sizes.highlight" = sizes.highlight,
                       "legend.ncol" = legend.ncol,
                       "legend.nrow" = legend.nrow,
                       "font.size" = font.size,
                       "legend.icon.size" = legend.icon.size,
                       "ncol" = ncol,
                       "raster.dpi" = raster.dpi,
                       "marginal.size" = marginal.size,
                       "border.size" = border.size,
                       "contour_expand_axes" = contour_expand_axes,
                       "label.size" = label.size,
                       "border.density" = border.density)
  check_type(parameters = numeric_list, required_type = "numeric", test_function = is.numeric)
  # Check character parameters.
  character_list <- list("legend.position" = legend.position,
                         "plot.title" = plot.title,
                         "plot.subtitle" = plot.subtitle,
                         "plot.caption" = plot.caption,
                         "legend.title" = legend.title,
                         "cells.highlight" = cells.highlight,
                         "idents.keep" = idents.keep,
                         "order" = order,
                         "na.value" = na.value,
                         "idents.highlight" = idents.highlight,
                         "legend.title.position" = legend.title.position,
                         "font.type" = font.type,
                         "marginal.type" = marginal.type,
                         "border.color" = border.color,
                         "contour.position" = contour.position,
                         "contour.color" = contour.color,
                         "contour.lineend" = contour.lineend,
                         "contour.linejoin" = contour.linejoin,
                         "plot.title.face" = plot.title.face,
                         "plot.subtitle.face" = plot.subtitle.face,
                         "plot.caption.face" = plot.caption.face,
                         "axis.title.face" = axis.title.face,
                         "axis.text.face" = axis.text.face,
                         "legend.title.face" = legend.title.face,
                         "legend.text.face" = legend.text.face)
  check_type(parameters = character_list, required_type = "character", test_function = is.character)

  # Checks to ensure proper function.
  group_by_and_split_by_used <- !(is.null(split.by)) & !(is.null(group.by))
  group_by_and_highlighting_cells <- (!(is.null(cells.highlight)) | !(is.null(idents.highlight))) & !(is.null(group.by))
  split_by_and_highlighting_cells <- (!(is.null(cells.highlight)) | !(is.null(idents.highlight))) & !(is.null(split.by))
  order_and_shuffle_used <- !(is.null(order)) & isTRUE(shuffle)

  assertthat::assert_that(!group_by_and_highlighting_cells,
                          msg = paste0(add_cross(), crayon_body("Either "),
                                       crayon_key("group.by"),
                                       crayon_body(" or "),
                                       crayon_key("cells.highlight | idents.hightlight"),
                                       crayon_body(" have to be set to "),
                                       crayon_key("NULL"),
                                       crayon_body(".")))

  assertthat::assert_that(!split_by_and_highlighting_cells,
                          msg = paste0(add_cross(), crayon_body("Either "),
                                       crayon_key("split.by"),
                                       crayon_body(" or "),
                                       crayon_key("cells.highlight | idents.hightlight"),
                                       crayon_body(" have to be set to "),
                                       crayon_key("NULL"),
                                       crayon_body(".")))

  if (order_and_shuffle_used){
    warning(paste0(add_warning(), crayon_body("Setting up a custom order with paramter "),
                   crayon_key("order"),
                   crayon_body(" when "),
                   crayon_key("shuffle = TRUE"),
                   crayon_body(" might result in unexpected behaviors.\nPlease, consider using it alongside "),
                   crayon_key("shuffle = FALSE"),
                   crayon_body(".")), call. = FALSE)
  }

  # Check for label.color.
  ## Check for the colors assigned to the labels if label = TRUE.
  check_colors(label.color, parameter_name = "label.color")
  ## Check the color assigned to NAs.
  check_colors(na.value, parameter_name = "na.value")
  ## Check the color assigned to border.color.
  check_colors(border.color, parameter_name = "border.color")
  ## Check the color assigned to contour.color.
  check_colors(contour.color, parameter_name = "contour.color")

  ## If the user provides more than one color to na.value, stop the function.
  assertthat::assert_that(length(na.value) == 1,
                          msg = paste0(add_cross(), crayon_body("Please, provide only "),
                                       crayon_key("one color"),
                                       crayon_body(" to parameter "),
                                       crayon_key("na.value"),
                                       crayon_body(".")))

  ## Check that the contour_expand_axes is between 0 and 1.
  assertthat::assert_that(contour_expand_axes <= 1,
                          msg = paste0(add_cross(), crayon_body("Please, provide a value "),
                                       crayon_key("lower or equal to 1"),
                                       crayon_body(" to parameter "),
                                       crayon_key("contour_expand_axes"),
                                       crayon_body(".")))

  assertthat::assert_that(contour_expand_axes >= 0,
                          msg = paste0(add_cross(), crayon_body("Please, provide a value "),
                                       crayon_key("lower or equal to 1"),
                                       crayon_body(" to parameter "),
                                       crayon_key("contour_expand_axes"),
                                       crayon_body(".")))

  # If the user provides raster = TRUE but the pt.size is less than 1, warn it.
  if (isTRUE(raster) & pt.size < 1){
    warning(paste0(add_warning(), crayon_body("Setting "),
                   crayon_key("raster = TRUE"),
                   crayon_body(" and "),
                   crayon_key("pt.size < 1"),
                   crayon_body("will result in the cells being plotted as a "),
                   crayon_key("cross"),
                   crayon_body(" instead of dots.\nThis behaviour can not be modified, but can be avoided by using "),
                   crayon_key("pt.size >= 1"),
                   crayon_body(".")), call. = FALSE)
  }

  check_parameters(parameter = font.type, parameter_name = "font.type")
  check_parameters(parameter = legend.position, parameter_name = "legend.position")
  check_parameters(parameter = marginal.type, parameter_name = "marginal.type")
  check_parameters(parameter = contour.lineend, parameter_name = "contour.lineend")
  check_parameters(parameter = contour.linejoin, parameter_name = "contour.linejoin")
  check_parameters(parameter = contour.position, parameter_name = "contour.position")
  check_parameters(parameter = border.density, parameter_name = "border.density")
  check_parameters(plot.title.face, parameter_name = "plot.title.face")
  check_parameters(plot.subtitle.face, parameter_name = "plot.subtitle.face")
  check_parameters(plot.caption.face, parameter_name = "plot.caption.face")
  check_parameters(axis.title.face, parameter_name = "axis.title.face")
  check_parameters(axis.text.face, parameter_name = "axis.text.face")
  check_parameters(legend.title.face, parameter_name = "legend.title.face")
  check_parameters(legend.text.face, parameter_name = "legend.text.face")
  
  # If the user has not provided colors.
  if (is.null(colors.use)){
    colors.use <- {
      # Default parameters.
      default_parameters <- is.null(group.by) & is.null(split.by) & is.null(cells.highlight) & is.null(idents.highlight)
      # Group.by was used.
      group_by_is_used <- !(is.null(group.by)) & is.null(split.by) & is.null(cells.highlight) & is.null(idents.highlight)
      # Split.by was used.
      split_by_is_used <- is.null(group.by) & !(is.null(split.by)) & is.null(cells.highlight) & is.null(idents.highlight)
      # Cells.highlight or idents.highlight was used.
      highlighting_cells <- is.null(group.by) & is.null(split.by) & (!(is.null(cells.highlight)) | !(is.null(idents.highlight)))
      if (isTRUE(default_parameters)){
        # Generate the color scale based on the levels assigned to the sample.
        generate_color_scale(levels(sample))
      } else if (isTRUE(group_by_is_used) | isTRUE(group_by_and_split_by_used)){
        # Retrieve the unique values in group.by metadata variable.
        data.use <- sample[[]][, group.by, drop = FALSE]
        # If the variable is a factor, use the levels as order. If not, order the values alphabetically.
        names.use <- if (is.factor(data.use[, 1])){levels(data.use[, 1])} else {sort(unique(data.use[, 1]))}
        # Generate the color scale to be used based on the unique values of group.by.
        generate_color_scale(names.use)
      } else if (isTRUE(split_by_is_used)){
        # Retrieve the unique values in split.by metadata variable.
        data.use <- sample[[]][, split.by, drop = FALSE]
        # If the variable is a factor, use the levels as order. If not, order the values alphabetically.
        names.use <- if (is.factor(data.use[, 1])){levels(data.use[, 1])} else {sort(unique(data.use[, 1]))}
        # Generate the color scale based on the unique values of split.by
        generate_color_scale(names.use)
      } else if (isTRUE(highlighting_cells)){
        # If the user wants to highlight some cells, use this color.
        colors.use <- "#0A305F"
      }
    }
    # For split.by + group.by cases.
    colors.use.original <- colors.use
    # But, if the user has provided a custom color palette.
  } else {
    # Check that the provided values are valid color representations.
    check_colors(colors.use, parameter_name = "colors.use")
    # If no further parameters are used.
    default_parameters <- is.null(group.by) & is.null(split.by) & is.null(cells.highlight) & is.null(idents.highlight)
    # Group.by was used.
    group_by_is_used <- !(is.null(group.by)) & is.null(split.by) & is.null(cells.highlight) & is.null(idents.highlight)
    # Split.by was used.
    split_by_is_used <- is.null(group.by) & !(is.null(split.by)) & is.null(cells.highlight) & is.null(idents.highlight)
    # When either cells.highlight or idents.highlight was used.
    highlighting_cells <- is.null(group.by) & is.null(split.by) & (!(is.null(cells.highlight)) | !(is.null(idents.highlight)))
    
    # For split.by + group.by cases.
    colors.use.original <- colors.use
    
    # When running under default parameters.
    if (isTRUE(default_parameters)){
      # Check that the color palette has the right amount of named colors with regards to the current identities.
      colors.use <- check_consistency_colors_and_names(sample = sample,
                                                       colors = colors.use,
                                                       idents.keep = idents.keep)
      # When using group.by or a combination of group.by and split.by.
    } else if (isTRUE(group_by_is_used) | isTRUE(group_by_and_split_by_used)){
      # Check that the color palette has the right amount of named colors with regards to group.by values.
      colors.use <- check_consistency_colors_and_names(sample = sample,
                                                       colors = colors.use,
                                                       grouping_variable = group.by,
                                                       idents.keep = idents.keep)
      # When using split.by.
    } else if (isTRUE(split_by_is_used)){
      # Check that the color palette has the right amount of named colors with regards to split.by values.
      colors.use <- check_consistency_colors_and_names(sample = sample,
                                                       colors = colors.use,
                                                       grouping_variable = split.by,
                                                       idents.keep = idents.keep)
    
      # When highlighting cells.
    } else if (isTRUE(highlighting_cells)){
      # Stop the execution if more than one color is provided to highlight the cells.
      assertthat::assert_that(length(colors.use) == 1,
                              msg = paste0(add_cross(), crayon_body("Please, provide only "),
                                           crayon_key("one color"),
                                           crayon_body(" to "),
                                           crayon_key("cells.highlight"),
                                           crayon_body(" or "),
                                           crayon_key("idents.highlight"),
                                           crayon_body(".")))
    }
  }
  
  # Compute the colors for the labels.
  if (isTRUE(label)){
    if (isTRUE(label.box)){
      if (is.null(label.fill)){
        colors.use.label.fill <- colors.use
      } else {
        # Check that only one color has been provided to label.fill.
        assertthat::assert_that(length(label.fill) == 1,
                                msg = paste0(add_cross(), crayon_body("Please, provide only "),
                                             crayon_key("one color"),
                                             crayon_body(" to "),
                                             crayon_key("label.fill"),
                                             crayon_body(" or "),
                                             crayon_key("NULL"),
                                             crayon_body(".")))
        
        # And check that is a valid color.
        check_colors(label.fill, parameter_name = "label.fill")
        
        colors.use.label.fill <- rep(label.fill, length(colors.use))
      }
    }
  }

  # Set cells to NA according to idents.keep.
  # If the user does not want to highlight cells or split by identities but wants to remove some identities.
  idents_keep_used <- is.null(cells.highlight) & is.null(idents.highlight) & !(is.null(idents.keep))
  if (isTRUE(idents_keep_used)){
    # CONDITION: both group.by and split.by are not used.
    group_by_and_split_by_are_null <- is.null(group.by) & is.null(split.by)
    # CONDITION: group.by is used.
    group_by_is_used <- !(is.null(group.by)) & is.null(split.by)
    # CONDITION: split.by is used.
    split_by_is_used <- is.null(group.by) & !(is.null(split.by))
    # When running under default parameters.
    if (isTRUE(group_by_and_split_by_are_null)){
      # Check that idents.keep matches the values and if not, stop the execution.
      assertthat::assert_that(isTRUE(length(idents.keep) == sum(idents.keep %in% levels(sample))),
                              msg = paste0(add_cross(), crayon_body("All the values in "),
                                           crayon_key("idents.keep"),
                                           crayon_body(" must be in "),
                                           crayon_key("levels(sample"),
                                           crayon_body(".")))
      # Set the identities that the user wants to exclude as NA.
      Seurat::Idents(sample)[!(Seurat::Idents(sample) %in% idents.keep)] <- NA

      colors.use <- check_consistency_colors_and_names(sample = sample, 
                                                       colors = colors.use,
                                                       idents.keep = idents.keep)
      # If split.by is used instead.
    } else if (group_by_and_split_by_used){
      # Check that the values in idents.keep are in the unique values of split.by.
      assertthat::assert_that(isTRUE(length(idents.keep) == sum(idents.keep %in% unique(sample@meta.data[, split.by]))),
                              msg = paste0(add_cross(), crayon_body("All the values in "),
                                           crayon_key("idents.keep"),
                                           crayon_body(" must be in the "),
                                           crayon_key("split.by"),
                                           crayon_body(" metadata provided.")))
      
      colors.use <- check_consistency_colors_and_names(sample = sample, 
                                                       colors = colors.use, 
                                                       grouping_variable = group.by)
    } else if (split_by_is_used){
      # Check that the values in idents.keep are in the unique values of split.by.
      assertthat::assert_that(isTRUE(length(idents.keep) == sum(idents.keep %in% unique(sample@meta.data[, split.by]))),
                              msg = paste0(add_cross(), crayon_body("All the values in "),
                                           crayon_key("idents.keep"),
                                           crayon_body(" must be in the "),
                                           crayon_key("split.by"),
                                           crayon_body(" metadata provided.")))
      
      colors.use <- check_consistency_colors_and_names(sample = sample, 
                                                       colors = colors.use, 
                                                       grouping_variable = split.by,
                                                       idents.keep = idents.keep)

      # When using group.by, check with the values in group.by.
    } else if (group_by_is_used) {
      # Check that idents.keep matches the values, if not, stop the execution.
      assertthat::assert_that(isTRUE(length(idents.keep) == sum(idents.keep %in% unique(sample@meta.data[, group.by]))),
                              msg = paste0(add_cross(), crayon_body("All the values in "),
                                           crayon_key("idents.keep"),
                                           crayon_body(" must be in the "),
                                           crayon_key("group.by"),
                                           crayon_body(" metadata variable provided.")))
      # Convert to NA values in group.by not included in the user's selected values.
      sample@meta.data[, group.by][!(sample@meta.data[, group.by] %in% idents.keep)] <- NA
      colors.use <- check_consistency_colors_and_names(sample = sample, 
                                                       colors = colors.use, 
                                                       grouping_variable = group.by,
                                                       idents.keep = idents.keep)
      
    }
  }


  # Generate base layer.
    out <- compute_umap_layer(sample = sample,
                              labels = colnames(sample@reductions[[reduction]][[]])[dims],
                              pt.size = pt.size,
                              border.density = border.density,
                              border.size = border.size,
                              border.color = border.color,
                              raster = raster,
                              raster.dpi = raster.dpi,
                              reduction = reduction,
                              group.by = group.by,
                              split.by = split.by,
                              na.value = na.value,
                              n = 100)
    base_layer <- out$base_layer
    na_layer <- out$na_layer

  # PLOTTING

  # If raster = TRUE, add 1 to pt.size to keep consistency between plots.

  # If the UMAP does not need to be split in multiple panes (default case).
  # CONDITION: Not highligting cells and not using split.by.
  not_highlighting_and_not_split_by <- is.null(cells.highlight) & is.null(idents.highlight) & is.null(split.by)
  # CONDITION: Using split.by.
  split_by_used <- is.null(cells.highlight) & is.null(idents.highlight) & !(is.null(split.by))
  # CONDITION: highlighting cells.
  highlighting_cells <- !(is.null(cells.highlight)) | !(is.null(idents.highlight))
  # When running under default parameters or using group.by
  if (not_highlighting_and_not_split_by){
    if (utils::packageVersion("Seurat") >= "4.1.0"){
      p <- Seurat::DimPlot(if (is.null(idents.keep)) {sample} else {if (is.null(group.by)) {sample[, Seurat::Idents(sample) %in% idents.keep]} else {sample[, sample@meta.data[, group.by] %in% idents.keep]}},
                           reduction = reduction,
                           label = label,
                           dims = dims,
                           repel = repel,
                           label.box = label.box,
                           label.color = label.color,
                           label.size = label.size,
                           na.value = na.value,
                           shuffle = shuffle,
                           order = order,
                           pt.size = pt.size,
                           group.by = group.by,
                           cols = colors.use,
                           raster = raster,
                           raster.dpi = c(raster.dpi, raster.dpi),
                           ncol = ncol)
    } else { # nocov start
      p <- Seurat::DimPlot(if (is.null(idents.keep)) {sample} else {if (is.null(group.by)) {sample[, Seurat::Idents(sample) %in% idents.keep]} else {sample[, sample@meta.data[, group.by] %in% idents.keep]}},
                           reduction = reduction,
                           label = label,
                           dims = dims,
                           repel = repel,
                           label.box = label.box,
                           label.color = label.color,
                           label.size = label.size,
                           na.value = na.value,
                           shuffle = shuffle,
                           order = order,
                           pt.size = pt.size,
                           group.by = group.by,
                           cols = colors.use,
                           raster = raster,
                           ncol = ncol)
    } # nocov end
    p <- p &
      ggplot2::guides(color = ggplot2::guide_legend(ncol = legend.ncol,
                                                    nrow = legend.nrow,
                                                    byrow = legend.byrow,
                                                    override.aes = list(size = legend.icon.size),
                                                    title = legend.title,
                                                    title.position = legend.title.position))

    if (isTRUE(label)){
      if (isTRUE(label.box)){
        if (is.null(label.fill)){
          colors.use.label.fill <- colors.use
        } else {
          colors.use.label.fill <- rep(label.fill, length(colors.use))
        }
        p <- add_scale(p = p,
                       function_use = ggplot2::scale_fill_manual(values = colors.use.label.fill),
                       scale = "fill")
        
      }
      p$layers[[length(p$layers)]]$aes_params$fontface <- "bold"
    }

    if (!(is.null(group.by))){
      # Remove automatic title inserted by Seurat.
      p <- p & ggplot2::ggtitle("")
    }
    
    
    # Add another layer of black dots to make the colored ones stand up.
    if (!is.null(idents.keep)){
      if (isTRUE(plot_cell_borders)){
        sample.use <- if (is.null(group.by)) {sample[, Seurat::Idents(sample) %in% idents.keep]} else {sample[, sample@meta.data[, group.by] %in% idents.keep]}
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
                                  na.value = na.value,
                                  n = 100)
        base_layer.subset <- out$base_layer
        p$layers <- append(base_layer.subset, p$layers)
      }
      # Add NA layer.
      p$layers <- append(na_layer, p$layers)
    }
    
    
    # Add cell borders.
    if (isTRUE(plot_cell_borders)){
      p$layers <- append(base_layer, p$layers)
    }

    if (isTRUE(plot_density_contour)){
      data <- ggplot2::ggplot_build(p)

      density_layer <- ggplot2::stat_density_2d(data = data$data[[1]],
                                                mapping = ggplot2::aes(x = .data$x,
                                                                       y = .data$y),
                                                color = contour.color,
                                                lineend = contour.lineend,
                                                linejoin = contour.linejoin)
      if (contour.position == "bottom"){
        p$layers <- append(density_layer, p$layers)
      } else if (contour.position == "top"){
        p$layers <- append(p$layers, density_layer)
      }

      min_x <- min(data$data[[1]]$x) * (1 + contour_expand_axes)
      max_x <- max(data$data[[1]]$x) * (1 + contour_expand_axes)
      min_y <- min(data$data[[1]]$y) * (1 + contour_expand_axes)
      max_y <- max(data$data[[1]]$y) * (1 + contour_expand_axes)
      # Expand axes limits to allocate the new contours.
      suppressMessages({
        p <- p +
             ggplot2::xlim(c(min_x, max_x)) +
             ggplot2::ylim(c(min_y, max_y))
      })
      
    }
    # Add theme settings to all plots.
    p <- p &
         ggplot2::theme_minimal(base_size = font.size) &
         ggplot2::theme(plot.title = ggplot2::element_text(face = plot.title.face, hjust = 0),
                        plot.subtitle = ggplot2::element_text(face = plot.subtitle.face, hjust = 0),
                        plot.caption = ggplot2::element_text(face = plot.caption.face, hjust = 1),
                        plot.title.position = "plot",
                        plot.caption.position = "plot",
                        text = ggplot2::element_text(family = font.type),
                        legend.justification = "center",
                        legend.text = ggplot2::element_text(face = legend.text.face),
                        legend.title = if (legend.position != "none") {ggplot2::element_text(face = legend.title.face)} else {ggplot2::element_blank()},
                        legend.position = legend.position,
                        panel.grid = ggplot2::element_blank(),
                        plot.margin = ggplot2::margin(t = 0, r = 0, b = 0, l = 0),
                        plot.background = ggplot2::element_rect(fill = "white", color = "white"),
                        panel.background = ggplot2::element_rect(fill = "white", color = "white"),
                        legend.background = ggplot2::element_rect(fill = "white", color = "white"))
  } else if (isTRUE(group_by_and_split_by_used) | isTRUE(split_by_used)){
    list.plots <- list()
    unique_values <- if(is.factor(sample@meta.data[, split.by])){levels(sample@meta.data[, split.by])} else {sort(unique(sample@meta.data[, split.by]))}
    if (!is.null(idents.keep)){
      unique_values <- unique_values[unique_values %in% idents.keep]
    }
    
    # If group.by and split.by are used, add a general view that will have a legend with all parameters.
    p.extra <- do_DimPlot(sample = sample,
                          reduction = reduction,
                          group.by = ifelse(isTRUE(group_by_and_split_by_used), group.by, split.by),
                          split.by = NULL,
                          colors.use = colors.use.original,
                          shuffle = shuffle,
                          order = order,
                          raster = raster,
                          pt.size = pt.size,
                          label = label,
                          label.color = label.color,
                          label.fill = label.fill,
                          label.size = label.size,
                          label.box = label.box,
                          repel = repel,
                          cells.highlight = NULL,
                          idents.highlight = NULL,
                          idents.keep = NULL,
                          sizes.highlight = sizes.highlight,
                          ncol = ncol,
                          plot.title = "Combined",
                          plot.subtitle = NULL,
                          plot.caption = NULL,
                          legend.title = legend.title,
                          legend.position = legend.position,
                          legend.title.position = legend.title.position,
                          legend.ncol = legend.ncol,
                          legend.nrow = legend.nrow,
                          legend.icon.size = legend.icon.size,
                          legend.byrow = legend.byrow,
                          raster.dpi = raster.dpi,
                          dims = dims,
                          font.size = font.size,
                          font.type = font.type,
                          na.value = na.value,
                          plot_cell_borders = plot_cell_borders,
                          border.size = border.size,
                          border.color = border.color,
                          border.density = border.density,
                          plot_marginal_distributions = plot_marginal_distributions,
                          marginal.type = marginal.type,
                          marginal.size = marginal.size,
                          marginal.group = marginal.group,
                          plot.axes = plot.axes,
                          plot_density_contour = plot_density_contour,
                          contour.position = contour.position,
                          contour.color = contour.color,
                          contour.lineend = contour.lineend,
                          contour.linejoin = contour.linejoin,
                          contour_expand_axes = contour_expand_axes,
                          plot.title.face = plot.title.face,
                          plot.subtitle.face = plot.subtitle.face,
                          plot.caption.face = plot.caption.face,
                          axis.title.face = axis.title.face,
                          axis.text.face = axis.text.face,
                          legend.title.face = legend.title.face,
                          legend.text.face = legend.text.face)
    p.extra <- p.extra + ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
    # Add plot to list.
    list.plots[[1]] <- p.extra

    num_values <- length(unique_values)
    for (i in seq_len(num_values)){
      value <- unique_values[i]
      # Generate a middle layer for the missing values after split.by.
      sample.use <- sample[, sample@meta.data[, split.by] == value]

      if (utils::packageVersion("Seurat") >= "4.1.0"){
        p.loop <- Seurat::DimPlot(sample.use,
                                  reduction = reduction,
                                  group.by = if (is.null(group.by)) {split.by} else {group.by},
                                  label = label,
                                  dims = dims,
                                  repel = repel,
                                  label.box = label.box,
                                  label.color = label.color,
                                  label.size = label.size,
                                  na.value = na.value,
                                  shuffle = shuffle,
                                  order = order,
                                  pt.size = pt.size,
                                  cols = colors.use,
                                  raster = raster,
                                  raster.dpi = c(raster.dpi, raster.dpi))
      } else { # nocov start
        p.loop <- Seurat::DimPlot(sample.use,
                                  reduction = reduction,
                                  group.by = if (is.null(group.by)) {split.by} else {group.by},
                                  label = label,
                                  dims = dims,
                                  repel = repel,
                                  label.box = label.box,
                                  label.color = label.color,
                                  label.size = label.size,
                                  na.value = na.value,
                                  shuffle = shuffle,
                                  order = order,
                                  pt.size = pt.size,
                                  cols = colors.use,
                                  raster = raster)
      } # nocov end
      p.loop <- p.loop +
                ggplot2::ggtitle(value) +
                ggplot2::guides(color = ggplot2::guide_legend(title = legend.title,
                                                              ncol = legend.ncol,
                                                              nrow = legend.nrow,
                                                              byrow = legend.byrow,
                                                              override.aes = list(size = legend.icon.size),
                                                              title.position = legend.title.position))
      if (isTRUE(label)){
        if (isTRUE(label.box)){
          p.loop <- add_scale(p = p.loop,
                              function_use = ggplot2::scale_fill_manual(values = colors.use.label.fill),
                              scale = "fill")
        }
        p.loop$layers[[length(p.loop$layers)]]$aes_params$fontface <- "bold"
      }

      # Add another layer of black dots to make the colored ones stand up.
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
                                  na.value = na.value,
                                  n = 100)
        base_layer.subset <- out$base_layer
        p.loop$layers <- append(base_layer.subset, p.loop$layers)
      }
      # Add NA layer.
      p.loop$layers <- append(na_layer, p.loop$layers)

      # Add cell borders.
      if (isTRUE(plot_cell_borders)){
        p.loop$layers <- append(base_layer, p.loop$layers)
        suppressMessages({
          p.loop <- p.loop + 
                    ggplot2::scale_x_continuous(limits = c(min(p.loop$layers[[1]]$data$x), 
                                                           max(p.loop$layers[[1]]$data$x))) + 
                    ggplot2::scale_y_continuous(limits = c(min(p.loop$layers[[1]]$data$y), 
                                                           max(p.loop$layers[[1]]$data$y)))
        })
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
      
      p.loop <- p.loop +
                ggplot2::theme_minimal(base_size = font.size) +
                ggplot2::theme(plot.title = ggplot2::element_text(face = plot.title.face, hjust = 0.5),
                               plot.subtitle = ggplot2::element_text(face = plot.subtitle.face, hjust = 0),
                               plot.caption = ggplot2::element_text(face = plot.caption.face, hjust = 1),
                               plot.title.position = "plot",
                               plot.caption.position = "plot",
                               text = ggplot2::element_text(family = font.type),
                               legend.justification = "center",
                               legend.text = ggplot2::element_text(face = legend.text.face),
                               legend.title = if (legend.position != "none") {ggplot2::element_text(face = legend.title.face)} else {ggplot2::element_blank()},
                               panel.grid = ggplot2::element_blank(),
                               plot.margin = ggplot2::margin(t = 0, r = 0, b = 0, l = 0),
                               plot.background = ggplot2::element_rect(fill = "white", color = "white"),
                               panel.background = ggplot2::element_rect(fill = "white", color = "white"),
                               legend.background = ggplot2::element_rect(fill = "white", color = "white"),
                               legend.position = "none")

      list.plots[[value]] <- p.loop
    }
    
    
    p <- patchwork::wrap_plots(list.plots, ncol = ncol, guides = "collect") 
    
    p <- p + 
         patchwork::plot_annotation(theme = ggplot2::theme(legend.position = legend.position))
    
  }
   # If the user wants to highlight some of the cells.
  else if (highlighting_cells){
    # Compute the cells to highlight.
    if (is.null(idents.highlight) & !(is.null(cells.highlight))){
      # Only if cells.highlight parameters is used.
      cells.use <- cells.highlight
    } else if (!(is.null(idents.highlight)) & is.null(cells.highlight)){
      # Only if idents.highlight parameter is used.
      # Check if the provided identities are part of the active identities in the object.
      check_identity(sample = sample, identities = idents.highlight)
      cells.use <- names(Seurat::Idents(sample)[Seurat::Idents(sample) %in% idents.highlight])
    } else if (!(is.null(idents.highlight)) & !(is.null(cells.highlight))){
      # Check if the provided identities are part of the active identities in the object.
      check_identity(sample = sample, identities = idents.highlight)
      # Both idents.highlight and cells.highlight are used.
      cells.1 <- cells.highlight
      cells.2 <- names(Seurat::Idents(sample)[Seurat::Idents(sample) %in% idents.highlight])
      cells.use <- unique(c(cells.1, cells.2))
    }

    sample$selected_cells <- ifelse(colnames(sample) %in% cells.use, "Selected cells", NA)
    colors.use.highlight <- c("Selected cells" = colors.use)
    if (utils::packageVersion("Seurat") >= "4.1.0"){
      p <- Seurat::DimPlot(sample[, cells.use],
                           reduction = reduction,
                           group.by = "selected_cells",
                           dims = dims,
                           pt.size = pt.size,
                           raster = raster,
                           raster.dpi = c(raster.dpi, raster.dpi),
                           ncol = ncol,
                           cols = colors.use.highlight,
                           na.value = "#bfbfbf00")
    } else { # nocov start
      p <- Seurat::DimPlot(sample[, cells.use],
                           group.by = "selected_cells",
                           reduction = reduction,
                           dims = dims,
                           pt.size = pt.size,
                           raster = raster,
                           ncol = ncol,
                           cols = colors.use.highlight,
                           na.value = "#bfbfbf00")
    } # nocov end

    p <- p &
         ggplot2::ggtitle("") &
         ggplot2::guides(color = ggplot2::guide_legend(title = legend.title,
                                                       ncol = legend.ncol,
                                                       nrow = legend.nrow,
                                                       byrow = legend.byrow,
                                                       override.aes = list(size = legend.icon.size),
                                                       title.position = legend.title.position))
   
    # Add cell borders.
    if (isTRUE(plot_cell_borders)){
      # Compute extra layer for the highlighted cells.
      out <- compute_umap_layer(sample = sample[, cells.use],
                                labels = colnames(sample[, cells.use]@reductions[[reduction]][[]])[dims],
                                pt.size = sizes.highlight,
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
      p$layers <- append(base_layer_subset, p$layers)
      p$layers <- append(na_layer, p$layers)
      p$layers <- append(base_layer, p$layers)
      
      suppressMessages({
        p <- p + 
             ggplot2::scale_x_continuous(limits = c(min(p$layers[[1]]$data$x), 
                                                    max(p$layers[[1]]$data$x))) + 
             ggplot2::scale_y_continuous(limits = c(min(p$layers[[1]]$data$y), 
                                                    max(p$layers[[1]]$data$y)))
      })
    }

    if (isTRUE(plot_density_contour)){
      data <- ggplot2::ggplot_build(p)

      density_layer <- ggplot2::stat_density_2d(data = data$data[[1]],
                                                mapping = ggplot2::aes(x = .data$x,
                                                                       y = .data$y),
                                                color = contour.color,
                                                lineend = contour.lineend,
                                                linejoin = contour.linejoin)
      if (contour.position == "bottom"){
        p$layers <- append(density_layer, p$layers)
      } else if (contour.position == "top"){
        p$layers <- append(p$layers, density_layer)
      }

      min_x <- min(data$data[[1]]$x) * (1 + contour_expand_axes)
      max_x <- max(data$data[[1]]$x) * (1 + contour_expand_axes)
      min_y <- min(data$data[[1]]$y) * (1 + contour_expand_axes)
      max_y <- max(data$data[[1]]$y) * (1 + contour_expand_axes)
      # Expand axes limits to allocate the new contours.
      suppressMessages({
        p <- p +
             ggplot2::xlim(c(min_x, max_x)) +
             ggplot2::ylim(c(min_y, max_y))
      })
    }
    
    # Titles in split.by are centered by default.
    hjust_use <- if(split_by_used){0.5} else {0}
    # Add theme settings to all plots.
    p <- p &
         ggplot2::theme_minimal(base_size = font.size) &
         ggplot2::theme(plot.title = ggplot2::element_text(face = plot.title.face, hjust = hjust_use),
                        plot.subtitle = ggplot2::element_text(face = plot.subtitle.face, hjust = 0),
                        plot.caption = ggplot2::element_text(face = plot.caption.face, hjust = 1),
                        plot.title.position = "plot",
                        plot.caption.position = "plot",
                        text = ggplot2::element_text(family = font.type),
                        legend.justification = "center",
                        legend.text = ggplot2::element_text(face = legend.text.face),
                        legend.title = if (legend.position != "none") {ggplot2::element_text(face = legend.title.face)} else {ggplot2::element_blank()},
                        legend.position = legend.position,
                        panel.grid = ggplot2::element_blank(),
                        plot.margin = ggplot2::margin(t = 0, r = 0, b = 0, l = 0),
                        plot.background = ggplot2::element_rect(fill = "white", color = "white"),
                        panel.background = ggplot2::element_rect(fill = "white", color = "white"),
                        legend.background = ggplot2::element_rect(fill = "white", color = "white"))
  }


  # Add plot title to the plots.
  if (!is.null(plot.title)){
    if (!(is.null(split.by))){
      p <- p +
        patchwork::plot_annotation(title = plot.title)
    } else {
      p <- p &
        ggplot2::labs(title = plot.title)
    }
  }


  # Add plot subtitle to the plots.
  if (!is.null(plot.subtitle)){
    if (!(is.null(split.by))){
      p <- p +
        patchwork::plot_annotation(subtitle = plot.subtitle)
    } else {
      p <- p +
        ggplot2::labs(subtitle = plot.subtitle)
    }
  }

  # Add plot caption to the plots.
  if (!is.null(plot.caption)){
    if (!(is.null(split.by))){
      p <- p +
        patchwork::plot_annotation(caption = plot.caption)
    } else {
      p <- p +
        ggplot2::labs(caption = plot.caption)
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

  # Add marginal plots.
  if (not_highlighting_and_not_split_by & isTRUE(plot_marginal_distributions & base::isFALSE(plot_cell_borders))){
    # Remove annoying warnings when violin is used as marginal distribution.
    if (marginal.type == "violin"){
      p <- suppressWarnings({ggExtra::ggMarginal(p = p,
                                                 groupColour = ifelse(isTRUE(marginal.group), TRUE, FALSE),
                                                 groupFill = ifelse(isTRUE(marginal.group), TRUE, FALSE),
                                                 type = marginal.type,
                                                 size = marginal.size)})
    } else {
      p <- ggExtra::ggMarginal(p = p,
                               groupColour = ifelse(isTRUE(marginal.group), TRUE, FALSE),
                               groupFill = ifelse(isTRUE(marginal.group), TRUE, FALSE),
                               type = marginal.type,
                               size = marginal.size)
    }
    # Transform back to ggplot2 object.
    p <- ggplotify::as.ggplot(p)

    # Fix for the plot backgrounds after applying ggMarginal.
    p$theme$plot.background <- ggplot2::element_rect(fill = "white", color = "white")
    p$theme$legend.background <- ggplot2::element_rect(fill = "white", color = "white")
    p$theme$panel.background <- ggplot2::element_rect(fill = "white", color = "white")
  } else if (isTRUE(plot_marginal_distributions)) {
    stop(paste0(add_cross(), 
                crayon_body("Marginal distributions can not be used alongside when splitting by categories ("), 
                crayon_key("split.by"), 
                crayon_body("), highlighting cells ("), 
                crayon_key("cells.highlight/idents.highlight"), 
                crayon_body(") or plotting cell borders ("), 
                crayon_key("plot_cell_borders"), 
                crayon_body(").")), call. = FALSE)
  }


  # Return the final plot.
  return(p)
}
