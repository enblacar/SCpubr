#' Wrapper for \link[Seurat]{DimPlot}.
#'
#' @inheritParams doc_function
#' @param label \strong{\code{\link[base]{logical}}} | Whether to plot the cluster labels in the UMAP. The cluster labels will have the same color as the cluster colors.
#' @param idents.keep \strong{\code{\link[base]{character}}} | Vector of identities to keep. This will effectively set the rest of the cells that do not match the identities provided to NA, therefore coloring them according to na.value parameter.
#' @param shuffle \strong{\code{\link[base]{logical}}} | Whether to shuffle the cells or not, so that they are not plotted cluster-wise. Recommended.
#' @param order \strong{\code{\link[base]{character}}} | Vector of identities to be plotted. Either one with all identities or just some, which will be plotted last.
#' @param sizes.highlight \strong{\code{\link[base]{numeric}}} | Point size of highlighted cells using cells.highlight parameter.
#' @param repel \strong{\code{\link[base]{logical}}} | Whether to repel the labels if label is set to TRUE.
#' @param label.color \strong{\code{\link[base]{character}}} | HEX code for the color of the text in the labels if label is set to TRUE.
#' @return  A ggplot2 object containing a DimPlot.
#'
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
                       pt.size = 1,
                       label = FALSE,
                       label.color = "white",
                       repel = TRUE,
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
                       raster = FALSE,
                       raster.dpi = 1024,
                       dims = c(1, 2),
                       font.size = 14,
                       font.type = "sans",
                       na.value = "grey75",
                       plot_cell_borders = TRUE,
                       border.size = 2,
                       border.color = "black",
                       plot_marginal_distributions = FALSE,
                       marginal.type = "density",
                       marginal.size = 5,
                       marginal.group = TRUE,
                       plot.axes = FALSE){
  # Checks for packages.
  check_suggests(function_name = "do_DimPlot")
  # Check if the sample provided is a Seurat object.
  check_Seurat(sample = sample)
  # Check the reduction.
  reduction <- check_and_set_reduction(sample = sample, reduction = reduction)
  # Check the dimensions.
  dimensions <- check_and_set_dimensions(sample = sample, reduction = reduction, dims = dims)
  # Check logical parameters.
  logical_list <- list("label" = label,
                       "repel" = repel,
                       "shuffle" = shuffle,
                       "legend.byrow" = legend.byrow,
                       "raster" = raster,
                       "plot_marginal_distributions" = plot_marginal_distributions,
                       "marginal.group" = marginal.group,
                       "plot_cell_borders" = plot_cell_borders,
                       "plot.axes" = plot.axes)
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
                       "border.size" = border.size)
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
                         "border.color" = border.color)
  check_type(parameters = character_list, required_type = "character", test_function = is.character)

  # Checks to ensure proper function.
  group_by_and_split_by_used <- !(is.null(split.by)) & !(is.null(group.by))
  group_by_and_highlighting_cells <- (!(is.null(cells.highlight)) | !(is.null(idents.highlight))) & !(is.null(group.by))
  split_by_and_highlighting_cells <- (!(is.null(cells.highlight)) | !(is.null(idents.highlight))) & !(is.null(split.by))
  order_and_shuffle_used <- !(is.null(order)) & isTRUE(shuffle)

  assertthat::assert_that(!group_by_and_split_by_used,
                          msg = "Either group.by or split.by has to be NULL.")

  assertthat::assert_that(!group_by_and_highlighting_cells,
                          msg = "Either group.by or cells.highlight has to be NULL.")

  assertthat::assert_that(!split_by_and_highlighting_cells,
                          msg = "Either split.by or cells.highlight has to be NULL.")

  if (order_and_shuffle_used){warning("Setting up a custom order while 'shuffle = TRUE' might result in unexpected behaviours.\nPlease consider using it alongside 'shuffle = FALSE'.", call. = FALSE)}

  # Check for label.color.
  ## Check for the colors assigned to the labels if label = TRUE.
  check_colors(label.color, parameter_name = "label.color")
  ## Check the color assigned to NAs.
  check_colors(na.value, parameter_name = "na.value")
  ## Check the color assigned to border.color.
  check_colors(border.color, parameter_name = "border.color")
  ## If the user provides more than one color to na.value, stop the function.
  assertthat::assert_that(length(na.value) == 1,
                          msg = "Please provide only one color to na.value.")

  # If the user provides raster = TRUE but the pt.size is less than 1, warn it.
  if (isTRUE(raster) & pt.size < 1){
    warning("Setting raster = TRUE and pt.size < 1 will result in the cells being ploted as a cross. This behaviour can not be modified, but setting pt.size to 1 or higher solves it. For DimPlots, optimized values would be pt.size = 3 and raster.dpi = 2048.", call. = FALSE)
  }

  check_parameters(parameter = font.type, parameter_name = "font.type")
  check_parameters(parameter = legend.position, parameter_name = "legend.position")
  check_parameters(parameter = marginal.type, parameter_name = "marginal.type")

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
      if (default_parameters){
        # Generate the color scale based on the levels assigned to the sample.
        generate_color_scale(levels(sample))
      } else if (group_by_is_used){
        # Retrieve the unique values in group.by metadata variable.
        data.use <- sample[[]][, group.by, drop = FALSE]
        # If the variable is a factor, use the levels as order. If not, order the values alphabetically.
        names.use <- if (is.factor(data.use[, 1])){levels(data.use[, 1])} else {sort(unique(data.use[, 1]))}
        # Generate the color scale to be used based on the unique values of group.by.
        generate_color_scale(names.use)
      } else if (split_by_is_used){
        # Retrieve the unique values in split.by metadata variable.
        data.use <- sample[[]][, split.by, drop = FALSE]
        # If the variable is a factor, use the levels as order. If not, order the values alphabetically.
        names.use <- if (is.factor(data.use[, 1])){levels(data.use[, 1])} else {sort(unique(data.use[, 1]))}
        # Generate the color scale based on the unique values of split.by
        generate_color_scale(names.use)
      } else if (highlighting_cells){
        # If the user wants to highlight some cells, use this color.
        colors.use <- "#0A305F"
      }
    }
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
    # When running under default parameters.
    if (default_parameters){
      # Check that the color palette has the right amount of named colors with regards to the current identities.
      colors.use <- check_consistency_colors_and_names(sample = sample,
                                                       colors = colors.use)
      # When using group.by.
    } else if (group_by_is_used){
      # Check that the color palette has the right amount of named colors with regards to group.by values.
      colors.use <- check_consistency_colors_and_names(sample = sample,
                                                       colors = colors.use,
                                                       grouping_variable = group.by)
      # When using split.by.
    } else if (split_by_is_used){
      if (length(colors.use) != 1){
        # Check that the color palette has the right amount of named colors with regards to split.by values.
        colors.use <- check_consistency_colors_and_names(sample = sample,
                                                         colors = colors.use,
                                                         grouping_variable = split.by)
      }
      # When highlighting cells.
    } else if (highlighting_cells){
      # Stop the execution if more than one color is provided to highlight the cells.
      assertthat::assert_that(length(colors.use) == 1,
                              msg = "Provide only one color if cells.highlight or idents.highlight is used.")
    }
  }

  # Set cells to NA according to idents.keep.
  # If the user does not want to highlight cells or split by identities but wants to remove some identities.
  idents_keep_used <- is.null(cells.highlight) & is.null(idents.highlight) & !(is.null(idents.keep))
  if (idents_keep_used){
    # CONDITION: both group.by and split.by are not used.
    group_by_and_split_by_are_null <- is.null(group.by) & is.null(split.by)
    # CONDITION: group.by is used.
    group_by_is_used <- !(is.null(group.by)) & is.null(split.by)
    # CONDITION: split.by is used.
    split_by_is_used <- is.null(group.by) & !(is.null(split.by))
    # When running under default parameters.
    if (group_by_and_split_by_are_null){
      # Check that idents.keep matches the values and if not, stop the execution.
      assertthat::assert_that(isTRUE(length(idents.keep) == sum(idents.keep %in% levels(sample))),
                              msg = "All the values in idents.keep must be in levels(sample).")
      # Set the identities that the user wants to exclude as NA.
      Seurat::Idents(sample)[!(Seurat::Idents(sample) %in% idents.keep)] <- NA

      colors.use <- check_consistency_colors_and_names(sample = sample, colors = colors.use)
      # When using group.by, check with the values in group.by.
    } else if (group_by_is_used) {
      # Check that idents.keep matches the values, if not, stop the execution.
      assertthat::assert_that(isTRUE(length(idents.keep) == sum(idents.keep %in% unique(sample@meta.data[, group.by]))),
                              msg = "All the values in idents.keep must be in the group.by variable provided.")
      # Convert to NA values in group.by not included in the user's selected values.
      sample@meta.data[, group.by][!(sample@meta.data[, group.by] %in% idents.keep)] <- NA
      colors.use <- check_consistency_colors_and_names(sample = sample, colors = colors.use, grouping_variable = group.by)
      # If split.by is used instead.
    } else if (split_by_is_used){
      # Check that the values in idents.keep are in the unique values of split.by.
      assertthat::assert_that(isTRUE(length(idents.keep) == sum(idents.keep %in% unique(sample@meta.data[, split.by]))),
                              msg = "All the values in idents.keep must be in the split.by variable provided.")
      colors.use <- check_consistency_colors_and_names(sample = sample, colors = colors.use, grouping_variable = split.by)
    }
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
    p <- Seurat::DimPlot(sample,
                         reduction = reduction,
                         label = label,
                         dims = dims,
                         repel = ifelse(is.null(label) == TRUE, NULL, TRUE),
                         label.box = ifelse(is.null(label) == TRUE, NULL, TRUE),
                         label.color = ifelse(is.null(label) == TRUE, NULL, label.color),
                         na.value = na.value,
                         shuffle = shuffle,
                         order = order,
                         pt.size = pt.size,
                         group.by = group.by,
                         cols = colors.use,
                         raster = raster,
                         raster.dpi = c(raster.dpi, raster.dpi),
                         ncol = ncol) &
      ggplot2::guides(color = ggplot2::guide_legend(ncol = legend.ncol,
                                                    nrow = legend.nrow,
                                                    byrow = legend.byrow,
                                                    override.aes = list(size = legend.icon.size),
                                                    title.position = legend.title.position))
    if (isTRUE(label)){
      p <- add_scale(p = p,
                     function_use = ggplot2::scale_fill_manual(values = colors.use),
                     scale = "fill")
    }
    if (!(is.null(group.by))){
      # Remove automatic title inserted by Seurat.
      p <- p & ggplot2::ggtitle("")
    }

    # Add cell borders.
    if (isTRUE(plot_cell_borders)){
      p$layers <- append(base_layer, p$layers)
    }
  }
  # If split.by is used, the UMAP has to be split in multiple panes.
  else if (split_by_used){
    # If the user provided multiple highlighting colors.
    multiple_colors <- ifelse(length(colors.use) > 1, TRUE, FALSE)
    # List to store each individual plots.
    list.plots <- list()
    # Recover metadata values associated with split.by.
    data.use <- sample@meta.data[, split.by, drop = FALSE]
    # Retrieve the plotting order, keep factor levels if the column is a factor.
    plot_order <- if (is.factor(data.use[, 1])){levels(data.use[, 1])} else {sort(unique(data.use[, 1]))}
    # If idents.keep is used, subset to only these values.
    if (!(is.null(idents.keep))){
      plot_order <- if (is.factor(data.use[, 1])){levels(data.use[, 1])[levels(data.use[, 1]) %in% idents.keep]} else {sort(unique(data.use[, 1])[unique(data.use[, 1]) %in% idents.keep])}
      # If the user wants more than one color.
      if (isTRUE(multiple_colors)){
        colors.use <- colors.use[names(colors.use) %in% idents.keep]
      }
    }
    # Iterate over each unique value in split.by parameter.
    for (iteration in plot_order){
      # Retrieve the cells that do belong to the iteration's split.by value.
      cells.highlight <- rownames(data.use)[which(data.use == iteration)]
      p <- Seurat::DimPlot(sample,
                           reduction = reduction,
                           dims = dims,
                           cells.highlight = cells.highlight,
                           sizes.highlight = sizes.highlight,
                           pt.size = pt.size,
                           raster = raster,
                           raster.dpi = c(raster.dpi, raster.dpi),
                           ncol = ncol) &
        ggplot2::labs(title = iteration)
      p <- add_scale(p = p,
                     function_use = ggplot2::scale_color_manual(labels = c("Not selected", iteration),
                                                                values = c(na.value, ifelse(multiple_colors == TRUE, colors.use[[iteration]], colors.use)),
                                                                na.value = na.value),
                     scale = "color") &
        ggplot2::guides(color = ggplot2::guide_legend(title = legend.title,
                                                      ncol = legend.ncol,
                                                      nrwo = legend.nrow,
                                                      byrow = legend.byrow,
                                                      override.aes = list(size = legend.icon.size),
                                                      title.position = legend.title.position))
      # Add cell borders.
      if (isTRUE(plot_cell_borders)){
        p$layers <- append(base_layer, p$layers)
      }
      list.plots[[iteration]] <- p
    }
    # Assemble individual plots as a patch.
    p <- patchwork::wrap_plots(list.plots, ncol = ncol)
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
    p <- Seurat::DimPlot(sample,
                         reduction = reduction,
                         cells.highlight = cells.use,
                         sizes.highlight = sizes.highlight,
                         dims = dims,
                         pt.size = pt.size,
                         raster = raster,
                         raster.dpi = c(raster.dpi, raster.dpi),
                         ncol = ncol)
    p <- add_scale(p = p,
                   function_use = ggplot2::scale_color_manual(labels = c("Not selected", "Selected"),
                                                              values = c(na.value, colors.use),
                                                              na.value = na.value),
                   scale = "color") &
      ggplot2::guides(color = ggplot2::guide_legend(title = legend.title,
                                                    ncol = legend.ncol,
                                                    nrow = legend.nrow,
                                                    byrow = legend.byrow,
                                                    override.aes = list(size = legend.icon.size),
                                                    title.position = legend.title.position))

    # Add cell borders.
    if (isTRUE(plot_cell_borders)){
      p$layers <- append(base_layer, p$layers)
    }

  }
  # Titles in split.by are centered by default.
  hjust_use <- if(split_by_used){0.5} else {0}
  # Add theme settings to all plots.
  p <- p &
    ggplot2::theme_minimal(base_size = font.size) &
    ggplot2::theme(plot.title = ggplot2::element_text(face = "bold", hjust = hjust_use),
                   plot.subtitle = ggplot2::element_text(hjust = 0),
                   plot.caption = ggplot2::element_text(hjust = 1),
                   plot.title.position = "plot",
                   plot.caption.position = "plot",
                   text = ggplot2::element_text(family = font.type),
                   legend.justification = "center",
                   legend.text = ggplot2::element_text(face = "bold"),
                   legend.title = if (legend.position != "none"){ggplot2::element_text(face = "bold")} else {ggplot2::element_blank()},
                   legend.position = legend.position,
                   panel.grid = ggplot2::element_blank(),
                   plot.margin = ggplot2::margin(t = 10, r = 10, b = 10, l = 10),
                   plot.background = ggplot2::element_rect(fill = "white", color = "white"),
                   panel.background = ggplot2::element_rect(fill = "white", color = "white"),
                   legend.background = ggplot2::element_rect(fill = "white", color = "white"))

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

  # For embeddings that are umap of tsne, we remove all axes.
  if (reduction %in% c("umap", "tsne")){
    # If dims is first and then second (most of the cases).
    if (sum(dims == c(1, 2)) == 2){
      # Remove axes completely.
      p <- p &
        ggplot2::theme(axis.title = if (isFALSE(plot.axes)){ggplot2::element_blank()} else {ggplot2::element_text(color = "black", face = "bold", hjust = 0.5)},
                       axis.text = if (isFALSE(plot.axes)){ggplot2::element_blank()} else {ggplot2::element_text(color = "black")},
                       axis.ticks = if (isFALSE(plot.axes)){ggplot2::element_blank()} else {ggplot2::element_line(color = "black")},
                       axis.line = if (isFALSE(plot.axes)){ggplot2::element_blank()} else {ggplot2::element_line(color = "black")})
      # If dims do not follow the usual order.
    } else {
      # Get the name of the selected dims.
      labels <- colnames(sample@reductions[[reduction]][[]])[dims]
      # Remove everything in the axes but the axis titles.
      p <- p &
        ggplot2::theme(axis.text = if (isFALSE(plot.axes)){ggplot2::element_blank()} else {ggplot2::element_text(color = "black")},
                       axis.ticks = if (isFALSE(plot.axes)){ggplot2::element_blank()} else {ggplot2::element_line(color = "black")},
                       axis.line =if (isFALSE(plot.axes)){ggplot2::element_blank()} else {ggplot2::element_line(color = "black")},
                       axis.title = ggplot2::element_text(face = "bold", hjust = 0.5, color = "black")) &
        ggplot2::xlab(labels[1]) &
        ggplot2::ylab(labels[2])
    }
    # For diffusion maps, we do want to keep at least the axis titles so that we know which DC are we plotting.
  } else {
    # Get the name of the selected dims.
    labels <- colnames(sample@reductions[[reduction]][[]])[dims]
    # Remove everything in the axes but not the axis titles.
    p <- p &
      ggplot2::theme(axis.text = if (isFALSE(plot.axes)){ggplot2::element_blank()} else {ggplot2::element_text(color = "black")},
                     axis.ticks = if (isFALSE(plot.axes)){ggplot2::element_blank()} else {ggplot2::element_line(color = "black")},
                     axis.line =if (isFALSE(plot.axes)){ggplot2::element_blank()} else {ggplot2::element_line(color = "black")},
                     axis.title = ggplot2::element_text(face = "bold", hjust = 0.5, color = "black")) &
      ggplot2::xlab(labels[1]) &
      ggplot2::ylab(labels[2])
  }

  # Turn the labels to bold, when label is set to TRUE.
  if (label == TRUE && is.null(cells.highlight)){
    p$layers[[2]]$aes_params$fontface <- "bold"
  }

  # Add marginal plots.
  if (not_highlighting_and_not_split_by & isTRUE(plot_marginal_distributions & isFALSE(plot_cell_borders))){
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
    stop("Marginal distributions can not be used alongside when splitting by categories or highlighting cells or plotting cell borders .", call. = FALSE)
  }


  # Return the final plot.
  return(p)
}
