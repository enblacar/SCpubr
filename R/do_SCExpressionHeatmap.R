#' Perform a single-cell-based heatmap showing the expression of genes.
#'
#' This function is heavily inspired by \strong{\code{\link[Seurat]{DoHeatmap}}}.
#'
#' @inheritParams doc_function
#' @param proportional.size \strong{\code{\link[base]{logical}}} | Whether the groups should take the same space in the plot or not.
#' @param main.heatmap.size \strong{\code{\link[base]{numeric}}} | Controls the size of the main heatmap (proportion-wise, defaults to 0.95).
#' @param metadata \strong{\code{\link[base]{character}}} | Categorical metadata variables to plot alongside the main heatmap.
#' @param metadata.colors \strong{\code{\link[SCpubr]{named_list}}} | Named list of valid colors for each of the variables defined in \strong{\code{metadata}}.
#' @return A ggplot2 object.
#' @export
#'
#' @example /man/examples/examples_do_SCExpressionHeatmap.R
do_SCExpressionHeatmap <- function(sample,
                                   features,
                                   assay = NULL,
                                   slot = NULL,
                                   group.by = NULL,
                                   features.order = NULL,
                                   metadata = NULL,
                                   metadata.colors = NULL,
                                   colorblind = FALSE,
                                   subsample = NA,
                                   cluster = TRUE,
                                   interpolate = FALSE,
                                   xlab = "Cells",
                                   ylab = "Genes",
                                   font.size = 14,
                                   font.type = "sans",
                                   plot.title = NULL,
                                   plot.subtitle = NULL,
                                   plot.caption = NULL,
                                   legend.position = "bottom",
                                   legend.title = "Expression",
                                   legend.type = "colorbar",
                                   legend.framewidth = 0.5,
                                   legend.tickwidth = 0.5,
                                   legend.length = 20,
                                   legend.width = 1,
                                   legend.framecolor = "grey50",
                                   legend.tickcolor = "white",
                                   strip.text.color = "black",
                                   strip.text.angle = 0,
                                   strip.spacing = 10,
                                   legend.ncol = NULL,
                                   legend.nrow = NULL,
                                   legend.byrow = FALSE,
                                   min.cutoff = NA,
                                   max.cutoff = NA,
                                   number.breaks = 5,
                                   main.heatmap.size = 0.95,
                                   enforce_symmetry = FALSE,
                                   use_viridis = FALSE,
                                   viridis.palette = "G",
                                   viridis.direction = -1,
                                   na.value = "grey75",
                                   diverging.palette = "RdBu",
                                   diverging.direction = -1,
                                   sequential.palette = "YlGnBu",
                                   sequential.direction = 1,
                                   proportional.size = TRUE,
                                   verbose = TRUE,
                                   border.color = "black",
                                   plot.title.face = "bold",
                                   plot.subtitle.face = "plain",
                                   plot.caption.face = "italic",
                                   axis.title.face = "bold",
                                   axis.text.face = "plain",
                                   legend.title.face = "bold",
                                   legend.text.face = "plain"){
  # Add lengthy error messages.
  withr::local_options(.new = list("warning.length" = 8170))

  check_suggests(function_name = "do_SCExpressionHeatmap")
  check_Seurat(sample)

  if (is.null(assay)){assay <- check_and_set_assay(sample)$assay}
  slot <- if(is.null(slot)){"data"} else {slot}

  # Check logical parameters.
  logical_list <- list("enforce_symmetry" = enforce_symmetry,
                       "proportional.size" = proportional.size,
                       "verbose" = verbose,
                       "legend.byrow" = legend.byrow,
                       "use_viridis" = use_viridis,
                       "cluster" = cluster,
                       "interpolate" = interpolate,
                       "colorblind" = colorblind)
  check_type(parameters = logical_list, required_type = "logical", test_function = is.logical)
  # Check numeric parameters.
  numeric_list <- list("font.size" = font.size,
                       "legend.framewidth" = legend.framewidth,
                       "legend.tickwidth" = legend.tickwidth,
                       "legend.length" = legend.length,
                       "legend.width" = legend.width,
                       "min.cutoff" = min.cutoff,
                       "max.cutoff" = max.cutoff,
                       "number.breaks" = number.breaks,
                       "viridis.direction" = viridis.direction,
                       "legend.ncol" = legend.ncol,
                       "legend.nrow" = legend.ncol,
                       "strip.spacing" = strip.spacing,
                       "strip.text.angle" = strip.text.angle,
                       "main.heatmap.size" = main.heatmap.size,
                       "sequential.direction" = sequential.direction,
                       "diverging.direction" = diverging.direction)
  check_type(parameters = numeric_list, required_type = "numeric", test_function = is.numeric)
  # Check character parameters.
  character_list <- list("features" = features,
                         "assay" = assay,
                         "slot" = slot,
                         "group.by" = group.by,
                         "xlab" = xlab,
                         "ylab" = ylab,
                         "font.type" = font.type,
                         "plot.title" = plot.title,
                         "plot.subtitle" = plot.subtitle,
                         "plot.caption" = plot.caption,
                         "legend.position" = legend.position,
                         "legend.title" = legend.title,
                         "legend.type" = legend.type,
                         "legend.framecolor" = legend.framecolor,
                         "legend.tickcolor" = legend.tickcolor,
                         "strip.text.color" = strip.text.color,
                         "viridis.palette" = viridis.palette,
                         "na.value" = na.value,
                         "metadata" = metadata,
                         "metadata.colors" = metadata.colors,
                         "diverging.palette" = diverging.palette,
                         "sequential.palette" = sequential.palette,
                         "border.color" = border.color,
                         "plot.title.face" = plot.title.face,
                         "plot.subtitle.face" = plot.subtitle.face,
                         "plot.caption.face" = plot.caption.face,
                         "axis.title.face" = axis.title.face,
                         "axis.text.face" = axis.text.face,
                         "legend.title.face" = legend.title.face,
                         "legend.text.face" = legend.text.face)
  check_type(parameters = character_list, required_type = "character", test_function = is.character)

  check_colors(na.value, parameter_name = "na.value")
  check_colors(legend.framecolor, parameter_name = "legend.framecolor")
  check_colors(legend.tickcolor, parameter_name = "legend.tickcolor")
  check_colors(border.color, parameter_name = "border.color")

  check_parameters(parameter = font.type, parameter_name = "font.type")
  check_parameters(parameter = legend.type, parameter_name = "legend.type")
  check_parameters(parameter = legend.position, parameter_name = "legend.position")
  check_parameters(parameter = viridis.palette, parameter_name = "viridis.palette")
  check_parameters(parameter = number.breaks, parameter_name = "number.breaks")
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

  `%>%` <- magrittr::`%>%`
  
  if (utils::packageVersion("Seurat") < "5.0.0"){
    genes.avail <- rownames(SeuratObject::GetAssayData(sample, slot = slot, assay = assay))
  } else {
    genes.avail <- rownames(SeuratObject::GetAssayData(sample, layer = slot, assay = assay))
  }
  
  assertthat::assert_that(sum(features %in% genes.avail) > 0,
                          msg = paste0(add_cross(), crayon_body("None of the features are present in the row names of the assay "),
                                       crayon_key(assay),
                                       crayon_body(" using the slot "),
                                       crayon_key(slot),
                                       crayon_body(".\nPlease make sure that you only provide "),
                                       crayon_key("genes"),
                                       crayon_body(" as input.\nIf you select the slot "),
                                       crayon_key("scale.data"),
                                       crayon_body(", sometimes some of the features are missing.")))


  missing_features <- features[!(features %in% genes.avail)]
  if (length(missing_features) > 0){
    if (isTRUE(verbose)){
      warning(paste0(add_warning(), crayon_body("Some features are missing in the following assay "),
                     crayon_key(assay),
                     crayon_body(" using the slot "),
                     crayon_key(slot),
                     crayon_body(":\n"),
                     paste(vapply(missing_features, crayon_key, FUN.VALUE = character(1)), collapse = crayon_body(", "))), call. = FALSE)
    }
  }

  features <- features[features %in% genes.avail]

  if (!is.null(features.order)){
    features.order <- features.order[features.order %in% genes.avail]
    assertthat::assert_that(sum(features.order %in% features) == length(features),
                            msg = paste0(add_cross(), crayon_body("The names provided to "),
                                         crayon_key("features.order"),
                                         crayon_body(" do not match the names of the gene sets in "),
                                         crayon_key("input_gene_list"),
                                         crayon_body(".")))
  }
  
  if (utils::packageVersion("Seurat") < "5.0.0"){
    matrix <- SeuratObject::GetAssayData(sample,
                                         assay = assay,
                                         slot = slot)[features, , drop = FALSE] %>%
              as.matrix()
  } else {
    matrix <- SeuratObject::GetAssayData(sample,
                                         assay = assay,
                                         layer = slot)[features, , drop = FALSE] %>%
              as.matrix()
  }
  
  # Check group.by.
  out <- check_group_by(sample = sample,
                        group.by = group.by,
                        is.heatmap = TRUE)
  sample <- out[["sample"]]
  group.by <- out[["group.by"]]

  assertthat::assert_that(length(group.by) == 1,
                          msg = paste0(add_cross(), crayon_body("Please provide only a single value to "),
                                       crayon_key("group.by"),
                                       crayon_body(".")))




  # nocov start
  # Perform hierarchical clustering cluster-wise
  order.use <- if (is.factor(sample@meta.data[, group.by])){levels(sample@meta.data[, group.by])} else {sort(unique(sample@meta.data[, group.by]))}
  # nocov end

  matrix <-  matrix %>%
             t() %>%
             as.data.frame() %>%
             tibble::rownames_to_column(var = "cell") %>%
             dplyr::left_join(y = {sample@meta.data %>%
                                   tibble::rownames_to_column(var = "cell") %>%
                                   dplyr::select(dplyr::all_of(c("cell", group.by)))},
                                   by = "cell") %>%
             dplyr::group_by(.data[[group.by]])
  if (!is.na(subsample)){
   matrix <- matrix %>%
             dplyr::slice_sample(n = subsample)
  }
  # Retrieve the order median-wise to cluster heatmap bodies.
  if (isTRUE(cluster)){
    median.matrix <- matrix %>%
                     dplyr::summarise(dplyr::across(dplyr::all_of(features), function(x){stats::median(x, na.rm = TRUE)})) %>%
                     dplyr::mutate("group.by" = as.character(.data[[group.by]])) %>%
                     dplyr::select(-dplyr::all_of(group.by)) %>%
                     as.data.frame() %>%
                     tibble::column_to_rownames(var = "group.by") %>%
                     as.matrix() %>%
                     t()
    group_order <- stats::hclust(stats::dist(t(median.matrix), method = "euclidean"), method = "ward.D")$order
    order.use <- order.use[group_order]
  }


  # Retrieve the order median-wise for the genes.
  if (length(features) == 1) {
    row_order <- features[1]
  } else {
    if (isTRUE(cluster)){
      row_order <- features[stats::hclust(stats::dist(median.matrix, method = "euclidean"), method = "ward.D")$order]
    } else {
      row_order <- features
    }
  }


  # Compute cell order to group cells withing heatmap bodies.
  # nocov start
  if (isTRUE(cluster)){
    if (sum(matrix %>% dplyr::pull(dplyr::all_of(c(group.by))) %>% table() > 65536)){
      warning(paste0(add_warning(), crayon_body("A given group in "),
                     crayon_key("group.by"),
                     crayon_body(" has more than "),
                     crayon_key("65536"),
                     crayon_body(" cells. Disabling clustering of the cells.")), call. = FALSE)
      cluster <- FALSE
    }
  }
  # nocov end

  if (isTRUE(cluster)){
    col_order <- list()
    for (item in order.use){
      cells.use <- matrix %>%
        dplyr::filter(.data[[group.by]] == item) %>%
        dplyr::pull(dplyr::all_of("cell"))

      matrix.subset <- matrix %>%
                       dplyr::ungroup() %>%
                       dplyr::select(-dplyr::all_of(c(group.by))) %>%
                       tibble::column_to_rownames(var = "cell") %>%
                       as.data.frame() %>%
                       as.matrix() %>%
                       t()
      matrix.subset <- matrix.subset[, cells.use]
      # nocov start
      if (sum(is.na(matrix.subset)) > 0){
        warning(paste0(add_warning(), crayon_key("NA"), crayon_body("found in the "),
                       crayon_key("expression matrix"),
                       crayon_body(". Replacing them with "),
                       crayon_key("0"),
                       crayon_body(".")), call. = FALSE)
        matrix.subset[is.na(matrix.subset)] <- 0
      }
      # nocov end
      if (length(features) == 1){
        matrix.use <- as.matrix(matrix.subset)
      } else {
        matrix.use <- t(matrix.subset)
      }
      col_order.use <- stats::hclust(stats::dist(matrix.use, method = "euclidean"), method = "ward.D")$order

      col_order[[item]] <- cells.use[col_order.use]
    }
    col_order <- unlist(unname(col_order))
  } else {
    col_order <- matrix %>% dplyr::pull("cell")
  }





  # Retrieve metadata matrix.
  metadata_plots <- list()
  if (!is.null(metadata)){
    metadata.matrix <- sample@meta.data %>%
                       dplyr::select(dplyr::all_of(c(metadata, group.by))) %>%
                       dplyr::mutate("group.by" = .data[[group.by]]) %>%
                       as.matrix() %>%
                       t()
    metadata.matrix <- metadata.matrix[, col_order]

    counter <- 0
    for (name in metadata){
      counter <- counter + 1
      if (counter == 1){
        name_labels <- name
      }
      plot_data <- metadata.matrix[c(name, "group.by"), ] %>%
                   t() %>%
                   as.data.frame() %>%
                   tibble::rownames_to_column(var = "cell") %>%
                   dplyr::mutate("group.by" = factor(.data$group.by, levels = order.use),
                                 "y" = .data[[name]],
                                 "y_row" = name,
                                 "cell" = factor(.data$cell, levels = col_order)) %>%
                   dplyr::select(-dplyr::all_of(name)) %>%
                   tibble::as_tibble()

      if (name %in% names(metadata.colors)){
        colors.use <- metadata.colors[[name]]
      } else {
        names.use <- if(is.factor(sample@meta.data[, name])){levels(sample@meta.data[, name])} else {sort(unique(sample@meta.data[, name]))}
        colors.use <- generate_color_scale(names_use = names.use, colorblind = colorblind)
      }
      p <- plot_data %>%
           ggplot2::ggplot(mapping = ggplot2::aes(x = .data$cell,
                                                  y = .data$y_row,
                                                  fill = .data$y)) +
           ggplot2::geom_tile() +
           ggplot2::facet_grid(~ .data$group.by,
                               scales = "free_x",
                               space = if(isTRUE(proportional.size)) {"fixed"} else {"free"}) +
           ggplot2::scale_fill_manual(values = colors.use) +
           ggplot2::guides(fill = ggplot2::guide_legend(title = name,
                                                       title.position = "top",
                                                       title.hjust = 0.5,
                                                       override.aes = list(color = "black",
                                                                           shape = 22),
                                                       ncol = legend.ncol,
                                                       nrow = legend.nrow,
                                                       byrow = legend.byrow)) +
           ggplot2::xlab(NULL) +
           ggplot2::ylab(NULL)

      metadata_plots[[name]] <- p
    }
  }

  # Generate the plotting data.
  plot_data <- matrix %>%
               dplyr::ungroup() %>%
               as.data.frame() %>%
               tidyr::pivot_longer(cols = -dplyr::all_of(c(group.by, "cell")),
                                   names_to = "gene",
                                   values_to = "expression") %>%
               dplyr::rename("group.by" = dplyr::all_of(c(group.by))) %>%
               dplyr::mutate("group.by" = factor(.data$group.by, levels = order.use),
                             "gene" = factor(.data$gene, levels = if (is.null(features.order)){rev(row_order)} else {features.order}),
                             "cell" = factor(.data$cell, levels = col_order))


  # Modify data to fit the cutoffs selected.
  plot_data_limits <- plot_data
  if (!is.na(min.cutoff)){
    plot_data <- plot_data %>%
                 dplyr::mutate("expression" = ifelse(.data$expression < min.cutoff, min.cutoff, .data$expression))
  }

  if (!is.na(max.cutoff)){
    plot_data <- plot_data %>%
                 dplyr::mutate("expression" = ifelse(.data$expression > max.cutoff, max.cutoff, .data$expression))
  }

  p <- plot_data %>%
       ggplot2::ggplot(mapping = ggplot2::aes(x = .data$cell,
                                              y = .data$gene,
                                              fill = .data$expression)) +
       ggplot2::geom_raster(interpolate = interpolate)


  p <- p + ggplot2::facet_grid(~ .data$group.by,
                               scales = "free_x",
                               space = if(isTRUE(proportional.size)) {"fixed"} else {"free"})

  limits.use <- c(min(plot_data_limits$expression, na.rm = TRUE),
                  max(plot_data_limits$expression, na.rm = TRUE))

  scale.setup <- compute_scales(sample = sample,
                                feature = NULL,
                                assay = assay,
                                reduction = NULL,
                                slot = slot,
                                number.breaks = number.breaks,
                                min.cutoff = min.cutoff,
                                max.cutoff = max.cutoff,
                                flavor = "Seurat",
                                enforce_symmetry = enforce_symmetry,
                                from_data = TRUE,
                                limits.use = limits.use)

  p <- p +
       ggplot2::ylab(ylab) +
       ggplot2::xlab(xlab) +
       ggplot2::scale_fill_gradientn(colors = colors.gradient,
                                     na.value = na.value,
                                     name = legend.title,
                                     breaks = scale.setup$breaks,
                                     labels = scale.setup$labels,
                                     limits = scale.setup$limits)

  p <- modify_continuous_legend(p = p,
                                legend.title = legend.title,
                                legend.aes = "fill",
                                legend.type = legend.type,
                                legend.position = legend.position,
                                legend.length = legend.length,
                                legend.width = legend.width,
                                legend.framecolor = legend.framecolor,
                                legend.tickcolor = legend.tickcolor,
                                legend.framewidth = legend.framewidth,
                                legend.tickwidth = legend.tickwidth)



  # Theme setup.
  metadata_plots[["main"]] <- p


  # Configure plot margins.

  for (name in names(metadata_plots)){

    metadata_plots[[name]] <- metadata_plots[[name]] +
                              ggplot2::scale_x_discrete(expand = c(0, 0)) +
                              ggplot2::scale_y_discrete(expand = c(0, 0)) +
                              ggplot2::theme_minimal(base_size = font.size) +
                              ggplot2::theme(axis.text.x = ggplot2::element_blank(),
                                             axis.text.y = ggplot2::element_text(face = axis.text.face,
                                                                                 color = "black"),
                                             axis.ticks.y = ggplot2::element_line(color = "black"),
                                             axis.ticks.x = ggplot2::element_blank(),
                                             axis.line = ggplot2::element_blank(),
                                             axis.title = ggplot2::element_text(face = axis.title.face, color = "black"),
                                             plot.title = ggplot2::element_text(face = plot.title.face, hjust = 0),
                                             plot.subtitle = ggplot2::element_text(face = plot.subtitle.face, hjust = 0),
                                             plot.caption = ggplot2::element_text(face = plot.caption.face, hjust = 1),
                                             legend.text = ggplot2::element_text(face = legend.text.face),
                                             legend.title = ggplot2::element_text(face = legend.title.face),
                                             plot.title.position = "plot",
                                             panel.grid = ggplot2::element_blank(),
                                             panel.grid.minor.y = ggplot2::element_line(color = "white"),
                                             strip.background = ggplot2::element_blank(),
                                             strip.clip = "off",
                                             panel.spacing = ggplot2::unit(strip.spacing, units = "pt"),
                                             text = ggplot2::element_text(family = font.type),
                                             plot.caption.position = "plot",
                                             legend.position = legend.position,
                                             legend.justification = "center",
                                             plot.margin = ggplot2::margin(t = 0, r = 10, b = 0, l = 10),
                                             panel.border = ggplot2::element_rect(color = border.color, fill = NA),
                                             panel.grid.major = ggplot2::element_blank(),
                                             plot.background = ggplot2::element_rect(fill = "white", color = "white"),
                                             panel.background = ggplot2::element_rect(fill = "white", color = "white"),
                                             legend.background = ggplot2::element_rect(fill = "white", color = "white"))

    if (!is.null(metadata)){
      if (name == name_labels){
        metadata_plots[[name]] <- metadata_plots[[name]] + ggplot2::theme(strip.text.x = ggplot2::element_text(family = font.type,
                                                                                                               face = "bold",
                                                                                                               color = strip.text.color,
                                                                                                               angle = strip.text.angle))
      } else {
        metadata_plots[[name]] <- metadata_plots[[name]] + ggplot2::theme(strip.text.x = ggplot2::element_blank())
      }
    } else {
      metadata_plots[[name]] <- metadata_plots[[name]] + ggplot2::theme(strip.text.x = ggplot2::element_text(family = font.type,
                                                                                                             face = "bold",
                                                                                                             color = strip.text.color,
                                                                                                             angle = strip.text.angle))
    }
  }

  if (!is.null(metadata)){
    plots_wrap <- c(metadata_plots[c(metadata, "main")])
    main_body_size <- main.heatmap.size
    height_unit <- c(rep((1 - main_body_size) / length(metadata), length(metadata)), main_body_size)

    out <- patchwork::wrap_plots(plots_wrap,
                                 ncol = 1,
                                 guides = "collect",
                                 heights = height_unit) +
           patchwork::plot_annotation(title = plot.title,
                                      subtitle = plot.subtitle,
                                      caption = plot.caption,
                                      theme = ggplot2::theme(legend.position = legend.position,
                                                             plot.title = ggplot2::element_text(family = font.type,
                                                                                                color = "black",
                                                                                                face = plot.title.face,
                                                                                                hjust = 0),
                                                             plot.subtitle = ggplot2::element_text(family = font.type,
                                                                                                   color = "black",
                                                                                                   face = plot.subtitle.face,
                                                                                                   hjust = 0),
                                                             plot.caption = ggplot2::element_text(family = font.type,
                                                                                                  color = "black",
                                                                                                  face = plot.caption.face,
                                                                                                  hjust = 1),
                                                             plot.caption.position = "plot"))

  } else {
    out <- metadata_plots[["main"]]
  }


  return(out)
}
