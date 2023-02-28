#' Perform a single-cell-based heatmap showing the expression of genes.
#'
#' This function is heavily inspired by \strong{\code{\link[Seurat]{DoHeatmap}}}.
#'
#' @inheritParams doc_function
#' @param make_size_proportional \strong{\code{\link[base]{logical}}} | Whether the groups should take the same space in the plot or not.
#' @param main.heatmap.size \strong{\code{\link[base]{numeric}}} | Controls the size of the main heatmap (proportion-wise, defaults to 0.95).
#' @param metadata \strong{\code{\link[base]{character}}} | Categorical metadata variables to plot alongside the main heatmap.
#' @param metadata.colors \strong{\code{\link[base]{list}}} | Named list of valid colors for each of the variables defined in \strong{\code{metadata}}.
#' @param metadata.location \strong{\code{\link[base]{character}}} | Location of the metadata rows. Either top or bottom.
#' @return A ggplot2 object.
#' @export

#'
#' @examples
#' \donttest{
#' TBD
#' }
do_SCExpressionHeatmap <- function(sample,
                                   features,
                                   assay = "SCT",
                                   slot = "data",
                                   group.by = NULL,
                                   metadata = NULL,
                                   metadata.colors = NULL,
                                   metadata.location = "top",
                                   subsample = NA,
                                   cluster_cells = TRUE,
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
                                   rotate_strip_labels = 0,
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
                                   viridis_color_map = "G",
                                   viridis_direction = -1,
                                   na.value = "grey75",
                                   diverging.palette = "RdBu",
                                   sequential.palette = "YlGnBu",
                                   sequential_direction = 1,
                                   make_size_proportional = TRUE,
                                   verbose = TRUE){

  check_suggests(function_name = "do_SCExpressionHeatmap")
  check_Seurat(sample)

  # Check the assay.
  out <- check_and_set_assay(sample = sample, assay = assay)
  sample <- out[["sample"]]
  assay <- out[["assay"]]

  # Check slot.
  slot <- check_and_set_slot(slot = slot)

  # Check logical parameters.
  logical_list <- list("enforce_symmetry" = enforce_symmetry,
                       "make_size_proportional" = make_size_proportional,
                       "verbose" = verbose,
                       "legend.byrow" = legend.byrow,
                       "use_viridis" = use_viridis,
                       "cluster_cells" = cluster_cells)
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
                       "viridis_direction" = viridis_direction,
                       "legend.ncol" = legend.ncol,
                       "legend.nrow" = legend.ncol,
                       "strip.spacing" = strip.spacing,
                       "rotate_strip_labels" = rotate_strip_labels,
                       "main.heatmap.size" = main.heatmap.size,
                       "sequential_direction" = sequential_direction)
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
                         "viridis_color_map" = viridis_color_map,
                         "na.value" = na.value,
                         "metadata" = metadata,
                         "metadata.colors" = metadata.colors,
                         "metadata.location" = metadata.location,
                         "diverging.palette" = diverging.palette,
                         "sequential.palette" = sequential.palette)


  check_colors(na.value, parameter_name = "na.value")
  check_colors(legend.framecolor, parameter_name = "legend.framecolor")
  check_colors(legend.tickcolor, parameter_name = "legend.tickcolor")

  check_parameters(parameter = font.type, parameter_name = "font.type")
  check_parameters(parameter = legend.type, parameter_name = "legend.type")
  check_parameters(parameter = legend.position, parameter_name = "legend.position")
  check_parameters(parameter = viridis_direction, parameter_name = "viridis_direction")
  check_parameters(parameter = viridis_color_map, parameter_name = "viridis_color_map")
  check_parameters(parameter = number.breaks, parameter_name = "number.breaks")
  check_parameters(parameter = diverging.palette, parameter_name = "diverging.palette")
  check_parameters(parameter = sequential.palette, parameter_name = "sequential.palette")
  check_parameters(parameter = sequential_direction, parameter_name = "sequential_direction")


  `%>%` <- magrittr::`%>%`
  genes.avail <- rownames(Seurat::GetAssayData(sample, slot = slot, assay = assay))

  assertthat::assert_that(sum(features %in% genes.avail) > 0,
                          msg = paste0(crayon_body("None of the features are present in the row names of the assay "),
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
      warning(paste0(crayon_body("Some features are missing in the following assay "),
                     crayon_key(assay),
                     crayon_body(" using the slot "),
                     crayon_key(slot),
                     crayon_body(":\n"),
                     paste(sapply(missing_features, crayon_key), collapse = crayon_body(", "))), call. = FALSE)
    }
  }

  features <- features[features %in% genes.avail]


  matrix <- Seurat::GetAssayData(sample,
                                 assay = assay,
                                 slot = slot)[features, , drop = FALSE] %>%
            as.matrix()
  assertthat::assert_that(length(group.by) == 1,
                          msg = paste0(crayon_body("Please provide only a single value to "),
                                       crayon_key("group.by"),
                                       crayon_body(".")))
  
  
  if (is.null(group.by)){
    sample$Groups <- Seurat::Idents(sample)
    group.by <- "Groups"
  }


  # Perform hierarchical clustering cluster-wise
  order.use <- if (is.factor(sample@meta.data[, group.by])){levels(sample@meta.data[, group.by])} else {sort(unique(sample@meta.data[, group.by]))}
  
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
  median.matrix <- matrix %>%
                   dplyr::summarise(dplyr::across(dplyr::all_of(features), stats::median)) %>%
                   dplyr::mutate("group.by" = as.character(.data[[group.by]])) %>%
                   dplyr::select(-dplyr::all_of(group.by)) %>%
                   as.data.frame() %>%
                   tibble::column_to_rownames(var = "group.by") %>%
                   as.matrix() %>%
                   t()
  group_order <- stats::hclust(stats::dist(t(median.matrix), method = "euclidean"), method = "ward.D")$order
  order.use <- order.use[group_order]

  # Retrieve the order median-wise for the genes.
  if (length(features) == 1) {
    row_order <- c(1)
  } else {
    row_order <- stats::hclust(stats::dist(median.matrix, method = "euclidean"), method = "ward.D")$order
    row_order <- features[row_order]
  }


  # Compute cell order to group cells withing heatmap bodies.
  if (isTRUE(cluster_cells)){
    if (sum(matrix %>% dplyr::pull(.data[[group.by]]) %>% table() > 65536)){
      warning(paste0(crayon_body("A given group in "),
                     crayon_key("group.by"),
                     crayon_body(" has more than "),
                     crayon_key("65536"),
                     crayon_body(" cells. Disabling clustering of the cells.")))
      cluster_cells <- FALSE
    }
  }
  
  if (isTRUE(cluster_cells)){
    col_order <- list()
    for (item in order.use){
      cells.use <- matrix %>%
        dplyr::filter(.data[[group.by]] == item) %>%
        dplyr::pull(.data$cell)
      
      matrix.subset <- matrix %>% 
                       dplyr::ungroup() %>% 
                       dplyr::select(-dplyr::all_of(c(group.by))) %>% 
                       tibble::column_to_rownames(var = "cell") %>% 
                       as.data.frame() %>% 
                       as.matrix() %>% 
                       t()
      matrix.subset <- matrix.subset[, cells.use]
      col_order.use <- stats::hclust(stats::dist(t(matrix.subset), method = "euclidean"), method = "ward.D")$order
      
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
                       dplyr::select(-dplyr::all_of(group.by)) %>%
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
        colors.use <- SCpubr:::generate_color_scale(names_use = names.use)
      }
      p <- plot_data %>%
           ggplot2::ggplot(mapping = ggplot2::aes(x = .data$cell,
                                                  y = .data$y_row,
                                                  fill = .data$y)) +
           ggplot2::geom_tile() +
           ggplot2::facet_grid(~ .data$group.by,
                               scales = "free_x",
                               space = if(isTRUE(make_size_proportional)) {"fixed"} else {"free"}) +
           ggplot2::scale_fill_manual(values = colors.use) + 
           ggplot2::guides(fill = ggplot2::guide_legend(title = name,
                                                       title.position = "top",
                                                       title.hjust = 0.5,
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
               dplyr::rename("group.by" = .data[[group.by]]) %>% 
               dplyr::mutate("group.by" = factor(.data$group.by, levels = order.use),
                             "gene" = factor(.data$gene, levels = rev(row_order)),
                             "cell" = factor(.data$cell, levels = col_order))


  # Modify data to fit the cutoffs selected.
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
       ggplot2::geom_raster()


  p <- p + ggplot2::facet_grid(~ .data$group.by,
                               scales = "free_x",
                               space = if(isTRUE(make_size_proportional)) {"fixed"} else {"free"})

  limits.use <- c(min(plot_data$expression, na.rm = TRUE),
                  max(plot_data$expression, na.rm = TRUE))

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
  p <- p + ggplot2::ylab(ylab)
  if (!is.null(metadata)){
    if (metadata.location == "top"){
      p <- p + ggplot2::xlab(xlab)
    } else {
      p <- p + ggplot2::xlab(NULL)
    }
  } else {
    p <- p + ggplot2::xlab(xlab)
  }

  if (isFALSE(enforce_symmetry)){
    if (isTRUE(use_viridis)){
      p <- p +
           ggplot2::scale_fill_viridis_c(na.value = na.value,
                                         option = viridis_color_map,
                                         direction = viridis_direction,
                                         name = legend.title,
                                         breaks = scale.setup$breaks,
                                         labels = scale.setup$labels,
                                         limits = scale.setup$limits)
    } else {
      p <- p +
           ggplot2::scale_fill_gradientn(colors = if(sequential_direction == 1){RColorBrewer::brewer.pal(n = 9, name = sequential.palette)[2:9]} else {rev(RColorBrewer::brewer.pal(n = 9, name = sequential.palette)[2:9])},
                                         na.value = na.value,
                                         name = legend.title,
                                         breaks = scale.setup$breaks,
                                         labels = scale.setup$labels,
                                         limits = scale.setup$limits)
    }



  } else if (isTRUE(enforce_symmetry)){
    p <- add_scale(p = p,
                   function_use = ggplot2::scale_fill_gradientn(colors = RColorBrewer::brewer.pal(n = 11, name = diverging.palette) %>% rev(),
                                                                na.value = na.value,
                                                                name = "Regulon Score",
                                                                breaks = scale.setup$breaks,
                                                                labels = scale.setup$labels,
                                                                limits = scale.setup$limits),
                   scale = "fill")
  }


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
                              ggplot2::labs(title = plot.title,
                                            subtitle = plot.subtitle,
                                            caption = plot.caption) +
                              ggplot2::theme_minimal(base_size = font.size) +
                              ggplot2::theme(axis.text.x = ggplot2::element_blank(),
                                             axis.text.y = ggplot2::element_text(face = "bold",
                                                                                 color = "black"),
                                             axis.ticks.y = ggplot2::element_line(color = "black"),
                                             axis.ticks.x = ggplot2::element_blank(),
                                             axis.line = ggplot2::element_blank(),
                                             axis.title = ggplot2::element_text(face = "bold", color = "black"),
                                             plot.title = ggplot2::element_text(face = "bold", hjust = 0),
                                             plot.subtitle = ggplot2::element_text(hjust = 0),
                                             plot.caption = ggplot2::element_text(hjust = 1),
                                             plot.title.position = "plot",
                                             panel.grid = ggplot2::element_blank(),
                                             panel.grid.minor.y = ggplot2::element_line(color = "white"),
                                             strip.background = ggplot2::element_blank(),
                                             strip.clip = "off",
                                             panel.spacing = ggplot2::unit(strip.spacing, units = "pt"),
                                             text = ggplot2::element_text(family = font.type),
                                             plot.caption.position = "plot",
                                             legend.text = ggplot2::element_text(face = "bold"),
                                             legend.position = legend.position,
                                             legend.title = ggplot2::element_text(face = "bold"),
                                             legend.justification = "center",
                                             plot.margin = ggplot2::margin(t = 0, r = 10, b = 0, l = 10),
                                             panel.border = ggplot2::element_rect(color = "black", fill = NA),
                                             panel.grid.major = ggplot2::element_blank(),
                                             plot.background = ggplot2::element_rect(fill = "white", color = "white"),
                                             panel.background = ggplot2::element_rect(fill = "white", color = "white"),
                                             legend.background = ggplot2::element_rect(fill = "white", color = "white"))

    if (!is.null(metadata)){
      if (name == name_labels){
        metadata_plots[[name]] <- metadata_plots[[name]] + ggplot2::theme(strip.text.x = ggplot2::element_text(family = font.type,
                                                                                                               face = "bold",
                                                                                                               color = strip.text.color,
                                                                                                               angle = rotate_strip_labels))
      } else {
        metadata_plots[[name]] <- metadata_plots[[name]] + ggplot2::theme(strip.text.x = ggplot2::element_blank())
      }
    } else {
      metadata_plots[[name]] <- metadata_plots[[name]] + ggplot2::theme(strip.text.x = ggplot2::element_text(family = font.type,
                                                                                                             face = "bold",
                                                                                                             color = strip.text.color,
                                                                                                             angle = rotate_strip_labels))
    }
  }

  if (!is.null(metadata)){
    plots_wrap <- {
      if (metadata.location == "bottom") {
        c(metadata_plots[c("main", metadata)])
      } else {
        c(metadata_plots[c(metadata, "main")])
      }
    }
    main_body_size <- main.heatmap.size
    height_unit <- {
      if(metadata.location == "bottom"){
        c(main_body_size, rep((1 - main_body_size) / length(metadata), length(metadata)))
      } else {
        c(rep((1 - main_body_size) / length(metadata), length(metadata)), main_body_size)
      }
    }
    out <- patchwork::wrap_plots(plots_wrap,
                                 ncol = 1,
                                 guides = "collect",
                                 heights = height_unit) +
      patchwork::plot_annotation(title = plot.title,
                                 subtitle = plot.subtitle,
                                 caption = plot.caption,
                                 theme = ggplot2::theme(legend.position = legend.position,
                                                        plot.title = ggplot2::element_text(size = font.size,
                                                                                           family = font.type,
                                                                                           color = "black",
                                                                                           face = "bold",
                                                                                           hjust = 0),
                                                        plot.subtitle = ggplot2::element_text(size = font.size,
                                                                                              family = font.type,
                                                                                              color = "black",
                                                                                              hjust = 0),
                                                        plot.caption = ggplot2::element_text(size = font.size,
                                                                                             family = font.type,
                                                                                             color = "black",
                                                                                             hjust = 1),
                                                        plot.caption.position = "plot"))
  } else {
    out <- metadata_plots[["main"]]
  }


  return(out)
}
