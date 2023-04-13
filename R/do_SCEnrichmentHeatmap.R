#' Perform a single-cell-based heatmap showing the enrichment in a list of gene sets.
#'
#' This function is heavily inspired by \strong{\code{\link[Seurat]{DoHeatmap}}}.
#'
#' @inheritParams doc_function
#' @param proportional.size \strong{\code{\link[base]{logical}}} | Whether the groups should take the same space in the plot or not.
#' @param main.heatmap.size \strong{\code{\link[base]{numeric}}} | Controls the size of the main heatmap (proportion-wise, defaults to 0.95).
#' @param metadata \strong{\code{\link[base]{character}}} | Categorical metadata variables to plot alongside the main heatmap.
#' @param metadata.colors \strong{\code{\link[SCpubr]{named_list}}} | Named list of valid colors for each of the variables defined in \strong{\code{metadata}}.
#' @param flavor \strong{\code{\link[base]{character}}} | One of: Seurat, UCell. Compute the enrichment scores using \link[Seurat]{AddModuleScore} or \link[UCell]{AddModuleScore_UCell}.
#' @param ncores \strong{\code{\link[base]{numeric}}} | Number of cores used to run UCell scoring.
#' @param storeRanks \strong{\code{\link[base]{logical}}} | Whether to store the ranks for faster UCell scoring computations. Might require large amounts of RAM.
#' @return A ggplot2 object.
#' @export
#'
#' @example /man/examples/examples_do_SCEnrichmentHeatmap.R
do_SCEnrichmentHeatmap <- function(sample,
                                   input_gene_list,
                                   assay = NULL,
                                   slot = NULL,
                                   group.by = NULL,
                                   metadata = NULL,
                                   metadata.colors = NULL,
                                   subsample = NA,
                                   cluster_cells = TRUE,
                                   flavor = "Seurat",
                                   return_object = FALSE,
                                   ncores = 1,
                                   storeRanks = TRUE,
                                   nbin = 24,
                                   ctrl = 100,
                                   xlab = "Cells",
                                   ylab = "Genes",
                                   font.size = 14,
                                   font.type = "sans",
                                   plot.title = NULL,
                                   plot.subtitle = NULL,
                                   plot.caption = NULL,
                                   legend.position = "bottom",
                                   legend.title = if (flavor != "AUCell") {"Enrichment"} else {"AUC"},
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
                                   sequential.palette = "YlGnBu",
                                   sequential.direction = 1,
                                   proportional.size = TRUE,
                                   verbose = FALSE,
                                   border.color = "black",
                                   plot.title.face = "bold",
                                   plot.subtitle.face = "plain",
                                   plot.caption.face = "italic",
                                   axis.title.face = "bold",
                                   axis.text.face = "bold",
                                   legend.title.face = "bold",
                                   legend.text.face = "plain"){
  # Get defaults user warning length.
  length.use <- getOption("warning.length")
  
  # Restore the warning length on exit.
  on.exit(options(warning.length = length.use))
  
  # Set warning length to maximum.
  options(warning.length = 8170)
  
  check_suggests(function_name = "do_SCEnrichmentHeatmap")
  check_Seurat(sample)
  
  # Check logical parameters.
  logical_list <- list("enforce_symmetry" = enforce_symmetry,
                       "proportional.size" = proportional.size,
                       "verbose" = verbose,
                       "legend.byrow" = legend.byrow,
                       "use_viridis" = use_viridis,
                       "cluster_cells" = cluster_cells,
                       "storeRanks" = storeRanks,
                       "return_object" = return_object) 
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
                       "nbin" = nbin,
                       "ctrl" = ctrl,
                       "ncores" = ncores)
  check_type(parameters = numeric_list, required_type = "numeric", test_function = is.numeric)
  # Check character parameters.
  character_list <- list("input_gene_list" = input_gene_list,
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
                         "flavor" = flavor,
                         "border.color" = border.color,
                         "plot.title.face" = plot.title.face,
                         "plot.subtitle.face" = plot.subtitle.face,
                         "plot.caption.face" = plot.caption.face,
                         "axis.title.face" = axis.title.face,
                         "axis.text.face" = axis.text.face,
                         "legend.title.face" = legend.title.face,
                         "legend.text.face" = legend.text.face) 
  
  
  check_colors(na.value, parameter_name = "na.value")
  check_colors(legend.framecolor, parameter_name = "legend.framecolor")
  check_colors(legend.tickcolor, parameter_name = "legend.tickcolor")
  check_colors(border.color, parameter_name = "border.color")
  
  check_parameters(parameter = font.type, parameter_name = "font.type")
  check_parameters(parameter = legend.type, parameter_name = "legend.type")
  check_parameters(parameter = legend.position, parameter_name = "legend.position")
  check_parameters(parameter = viridis.direction, parameter_name = "viridis.direction")
  check_parameters(parameter = viridis.palette, parameter_name = "viridis.palette")
  check_parameters(parameter = number.breaks, parameter_name = "number.breaks")
  check_parameters(parameter = diverging.palette, parameter_name = "diverging.palette")
  check_parameters(parameter = sequential.palette, parameter_name = "sequential.palette")
  check_parameters(parameter = sequential.direction, parameter_name = "sequential.direction")
  check_parameters(plot.title.face, parameter_name = "plot.title.face")
  check_parameters(plot.subtitle.face, parameter_name = "plot.subtitle.face")
  check_parameters(plot.caption.face, parameter_name = "plot.caption.face")
  check_parameters(axis.title.face, parameter_name = "axis.title.face")
  check_parameters(axis.text.face, parameter_name = "axis.text.face")
  check_parameters(legend.title.face, parameter_name = "legend.title.face")
  check_parameters(legend.text.face, parameter_name = "legend.text.face")
  
  `%>%` <- magrittr::`%>%`
  
  
  
  if (!(is.null(assay)) & flavor == "UCell"){
    warning(paste0(add_warning(), crayon_body("When using "),
                   crayon_key("flavor = UCell"),
                   crayon_body(" do not use the "),
                   crayon_key("assay"),
                   crayon_body(" parameter.\nInstead, make sure that the "),
                   crayon_key("assay"),
                   crayon_body(" you want to compute the scores with is set as the "),
                   crayon_key("default"),
                   crayon_body(" assay. Setting it to "),
                   crayon_key("NULL"),
                   crayon_body(".")), call. = FALSE)
  }
  
  if (!(is.null(slot)) & flavor == "Seurat"){
    warning(paste0(add_warning(), crayon_body("When using "),
                   crayon_key("flavor = Seurat"),
                   crayon_body(" do not use the "),
                   crayon_key("slot"),
                   crayon_body(" parameter.\nThis is determiend by default in "),
                   crayon_key("Seurat"),
                   crayon_body(". Setting it to "),
                   crayon_key("NULL"),
                   crayon_body(".")), call. = FALSE)
  }
  
  if (is.null(assay)){assay <- check_and_set_assay(sample)$assay}
  if (is.null(slot)){slot <- check_and_set_slot(slot)}
  
  if (is.character(input_gene_list)){
    stop(paste0(add_cross(),
                crayon_body("You have provided a string of genes to "),
                crayon_key("input_gene_list"),
                crayon_body(". Please provide a "),
                crayon_key("named list"),
                crayon_body(" instead.")), call. = FALSE)
    
  }
  
  input_list <- input_gene_list
  assertthat::assert_that(!is.null(names(input_list)),
                          msg = paste0(add_cross(), crayon_body("Please provide a "),
                                       crayon_key("named list"),
                                       crayon_body(" to "),
                                       crayon_key("input_gene_list"),
                                       crayon_body(".")))
  if (length(unlist(stringr::str_match_all(names(input_list), "_"))) > 0){
    warning(paste0(add_warning(), crayon_body("Found "),
                   crayon_key("underscores (_)"),
                   crayon_body(" in the name of the gene sets provided. Replacing them with "),
                   crayon_key("dots (.)"),
                   crayon_body(" to avoid conflicts when generating the Seurat assay.")), call. = FALSE)
    names.use <- stringr::str_replace_all(names(input_list), "_", ".")
    names(input_list) <- names.use
  }
  
  
  assertthat::assert_that(sum(names(input_list) %in% colnames(sample@meta.data)) == 0,
                          msg = paste0(add_cross(), crayon_body("Please make sure you do not provide a list of gene sets whose "),
                                       crayon_key("names"),
                                       crayon_body(" match any of the "),
                                       crayon_key("metadata columns"), 
                                       crayon_body(" of the Seurat object.")))
  # Compute the enrichment scores.
  sample <- compute_enrichment_scores(sample = sample,
                                      input_gene_list = input_list,
                                      verbose = verbose,
                                      nbin = nbin,
                                      ctrl = ctrl,
                                      flavor = flavor,
                                      ncores = ncores,
                                      storeRanks = storeRanks,
                                      # nocov start
                                      assay = if (flavor == "UCell"){NULL} else {assay},
                                      slot = if (flavor == "Seurat"){NULL} else {slot})
                                      # nocov end
  
  if (is.null(group.by)){
    sample$Groups <- Seurat::Idents(sample)
    group.by <- "Groups"
  }
  
  
  assertthat::assert_that(length(group.by) == 1,
                          msg = paste0(add_cross(), crayon_body("Please provide only a single value to "),
                                       crayon_key("group.by"),
                                       crayon_body(".")))
  
  
 
  
  # nocov start
  # Perform hierarchical clustering cluster-wise
  order.use <- if (is.factor(sample@meta.data[, group.by])){levels(sample@meta.data[, group.by])} else {sort(unique(sample@meta.data[, group.by]))}
  # nocov end
  
  matrix <-  sample@meta.data[, c(names(input_list), group.by)] %>% 
             tibble::rownames_to_column(var = "cell") %>%
             dplyr::group_by(.data[[group.by]])
  
  if (!is.na(subsample)){
    matrix <- matrix %>% 
              dplyr::slice_sample(n = subsample)
  }
  # Retrieve the order median-wise to cluster heatmap bodies.
  median.matrix <- matrix %>%
                   dplyr::summarise(dplyr::across(dplyr::all_of(names(input_list)), function(x){stats::median(x, na.rm = TRUE)})) %>%
                   dplyr::mutate("group.by" = as.character(.data[[group.by]])) %>%
                   dplyr::select(-dplyr::all_of(group.by)) %>%
                   as.data.frame() %>%
                   tibble::column_to_rownames(var = "group.by") %>%
                   as.matrix() %>%
                   t()
  group_order <- stats::hclust(stats::dist(t(median.matrix), method = "euclidean"), method = "ward.D")$order
  order.use <- order.use[group_order]
  
  # Retrieve the order median-wise for the genes.
  if (length(names(input_list)) == 1) {
    row_order <- names(input_list)[1]
  } else {
    row_order <- stats::hclust(stats::dist(median.matrix, method = "euclidean"), method = "ward.D")$order
    row_order <- names(input_list)[row_order]
  }
  
  
  # Compute cell order to group cells withing heatmap bodies.
  # nocov start
  if (isTRUE(cluster_cells)){
    if (sum(matrix %>% dplyr::pull(.data[[group.by]]) %>% table() > 65536)){
      warning(paste0(add_warning(), crayon_body("A given group in "),
                     crayon_key("group.by"),
                     crayon_body(" has more than "),
                     crayon_key("65536"),
                     crayon_body(" cells. Disabling clustering of the cells.")))
      cluster_cells <- FALSE
    }
  }
  # nocov end
  
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
      if (length(names(input_list)) == 1){
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
        colors.use <- generate_color_scale(names_use = names.use)
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
                             "gene" = factor(.data$gene, levels = rev(row_order)),
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
       ggplot2::geom_raster()
  
  
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
  
  p <- p + ggplot2::ylab(ylab) + 
           ggplot2::xlab(xlab)
  
  if (isFALSE(enforce_symmetry)){
    if (isTRUE(use_viridis)){
      p <- p +
           ggplot2::scale_fill_viridis_c(na.value = na.value,
                                         option = viridis.palette,
                                         direction = viridis.direction,
                                         name = legend.title,
                                         breaks = scale.setup$breaks,
                                         labels = scale.setup$labels,
                                         limits = scale.setup$limits)
    } else {
      p <- p +
           # nocov start
           ggplot2::scale_fill_gradientn(colors = if(sequential.direction == 1){RColorBrewer::brewer.pal(n = 9, name = sequential.palette)[2:9]} else {rev(RColorBrewer::brewer.pal(n = 9, name = sequential.palette)[2:9])},
                                         na.value = na.value,
                                         name = legend.title,
                                         breaks = scale.setup$breaks,
                                         labels = scale.setup$labels,
                                         limits = scale.setup$limits)
           # nocov end
    }
    
    
    
  } else if (isTRUE(enforce_symmetry)){
    p <- p + 
         ggplot2::scale_fill_gradientn(colors = RColorBrewer::brewer.pal(n = 11, name = diverging.palette) %>% rev(),
                                       na.value = na.value,
                                       name = "Regulon Score",
                                       breaks = scale.setup$breaks,
                                       labels = scale.setup$labels,
                                       limits = scale.setup$limits)
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
                                                                                                   face = plot.subtitle.face,
                                                                                                   color = "black",
                                                                                                   hjust = 0),
                                                             plot.caption = ggplot2::element_text(family = font.type,
                                                                                                  face = plot.caption.face,
                                                                                                  color = "black",
                                                                                                  hjust = 1),
                                                             plot.caption.position = "plot"))

  } else {
    out <- metadata_plots[["main"]]
  }
  out.list <- list()
  out.list[["Heatmap"]] <- out
  
  if (isTRUE(return_object)){
    sample[["Enrichment"]] <- sample@meta.data %>% 
                              dplyr::select(dplyr::all_of(names(input_list))) %>% 
                              t() %>% 
                              as.data.frame() %>% 
                              Seurat::CreateAssayObject(.)
    
    sample@meta.data <- sample@meta.data %>% 
                        dplyr::select(-dplyr::all_of(names(input_list)))
    
    sample@assays$Enrichment@key <- "Enrichment_"
    
    out.list[["Object"]] <- sample
    
    return(out.list)
  } else {
    return(out.list[["Heatmap"]])
  }
}
