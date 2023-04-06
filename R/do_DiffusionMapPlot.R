#' Compute a heatmap of enrichment of gene sets on the context of a diffusion component.
#'
#' @inheritParams doc_function
#' @param colors.use \strong{\code{\link[base]{list}}} | A named list of named vectors. The names of the list correspond to the names of the values provided to metadata and the names of the items in the named vectors correspond to the unique values of that specific metadata variable. The values are the desired colors in HEX code for the values to plot. The used are pre-defined by the pacakge but, in order to get the most out of the plot, please provide your custom set of colors for each metadata column! 
#'
#' @return A list of ggplot2 objects and a Seurat object if desired.
#' @export
#'
#' @examples
#' \donttest{
#' TBD
#' }
do_DiffusionMapPlot <- function(sample,
                                input_gene_list,
                                assay = "SCT",
                                slot = "data",
                                scale.enrichment = TRUE,
                                dims = 1:5,
                                subsample = 2500,
                                reduction = "diffusion",
                                group.by = NULL,
                                colors.use = NULL,
                                nbin = 24,
                                ctrl = 100,
                                flavor = "Seurat",
                                main.heatmap.size = 0.95,
                                enforce_symmetry = ifelse(isTRUE(scale.enrichment), TRUE, FALSE),
                                use_viridis = FALSE,
                                viridis.palette = "G",
                                viridis.direction = -1,
                                sequential.palette = "YlGnBu",
                                sequential.direction = 1,
                                font.size = 14,
                                font.type = "sans",
                                na.value = "grey75",
                                legend.width = 1,
                                legend.length = 20,
                                legend.framewidth = 0.5,
                                legend.tickwidth = 0.5,
                                legend.framecolor = "grey50",
                                legend.tickcolor = "white",
                                legend.type = "colorbar",
                                legend.position = "bottom",
                                legend.nrow = NULL,
                                legend.ncol = NULL,
                                legend.byrow = FALSE,
                                number.breaks = 5,
                                diverging.palette = "RdBu",
                                axis.text.x.angle = 45,
                                border.color = "black",
                                return_object = FALSE,
                                verbose = TRUE,
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
  
  check_suggests("do_DiffusionMapPlot")
  check_Seurat(sample = sample)
  
  # Check logical parameters.
  logical_list <- list("enforce_symmetry" = enforce_symmetry,
                       "legend.byrow" = legend.byrow,
                       "return_object" = return_object,
                       "scale.enrichment" = scale.enrichment,
                       "use_viridis" = use_viridis,
                       "verbose" = verbose)
  check_type(parameters = logical_list, required_type = "logical", test_function = is.logical)
  
  # Check numeric parameters.
  numeric_list <- list("dims" = dims,
                       "subsample" = subsample,
                       "nbin" = nbin,
                       "ctrl" = ctrl,
                       "font.size" = font.size,
                       "legend.width" = legend.width,
                       "legend.length" = legend.length,
                       "legend.framewidth" = legend.framewidth,
                       "legend.tickwidth" = legend.tickwidth,
                       "number.breaks" = number.breaks,
                       "axis.text.x.angle" = axis.text.x.angle,
                       "legend.nrow" = legend.nrow,
                       "legend.ncol" = legend.ncol,
                       "main.heatmap.size" = main.heatmap.size,
                       "viridis.direction" = viridis.direction,
                       "sequential.direction" = sequential.direction)
  check_type(parameters = numeric_list, required_type = "numeric", test_function = is.numeric)
  
  # Check character parameters.
  character_list <- list("assay" = assay,
                         "reduction" = reduction,
                         "slot" = slot,
                         "group.by" = group.by,
                         "flavor" = flavor,
                         "font.type" = font.type,
                         "na.value" = na.value,
                         "legend.framecolor" = legend.framecolor,
                         "legend.tickcolor" = legend.tickcolor,
                         "legend.type" = legend.type,
                         "legend.position" = legend.position,
                         "viridis.palette" = viridis.palette,
                         "sequential.palette" = sequential.palette,
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
  
  check_parameters(parameter = legend.position, parameter_name = "legend.position")
  check_parameters(parameter = font.type, parameter_name = "font.type")
  check_parameters(parameter = legend.type, parameter_name = "legend.type")
  check_parameters(parameter = number.breaks, parameter_name = "number.breaks")
  check_parameters(parameter = diverging.palette, parameter_name = "diverging.palette")
  check_parameters(parameter = sequential.palette, parameter_name = "sequential.palette")
  check_parameters(parameter = sequential.direction, parameter_name = "sequential.direction")
  check_parameters(parameter = viridis.direction, parameter_name = "viridis.direction")
  check_parameters(parameter = viridis.palette, parameter_name = "viridis.palette")
  check_parameters(parameter = flavor, parameter_name = "flavor")
  check_parameters(plot.title.face, parameter_name = "plot.title.face")
  check_parameters(plot.subtitle.face, parameter_name = "plot.subtitle.face")
  check_parameters(plot.caption.face, parameter_name = "plot.caption.face")
  check_parameters(axis.title.face, parameter_name = "axis.title.face")
  check_parameters(axis.text.face, parameter_name = "axis.text.face")
  check_parameters(legend.title.face, parameter_name = "legend.title.face")
  check_parameters(legend.text.face, parameter_name = "legend.text.face")
  
  `%>%` <- magrittr::`%>%`
  
  
  if (!is.na(subsample)){
    sample <- sample[, sample(colnames(sample, subsample))]
  }
  
  if (is.null(group.by)){
    sample$Groups <- Seurat::Idents(sample)
    group.by <- "Groups"
  }
  
  genes.use <- unlist(input_gene_list) %>% unname() %>% unique()
  genes.use <- genes.use[genes.use %in% rownames(sample)]
  
  if (isTRUE(verbose)){message(paste0(add_info(), crayon_body("Computing "), crayon_key("enrichment scores"), crayon_body("...")))}
  
  sample <- SCpubr:::compute_enrichment_scores(sample, 
                                               input_gene_list = input_gene_list,
                                               nbin = nbin,
                                               ctrl = ctrl,
                                               flavor = flavor,
                                               assay = if (flavor == "UCell"){NULL} else {assay},
                                               slot = if (flavor == "Seurat"){NULL} else {slot})
  
  
  if (isTRUE(verbose)){message(paste0(add_info(), crayon_body("Plotting "), crayon_key("heatmaps"), crayon_body("...")))}
  
  # Obtain the DC embeddings, together with the enrichment scores.
  data.use <- sample@reductions[[reduction]]@cell.embeddings %>% 
              as.data.frame() %>% 
              tibble::rownames_to_column(var = "Cell") %>% 
              as.data.frame() %>% 
              tibble::as_tibble() %>% 
              tidyr::pivot_longer(cols = -dplyr::all_of(c("Cell")),
                                  names_to = "DC",
                                  values_to = "Score") %>% 
              dplyr::filter(.data$DC %in% sapply(dims, function(x){paste0("DC_", x)})) %>% 
              dplyr::group_by(.data$DC) %>% 
              dplyr::summarise("rank" = rank(.data$Score),
                               "Cell" = .data$Cell,
                               "Score" = .data$Score) %>% 
              dplyr::mutate("DC" = factor(.data$DC, levels = rev(sapply(dims, function(x){paste0("DC_", x)})))) %>% 
              dplyr::left_join(y = {sample@meta.data %>% 
                                    tibble::rownames_to_column(var = "Cell") %>% 
                                    tibble::as_tibble() %>% 
                                    dplyr::select(dplyr::all_of(c("Cell", group.by, names(input_gene_list))))},
                                    by = "Cell")
  
  if (isTRUE(scale.enrichment)){
    # Scale the enrichment scores as we are just interested in where they are enriched the most and not to compare across them.
    for (name in names(input_gene_list)){
      data.use[, name] <- scale(data.use[, name])[, 1]
    }
  }
  
  
  # Prepare the data to plot.
  data.use <- data.use %>% 
              tidyr::pivot_longer(cols = dplyr::all_of(c(names(input_gene_list))),
                                  names_to = "Gene_Set",
                                  values_to = "Enrichment")
  
  
  # Generate DC-based heatmaps.
  list.out <- list()
  for (dc.use in sapply(dims, function(x){paste0("DC_", x)})){
    # Filter for the DC.
    data.plot <- data.use %>% 
                  dplyr::filter(.data$DC == dc.use)
    
    # Limit the scale to quantiles 0.1 and 0.9 to avoid extreme outliers.
    limits <- c(stats::quantile(data.plot$Enrichment, 0.1),
                stats::quantile(data.plot$Enrichment, 0.9))
    
    # Bring extreme values to the cutoffs.
    data.plot <- data.plot %>% 
                 dplyr::mutate("Enrichment" = ifelse(.data$Enrichment <= limits[1], limits[1], .data$Enrichment)) %>% 
                 dplyr::mutate("Enrichment" = ifelse(.data$Enrichment >= limits[2], limits[2], .data$Enrichment))
    
    # Compute scale limits, breaks etc.
    scale.setup <- SCpubr:::compute_scales(sample = NULL,
                                           feature = NULL,
                                           assay = NULL,
                                           reduction = NULL,
                                           slot = NULL,
                                           number.breaks = 5,
                                           min.cutoff = NA,
                                           max.cutoff = NA,
                                           flavor = "Seurat",
                                           enforce_symmetry = enforce_symmetry,
                                           from_data = TRUE,
                                           limits.use = limits)
    
    # Generate the plot.
    p <- data.plot %>% 
         ggplot2::ggplot(mapping = ggplot2::aes(x = .data$rank,
                                                y = .data$Gene_Set,
                                                fill = .data$Enrichment)) + 
         ggplot2::geom_raster()
    
    legend.name <- if (flavor == "Seurat"){"Enrichment"} else if (flavor == "UCell"){"UCell score"} else if (flavor == "AUCell") {"AUC"}
    legend.name.use <- ifelse(isTRUE(scale.enrichment), paste0("Scaled + centered | ", legend.name), legend.name)
    if (isTRUE(enforce_symmetry)){
      p <- p + 
           ggplot2::scale_fill_gradientn(colors = rev(RColorBrewer::brewer.pal(n = 11, name = diverging.palette)),
                                         breaks = scale.setup$breaks,
                                         labels = scale.setup$labels,
                                         limits = scale.setup$limits,
                                         name = legend.name.use)
    } else {
      if (isTRUE(use_viridis)){
        p <- p + 
             ggplot2::scale_fill_viridis_c(option = viridis.palette,
                                           direction = viridis.direction,
                                           breaks = scale.setup$breaks,
                                           labels = scale.setup$labels,
                                           limits = scale.setup$limits,
                                           name = legend.name.use)
      } else {
        colors <- RColorBrewer::brewer.pal(n = 9, name = sequential.palette)[2:9]
        if (sequential.direction == -1){colors <- rev(colors)}
        p <- p + 
             ggplot2::scale_fill_gradientn(colors = colors,
                                           breaks = scale.setup$breaks,
                                           labels = scale.setup$labels,
                                           limits = scale.setup$limits,
                                           name = legend.name.use)
      }
    }
    p <- p + 
         ggplot2::xlab(paste0("Ordering of cells along ", dc.use)) + 
         ggplot2::ylab("Gene sets") +
         ggplot2::guides(y.sec = SCpubr:::guide_axis_label_trans(~paste0(levels(.data$Gene_Set)))) 
    
    # Modify the appearance of the plot.
    p <- SCpubr:::modify_continuous_legend(p = p,
                                           legend.title = legend.name.use,
                                           legend.aes = "fill",
                                           legend.type = legend.type,
                                           legend.position = legend.position,
                                           legend.length = legend.length,
                                           legend.width = legend.width,
                                           legend.framecolor = legend.framecolor,
                                           legend.tickcolor = legend.tickcolor,
                                           legend.framewidth = legend.framewidth,
                                           legend.tickwidth = legend.tickwidth)
    
    # Generate metadata plots to use on top of the main heatmap.
    counter <- 0
    list.plots <- list()
    list.plots[["main"]] <- p
    for (name in group.by){
      counter <- counter + 1
      if (counter == 1){
        name_labels <- name
      }
      
      # Select color palette for metadata.
      if (name %in% names(colors.use)){
        colors.use.iteration <- colors.use[[name]]
      } else {
        names.use <- if(is.factor(sample@meta.data[, name])){levels(sample@meta.data[, name])} else {sort(unique(sample@meta.data[, name]))}
        colors.use.iteration <- SCpubr:::generate_color_scale(names_use = names.use)
      }
      
      # Generate the metadata heatmap.
      p <- data.use %>% 
           dplyr::filter(.data$DC == dc.use) %>% 
           dplyr::mutate("grouped.var" = .env$group.by) %>% 
           ggplot2::ggplot(mapping = ggplot2::aes(x = .data$rank,
                                                  y = .data$grouped.var,
                                                  fill = .data[[group.by]])) + 
           ggplot2::geom_raster() + 
           ggplot2::scale_fill_manual(values = colors.use.iteration) + 
           ggplot2::guides(fill = ggplot2::guide_legend(title = name,
                                                        title.position = "top",
                                                        title.hjust = 0.5,
                                                        ncol = legend.ncol,
                                                        nrow = legend.nrow,
                                                        byrow = legend.byrow)) +
           ggplot2::xlab(NULL) +
           ggplot2::ylab(NULL) +
           ggplot2::guides(y.sec = SCpubr:::guide_axis_label_trans(~paste0(levels(.data$grouped.var)))) 
         
      list.plots[[name]] <- p
    }
    
    # Add theme to all plots.
    for (name in names(list.plots)){
      
      list.plots[[name]] <- list.plots[[name]] +
                            ggplot2::scale_x_discrete(expand = c(0, 0)) +
                            ggplot2::scale_y_discrete(expand = c(0, 0)) +
                            ggplot2::theme_minimal(base_size = font.size) +
                            ggplot2::theme(axis.text.x = ggplot2::element_blank(),
                                           axis.text.y.right = ggplot2::element_text(face = axis.text.face,
                                                                                     color = "black"),
                                           axis.text.y.left = ggplot2::element_blank(),
                                           axis.ticks.y.right = ggplot2::element_line(color = "black"),
                                           axis.ticks.y.left = ggplot2::element_blank(),
                                           axis.ticks.x = ggplot2::element_blank(),
                                           axis.line = ggplot2::element_blank(),
                                           axis.title.y = ggplot2::element_text(face = axis.title.face, color = "black", angle = 0, hjust = 0.5, vjust = 0.5),
                                           axis.title.x = ggplot2::element_text(face = axis.title.face, color = "black", angle = 0),
                                           plot.title = ggplot2::element_text(face = plot.title.face, hjust = 0),
                                           plot.subtitle = ggplot2::element_text(face = plot.subtitle.face, hjust = 0),
                                           plot.caption = ggplot2::element_text(face = plot.caption.face, hjust = 1),
                                           plot.title.position = "plot",
                                           panel.grid = ggplot2::element_blank(),
                                           panel.grid.minor.y = ggplot2::element_line(color = "white"),
                                           text = ggplot2::element_text(family = font.type),
                                           plot.caption.position = "plot",
                                           legend.text = ggplot2::element_text(face = legend.text.face),
                                           legend.position = legend.position,
                                           legend.title = ggplot2::element_text(face = legend.title.face),
                                           legend.justification = "center",
                                           plot.margin = ggplot2::margin(t = 10, r = 10, b = 0, l = 10),
                                           panel.border = ggplot2::element_rect(color = border.color, fill = NA),
                                           panel.grid.major = ggplot2::element_blank(),
                                           plot.background = ggplot2::element_rect(fill = "white", color = "white"),
                                           panel.background = ggplot2::element_rect(fill = "white", color = "white"),
                                           legend.background = ggplot2::element_rect(fill = "white", color = "white"))
    }
    
    # Reorder heatmaps for correct plotting.
    list.plots <- list.plots[c(group.by, "main")]
    height_unit <- c(rep((1 - main.heatmap.size) / length(group.by), length(group.by)), main.heatmap.size)
    
    
    # Assemble the final heatmap.
    p <- patchwork::wrap_plots(list.plots,
                               ncol = 1,
                               guides = "collect",
                               heights = height_unit) +
         patchwork::plot_annotation(theme = ggplot2::theme(legend.position = legend.position))
    
    list.out[[dc.use]] <- p
  }
  
  # Return the object.
  if (isTRUE(return_object)){
    list.out[["Object"]] <- sample
  }
  
  return(list.out)
}