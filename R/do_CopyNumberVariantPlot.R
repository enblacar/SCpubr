#' Display CNV scores from inferCNV as Feature Plots.
#'
#'
#' @inheritParams doc_function
#' @param infercnv_object \strong{\code{\link[infercnv]{infercnv}}} | Output inferCNV object run on the same Seurat object.
#' @param using_metacells \strong{\code{\link[base]{logical}}} | Whether inferCNV was run using metacells or not.
#' @param metacell_mapping \strong{\code{\link[SCpubr]{named_vector}}} | Vector or cell - metacell mapping.
#' @param chromosome_focus \strong{\code{\link[base]{character}}} | Region stating which chromosome to plot. Eg: 1p, 19q. NULL will plot all regions.
#' @param chromosome_locations \strong{\code{\link[tibble]{tibble}}} | Tibble containing the chromosome regions to use. Can be obtained using \code{utils::data("human_chr_locations", package = "SCpubr")}.
#'
#' @return A list containing Feature Plots for different chromosome regions and corresponding dot plots by groups..
#' @export
#'
#' @example man/examples/examples_do_CopyNumberVariantPlot.R
do_CopyNumberVariantPlot <- function(sample,
                                     infercnv_object,
                                     chromosome_locations,
                                     group.by = NULL,
                                     using_metacells = FALSE,
                                     metacell_mapping = NULL,
                                     chromosome_focus = NULL,
                                     legend.type = "colorbar",
                                     legend.position = "bottom",
                                     legend.length = 20,
                                     legend.width = 1,
                                     legend.framewidth = 0.5,
                                     legend.tickwidth = 0.5,
                                     legend.framecolor = "grey50",
                                     legend.tickcolor = "white",
                                     font.size = 14,
                                     pt.size = 1,
                                     font.type = "sans",
                                     rotate_x_axis_labels = 45,
                                     enforce_symmetry = TRUE,
                                     legend.title = NULL,
                                     na.value = "grey75",
                                     viridis.palette = "G",
                                     viridis.direction = 1,
                                     verbose = FALSE,
                                     min.cutoff = NA,
                                     max.cutoff = NA,
                                     number.breaks = 5,
                                     diverging.palette = "RdBu",
                                     sequential.palette = "YlGnBu",
                                     sequential.direction = -1,
                                     use_viridis = TRUE,
                                     return_object = FALSE,
                                     grid.color = "white",
                                     border.color = "black",
                                     flip = FALSE){
  check_suggests("do_CopyNumberVariantPlot")
  
  # Check logical parameters.
  logical_list <- list("using_metacells" = using_metacells,
                       "enforce_symmetry" = enforce_symmetry,
                       "use_viridis" = use_viridis)
  check_type(parameters = logical_list, required_type = "logical", test_function = is.logical)
  # Check numeric parameters.
  numeric_list <- list("font.size" = font.size,
                       "legend.length" = legend.length,
                       "legend.width" = legend.width,
                       "legend.framewidth" = legend.framewidth,
                       "legend.tickwidth" = legend.tickwidth,
                       "pt.size" = pt.size,
                       "viridis.direction" = viridis.direction,
                       "rotate_x_axis_labels" = rotate_x_axis_labels,
                       "min.cutoff" = min.cutoff,
                       "max.cutoff" = max.cutoff,
                       "number.breaks" = number.breaks,
                       "sequential.direction" = sequential.direction)
  check_type(parameters = numeric_list, required_type = "numeric", test_function = is.numeric)
  # Check character parameters.
  character_list <- list("group.by" = group.by,
                         "legend.type" = legend.type,
                         "legend.position" = legend.position,
                         "legend.framecolor" = legend.framecolor,
                         "legend.tickcolor" = legend.tickcolor,
                         "font.type" = font.type,
                         "legend.title" = legend.title,
                         "viridis.palette" = viridis.palette,
                         "diverging.palette" = diverging.palette,
                         "sequential.palette" = sequential.palette,
                         "grid.color" = grid.color,
                         "border.color" = border.color)
  check_type(parameters = character_list, required_type = "character", test_function = is.character)
  
  
  `:=` <- rlang::`:=`
  `%>%` <- magrittr::`%>%`
  
  check_colors(legend.framecolor, parameter_name = "legend.framecolor")
  check_colors(legend.tickcolor, parameter_name = "legend.tickcolor")
  check_colors(na.value, parameter_name = "na.value")
  check_colors(grid.color, parameter_name = "grid.color")
  check_colors(border.color, parameter_name = "border.color")
  
  check_parameters(parameter = font.type, parameter_name = "font.type")
  check_parameters(parameter = legend.type, parameter_name = "legend.type")
  check_parameters(parameter = legend.position, parameter_name = "legend.position")
  check_parameters(parameter = viridis.palette, parameter_name = "viridis.palette")
  check_parameters(parameter = viridis.direction, parameter_name = "viridis.direction")
  check_parameters(parameter = rotate_x_axis_labels, parameter_name = "rotate_x_axis_labels")
  check_parameters(parameter = number.breaks, parameter_name = "number.breaks")
  check_parameters(parameter = diverging.palette, parameter_name = "diverging.palette")
  check_parameters(parameter = sequential.palette, parameter_name = "sequential.palette")
  check_parameters(parameter = sequential.direction, parameter_name = "sequential.direction")
  
  chromosome_list <- c(as.character(seq(1, 22)))
  
  # Fix for group.by.
  if (is.null(group.by)){
    aggr_entities <- levels(sample)
    sample@meta.data[, "Groups"] <- sample@active.ident
    group.by <- "Groups"
  }
  
  
  # Retrieve the genes.
  genes <- infercnv_object@gene_order
  
  # Retrieve chr 1p start and end coordinates.
  chr_locations <- chromosome_locations
  
  # This list will contain all the outputs.
  return_list <- list()
  scores.assay <- data.frame(row.names = colnames(sample))
  for (chromosome in chromosome_list){
    # Retrieve chr locations of the chromosome.
    locations <- chr_locations %>%
                 dplyr::filter(.data[["chr"]] == chromosome)
    
    for (chr_arm in c("p", "q", "whole")){
      if (chr_arm != "whole"){
        # Retrieve the start.
        start <- locations %>%
                 dplyr::filter(.data[["arm"]] == chr_arm) %>%
                 dplyr::pull(start)
        # Retrieve the end.
        end <- locations %>%
               dplyr::filter(.data[["arm"]] == chr_arm) %>%
               dplyr::pull(end)
        
        # Retrieve the genes present in the chromosome arm.
        genes_use <- rownames(genes %>%
                              dplyr::filter(.data[["chr"]] == paste0("chr", chromosome),
                                            stop <= end))
      } else {
        # Retrieve the start.
        start <- locations %>%
                 dplyr::filter(.data[["arm"]] == "p") %>%
                 dplyr::pull(start)
        # Retrieve the end.
        end <- locations %>%
               dplyr::filter(.data[["arm"]] == "q") %>%
               dplyr::pull(end)
        
        # Retrieve the genes present in the chromosome arm.
        genes_use <- rownames(genes %>%
                              dplyr::filter(.data[["chr"]] == paste0("chr", chromosome)))
      }
      
      # Retrieve the CNV scores from the inferCNV object.
      CNV_scores <- infercnv_object@expr.data
      
      # Make it at least 2 genes in the object (this will only be applicable in VERY LOW QUALITY DATASETS)
      if (sum(genes_use %in% rownames(CNV_scores)) > 1){
        # Filter the scores for only the genes in the chromosome arm.
        CNV_scores <- CNV_scores[genes_use[genes_use %in% rownames(CNV_scores)], ]
        
        scores_name <- if (chr_arm != "whole"){paste0(chromosome, chr_arm)} else {chromosome}
        CNV_scores_final <- tibble::tibble(!!scores_name := colMeans(CNV_scores),
                                           "cells" = colnames(CNV_scores))
        
        
        # If metacells were used.
        if (isTRUE(using_metacells)){
          sample[["metacell_mapping"]] <- metacell_mapping
          scores.assay[[scores_name]] <- sample@meta.data %>%
                                         dplyr::mutate("cells" = colnames(sample)) %>%
                                         dplyr::left_join(y = {sample@meta.data %>%
                                                               dplyr::select(dplyr::all_of(c("metacell_mapping"))) %>%
                                                               tibble::rownames_to_column(var = "cells") %>%
                                                               dplyr::left_join(y = {CNV_scores_final %>% dplyr::rename("metacell_mapping" = dplyr::all_of(c("cells")))},
                                                                                by = "metacell_mapping") %>%
                                                               dplyr::select(-dplyr::all_of(c("metacell_mapping")))},
                                                          by = "cells") %>%
                                         tibble::column_to_rownames(var = "cells") %>% 
                                         dplyr::select(.env$scores_name)
          
          # If no metacells were used.
        } else if (isFALSE(using_metacells)){
          scores.assay[[scores_name]] <- CNV_scores_final[, scores_name]
        }
      } else {
        #nocov start
        if(isTRUE(verbose)){message(paste0("Your sample has only one gene in ", chromosome, chr_arm, ". Skipping this chromosome arm."))}
        #nocov end
      }
    }
  }
  
  # Generate an assay.
  assay <- scores.assay %>% 
           t() %>% 
           Seurat::CreateAssayObject(.)
  
  sample@assays$CNV_scores <- assay
  Seurat::DefaultAssay(sample) <- "CNV_scores"
  sample@assays$CNV_scores@key <- "CNV_scores_"
  
  
  list.data <- list()
  for (group in group.by){
    data.use  <- sample@assays$CNV_scores@data %>% 
                 as.matrix() %>% 
                 t() %>% 
                 as.data.frame() %>% 
                 tibble::rownames_to_column(var = "Cell") %>% 
                 dplyr::left_join(y = {sample@meta.data %>% 
                                       tibble::rownames_to_column(var = "Cell") %>% 
                                       dplyr::select(dplyr::all_of(c("Cell", group)))},
                                  by = "Cell") %>% 
                 tidyr::pivot_longer(cols = -dplyr::all_of(c("Cell", group)),
                                     values_to = "CNV_score",
                                     names_to = "Event") %>% 
                 
                 dplyr::group_by(.data[[group]], .data$Event) %>% 
                 dplyr::summarise("mean" = mean(.data$CNV_score, na.rm = TRUE))
    
    events <- c(as.character(seq(1,22)), sapply(seq(1,22), function(x){return(c(paste0(x, "p"), paste0(x, "q")))}))
    
    if (isFALSE(flip)){
      factor.levels <- events[events %in% unique(data.use$Event)]
    } else {
      factor.levels <- rev(events[events %in% unique(data.use$Event)])
    }
    data.use <- data.use %>% 
                dplyr::mutate("Event" = factor(.data$Event, levels = factor.levels))
    
    list.data[[group]][["data"]] <- data.use
  }
  
  # Compute limits.
  min.vector <- c()
  max.vector <- c()
  
  for (group in group.by){
    data.limits <- list.data[[group]][["data"]]
    
    min.vector <- append(min.vector, min(data.limits$mean, na.rm = TRUE))
    max.vector <- append(max.vector, max(data.limits$mean, na.rm = TRUE))
  }
  
  # Get the absolute limits of the datasets.
  limits <- c(min(min.vector, na.rm = TRUE),
              max(max.vector, na.rm = TRUE))
  
  # Compute overarching scales for all heatmaps.
  scale.setup <- compute_scales(sample = sample,
                                feature = " ",
                                assay = assay,
                                reduction = NULL,
                                slot = slot,
                                number.breaks = number.breaks,
                                min.cutoff = min.cutoff,
                                max.cutoff = max.cutoff,
                                flavor = "Seurat",
                                enforce_symmetry = TRUE,
                                from_data = TRUE,
                                limits.use = limits,
                                center_on_value = TRUE,
                                value_center = 1)
  list.plots <- list()
  if (isFALSE(flip)){
    values.use <- rev(group.by)
  } else {
    values.use <- group.by
  }
  for (group in values.use){
    data <- list.data[[group]][["data"]]
    p <- data %>% 
         ggplot2::ggplot(mapping = ggplot2::aes(x = if(isFALSE(flip)){.data$Event} else {.data[[group]]},
                                                y = if(isFALSE(flip)){.data[[group]]} else {.data$Event},
                                                fill = .data$mean)) + 
         ggplot2::geom_tile(color = grid.color, linewidth = 0.5) +
         ggplot2::scale_y_discrete(expand = c(0, 0)) +
         ggplot2::scale_x_discrete(expand = c(0, 0),
                                   position = "top") +
         ggplot2::guides(x.sec = guide_axis_label_trans(~paste0(levels(if(isFALSE(flip)){.data[[group]]} else {.data$Event}))),
                         y.sec = guide_axis_label_trans(~paste0(levels(if(isFALSE(flip)){.data[[group]]} else {.data$Event})))) + 
         ggplot2::coord_equal() + 
         ggplot2::scale_fill_gradientn(colors = RColorBrewer::brewer.pal(n = 11, name = diverging.palette),
                                       na.value = na.value,
                                       name =  "CNV scores",
                                       breaks = scale.setup$breaks,
                                       labels = scale.setup$labels,
                                       limits = scale.setup$limits)

    
    list.plots[[group]] <- p
  }
  
  # Modify legends.
  for (name in names(list.plots)){
    p <- list.plots[[name]]
    p <- modify_continuous_legend(p = p,
                                  legend.aes = "fill",
                                  legend.type = legend.type,
                                  legend.position = legend.position,
                                  legend.length = legend.length,
                                  legend.width = legend.width,
                                  legend.framecolor = legend.framecolor,
                                  legend.tickcolor = legend.tickcolor,
                                  legend.framewidth = legend.framewidth,
                                  legend.tickwidth = legend.tickwidth)
    list.plots[[name]] <- p
  }
  
  # Add theme
  counter <- 0
  for (name in rev(names(list.plots))){
    counter <- counter + 1
    
    if (isTRUE(flip)){
      if (counter == 1){
        xlab <- name
        ylab <- "CNV event"
      } else if (counter == length(names(list.plots))){
        xlab <- name
        ylab <- NULL
      } else {
        xlab <- name
        ylab <- NULL
      }
    } else {
      if (counter == 1){
        xlab <- NULL
        ylab <- name
      } else if (counter == length(names(list.plots))){
        xlab <- "CNV event"
        ylab <- name
      } else {
        xlab <- NULL
        ylab <- name
      }
    }
    
    p <- list.plots[[name]]
    
    axis.parameters <- handle_axis(flip = flip,
                                   group.by = rep("A", length(names(list.plots))),
                                   group = name,
                                   counter = counter,
                                   rotate_x_axis_labels = rotate_x_axis_labels)
    
    p <- p +
         ggplot2::xlab(xlab) +
         ggplot2::ylab(ylab) +
         ggplot2::theme_minimal(base_size = font.size) +
         ggplot2::theme(axis.ticks.x.bottom = axis.parameters$axis.ticks.x.bottom,
                        axis.ticks.x.top = axis.parameters$axis.ticks.x.top,
                        axis.ticks.y.left = axis.parameters$axis.ticks.y.left,
                        axis.ticks.y.right = axis.parameters$axis.ticks.y.right,
                        axis.text.y.left = axis.parameters$axis.text.y.left,
                        axis.text.y.right = axis.parameters$axis.text.y.right,
                        axis.text.x.top = axis.parameters$axis.text.x.top,
                        axis.text.x.bottom = axis.parameters$axis.text.x.bottom,
                        axis.title.x.bottom = axis.parameters$axis.title.x.bottom,
                        axis.title.x.top = axis.parameters$axis.title.x.top,
                        axis.title.y.right = axis.parameters$axis.title.y.right,
                        axis.title.y.left = axis.parameters$axis.title.y.left,
                        strip.background = axis.parameters$strip.background,
                        strip.clip = axis.parameters$strip.clip,
                        strip.text = axis.parameters$strip.text,
                        legend.position = legend.position,
                        axis.line = ggplot2::element_blank(),
                        plot.title = ggplot2::element_text(face = "bold", hjust = 0),
                        plot.subtitle = ggplot2::element_text(hjust = 0),
                        plot.caption = ggplot2::element_text(hjust = 1),
                        plot.title.position = "plot",
                        panel.grid = ggplot2::element_blank(),
                        panel.grid.minor.y = ggplot2::element_line(color = "white", linewidth = 1),
                        text = ggplot2::element_text(family = font.type),
                        plot.caption.position = "plot",
                        legend.text = ggplot2::element_text(face = "bold"),
                        legend.title = ggplot2::element_text(face = "bold"),
                        legend.justification = "center",
                        plot.margin = ggplot2::margin(t = 5, r = 0, b = 0, l = 5),
                        panel.border = ggplot2::element_rect(fill = NA, color = border.color, linewidth = 1),
                        panel.grid.major = ggplot2::element_blank(),
                        plot.background = ggplot2::element_rect(fill = "white", color = "white"),
                        panel.background = ggplot2::element_rect(fill = "white", color = "white"),
                        legend.background = ggplot2::element_rect(fill = "white", color = "white"),
                        panel.spacing.x = ggplot2::unit(0, "cm"))
       
    list.plots[[name]] <- p
  }
  
  # Plot the combined plot
  p <- patchwork::wrap_plots(list.plots[rev(group.by)],
                             ncol = if (isFALSE(flip)){1} else {NULL},
                             nrow = if(isTRUE(flip)) {1} else {NULL},
                             guides = "collect")
  p <- p +
       patchwork::plot_annotation(theme = ggplot2::theme(legend.position = legend.position,
                                                         plot.title = ggplot2::element_text(family = font.type,
                                                                                            color = "black",
                                                                                            face = "bold",
                                                                                            hjust = 0),
                                                         plot.subtitle = ggplot2::element_text(family = font.type,
                                                                                               color = "black",
                                                                                               hjust = 0),
                                                         plot.caption = ggplot2::element_text(family = "mono",
                                                                                              color = "black",
                                                                                              hjust = 1),
                                                         plot.caption.position = "plot"))
  
  if (isTRUE(return_object)){
    return_list <- list("Plot" = p,
                        "Object" = sample)
  } else {
    return_list <- p
  }
  return(return_list)
}
