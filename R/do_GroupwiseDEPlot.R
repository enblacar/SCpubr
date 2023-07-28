#' Compute a heatmap with the results of a group-wise DE analysis.
#'
#' @inheritParams doc_function
#' @param de_genes \strong{\code{\link[tibble]{tibble}}} | DE genes matrix resulting of running `Seurat::FindAllMarkers()`.
#' @param top_genes \strong{\code{\link[base]{numeric}}} | Top N differentially expressed (DE) genes by group to retrieve.
#'
#' @return A heatmap composed of 3 main panels: -log10(adjusted p-value), log2(FC) and mean expression by cluster.
#' @export
#'
#' @example /man/examples/examples_do_GroupwiseDEPlot.R
do_GroupwiseDEPlot <- function(sample,
                               de_genes,
                               group.by = NULL,
                               number.breaks = 5,
                               top_genes = 5,
                               use_viridis = FALSE,
                               viridis.direction = -1,
                               viridis.palette.pvalue = "C",
                               viridis.palette.logfc = "E",
                               viridis.palette.expression = "G",
                               sequential.direction = 1,
                               sequential.palette.pvalue = "YlGn",
                               sequential.palette.logfc = "YlOrRd",
                               sequential.palette.expression = "YlGnBu",
                               assay = NULL,
                               slot = "data",
                               legend.position = "bottom",
                               legend.width = 1,
                               legend.length = 20,
                               legend.framewidth = 0.5,
                               legend.tickwidth = 0.5,
                               legend.framecolor = "grey50",
                               legend.tickcolor = "white",
                               legend.type = "colorbar",
                               font.size = 14,
                               font.type = "sans",
                               axis.text.x.angle = 45,
                               min.cutoff = NA,
                               max.cutoff = NA,
                               na.value = "grey75",
                               grid.color = "white",
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
  
  check_suggests(function_name = "do_GroupwiseDEPlot")
  # Check if the sample provided is a Seurat object.
  check_Seurat(sample = sample)

  logical_list <- list("use_viridis" = use_viridis)
  check_type(parameters = logical_list, required_type = "logical", test_function = is.logical)
  # Check numeric parameters.
  numeric_list <- list("number.breaks" = number.breaks,
                       "top_genes" = top_genes,
                       "viridis.direction" = viridis.direction,
                       "legend.width" = legend.width,
                       "legend.length" = legend.length,
                       "legend.framewidth" = legend.framewidth,
                       "legend.tickwidth" = legend.tickwidth,
                       "font.size" = font.size,
                       "axis.text.x.angle" = axis.text.x.angle,
                       "min.cutoff" = min.cutoff,
                       "max.cutoff" = max.cutoff)
  check_type(parameters = numeric_list, required_type = "numeric", test_function = is.numeric)
  # Check character parameters.
  character_list <- list("group.by" = group.by,
                         "slot" = slot,
                         "legend.position" = legend.position,
                         "legend.framecolor" = legend.framecolor,
                         "legend.tickcolor" = legend.tickcolor,
                         "legend.type" = legend.type,
                         "font.type" = font.type,
                         "viridis.palette.pvalue" = viridis.palette.pvalue,
                         "viridis.palette.logfc" = viridis.palette.logfc,
                         "viridis.palette.expression" = viridis.palette.expression,
                         "sequential.palette.pvalue" = sequential.palette.pvalue,
                         "sequential.palette.logfc" = sequential.palette.logfc,
                         "sequential.palette.expression" = sequential.palette.expression,
                         "na.value" = na.value,
                         "grid.color" = grid.color,
                         "border.color" = border.color,
                         "plot.title.face" = plot.title.face,
                         "plot.subtitle.face" = plot.subtitle.face,
                         "plot.caption.face" = plot.caption.face,
                         "axis.title.face" = axis.title.face,
                         "axis.text.face" = axis.text.face,
                         "legend.title.face" = legend.title.face,
                         "legend.text.face" = legend.text.face)
  check_type(parameters = character_list, required_type = "character", test_function = is.character)


  `%>%` <- magrittr::`%>%`
  `:=` <- rlang::`:=`
  
  check_colors(legend.framecolor, parameter_name = "legend.framecolor")
  check_colors(legend.tickcolor, parameter_name = "legend.tickcolor")
  check_colors(na.value, parameter_name = "na.value")
  check_colors(grid.color, parameter_name = "grid.color")
  check_colors(border.color, parameter_name = "border.color")
  
  check_parameters(parameter = legend.position, parameter_name = "legend.position")
  check_parameters(parameter = viridis.palette.pvalue, parameter_name = "viridis_color_map")
  check_parameters(parameter = viridis.palette.logfc, parameter_name = "viridis_color_map")
  check_parameters(parameter = viridis.palette.expression, parameter_name = "viridis_color_map")
  check_parameters(parameter = sequential.palette.pvalue, parameter_name = "sequential.palette")
  check_parameters(parameter = sequential.palette.logfc, parameter_name = "sequential.palette")
  check_parameters(parameter = sequential.palette.expression, parameter_name = "sequential.palette")
  check_parameters(plot.title.face, parameter_name = "plot.title.face")
  check_parameters(plot.subtitle.face, parameter_name = "plot.subtitle.face")
  check_parameters(plot.caption.face, parameter_name = "plot.caption.face")
  check_parameters(axis.title.face, parameter_name = "axis.title.face")
  check_parameters(axis.text.face, parameter_name = "axis.text.face")
  check_parameters(legend.title.face, parameter_name = "legend.title.face")
  check_parameters(legend.text.face, parameter_name = "legend.text.face")
  check_parameters(viridis.direction, parameter_name = "viridis.direction")
  check_parameters(sequential.direction, parameter_name = "sequential.direction")
  
  # Check the assay.
  out <- check_and_set_assay(sample = sample, assay = assay)
  sample <- out[["sample"]]
  assay <- out[["assay"]]

  # Check group.by.
  out <- check_group_by(sample = sample,
                        group.by = group.by,
                        is.heatmap = TRUE)
  sample <- out[["sample"]]
  group.by <- out[["group.by"]]
  
  
  colors.gradient.pvalue <- compute_continuous_palette(name = ifelse(isTRUE(use_viridis), viridis.palette.pvalue, sequential.palette.pvalue),
                                                       use_viridis = use_viridis,
                                                       direction = ifelse(isTRUE(use_viridis), viridis.direction, sequential.direction),
                                                       enforce_symmetry = FALSE)
  
  colors.gradient.expression <- compute_continuous_palette(name = ifelse(isTRUE(use_viridis), viridis.palette.expression, sequential.palette.expression),
                                                           use_viridis = use_viridis,
                                                           direction = ifelse(isTRUE(use_viridis), viridis.direction, sequential.direction),
                                                           enforce_symmetry = FALSE)
  
  colors.gradient.logfc <- compute_continuous_palette(name = ifelse(isTRUE(use_viridis), viridis.palette.logfc, sequential.palette.logfc),
                                                      use_viridis = use_viridis,
                                                      direction = ifelse(isTRUE(use_viridis), viridis.direction, sequential.direction),
                                                      enforce_symmetry = FALSE)
  

  magnitude <- ifelse(slot == "data", "avg_log2FC", "avg_diff")
  specificity <- "p_val_adj"
  
  # Compute the top N genes per cluster.
  genes.use <- de_genes %>%
               dplyr::arrange(.data$p_val_adj, dplyr::desc(.data[[magnitude]])) %>%
               dplyr::group_by(.data$cluster) %>%
               dplyr::slice_head(n = top_genes) %>%
               dplyr::pull("gene") %>%
               unique()
  
  # Compute heatmap of log2FC.
  data.use <- de_genes %>%
              dplyr::arrange(.data[[specificity]], dplyr::desc(.data[[magnitude]])) %>%
              dplyr::group_by(.data$cluster) %>%
              dplyr::slice_head(n = top_genes) %>%
              dplyr::select(dplyr::all_of(c("gene", "cluster", magnitude, specificity)))
  max.cutoff.pval <- ceiling(-1 * (data.use %>% dplyr::filter(.data$p_val_adj != 0) %>% dplyr::pull(.data$p_val_adj) %>% min(na.rm = TRUE) %>% log10()))
  data.use <- data.use %>% 
              dplyr::mutate("-log10_padj" = -1 * log10(.data$p_val_adj),
                            "-log10_padj" = ifelse(is.infinite(.data$`-log10_padj`), max.cutoff.pval, .data$`-log10_padj`)) %>% 
              dplyr::mutate("specificity" = .data$`-log10_padj`,
                            "magnitude" = .data[[magnitude]])

  
  # Add missing data.
  data.use <- data.use %>% 
              dplyr::mutate("combination" = paste0(.data$gene, "_", .data$cluster))
  for (cluster in unique(data.use$cluster)){
    for (gene in unique(data.use$gene)){
      combination <- paste0(gene, "_", cluster)
      if (!combination %in% data.use$combination){
        row.add <- tibble::tibble("gene" = gene,
                                  "cluster" = cluster,
                                  "{magnitude}" := NA,
                                  "p_val_adj" = NA,
                                  "-log10_padj" = NA,
                                  "specificity" = NA,
                                  "magnitude" = NA,
                                  "combination" = combination)
        data.use <- rbind(data.use, row.add)
      }
    }
  }
  
  data.use <- data.use %>% 
              dplyr::select(-dplyr::all_of("combination")) %>% 
              dplyr::mutate("gene" = factor(.data$gene, levels = genes.use),
                            "cluster" = factor(.data$cluster, levels = rev(unique(data.use$cluster))))
  
  
  limits <- c(min(data.use$specificity, na.rm = TRUE),
              max.cutoff.pval)
  scale.setup <- compute_scales(sample = sample,
                                feature = " ",
                                assay = assay,
                                reduction = NULL,
                                slot = slot,
                                number.breaks = number.breaks,
                                min.cutoff = NA,
                                max.cutoff = max.cutoff,
                                flavor = "Seurat",
                                enforce_symmetry = FALSE,
                                from_data = TRUE,
                                limits.use = limits)
  
  list.plots <- list()
  
  p <- data.use %>% 
       ggplot2::ggplot(mapping = ggplot2::aes(x = .data$gene,
                                              y = .data$cluster,
                                              fill = .data$specificity)) + 
       ggplot2::geom_tile(color = grid.color, linewidth = 0.5) +
       ggplot2::scale_y_discrete(expand = c(0, 0)) +
       ggplot2::scale_x_discrete(expand = c(0, 0),
                                 position = "top") +
       ggplot2::guides(y.sec = guide_axis_label_trans(~paste0(levels(.data$cluster))),
                       x.sec = guide_axis_label_trans(~paste0(levels(.data$gene)))) + 
       ggplot2::coord_equal() + 
       ggplot2::scale_fill_gradientn(colors = colors.gradient.pvalue,
                                     na.value = na.value,
                                     name = expression(bold(paste("-", log["10"], "(p.adjust)"))),
                                     breaks = scale.setup$breaks,
                                     labels = scale.setup$labels,
                                     limits = scale.setup$limits)
          
  list.plots[["pval"]] <- p
  
  
  
  # Heatmap of avg_log2FC.
  limits <- c(min(data.use$magnitude, na.rm = TRUE),
              max(data.use$magnitude, na.rm = TRUE))
  scale.setup <- compute_scales(sample = sample,
                                feature = " ",
                                assay = assay,
                                reduction = NULL,
                                slot = slot,
                                number.breaks = number.breaks,
                                min.cutoff = NA,
                                max.cutoff = NA,
                                flavor = "Seurat",
                                enforce_symmetry = FALSE,
                                from_data = TRUE,
                                limits.use = limits)
  

  p <- data.use %>% 
       ggplot2::ggplot(mapping = ggplot2::aes(x = .data$gene,
                                              y = .data$cluster,
                                              fill = .data$magnitude)) + 
       ggplot2::geom_tile(color = grid.color, linewidth = 0.5) +
       ggplot2::scale_y_discrete(expand = c(0, 0)) +
       ggplot2::scale_x_discrete(expand = c(0, 0),
                                 position = "top") +
       ggplot2::guides(y.sec = guide_axis_label_trans(~paste0(levels(.data$cluster))),
                       x.sec = guide_axis_label_trans(~paste0(levels(.data$gene)))) + 
       ggplot2::coord_equal() + 
       ggplot2::scale_fill_gradientn(colors = colors.gradient.logfc,
                                     na.value = na.value,
                                     name = expression(bold(paste("Avg. ", log["2"], "(FC)"))),
                                     breaks = scale.setup$breaks,
                                     labels = scale.setup$labels,
                                     limits = scale.setup$limits)
  list.plots[["FC"]] <- p
  
  
  
  # Add averaged expression data.
  list.exp <- list()
  for (group in group.by){
    order.use <- if (is.factor(sample@meta.data[, group])){levels(sample@meta.data[, group])} else {sort(unique(sample@meta.data[, group]))}
    data <- .GetAssayData(sample,
                                 assay = assay,
                                 slot = slot)[genes.use, ] %>% 
            as.data.frame() %>% 
            tibble::rownames_to_column(var = "gene") %>% 
            tidyr::pivot_longer(cols = -dplyr::all_of("gene"),
                                names_to = "cell",
                                values_to = "expression") %>% 
            dplyr::left_join(y = {sample@meta.data %>% 
                                  tibble::rownames_to_column(var = "cell") %>% 
                                  dplyr::select(dplyr::all_of(c("cell", group)))},
                             by = "cell") %>% 
            dplyr::group_by(.data$gene, .data[[group]]) %>% 
            dplyr::summarize("Avg.Exp" = mean(.data$expression)) %>% 
            dplyr::mutate("gene" = factor(.data$gene, levels = genes.use),
                          "Group" = factor(.data[[group]], levels = rev(order.use)))
    
    list.exp[[group]] <- data
  }
  
  
  
  # Compute limits.
  min.vector <- NULL
  max.vector <- NULL
  
  for (group in group.by){
    data.limits <- list.exp[[group]]
    
    min.vector <- append(min.vector, min(data.limits$Avg.Exp, na.rm = TRUE))
    max.vector <- append(max.vector, max(data.limits$Avg.Exp, na.rm = TRUE))
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
                                enforce_symmetry = FALSE,
                                from_data = TRUE,
                                limits.use = limits)
  
  for (group in group.by){
    data <- list.exp[[group]]
    
    if (!is.na(min.cutoff)){
      data <- data %>% 
              dplyr::mutate("Avg.Exp" = ifelse(.data$Avg.Exp < min.cutoff, min.cutoff, .data$Avg.Exp))
    }
    
    if (!is.na(max.cutoff)){
      data <- data %>% 
              dplyr::mutate("Avg.Exp" = ifelse(.data$Avg.Exp > max.cutoff, max.cutoff, .data$Avg.Exp))
    }
    
    p <- data %>% 
         ggplot2::ggplot(mapping = ggplot2::aes(x = .data$gene,
                                                y = .data$Group,
                                                fill = .data$Avg.Exp)) + 
         ggplot2::geom_tile(color = grid.color, linewidth = 0.5) +
         ggplot2::scale_y_discrete(expand = c(0, 0)) +
         ggplot2::scale_x_discrete(expand = c(0, 0),
                                   position = "top") +
         ggplot2::guides(y.sec = guide_axis_label_trans(~paste0(levels(.data$Group))),
                         x.sec = guide_axis_label_trans(~paste0(levels(.data$gene)))) + 
         ggplot2::coord_equal() + 
         ggplot2::scale_fill_gradientn(colors = colors.gradient.expression,
                                       na.value = na.value,
                                       name = "Avg. Expression",
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
    if (name == "pval"){
      xlab <- "Genes"
      ylab <- expression(bold(paste("-", log["10"], "(p.adjust)")))
    } else if (name == "FC"){
      xlab <- NULL
      ylab <- expression(bold(paste("Avg. ", log["2"], "(FC)")))
    } else {
      xlab <- NULL
      ylab <- paste0("Avg. Exp | ", name)
    }
    
    counter <- counter + 1
    p <- list.plots[[name]]
    
    axis.parameters <- handle_axis(flip = FALSE,
                                   group.by = rep("A", length(names(list.plots))),
                                   group = name,
                                   counter = counter,
                                   axis.text.x.angle = axis.text.x.angle,
                                   plot.title.face = plot.title.face,
                                   plot.subtitle.face = plot.subtitle.face,
                                   plot.caption.face = plot.caption.face,
                                   axis.title.face = axis.title.face,
                                   axis.text.face = axis.text.face,
                                   legend.title.face = legend.title.face,
                                   legend.text.face = legend.text.face)
    
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
                        plot.title = ggplot2::element_text(face = plot.title.face, hjust = 0),
                        plot.subtitle = ggplot2::element_text(face = plot.subtitle.face, hjust = 0),
                        plot.caption = ggplot2::element_text(face = plot.caption.face, hjust = 1),
                        legend.text = ggplot2::element_text(face = legend.text.face),
                        legend.title = ggplot2::element_text(face = legend.title.face),
                        plot.title.position = "plot",
                        panel.grid = ggplot2::element_blank(),
                        panel.grid.minor.y = ggplot2::element_line(color = "white", linewidth = 1),
                        text = ggplot2::element_text(family = font.type),
                        plot.caption.position = "plot",
                        legend.justification = "center",
                        plot.margin = ggplot2::margin(t = 5, r = 0, b = 0, l = 0),
                        panel.border = ggplot2::element_rect(fill = NA, color = border.color, linewidth = 1),
                        panel.grid.major = ggplot2::element_blank(),
                        plot.background = ggplot2::element_rect(fill = "white", color = "white"),
                        panel.background = ggplot2::element_rect(fill = "white", color = "white"),
                        legend.background = ggplot2::element_rect(fill = "white", color = "white"),
                        panel.spacing.x = ggplot2::unit(0, "cm"))
    
    list.plots[[name]] <- p
  }
  
  p <- patchwork::wrap_plots(list.plots,
                             ncol = 1,
                             guides = "collect")
  p <- p +
       patchwork::plot_annotation(theme = ggplot2::theme(legend.position = legend.position,
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
  # Return the final heatmap.
  return(p)
}
