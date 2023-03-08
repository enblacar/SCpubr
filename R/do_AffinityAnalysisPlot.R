do_AffinityAnalysisPlot <- function(sample,
                                    input_gene_list,
                                    subsample = 100,
                                    group.by = NULL,
                                    assay = "SCT",
                                    slot = "data",
                                    statistic = "norm_wmean",
                                    number.breaks = 5,
                                    use_viridis = FALSE,
                                    viridis.palette = "G",
                                    viridis.direction = -1,
                                    sequential.palette = "YlGnBu",
                                    sequential.direction = 1,
                                    diverging.palette = "RdBu",
                                    enforce_symmetry = TRUE,
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
                                    rotate_x_axis_labels = 45,
                                    flip = FALSE){
  `%>%` <- magrittr::`%>%`
  
  # Subsample input object based on provided identities.
  sample$subset.me <- Seurat::Idents(sample)
  
  subset.vector <- sample@meta.data %>% 
                   tibble::rownames_to_column(var = "cell") %>% 
                   dplyr::select(dplyr::all_of(c("cell", "subset.me"))) %>% 
                   dplyr::group_by(.data$subset.me) %>% 
                   dplyr::slice_sample(n = subsample) %>% 
                   dplyr::pull("cell")
  
  sample <- sample[, subset.vector]
  
  
  # Generate a network with the names of the list of genes as source and the gene sets as targets with 1 of mode of regulation.
  network <- input_gene_list %>% 
             tibble::as_tibble() %>% 
             tidyr::pivot_longer(cols = dplyr::everything(),
                                 names_to = "source",
                                 values_to = "target") %>% 
             dplyr::mutate("mor" = 1)
  
  # Get expression data.
  mat <- Seurat::GetAssayData(sample,
                              assay = assay,
                              slot = slot)
  
  # Compute activities.
  acts <- decoupleR::run_wmean(mat = mat, 
                               network = network)
  
  # Turn them into a matrix compatible to turn into a Seurat assay.
  acts.matrix <- acts %>% 
                 dplyr::filter(.data$statistic == .env$statistic) %>% 
                 dplyr::mutate("score" = ifelse(.data$p_value <= 0.05, .data$score, NA)) %>% 
                 tidyr::pivot_wider(id_cols = dplyr::all_of(c("source")),
                                    names_from = "condition",
                                    values_from = "score") %>%
                 tibble::column_to_rownames('source')
  
  # Generate a Seurat assay.
  assay.add <- Seurat::CreateAssayObject(acts.matrix)
  
  # Add the assay to the Seurat object.
  sample@assays$affinity <- assay.add
  
  # Set it as default assay.
  Seurat::DefaultAssay(sample) <- "affinity"
  
  # Scale and center the activity data.
  scale.data <- Seurat::GetAssayData(sample,
                                     assay = "affinity",
                                     slot = "data") %>%
                as.matrix() %>%
                t() %>%
                as.data.frame() %>%
                scale() %>%
                t()
  
  # Set it to the scale.data slot.
  sample@assays$affinity@scale.data <- scale.data

  # Plotting.
  # Get the data frames per group.by value for plotting.
  list.data <- list()
  counter <- 0
  for (group in group.by){
    counter <- counter + 1
    data.use <- Seurat::GetAssayData(sample, 
                                     assay = "affinity", 
                                     slot = "scale.data") %>% 
                t() %>% 
                as.data.frame() %>% 
                tibble::rownames_to_column(var = "cell") %>% 
                dplyr::left_join(y = {sample@meta.data %>% 
                                      tibble::rownames_to_column(var = "cell") %>% 
                                      dplyr::select(dplyr::all_of(c("cell", group)))},
                                 by = "cell") %>% 
                tidyr::pivot_longer(cols = -dplyr::all_of(c("cell", group)),
                                    names_to = "source",
                                    values_to = "score")
    
    # Clustering based on the median across all cells.
    data.cluster <- data.use %>% 
                    tidyr::pivot_wider(id_cols = dplyr::all_of(c("cell", group)),
                                       names_from = "source",
                                       values_from = "score") %>% 
                    dplyr::group_by(.data[[group]]) %>% 
                    dplyr::summarise(dplyr::across(.cols = dplyr::all_of(c(names(input_gene_list))),
                                                   stats::median,
                                                   na.rm = TRUE)) %>% 
                    as.data.frame() %>% 
                    tibble::column_to_rownames(var = group) 
    
    list.data[[group]][["data"]] <- data.use
    list.data[[group]][["data.cluster"]] <- data.cluster
  }
  
  # Plot individual heatmaps.
  
  list.heatmaps <- list()
  counter <- 0
  for (group in group.by){
    counter <- counter + 1
    
    data.use <- list.data[[group]][["data"]]
    data.cluster <- list.data[[group]][["data.cluster"]]
    
    if (counter == 1){
      if (length(colnames(data.cluster)) == 1){
        col_order <- colnames(data.cluster)[1]
      } else {
        col_order <- colnames(data.cluster)[stats::hclust(stats::dist(t(data.cluster), method = "euclidean"), method = "ward.D")$order]
      }
    }
    if(length(rownames(data.cluster)) == 1){
      row_order <- rownames(data.cluster)[1]
    } else {
      row_order <- rownames(data.cluster)[stats::hclust(stats::dist(data.cluster, method = "euclidean"), method = "ward.D")$order]
    }
    
    
    p <- data.use %>% 
         dplyr::group_by(.data[[group]], .data$source) %>% 
         dplyr::summarise("mean" = mean(.data$score, na.rm = TRUE)) %>% 
         dplyr::mutate("source" = factor(.data$source, levels = col_order),
                       "target" = factor(.data[[group]], levels = row_order)) %>% 
         ggplot2::ggplot(mapping = ggplot2::aes(x = if (isTRUE(flip)){.data$source} else {.data$target},
                                                y = if (isTRUE(flip)){.data$target} else {.data$source},
                                                fill = .data$mean)) +
         ggplot2::geom_tile(color = "white", linewidth = 0.5, na.rm = TRUE) +
         ggplot2::scale_y_discrete(expand = c(0, 0)) +
         ggplot2::scale_x_discrete(expand = c(0, 0),
                                   position = "top") +
         ggplot2::guides(y.sec = SCpubr:::guide_axis_label_trans(~paste0(levels(if (isTRUE(flip)){.data$target} else {.data$source}))),
                         x.sec = SCpubr:::guide_axis_label_trans(~paste0(levels(if (isTRUE(flip)){.data$source} else {.data$target})))) + 
         ggplot2::coord_equal() 
    list.heatmaps[[group]] <- p
  }
  
  
  # Compute limits.
  min.vector <- c()
  max.vector <- c()
  
  for (group in group.by){
    data.limits <- list.data[[group]][["data"]]
    
    min.vector <- append(min.vector, min(data.limits$score, na.rm = TRUE))
    max.vector <- append(max.vector, max(data.limits$score, na.rm = TRUE))
  }
  
  # Get the absolute limits of the datasets.
  limits <- c(min(min.vector, na.rm = TRUE),
              max(max.vector, na.rm = TRUE))
  
  # Compute overarching scales for all heatmaps.
  scale.setup <- SCpubr:::compute_scales(sample = sample,
                                feature = " ",
                                assay = assay,
                                reduction = NULL,
                                slot = slot,
                                number.breaks = number.breaks,
                                min.cutoff = NA,
                                max.cutoff = NA,
                                flavor = "Seurat",
                                enforce_symmetry = enforce_symmetry,
                                from_data = TRUE,
                                limits.use = limits)
  
  for (group in group.by){
    p <- list.heatmaps[[group]]
    
    if (isFALSE(enforce_symmetry)){
      if (isTRUE(use_viridis)){
        p <- p + 
             ggplot2::scale_fill_viridis_c(direction = viridis.direction,
                                           option = viridis.palette,
                                           na.value = "grey75",
                                           breaks = scale.setup$breaks,
                                           labels = scale.setup$labels,
                                           limits = scale.setup$limits,
                                           name = statistic)
      } else {
        p <- p + 
             ggplot2::scale_fill_gradientn(colors = if(sequential.direction == 1){RColorBrewer::brewer.pal(n = 9, name = sequential.palette)[2:9]} else {rev(RColorBrewer::brewer.pal(n = 9, name = sequential.palette)[2:9])},
                                           na.value = "grey75",
                                           name =  statistic,
                                           breaks = scale.setup$breaks,
                                           labels = scale.setup$labels,
                                           limits = scale.setup$limits)
      }
    } else {
      p <- p + 
           ggplot2::scale_fill_gradientn(colors = RColorBrewer::brewer.pal(n = 11, name = diverging.palette) %>% rev(),
                                         na.value = "grey75",
                                         name = statistic,
                                         breaks = scale.setup$breaks,
                                         labels = scale.setup$labels,
                                         limits = scale.setup$limits)
    }
    
    list.heatmaps[[group]] <- p
  }
  
  # Modify legends.
  for (group in group.by){
    p <- list.heatmaps[[group]]
    
    p <- SCpubr:::modify_continuous_legend(p = p,
                                  legend.aes = "fill",
                                  legend.type = legend.type,
                                  legend.position = legend.position,
                                  legend.length = legend.length,
                                  legend.width = legend.width,
                                  legend.framecolor = legend.framecolor,
                                  legend.tickcolor = legend.tickcolor,
                                  legend.framewidth = legend.framewidth,
                                  legend.tickwidth = legend.tickwidth)
    list.heatmaps[[group]] <- p
  }
  
  # Add theme
  counter <- 0
  for (group in group.by){
    counter <- counter + 1
    
    p <- list.heatmaps[[group]]
    
    # Set axis titles.
    if (isTRUE(flip)){
      if (counter == 1){
        ylab <- group
        xlab <- NULL
      } else {
        xlab <- "Gene set"
        ylab <- group
      }
    } else {
      if (counter == 1){
        ylab <- "Gene set"
        xlab <- group
      } else {
        ylab <- NULL
        xlab <- group
      }
    }
    
    
    p <- list.heatmaps[[group]]
    
    axis.parameters <- SCpubr:::handle_axis(flip = !flip,
                                   group.by = rep("A", length(group.by)),
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
                        panel.border = ggplot2::element_rect(fill = NA, color = "black", linewidth = 1),
                        panel.grid.major = ggplot2::element_blank(),
                        plot.background = ggplot2::element_rect(fill = "white", color = "white"),
                        panel.background = ggplot2::element_rect(fill = "white", color = "white"),
                        legend.background = ggplot2::element_rect(fill = "white", color = "white"),
                        panel.spacing.x = ggplot2::unit(0, "cm"))
    
    list.heatmaps[[group]] <- p
  }
  
  if (isTRUE(flip)){
    list.heatmaps <- list.heatmaps[rev(group.by)]
  }
  p <- patchwork::wrap_plots(list.heatmaps,
                             ncol = if (isFALSE(flip)){NULL} else {1},
                             nrow = if(isFALSE(flip)){1} else {NULL},
                             guides = "collect")
  p <- p +
       patchwork::plot_annotation(theme = ggplot2::theme(legend.position = legend.position,
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
  
  list.output <- list()
  
  list.output[["Heatmap"]] <- p
  
  for (group in group.by){
    for (gene_set in names(input_gene_list)){
      p <- SCpubr::do_BoxPlot(sample, 
                              feature = gene_set,
                              group.by = group,
                              order = TRUE, 
                              flip = !flip)
      list.output[["Box plots"]][[gene_set]][[group]] <- p
    }
  }

  return(p)
}