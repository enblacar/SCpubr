#' Compute a heatmap summary of the top and bottom genes in the PCA loadings for the desired PCs in a Seurat object.
#'
#' @inheritParams doc_function
#' @param subsample \strong{\code{\link[base]{numeric}}} | Number of cells to subsample the Seurat object to increase computational speed. Use NA to include the Seurat object as is. 
#' @param dims \strong{\code{\link[base]{numeric}}} | PCs to include in the analysis.
#' @param top_loadings \strong{\code{\link[base]{numeric}}} | Number of top and bottom scored genes in the PCA Loadings for each PC.
#' @param min.cutoff.loadings,max.cutoff.loadings \strong{\code{\link[base]{numeric}}} | Cutoff to subset the scale of the Loading score heatmap. NA will use quantiles 0.05 and 0.95.
#' @param min.cutoff.expression,max.cutoff.expression \strong{\code{\link[base]{numeric}}} | Cutoff to subset the scale of the expression heatmap. NA will use 0 (no quantile) and quantile 0.95.
#'
#' @return A ggplot2 object.
#' @export
#'
#' @examples NULL
do_LoadingsPlot <- function(sample,
                            group.by = NULL,
                            subsample = NA,
                            dims = 1:10,
                            top_loadings = 5,
                            assay = "SCT",
                            slot = "data",
                            grid.color = "white",
                            border.color = "black",
                            number.breaks = 5,
                            na.value = "grey75",
                            legend.position = "bottom",
                            legend.title = "Expression",
                            legend.type = "colorbar",
                            legend.framewidth = 0.5,
                            legend.tickwidth = 0.5,
                            legend.length = 20,
                            legend.width = 1,
                            legend.framecolor = "grey50",
                            legend.tickcolor = "white",
                            font.size = 14,
                            font.type = "sans",
                            axis.text.x.angle = 45,
                            use_viridis = FALSE,
                            sequential.direction = 1,
                            sequential.palette = "YlGnBu",
                            viridis.palette = "G",
                            viridis.direction = -1,
                            diverging.palette = "RdBu",
                            flip = FALSE,
                            min.cutoff.loadings = NA,
                            max.cutoff.loadings = NA,
                            min.cutoff.expression = NA,
                            max.cutoff.expression = NA,
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
  
  check_suggests("do_LoadingsPlot")
  
  # Check logical parameters.
  logical_list <- list("use_viridis" = use_viridis,
                       "flip" = flip)
  check_type(parameters = logical_list, required_type = "logical", test_function = is.logical)
  # Check numeric parameters.
  numeric_list <- list("axis.text.x.angle" = axis.text.x.angle,
                       "legend.width" = legend.width,
                       "legend.length" = legend.length,
                       "legend.framewidth" = legend.framewidth,
                       "legend.tickwidth" = legend.tickwidth,
                       "font.size" = font.size,
                       "number.breaks" = number.breaks,
                       "viridis.direction" = viridis.direction,
                       "sequential.direction" = sequential.direction,
                       "min.cutoff.loadings" = min.cutoff.loadings,
                       "max.cutoff.loadings" = max.cutoff.loadings,
                       "min.cutoff.expression" = min.cutoff.expression,
                       "max.cutoff.expression" = max.cutoff.expression)
  check_type(parameters = numeric_list, required_type = "numeric", test_function = is.numeric)
  
  # Check character parameters.
  character_list <- list("legend.type" = legend.type,
                         "font.type" = font.type,
                         "legend.position" = legend.position,
                         "legend.framecolor" = legend.framecolor,
                         "legend.tickcolor" = legend.tickcolor,
                         "group.by" = group.by,
                         "na.value" = na.value,
                         "slot" = slot,
                         "assay" = assay,
                         "group.by" = group.by,
                         "diverging.palette" = diverging.palette,
                         "sequential.palette" = sequential.palette,
                         "viridis.palette" = viridis.palette,
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
  
  check_colors(na.value)
  check_colors(legend.framecolor)
  check_colors(legend.tickcolor)
  check_colors(grid.color)
  check_colors(border.color)
  
  check_parameters(parameter = legend.position, parameter_name = "legend.position")
  check_parameters(parameter = font.type, parameter_name = "font.type")
  check_parameters(parameter = legend.type, parameter_name = "legend.type")
  check_parameters(parameter = number.breaks, parameter_name = "number.breaks")
  check_parameters(parameter = diverging.palette, parameter_name = "diverging.palette")
  check_parameters(parameter = sequential.palette, parameter_name = "sequential.palette")
  check_parameters(parameter = viridis.palette, parameter_name = "viridis.palette")
  check_parameters(parameter = viridis.direction, parameter_name = "viridis.direction")
  check_parameters(parameter = sequential.direction, parameter_name = "sequential.direction")
  check_parameters(plot.title.face, parameter_name = "plot.title.face")
  check_parameters(plot.subtitle.face, parameter_name = "plot.subtitle.face")
  check_parameters(plot.caption.face, parameter_name = "plot.caption.face")
  check_parameters(axis.title.face, parameter_name = "axis.title.face")
  check_parameters(axis.text.face, parameter_name = "axis.text.face")
  check_parameters(legend.title.face, parameter_name = "legend.title.face")
  check_parameters(legend.text.face, parameter_name = "legend.text.face")
  
  
  `%>%` <- magrittr::`%>%`
  `:=` <- rlang::`:=`
  
  if (is.null(group.by)){
    sample$Groups <- Seurat::Idents(sample)
    group.by <- "Groups"
  }
  
  if (!is.na(subsample)){
    sample <- sample[, sample(colnames(sample), 2500, replace = FALSE)]
  }
  
  loadings <- Seurat::Loadings(sample)[, dims] %>% 
              as.data.frame() %>% 
              tibble::rownames_to_column(var = "Gene") %>% 
              tidyr::pivot_longer(cols = -dplyr::all_of(dplyr::all_of(c("Gene"))),
                                  values_to = "Loading_Score",
                                  names_to = "PC")
  
  top_loadings.up <- loadings %>% 
                     dplyr::group_by(.data$PC) %>% 
                     dplyr::arrange(dplyr::desc(.data$Loading_Score)) %>% 
                     dplyr::slice_head(n = top_loadings) %>% 
                     dplyr::pull(.data$Gene)
  
  top_loadings.down <- loadings %>% 
                       dplyr::group_by(.data$PC) %>% 
                       dplyr::arrange(.data$Loading_Score) %>% 
                       dplyr::slice_head(n = top_loadings) %>% 
                       dplyr::pull(.data$Gene)
  
  genes.use <- c()
  
  for (i in seq(1, length(dims) * top_loadings, by = top_loadings)){
    range <- seq(i, i + (top_loadings  - 1))
    genes.add <- c(top_loadings.up[range], top_loadings.down[range])
    genes.add <- genes.add[!(genes.add %in% genes.use)]
    genes.use <- append(genes.use, genes.add)
  }
  
  loadings <- loadings %>% 
              dplyr::filter(.data$Gene %in% genes.use)
  
  embeddings <- Seurat::Embeddings(sample, reduction = "pca")[, dims] %>% 
                as.data.frame() %>% 
                tibble::rownames_to_column(var = "Cell") %>% 
                tidyr::pivot_longer(cols = -dplyr::all_of(dplyr::all_of(c("Cell"))),
                                    values_to = "Embedding_Score",
                                    names_to = "PC")
  
  metadata <- sample@meta.data %>% 
              as.data.frame() %>% 
              tibble::rownames_to_column(var = "Cell") %>% 
              dplyr::select(dplyr::all_of(c("Cell", group.by))) %>% 
              tibble::as_tibble()
  
  data.use <- metadata %>% 
              dplyr::left_join(y = embeddings,
                               by = "Cell") %>% 
              dplyr::left_join(y = loadings,
                               by = "PC")
  
  data.use <- data.use %>% 
              dplyr::left_join(y = {Seurat::GetAssayData(sample,
                                                         assay = assay,
                                                         slot = slot)[unique(data.use$Gene), ] %>%
                                    as.matrix() %>% 
                                    t() %>% 
                                    as.data.frame() %>% 
                                    tibble::rownames_to_column(var = "Cell") %>% 
                                    tidyr::pivot_longer(cols = -dplyr::all_of(c("Cell")),
                                                        names_to = "Gene",
                                                        values_to = "Expression")},
                               by = c("Gene", "Cell")) %>% 
              dplyr::mutate("Gene" = factor(.data$Gene, levels = genes.use))
  
  data.loading <- data.use %>% 
                  dplyr::group_by(.data$Gene, .data$PC) %>% 
                  dplyr::summarise("mean_Loading_Score" = mean(.data$Loading_Score, na.rm = TRUE))
  
  data.expression <- data.use %>% 
                     dplyr::group_by(.data[[group.by]], .data$Gene) %>% 
                     dplyr::summarise("mean_Expression" = mean(.data$Expression, na.rm = TRUE))
  
  data.expression.wide <- data.expression %>% 
                          tidyr::pivot_wider(names_from = "Gene",
                                             values_from = "mean_Expression") %>% 
                          as.data.frame() %>% 
                          tibble::column_to_rownames(var = group.by)
  
  data.loadings.wide <- data.loading %>% 
                        tidyr::pivot_wider(names_from = "Gene",
                                           values_from = "mean_Loading_Score") %>% 
                        as.data.frame() %>% 
                        tibble::column_to_rownames(var = "PC")
  
  # Cluster items.
  gene.order <- genes.use[stats::hclust(stats::dist(t(data.expression.wide), method = "euclidean"), method = "ward.D")$order]
  group.order <- if(is.factor(data.expression[[group.by]])){levels(data.expression[[group.by]])} else {sort(unique(data.expression[[group.by]]))}
  group.order <- group.order[stats::hclust(stats::dist(data.expression.wide, method = "euclidean"), method = "ward.D")$order]
  pc.order <- if(is.factor(data.loading[["PC"]])){levels(data.loading[["PC"]])} else {sort(unique(data.loading[["PC"]]))}
  pc.order <- pc.order[stats::hclust(stats::dist(data.loadings.wide, method = "euclidean"), method = "ward.D")$order]
  
  # Reorder items.
  data.loading <- data.loading %>% 
                  dplyr::mutate("PC" = factor(.data$PC, levels = pc.order),
                                "Gene" = factor(.data$Gene, levels = gene.order))
  
  data.expression <- data.expression %>% 
                     dplyr::mutate("{group.by}" := factor(.data[[group.by]], levels = group.order),
                                   "Gene" = factor(.data$Gene, levels = gene.order))
  
 
  
  # Apply cutoffs.
  if (!is.na(min.cutoff.loadings)){
    data.loading <- data.loading %>% 
                    dplyr::mutate("mean_Loading_Score" = ifelse(.data$mean_Loading_Score < min.cutoff.loadings, min.cutoff.loadings, .data$mean_Loading_Score))
  } else {
    data.loading <- data.loading %>% 
                    dplyr::mutate("mean_Loading_Score" = ifelse(.data$mean_Loading_Score < round(stats::quantile(.data$mean_Loading_Score, 0.05), 1), round(stats::quantile(.data$mean_Loading_Score, 0.05), 1), .data$mean_Loading_Score))
  }
  
  if (!is.na(max.cutoff.loadings)){
    data.loading <- data.loading %>% 
                    dplyr::mutate("mean_Loading_Score" = ifelse(.data$mean_Loading_Score > max.cutoff.loadings, max.cutoff.loadings, .data$mean_Loading_Score))
  } else {
    data.loading <- data.loading %>% 
                    dplyr::mutate("mean_Loading_Score" = ifelse(.data$mean_Loading_Score > round(stats::quantile(.data$mean_Loading_Score, 0.95), 1), round(stats::quantile(.data$mean_Loading_Score, 0.95), 1), .data$mean_Loading_Score))
  }
  
  
  if (!is.na(min.cutoff.expression)){
    data.expression <- data.expression %>% 
                       dplyr::mutate("mean_Expression" = ifelse(.data$mean_Expression < min.cutoff.loadings, min.cutoff.loadings, .data$mean_Expression))
  }
  
  if (!is.na(max.cutoff.expression)){
    data.expression <- data.expression %>% 
                       dplyr::mutate("mean_Expression" = ifelse(.data$mean_Expression > max.cutoff.loadings, max.cutoff.loadings, .data$mean_Expression))
  } else {
    data.expression <- data.expression %>% 
                       dplyr::mutate("mean_Expression" = ifelse(.data$mean_Expression > round(stats::quantile(.data$mean_Expression, 0.95), 1), round(stats::quantile(.data$mean_Expression, 0.95), 1), .data$mean_Expression))
  }
  
  # Compute scales.
  limits <- c(min(data.loading$mean_Loading_Score, na.rm = TRUE),
              max(data.loading$mean_Loading_Score, na.rm = TRUE))
  
  
  scale.setup <- compute_scales(sample = sample,
                                feature = " ",
                                assay = "SCT",
                                reduction = NULL,
                                slot = "scale.data",
                                number.breaks = number.breaks,
                                min.cutoff = NA,
                                max.cutoff = NA,
                                flavor = "Seurat",
                                enforce_symmetry = TRUE,
                                from_data = TRUE,
                                limits.use = limits)
  
  p.loading <- data.loading %>% 
               ggplot2::ggplot(mapping = ggplot2::aes(x = .data$Gene,
                                                      y = .data$PC,
                                                      fill = .data$mean_Loading_Score)) +
               ggplot2::geom_tile(color = grid.color, linewidth = 0.5) +
               ggplot2::scale_y_discrete(expand = c(0, 0)) +
               ggplot2::scale_x_discrete(expand = c(0, 0),
                                         position = "top") +
               ggplot2::guides(y.sec = guide_axis_label_trans(~paste0(levels(.data$PC))),
                               x.sec = guide_axis_label_trans(~paste0(levels(.data$Gene)))) +
               ggplot2::scale_fill_gradientn(colors = RColorBrewer::brewer.pal(n = 11, name = diverging.palette) %>% rev(),
                                             na.value = na.value,
                                             name = "Avg. Loading score",
                                             breaks = scale.setup$breaks,
                                             labels = scale.setup$labels,
                                             limits = scale.setup$limits) +
               ggplot2::coord_equal() +
               ggplot2::xlab("Top genes") +
               ggplot2::ylab("PCs")
  
  limits <- c(min(data.expression$mean_Expression, na.rm = TRUE),
              max(data.expression$mean_Expression, na.rm = TRUE))
  scale.setup <- compute_scales(sample = sample,
                                         feature = " ",
                                         assay = "SCT",
                                         reduction = NULL,
                                         slot = "scale.data",
                                         number.breaks = number.breaks,
                                         min.cutoff = NA,
                                         max.cutoff = NA,
                                         flavor = "Seurat",
                                         enforce_symmetry = FALSE,
                                         from_data = TRUE,
                                         limits.use = limits)
  
  p.expression <- data.expression %>% 
                  ggplot2::ggplot(mapping = ggplot2::aes(x = .data$Gene,
                                                         y = .data[[group.by]],
                                                         fill = .data$mean_Expression)) +
                  ggplot2::geom_tile(color = grid.color, linewidth = 0.5) +
                  ggplot2::scale_y_discrete(expand = c(0, 0)) +
                  ggplot2::scale_x_discrete(expand = c(0, 0),
                                            position = "top") +
                  ggplot2::guides(y.sec = guide_axis_label_trans(~paste0(levels(.data[[group.by]]))),
                                  x.sec = guide_axis_label_trans(~paste0(levels(.data$Gene)))) +
                  ggplot2::coord_equal() +
                  ggplot2::xlab(NULL) +
                  ggplot2::ylab(group.by)
  
  if (isFALSE(use_viridis)){
    colors <- RColorBrewer::brewer.pal(n = 9, name = "YlGnBu")[2:9]
    if (sequential.direction == -1){
      colors <- rev(colors)
    }
    p.expression <- p.expression + 
                    ggplot2::scale_fill_gradientn(colors = colors,
                                                  na.value = na.value,
                                                  name = "Avg. Expression",
                                                  breaks = scale.setup$breaks,
                                                  labels = scale.setup$labels,
                                                  limits = scale.setup$limits)
  } else {
    p.expression <- p.expression + 
                    ggplot2::scale_fill_viridis_c(option = viridis.palette,
                                                  direction = viridis.direction,
                                                  name = "Avg. Expression",
                                                  breaks = scale.setup$breaks,
                                                  labels = scale.setup$labels,
                                                  limits = scale.setup$limits)
  }
                  
                  
  
  p.loading <- modify_continuous_legend(p = p.loading,
                                        legend.title = "Avg. Loading score",
                                        legend.aes = "fill",
                                        legend.type = legend.type,
                                        legend.position = legend.position,
                                        legend.length = legend.length,
                                        legend.width = legend.width,
                                        legend.framecolor = legend.framecolor,
                                        legend.tickcolor = legend.tickcolor,
                                        legend.framewidth = legend.framewidth,
                                        legend.tickwidth = legend.tickwidth)
  
  p.expression <- modify_continuous_legend(p = p.expression,
                                           legend.title = "Avg. Expression",
                                           legend.aes = "fill",
                                           legend.type = legend.type,
                                           legend.position = legend.position,
                                           legend.length = legend.length,
                                           legend.width = legend.width,
                                           legend.framecolor = legend.framecolor,
                                           legend.tickcolor = legend.tickcolor,
                                           legend.framewidth = legend.framewidth,
                                           legend.tickwidth = legend.tickwidth)
  
  list.plots <- list("Loadings" = p.loading,
                     "Expression" = p.expression)
  counter <- 0
  for (name in rev(names(list.plots))){
    counter <- counter + 1
    
    axis.parameters <- handle_axis(flip = FALSE,
                                   group.by = "A",
                                   group = "A",
                                   counter = counter,
                                   axis.text.x.angle = axis.text.x.angle,
                                   plot.title.face = plot.title.face,
                                   plot.subtitle.face = plot.subtitle.face,
                                   plot.caption.face = plot.caption.face,
                                   axis.title.face = axis.title.face,
                                   axis.text.face = axis.text.face,
                                   legend.title.face = legend.title.face,
                                   legend.text.face = legend.text.face)
    
    list.plots[[name]] <- list.plots[[name]] +
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
                                         plot.margin = ggplot2::margin(t = 0, r = 0, b = 0, l = 0),
                                         panel.border = ggplot2::element_rect(fill = NA, color = border.color, linewidth = 1),
                                         panel.grid.major = ggplot2::element_blank(),
                                         legend.position = legend.position,
                                         plot.background = ggplot2::element_rect(fill = "white", color = "white"),
                                         panel.background = ggplot2::element_rect(fill = "white", color = "white"),
                                         legend.background = ggplot2::element_rect(fill = "white", color = "white"))
  }
  list.plots[["Loadings"]] <- list.plots[["Loadings"]] + 
                              ggplot2::xlab(paste0("Top and bottom ", top_loadings, " genes in PCA loadings")) +
                              ggplot2::theme(axis.title.x.top = ggplot2::element_text(face = "bold", color = "black"))
  
  p <- patchwork::wrap_plots(A = list.plots$Loadings,
                             B = list.plots$Expression,
                             design = "A
                                       B",
                             guides = "collect") + 
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
                                                                                              caption = plot.caption.face,
                                                                                              color = "black",
                                                                                              hjust = 1),
                                                         plot.caption.position = "plot"))
  
  return(p)
}