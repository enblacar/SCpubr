
#' Compute Pseudotime analysis plots.
#'
#' @param sample Seurat object.
#' @param cds Cell Data Set of the same seurat object. Can be obtained using `SeuratWrappers::as.cell_data_set(sample)`.
#' @param compute_monocle_partitions Whether to tell monocle3 to compute different partitions. FALSE will treat all the UMAP as a single partition.
#' @param compute_monocle_clusters Whether to make monocle3 to re-compute clustering
#' @param trajectory_graph_color Character. Color of the trajectory graph plotted on top of the UMAP.
#' @param trajectory_graph_segment_size Integer. Size of the trajectory graph.
#' @param pseudotime_genes Character. List of genes that will be used to compute enrichment scores, that will be used for pseudotime.
#' @param is_max_score_the_start Logical. Do the cells with the highest enrichment scores depict the beginning of the trajectory (TRUE) or the end (FALSE)?
#' @param label_roots,label_branches,label_leaves Logical. Label roots, branches or leaves in the trajectory graph.
#' @param pt.size Numeric. Point size for the plots.
#' @param border.size Numeric. Point size for the border of the cells in the plots.
#' @param rotate_x_axis_labels Logical. Whether to rotate X axis labels in the
#' @param font.size Numeric. Overall fontsize for the plots.
#' @param font.type Character. Base font for the plot. One of mono, serif or sans.
#' @param legend.type Character. Type of legend to display. One of: normal, colorbar, colorsteps.
#' @param legend.position Position of the legend in the plot. Will only work if legend is set to TRUE.
#' @param legend.framewidth,legend.tickwidth Width of the lines of the box in the legend.
#' @param legend.framecolor,legend.tickcolor Color of the lines of the box in the legend.
#' @param legend.length,legend.width Length and width of the legend. Will adjust automatically depending on legend side.
#' @param viridis_color_map Character. Viridis color map to use in the FeaturePlot and dotplot for the enrichment scores.
#'
#' @export
#' @return A list containing a collection of ggplot2 objects.
#'
#' @examples
#' \dontrun{
#' TBD
#' }
do_PseudotimePlot <- function(sample,
                              cds,
                              compute_monocle_partitions = TRUE,
                              compute_monocle_clusters = FALSE,
                              trajectory_graph_color = "black",
                              trajectory_graph_segment_size = 1,
                              pseudotime_genes = NULL,
                              is_max_score_the_start = TRUE,
                              label_roots = FALSE,
                              label_branches = FALSE,
                              label_leaves = FALSE,
                              pt.size = 1,
                              border.size = 1.5,
                              rotate_x_axis_labels = F,
                              legend.position = "bottom",
                              legend.type = "colorbar",
                              font.size = 14,
                              font.type = "sans",
                              legend.length = 30,
                              legend.width = 1,
                              legend.framewidth = 1.5,
                              legend.tickwidth = 1.5,
                              legend.framecolor = "grey50",
                              legend.tickcolor = "white",
                              viridis_color_map = "D"){

  # Checks for packages.
  check_suggests(function_name = "do_PseudotimePlot")
  # Check if the sample provided is a Seurat object.
  check_Seurat(sample = sample)

  # Check logical parameters.
  logical_list <- list("compute_monocle_partitions" = compute_monocle_partitions,
                       "compute_monocle_clusters" = compute_monocle_clusters,
                       "is_max_score_the_start" = is_max_score_the_start,
                       "label_roots" = label_roots,
                       "label_branches" = label_branches,
                       "label_leaves" = label_leaves)
  check_type(parameters = logical_list, required_type = "logical", test_function = is.logical)
  # Check numeric parameters.
  numeric_list <- list("trajectory_graph_segment_size" = trajectory_graph_segment_size,
                       "pt.size" = pt.size,
                       "border.size" = border.size,
                       "font.size" = font.size,
                       "legend.length" = legend.length,
                       "legend.width" = legend.width,
                       "legend.framewidth" = legend.framewidth,
                       "legend.tickwidth" = legend.tickwidth)
  check_type(parameters = numeric_list, required_type = "numeric", test_function = is.numeric)
  # Check character parameters.
  character_list <- list("trajectory_graph_color" = trajectory_graph_color,
                         "pseudotime_genes" = pseudotime_genes,
                         "legend.position" = legend.position,
                         "legend.type" = legend.type,
                         "font.type" = font.type,
                         "legend.framecolor" = legend.framecolor,
                         "legend.tickcolor" = legend.tickcolor,
                         "viridis_color_map" = viridis_color_map)
  check_type(parameters = character_list, required_type = "character", test_function = is.character)

  check_colors(legend.framecolor, parameter_name = "legend.framecolor")
  check_colors(legend.tickcolor, parameter_name = "legend.tickcolor")
  check_colors(trajectory_graph_color, parameter_name = "trajectory_graph_color")


  # Define pipe operator internally.
  `%>%` <- purrr::`%>%`

  # Whether to tweak or not the output partitions and clusters.
  if (isFALSE(compute_monocle_clusters) & isFALSE(compute_monocle_partitions)){
    # Tweak to modify the partitions.
    sample[["monocle3_clusters"]] <- sample$seurat_clusters
    sample[["monocle3_partitions"]] <- 1
  } else if (isTRUE(compute_monocle_clusters) & isFALSE(compute_monocle_partitions)){
    cds <- monocle3::cluster_cells(cds)
    sample[["monocle3_partitions"]] <- 1
  } else if (isFALSE(compute_monocle_clusters) & isTRUE(compute_monocle_partitions)){
    cds <- monocle3::cluster_cells(cds)
    sample[["monocle3_clusters"]] <- sample$seurat_clusters
  } else if (isTRUE(compute_monocle_clusters) & isTRUE(compute_monocle_partitions)){
    cds <- monocle3::cluster_cells(cds)
  }

  cds <- suppressWarnings(monocle3::learn_graph(cds))

  list.out <- list()
  p <- monocle3::plot_cells(cds,
                            color_cells_by = "partition",
                            label_groups_by_cluster = FALSE,
                            label_leaves = label_leaves,
                            label_branch_points = label_branches,
                            label_cell_groups = FALSE,
                            label_roots = label_roots,
                            trajectory_graph_color = trajectory_graph_color,
                            trajectory_graph_segment_size = trajectory_graph_segment_size,
                            cell_size = pt.size)
  build <- ggplot2::ggplot_build(p)
  sample$partitions <-  build$plot$data$cell_color


  p.umap <- do_DimPlot(sample,
                       plot_cell_borders = TRUE,
                       border.size = border.size,
                       pt.size = pt.size)
  list.out[["umap"]] <- p.umap

  p.out <- do_DimPlot(sample,
                      group.by = "partitions",
                      plot_cell_borders = TRUE,
                      border.size = border.size,
                      pt.size = pt.size)
  p$layers[c(1, 2)] <- NULL
  p.out$layers <- append(p.out$layers, p$layers)

  list.out[["trajectory"]] <- p.out

  # Compute enrichment scores for the desired marker genes.
  sample <- compute_enrichment_scores(sample = sample, list_genes = list("Enrichment" = pseudotime_genes), verbose = F)

  # Order cells based on the cell with the highest value of the enrichment scores.
  if (isTRUE(is_max_score_the_start)){
    root_cells <- colnames(sample)[which.max(sample$Enrichment)]
  } else if (isFALSE(is_max_score_the_start)){
    root_cells <- colnames(sample)[which.min(sample$Enrichment)]
  }
  cds.ordered <- monocle3::order_cells(cds, root_cells = root_cells)

  #Plot Pseudotime.
  p.pseudotime <- monocle3::plot_cells(cds.ordered,
                                       color_cells_by = "pseudotime",
                                       group_cells_by = "partition",
                                       label_groups_by_cluster = FALSE,
                                       label_cell_groups = FALSE,
                                       label_leaves = label_leaves,
                                       label_branch_points = label_branches,
                                       label_roots = label_roots,
                                       trajectory_graph_color = trajectory_graph_color,
                                       trajectory_graph_segment_size = trajectory_graph_segment_size,
                                       cell_size = pt.size)
  build <- ggplot2::ggplot_build(p.pseudotime)
  sample$pseudotime <- build$plot$data$cell_color
  p.pseudotime.out <- do_FeaturePlot(sample = sample,
                                     features = "pseudotime",
                                     plot_cell_borders = TRUE,
                                     border.size = border.size,
                                     pt.size = pt.size,
                                     viridis_color_map = "C")
  p.pseudotime.out$layers <- append(p.pseudotime.out$layers, p.pseudotime$layers[seq(3, length(p.pseudotime$layers))])

  list.out[["pseudotime"]] <- p.pseudotime.out

  p.enrichment <- do_FeaturePlot(sample = sample,
                                 features = "Enrichment",
                                 pt.size = pt.size,
                                 border.size = border.size,
                                 plot_cell_borders = TRUE,
                                 viridis_color_map = viridis_color_map)
  list.out[["enrichment"]] <- p.enrichment


  # Define legend parameters. Width and height values will change depending on the legend orientation.
  if (legend.position %in% c("top", "bottom")){
    legend.barwidth <- legend.length
    legend.barheight <- legend.width
  } else if (legend.position %in% c("left", "right")){
    legend.barwidth <- legend.width
    legend.barheight <- legend.length
  }

  p.dot.enrichment <- sample@meta.data %>%
                      dplyr::select("seurat_clusters", "Enrichment", "pseudotime") %>%
                      dplyr::filter(.data$pseudotime != Inf) %>%
                      dplyr::mutate("seurat_clusters" = factor(.data$seurat_clusters, levels = {sample@meta.data %>%
                                                                                                dplyr::select("seurat_clusters", "Enrichment", "pseudotime") %>%
                                                                                                dplyr::filter(.data$pseudotime != Inf) %>%
                                                                                                dplyr::group_by(.data$seurat_clusters) %>%
                                                                                                dplyr::summarise("mean" = mean(.data$Enrichment)) %>%
                                                                                                dplyr::arrange(if (isTRUE(is_max_score_the_start)){dplyr::desc(.data$mean)} else {.data$mean}) %>%
                                                                                                dplyr::pull(.data$seurat_clusters)})) %>%
                      ggplot2::ggplot(mapping = ggplot2::aes(x = .data$seurat_clusters,
                                                             y = .data$Enrichment,
                                                             color = .data$Enrichment)) +
                      ggplot2::geom_point(position = ggplot2::position_jitter(width = 0.455,
                                                                              seed = 0),
                                          size = pt.size * border.size,
                                          color = "black",
                                          na.rm = TRUE) +
                      ggplot2::geom_point(position = ggplot2::position_jitter(width = 0.45,
                                                                              seed = 0),size = pt.size,
                                          na.rm = TRUE) +
                      ggdist::stat_pointinterval(position = ggplot2::position_dodge(width = 1),
                                                 na.rm = TRUE,
                                                 color = "black") +
                      ggplot2::scale_color_viridis_c(option = viridis_color_map) +
                      ggplot2::theme_minimal(base_size = font.size) +
                      ggplot2::theme(axis.title = ggplot2::element_blank(),
                                     axis.line.x = ggplot2::element_line(color = "black"),
                                     axis.text.x = ggplot2::element_text(color = "black",
                                                                         face = "bold",
                                                                         angle = ifelse(isTRUE(rotate_x_axis_labels), 90, 0),
                                                                         hjust = ifelse(isTRUE(rotate_x_axis_labels), 1, 0.5),
                                                                         vjust = ifelse(isTRUE(rotate_x_axis_labels), 0.5, 1)),
                                     axis.text.y = ggplot2::element_text(color = "black", face = "bold", hjust = 0),
                                     axis.ticks = ggplot2::element_line(color = "black"),
                                     panel.grid.major = ggplot2::element_blank(),
                                     plot.title.position = "plot",
                                     plot.title = ggtext::element_markdown(face = "bold", hjust = 0),
                                     plot.subtitle = ggtext::element_markdown(hjust = 0),
                                     plot.caption = ggtext::element_markdown(hjust = 1),
                                     panel.grid = ggplot2::element_blank(),
                                     text = ggplot2::element_text(family = font.type),
                                     plot.caption.position = "plot",
                                     legend.text = ggplot2::element_text(face = "bold"),
                                     legend.position = legend.position,
                                     legend.title = ggplot2::element_text(face = "bold"),
                                     legend.justification = "center",
                                     plot.margin = ggplot2::margin(t = 10, r = 40, b = 10, l = 40),
                                     plot.background = ggplot2::element_rect(fill = "white", color = "white"),
                                     panel.background = ggplot2::element_rect(fill = "white", color = "white"),
                                     legend.background = ggplot2::element_rect(fill = "white", color = "white"))

  p.dot.enrichment <- modify_continuous_legend(p = p.dot.enrichment,
                                               legend.aes = "color",
                                               legend.type = legend.type,
                                               legend.position = legend.position,
                                               legend.length = legend.length,
                                               legend.width = legend.width,
                                               legend.framecolor = legend.framecolor,
                                               legend.tickcolor = legend.tickcolor,
                                               legend.framewidth = legend.framewidth,
                                               legend.tickwidth = legend.tickwidth)
  list.out[["enrichment_dotplot"]] <- p.dot.enrichment

  p.dot.pseudotime <- sample@meta.data %>%
                      dplyr::select("seurat_clusters", "pseudotime") %>%
                      dplyr::filter(.data$pseudotime != Inf) %>%
                      dplyr::mutate("seurat_clusters" = factor(.data$seurat_clusters, levels = {sample@meta.data %>%
                                                                                                dplyr::select("seurat_clusters", "pseudotime") %>%
                                                                                                dplyr::filter(.data$pseudotime != Inf) %>%
                                                                                                dplyr::group_by(.data$seurat_clusters) %>%
                                                                                                dplyr::summarise("mean" = mean(.data$pseudotime)) %>%
                                                                                                dplyr::arrange(.data$mean) %>%
                                                                                                dplyr::pull(.data$seurat_clusters)})) %>%
                      ggplot2::ggplot(mapping = ggplot2::aes(x = .data$seurat_clusters,
                                                             y = .data$pseudotime,
                                                             color = .data$pseudotime)) +
                      ggplot2::geom_point(position = ggplot2::position_jitter(width = 0.455,
                                                                              seed = 0),
                                          size = pt.size * border.size,
                                          color = "black",
                                          na.rm = TRUE) +
                      ggplot2::geom_point(position = ggplot2::position_jitter(width = 0.45,
                                                                              seed = 0),size = pt.size,
                                          na.rm = TRUE) +
                      ggdist::stat_pointinterval(position = ggplot2::position_dodge(width = 1),
                                                 na.rm = TRUE,
                                                 color = "black") +
                      ggplot2::scale_color_viridis_c(option = "C") +
                      ggplot2::theme_minimal(base_size = font.size) +
                      ggplot2::theme(axis.title = ggplot2::element_blank(),
                                     axis.line.x = ggplot2::element_line(color = "black"),
                                     axis.text.x = ggplot2::element_text(color = "black",
                                                                         face = "bold",
                                                                         angle = ifelse(isTRUE(rotate_x_axis_labels), 90, 0),
                                                                         hjust = ifelse(isTRUE(rotate_x_axis_labels), 1, 0.5),
                                                                         vjust = ifelse(isTRUE(rotate_x_axis_labels), 0.5, 1)),
                                     axis.text.y = ggplot2::element_text(color = "black", face = "bold", hjust = 0),
                                     axis.ticks = ggplot2::element_line(color = "black"),
                                     panel.grid.major = ggplot2::element_blank(),
                                     plot.title.position = "plot",
                                     plot.title = ggtext::element_markdown(face = "bold", hjust = 0),
                                     plot.subtitle = ggtext::element_markdown(hjust = 0),
                                     plot.caption = ggtext::element_markdown(hjust = 1),
                                     panel.grid = ggplot2::element_blank(),
                                     text = ggplot2::element_text(family = font.type),
                                     plot.caption.position = "plot",
                                     legend.text = ggplot2::element_text(face = "bold"),
                                     legend.position = legend.position,
                                     legend.title = ggplot2::element_text(face = "bold"),
                                     legend.justification = "center",
                                     plot.margin = ggplot2::margin(t = 10, r = 40, b = 10, l = 40),
                                     plot.background = ggplot2::element_rect(fill = "white", color = "white"),
                                     panel.background = ggplot2::element_rect(fill = "white", color = "white"),
                                     legend.background = ggplot2::element_rect(fill = "white", color = "white"))




  p.dot.pseudotime <- modify_continuous_legend(p = p.dot.pseudotime,
                                               legend.aes = "color",
                                               legend.type = legend.type,
                                               legend.position = legend.position,
                                               legend.length = legend.length,
                                               legend.width = legend.width,
                                               legend.framecolor = legend.framecolor,
                                               legend.tickcolor = legend.tickcolor,
                                               legend.framewidth = legend.framewidth,
                                               legend.tickwidth = legend.tickwidth)

  list.out[["pseudotime_dotplot"]] <- p.dot.pseudotime

  layout <- "ABE
             CDF"
  patch <- patchwork::wrap_plots(A = list.out$umap,
                                 B = list.out$trajectory,
                                 C = list.out$enrichment,
                                 D = list.out$pseudotime,
                                 E = list.out$enrichment_dotplot,
                                 F = list.out$pseudotime_dotplot,
                                 design = layout)
  list.out[["combined_analysis"]] <- patch

  return(list.out)
}

