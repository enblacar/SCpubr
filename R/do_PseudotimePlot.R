
#' Compute Pseudotime analysis plots.
#'
#' @param sample Seurat object.
#' @param cds Cell Data Set of the same seurat object. Can be obtained using `SeuratWrappers::as.cell_data_set(sample)`.
#' @param group.by Character. Metadata variable to use as grouping if monocle3 clusters are not computed.
#' @param compute_monocle_partitions Whether to tell monocle3 to compute different partitions. FALSE will treat all the UMAP as a single partition.
#' @param compute_monocle_clusters Whether to make monocle3 to re-compute clustering
#' @param trajectory_graph_color Character. Color of the trajectory graph plotted on top of the UMAP.
#' @param trajectory_graph_segment_size Integer. Size of the trajectory graph.
#' @param pseudotime_genes Character. List of genes that will be used to compute enrichment scores, that will be used for pseudotime.
#' @param is_max_score_the_start Logical. Do the cells with the highest enrichment scores depict the beginning of the trajectory (TRUE) or the end (FALSE)?
#' @param label_roots,label_branches,label_leaves Logical. Label roots, branches or leaves in the trajectory graph.
#' @param pt.size Numeric. Point size for the plots.
#' @param border.size Numeric. Point size for the border of the cells in the plots.
#' @param font.size Numeric. Overall fontsize for the plots.
#' @param font.type Character. Base font for the plot. One of mono, serif or sans.
#' @param legend.type Character. Type of legend to display. One of: normal, colorbar, colorsteps.
#' @param legend.position Position of the legend in the plot. Will only work if legend is set to TRUE.
#' @param legend.framewidth,legend.tickwidth Width of the lines of the box in the legend.
#' @param legend.framecolor,legend.tickcolor Color of the lines of the box in the legend.
#' @param legend.length,legend.width Length and width of the legend. Will adjust automatically depending on legend side.
#' @param viridis_color_map Character. Viridis color map to use in the FeaturePlot and dotplot for the enrichment scores.
#' @param plot_cell_borders Logical. Whether to plot borders around the cells.
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
                              pseudotime_genes,
                              group.by = NULL,
                              compute_monocle_partitions = TRUE,
                              compute_monocle_clusters = FALSE,
                              trajectory_graph_color = "black",
                              trajectory_graph_segment_size = 1,
                              is_max_score_the_start = TRUE,
                              label_roots = FALSE,
                              label_branches = FALSE,
                              label_leaves = FALSE,
                              pt.size = 1,
                              border.size = 1.5,
                              legend.position = "bottom",
                              legend.type = "colorbar",
                              font.size = 14,
                              font.type = "sans",
                              legend.length = 20,
                              legend.width = 1,
                              legend.framewidth = 1.5,
                              legend.tickwidth = 1.5,
                              legend.framecolor = "grey50",
                              legend.tickcolor = "white",
                              viridis_color_map = "D",
                              plot_cell_borders = TRUE){
  sink(tempfile())
  on.exit(sink())
  invisible(force({
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
                         "label_leaves" = label_leaves,
                         "plot_cell_borders" = plot_cell_borders)
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
                           "viridis_color_map" = viridis_color_map,
                           "group.by" = group.by)
    check_type(parameters = character_list, required_type = "character", test_function = is.character)

    check_colors(legend.framecolor, parameter_name = "legend.framecolor")
    check_colors(legend.tickcolor, parameter_name = "legend.tickcolor")
    check_colors(trajectory_graph_color, parameter_name = "trajectory_graph_color")


    # Define pipe operator internally.
    `%>%` <- purrr::`%>%`

    if (is.null(group.by)){
      sample$group.by <- Seurat::Idents(sample)
    } else {
      if (!(group.by %in% colnames(sample@meta.data))){
        stop("group.by needs to be in the metadata.", call. = FALSE)
      }
      sample$group.by <- sample@meta.data[, group.by]
    }
    # Whether to tweak or not the output partitions and clusters.
    if (isFALSE(compute_monocle_clusters) & isFALSE(compute_monocle_partitions)){
      # Tweak to modify the partitions.
      sample[["monocle3_clusters"]] <- sample$group.by
      sample[["monocle3_partitions"]] <- 1
    } else if (isTRUE(compute_monocle_clusters) & isFALSE(compute_monocle_partitions)){
      cds <- monocle3::cluster_cells(cds)
      sample[["monocle3_partitions"]] <- 1
    } else if (isFALSE(compute_monocle_clusters) & isTRUE(compute_monocle_partitions)){
      cds <- monocle3::cluster_cells(cds)
      sample[["monocle3_clusters"]] <- sample$group.by
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

    p.out <- do_DimPlot(sample,
                        group.by = "partitions",
                        plot_cell_borders = plot_cell_borders,
                        border.size = border.size,
                        pt.size = pt.size,
                        font.size = font.size,
                        font.type = font.type)
    p$layers[c(1, 2)] <- NULL
    p.out$layers <- append(p.out$layers, p$layers)

    list.out[["trajectory_partitions"]] <- p.out

    p.out <- do_DimPlot(sample,
                        group.by = group.by,
                        plot_cell_borders = plot_cell_borders,
                        border.size = border.size,
                        pt.size = pt.size,
                        font.size = font.size,
                        font.type = font.type)
    p$layers[c(1, 2)] <- NULL
    p.out$layers <- append(p.out$layers, p$layers)
    list.out[["trajectory_groups"]] <- p.out

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
                                       plot_cell_borders = plot_cell_borders,
                                       border.size = border.size,
                                       pt.size = pt.size,
                                       viridis_color_map = "C",
                                       font.size = font.size,
                                       font.type = font.type,
                                       legend.length = legend.length,
                                       legend.width = legend.width,
                                       legend.framewidth = legend.framewidth,
                                       legend.tickwidth = legend.tickwidth,
                                       legend.framecolor = legend.framecolor,
                                       legend.tickcolor = legend.tickcolor)
    p.pseudotime.out$layers <- append(p.pseudotime.out$layers, p.pseudotime$layers[seq(3, length(p.pseudotime$layers))])

    list.out[["pseudotime"]] <- p.pseudotime.out

    p.enrichment <- do_FeaturePlot(sample = sample,
                                   features = "Enrichment",
                                   pt.size = pt.size,
                                   border.size = border.size,
                                   plot_cell_borders = plot_cell_borders,
                                   viridis_color_map = viridis_color_map,
                                   font.size = font.size,
                                   font.type = font.type,
                                   legend.length = legend.length,
                                   legend.width = legend.width,
                                   legend.framewidth = legend.framewidth,
                                   legend.tickwidth = legend.tickwidth,
                                   legend.framecolor = legend.framecolor,
                                   legend.tickcolor = legend.tickcolor)
    list.out[["enrichment"]] <- p.enrichment
  }))
  return(list.out)
}

