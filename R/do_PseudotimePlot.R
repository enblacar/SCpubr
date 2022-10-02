#' Compute Pseudotime analysis plots.
#'
#' @inheritParams doc_function
#' @param cds \strong{\code{\link[monocle3]{cell_data_set}}} | Cell Data Set of the same seurat object. Can be obtained using `SeuratWrappers::as.cell_data_set(sample)`.
#' @param group.by \strong{\code{\link[base]{character}}} | Metadata variable to use as grouping if monocle3 clusters are not computed.
#' @param compute_monocle_partitions \strong{\code{\link[base]{logical}}} | Whether to tell monocle3 to compute different partitions. FALSE will treat all the UMAP as a single partition.
#' @param compute_monocle_clusters \strong{\code{\link[base]{logical}}} | Whether to make monocle3 to re-compute clustering
#' @param trajectory_graph_color \strong{\code{\link[base]{character}}} | Color of the trajectory graph plotted on top of the UMAP.
#' @param trajectory_graph_segment_size \strong{\code{\link[base]{numeric}}} | Size of the trajectory graph.
#' @param pseudotime_genes \strong{\code{\link[base]{character}}} | List of genes that will be used to compute enrichment scores, that will be used for pseudotime.
#' @param is_max_score_the_start \strong{\code{\link[base]{logical}}} | Do the cells with the highest enrichment scores depict the beginning of the trajectory (TRUE) or the end (FALSE)?
#' @param label_roots,label_branches,label_leaves \strong{\code{\link[base]{logical}}} | Label roots, branches or leaves in the trajectory graph.
#' @param nbin \strong{\code{\link[base]{numeric}}} | Number of bins to use while computing enrichment scores.
#' @param ctrl \strong{\code{\link[base]{numeric}}} | Number of control genes per bin while computing enrichment scores.
#'
#' @export
#' @return A list containing a collection of ggplot2 objects.
#'
#' @example /man/examples/examples_do_PseudotimePlot.R
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
                              viridis_direction = 1,
                              plot_cell_borders = TRUE,
                              enforce_symmetry = FALSE,
                              nbin = 24,
                              ctrl = 100){
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
                         "plot_cell_borders" = plot_cell_borders,
                         "enforce_symmetry" = enforce_symmetry)
    check_type(parameters = logical_list, required_type = "logical", test_function = is.logical)
    # Check numeric parameters.
    numeric_list <- list("trajectory_graph_segment_size" = trajectory_graph_segment_size,
                         "pt.size" = pt.size,
                         "border.size" = border.size,
                         "font.size" = font.size,
                         "legend.length" = legend.length,
                         "legend.width" = legend.width,
                         "legend.framewidth" = legend.framewidth,
                         "legend.tickwidth" = legend.tickwidth,
                         "nbin" = nbin,
                         "ctrl" = ctrl,
                         "viridis_direction" = viridis_direction)
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

    check_parameters(parameter = font.type, parameter_name = "font.type")
    check_parameters(parameter = legend.type, parameter_name = "legend.type")
    check_parameters(parameter = legend.position, parameter_name = "legend.position")
    check_parameters(parameter = viridis_direction, parameter_name = "viridis_direction")
    check_parameters(parameter = viridis_color_map, parameter_name = "viridis_color_map")

    # Define pipe operator internally.
    `%>%` <- magrittr::`%>%`

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
      sample[["monocle3_clusters"]] <- as.factor(sample$group.by)
      sample[["monocle3_partitions"]] <- 1
      cds@clusters[["UMAP"]]$clusters <- sample$monocle3_clusters
      cds@clusters[["UMAP"]]$partitions <- sample$monocle3_partitions
    } else if (isTRUE(compute_monocle_clusters) & isFALSE(compute_monocle_partitions)){
      cds <- monocle3::cluster_cells(cds)
      sample[["monocle3_partitions"]] <- 1
      cds@clusters[["UMAP"]]$partitions <- sample$monocle3_partitions
      sample[["monocle_3_clusters"]] <- cds@clusters[["UMAP"]]$clusters
    } else if (isFALSE(compute_monocle_clusters) & isTRUE(compute_monocle_partitions)){
      cds <- monocle3::cluster_cells(cds)
      sample[["monocle3_clusters"]] <- factor(sample$group.by)
      cds@clusters[["UMAP"]]$clusters <- sample$monocle3_clusters
      sample[["monocle3_partitions"]] <- cds@clusters[["UMAP"]]$partitions
    } else if (isTRUE(compute_monocle_clusters) & isTRUE(compute_monocle_partitions)){
      cds <- monocle3::cluster_cells(cds)
      sample[["monocle_3_clusters"]] <- cds@clusters[["UMAP"]]$clusters
      sample[["monocle3_partitions"]] <- cds@clusters[["UMAP"]]$partitions
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

    p.out <- do_DimPlot(sample,
                        group.by = "monocle3_partitions",
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
    sample <- compute_enrichment_scores(sample = sample, input_gene_list = list("Enrichment" = pseudotime_genes), verbose = F, nbin = nbin, ctrl = ctrl)

    # Order cells based on the cell with the highest value of the enrichment scores.
    if (isTRUE(is_max_score_the_start)){
      cell_use <- sample@meta.data %>%
                  dplyr::select("Enrichment", "group.by") %>%
                  tibble::rownames_to_column(var = "cells") %>%
                  dplyr::filter(.data$group.by == {sample@meta.data %>%
                                                   dplyr::select("monocle3_partitions", "Enrichment", "group.by") %>%
                                                   dplyr::group_by(.data$group.by) %>%
                                                   dplyr::summarise(mean = mean(.data$Enrichment)) %>%
                                                   dplyr::arrange(dplyr::desc(.data$mean)) %>%
                                                   dplyr::slice_head(n = 1) %>%
                                                   dplyr::pull("group.by") %>%
                                                   as.character()}) %>%
                  dplyr::arrange(dplyr::desc(.data$Enrichment)) %>%
                  dplyr::slice_head(n = 1) %>%
                  dplyr::pull(.data$cells)
    } else if (isFALSE(is_max_score_the_start)){
      cell_use <- sample@meta.data %>%
                  dplyr::select("Enrichment", "group.by") %>%
                  tibble::rownames_to_column(var = "cells") %>%
                  dplyr::filter(.data$group.by == {sample@meta.data %>%
                                                   dplyr::select("monocle3_partitions", "Enrichment", "group.by") %>%
                                                   dplyr::group_by(.data$group.by) %>%
                                                   dplyr::summarise(mean = mean(.data$Enrichment)) %>%
                                                   dplyr::arrange(dplyr::desc(.data$mean)) %>%
                                                   dplyr::slice_tail(n = 1) %>%
                                                   dplyr::pull("group.by") %>%
                                                   as.character()}) %>%
                  dplyr::arrange(dplyr::desc(.data$Enrichment)) %>%
                  dplyr::slice_tail(n = 1) %>%
                  dplyr::pull(.data$cells)
    }
    cds.ordered <- monocle3::order_cells(cds, root_cells = cell_use)

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
                                       legend.tickcolor = legend.tickcolor,
                                       viridis_direction = viridis_direction)
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
                                   enforce_symmetry = enforce_symmetry,
                                   legend.length = legend.length,
                                   legend.width = legend.width,
                                   legend.framewidth = legend.framewidth,
                                   legend.tickwidth = legend.tickwidth,
                                   legend.framecolor = legend.framecolor,
                                   legend.tickcolor = legend.tickcolor,
                                   viridis_direction = viridis_direction)
    list.out[["enrichment"]] <- p.enrichment
  }))
  return(list.out)
}

