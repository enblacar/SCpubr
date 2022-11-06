#' Create enrichment scores heatmaps.
#'
#' This function computes the enrichment scores for the cells using \link[Seurat]{AddModuleScore} and then aggregates the scores by the metadata variables provided by the user and displays it as a heatmap, computed by \link[ComplexHeatmap]{Heatmap}.
#'
#' @inheritParams doc_function
#' @param enforce_symmetry \strong{\code{\link[base]{logical}}} | Whether the geyser and feature plot has a symmetrical color scale.
#' @param flavor \strong{\code{\link[base]{character}}} | One of: Seurat, UCell. Compute the enrichment scores using \link[Seurat]{AddModuleScore} or \link[UCell]{AddModuleScore_UCell}.
#' @param ncores \strong{\code{\link[base]{numeric}}} | Number of cores used to run UCell scoring.
#' @param storeRanks \strong{\code{\link[base]{logical}}} | Whether to store the ranks for faster UCell scoring computations. Might require large amounts of RAM.
#' @param geyser_order_by_mean \strong{\code{\link[base]{logical}}} | Whether to order the X axis by the mean of the values.
#' @param geyser_scale_type \strong{\code{\link[base]{character}}} | Type of scale to use. Either "continuous" or "categorical.
#' @param plot_FeaturePlots \strong{\code{\link[base]{logical}}} | Compute output FeaturePlots for each of the gene lists.
#' @param plot_GeyserPlots \strong{\code{\link[base]{logical}}} | Compute output GeyserPlots for each of the gene lists and group.by variable.
#' @return A ComplexHeatmap object.
#' @export
#'
#' @example /man/examples/examples_do_EnrichmentHeatmap.R
do_EnrichmentHeatmap <- function(sample,
                                 input_gene_list,
                                 assay = NULL,
                                 slot = NULL,
                                 reduction = NULL,
                                 group.by = NULL,
                                 verbose = FALSE,
                                 flip = FALSE,
                                 cluster_cols = FALSE,
                                 cluster_rows = FALSE,
                                 legend.title = "Enrichment",
                                 row_names_rot = 0,
                                 column_names_rot = 45,
                                 cell_size = 5,
                                 na.value = "grey75",
                                 legend.position = "bottom",
                                 use_viridis = TRUE,
                                 viridis_color_map = "G",
                                 viridis_direction = -1,
                                 heatmap.legend.length = 75,
                                 heatmap.legend.width = 5,
                                 heatmap.legend.framecolor = "black",
                                 legend.width = 1,
                                 legend.length = 20,
                                 legend.framewidth = 1.5,
                                 legend.tickwidth = 1.5,
                                 legend.framecolor = "grey50",
                                 legend.tickcolor = "white",
                                 legend.type = "colorbar",
                                 font.size = 14,
                                 font.type = "sans",
                                 rotate_x_axis_labels = 45,
                                 enforce_symmetry = FALSE,
                                 heatmap_gap = 0.5,
                                 row_names_side = "right",
                                 row_title_side = "left",
                                 row_title_rot = 90,
                                 column_title = if (isFALSE(flip)){"List of Genes"} else {"Groups"},
                                 row_title = if (isFALSE(flip)){"Groups"} else {"List of Genes"},
                                 nbin = 24,
                                 ctrl = 100,
                                 flavor = "Seurat",
                                 ncores = 1,
                                 storeRanks = TRUE,
                                 min.cutoff = NULL,
                                 max.cutoff = NULL,
                                 plot_FeaturePlots = FALSE,
                                 plot_GeyserPlots = FALSE,
                                 pt.size = 1,
                                 plot_cell_borders = TRUE,
                                 border.size = 2,
                                 geyser_order_by_mean = TRUE,
                                 geyser_scale_type = "continuous"){
  check_suggests(function_name = "do_EnrichmentHeatmap")
  # Check if the sample provided is a Seurat object.
  check_Seurat(sample = sample)

  # Check logical parameters.
  logical_list <- list("cluster_cols" = cluster_cols,
                       "cluster_rows" = cluster_rows,
                       "use_viridis" = use_viridis,
                       "enforce_symmetry" = enforce_symmetry,
                       "plot_FeaturePlots" = plot_FeaturePlots,
                       "plot_GeyserPlots" = plot_GeyserPlots,
                       "plot_cell_borders" = plot_cell_borders,
                       "geyser_order_by_mean" = geyser_order_by_mean)
  check_type(parameters = logical_list, required_type = "logical", test_function = is.logical)
  # Check numeric parameters.
  numeric_list <- list("row_names_rot" = row_names_rot,
                       "column_names_rot" = column_names_rot,
                       "cell_size" = cell_size,
                       "viridis_direction" = viridis_direction,
                       "row_title_rot" = row_title_rot,
                       "nbin" = nbin,
                       "ctrl" = ctrl,
                       "ncores" = ncores,
                       "pt.size" = pt.size,
                       "border.size" = border.size,
                       "font.size" = font.size,
                       "legend.width" = legend.width,
                       "legend.length" = legend.length,
                       "legend.framewidth" = legend.framewidth,
                       "legend.tickwidth" = legend.tickwidth,
                       "viridis_direction" = viridis_direction,
                       "rotate_x_axis_labels" = rotate_x_axis_labels,
                       "min.cutoff" = min.cutoff,
                       "max.cutoff" = max.cutoff)
  check_type(parameters = numeric_list, required_type = "numeric", test_function = is.numeric)
  # Check character parameters.
  character_list <- list("input_gene_list" = input_gene_list,
                         "column_title" = column_title,
                         "row_title" = row_title,
                         "legend.title" = legend.title,
                         "legend.position" = legend.position,
                         "legend.framecolor" = legend.framecolor,
                         "font.type" = font.type,
                         "group.by" = group.by,
                         "na.value" = na.value,
                         "legend.position" = legend.position,
                         "viridis_color_map" = viridis_color_map,
                         "row_names_side" = row_names_side,
                         "row_title_side" = row_title_side,
                         "flavor" = flavor,
                         "geyser_scale_type" = geyser_scale_type)
  check_type(parameters = character_list, required_type = "character", test_function = is.character)
  check_colors(na.value)
  check_colors(heatmap.legend.framecolor, parameter_name = "heatmap.legend.framecolor")

  check_parameters(parameter = viridis_direction, parameter_name = "viridis_direction")
  check_parameters(parameter = legend.position, parameter_name = "legend.position")
  check_parameters(parameter = flavor, parameter_name = "flavor")


  `%v%` <- ComplexHeatmap::`%v%`
  `%>%` <- magrittr::`%>%`

  if (!(is.null(assay)) & flavor == "UCell"){
    stop("If setting flavor = UCell, do not use assay parameter. Instead, make sure the assay you want is set as default.", .call = FALSE)
  }

  if (!(is.null(slot)) & flavor == "Seurat"){
    stop("If setting flavor = Seurat, do not use slot parameter. The slot is determined by default in Seurat.", .call = FALSE)
  }

  if (is.character(input_gene_list)){
    # If input_gene_list is a character of genes.
    input_list <- list("Input" = input_gene_list)
  } else if (is.list(input_gene_list)){
    input_list <- input_gene_list
    assertthat::assert_that(!is.null(names(input_list)),
                            msg = "Please provide a named list. This is, each gene list has to come with a name.")
  }

  # Compute the enrichment scores.
  sample <- compute_enrichment_scores(sample = sample,
                                      input_gene_list = input_gene_list,
                                      verbose = verbose,
                                      nbin = nbin,
                                      ctrl = ctrl,
                                      flavor = flavor,
                                      ncores = ncores,
                                      storeRanks = storeRanks,
                                      assay = assay,
                                      slot = slot)

  list.heatmaps <- list()
  list.legends <- list()
  if (is.null(group.by)){
    aggr_entities <- levels(sample)
    sample@meta.data[, "Groups"] <- sample@active.ident
    group.by <- "Groups"
  }

  for (g in group.by){
    assertthat::assert_that(g %in% colnames(sample@meta.data),
                            msg = paste0("Value ", g, " in group.by is not in the sample metadata."))

    assertthat::assert_that(class(sample@meta.data[, g]) %in% c("character", "factor"),
                            msg = paste0("Value ", g, " in group.by is not a character or factor column in the metadata."))
  }


  max_value_list <- c()
  min_value_list <- c()
  # Compute the max values for all heatmaps.
  for (variant in group.by){
    max_value <- sample@meta.data %>%
                 dplyr::select(dplyr::all_of(c(variant, names(input_list)))) %>%
                 dplyr::group_by(.data[[variant]]) %>%
                 dplyr::summarize(dplyr::across(.cols = names(input_list), mean)) %>%
                 tibble::column_to_rownames(var = variant) %>%
                 as.matrix() %>%
                 max()

    min_value <- sample@meta.data %>%
                 dplyr::select(dplyr::all_of(c(variant, names(input_list)))) %>%
                 dplyr::group_by(.data[[variant]]) %>%
                 dplyr::summarize(dplyr::across(.cols = names(input_list), mean)) %>%
                 tibble::column_to_rownames(var = variant) %>%
                 as.matrix() %>%
                 min()

    max_value_list <- c(max_value_list, max_value)
    min_value_list <- c(min_value_list, min_value)
  }
  range.data <- c(min(min_value_list), max(max_value_list))
  if (!is.null(min.cutoff) & !is.null(max.cutoff)){
    assertthat::assert_that(min.cutoff < max.cutoff,
                            msg = paste0("The value provided for min.cutoff (", min.cutoff, ") has to be lower than the value provided to max.cutoff (", max.cutoff, "). Please select another value."))

    assertthat::assert_that(max.cutoff > min.cutoff,
                            msg = paste0("The value provided for max.cutoff (", max.cutoff, ") has to be higher than the value provided to min.cutoff (", min.cutoff, "). Please select another value."))

    assertthat::assert_that(max.cutoff != min.cutoff,
                            msg = paste0("The value provided for max.cutoff (", max.cutoff, ") can not be the same than the value provided to min.cutoff (", min.cutoff, "). Please select another value."))

  }

  if (!is.null(min.cutoff)){
    assertthat::assert_that(min.cutoff >= range.data[1],
                            msg = paste0("The value provided for min.cutoff (", min.cutoff, ") is lower than the minimum value in the enrichment matrix (", range.data[1], "). Please select another value."))
    range.data <- c(min.cutoff, range.data[2])
  }

  if (!is.null(max.cutoff)){
    assertthat::assert_that(max.cutoff <= range.data[2],
                            msg = paste0("The value provided for max.cutoff (", max.cutoff, ") is lower than the maximum value in the enrichment matrix (", range.data[2], "). Please select another value."))
    range.data <- c(range.data[1], max.cutoff)

  }


  for (variant in group.by){
    data <- sample@meta.data %>%
            dplyr::select(dplyr::all_of(c(variant, names(input_list)))) %>%
            dplyr::group_by(.data[[variant]]) %>%
            dplyr::summarize(dplyr::across(.cols = names(input_list), mean)) %>%
            tibble::column_to_rownames(var = variant) %>%
            as.matrix()

    if (isTRUE(flip)){
      data <- t(data)
      row_title_use <- column_title
      column_title_use <- row_title
    } else {
      row_title_use <- row_title
      column_title_use <- column_title
    }


    out <- heatmap_inner(data,
                         legend.title = legend.title,
                         column_title = column_title,
                         row_title = row_title,
                         cluster_columns = cluster_cols,
                         cluster_rows = cluster_rows,
                         column_names_rot = column_names_rot,
                         row_names_rot = row_names_rot,
                         row_names_side = row_names_side,
                         row_title_side = row_title_side,
                         row_title_rotation = row_title_rot,
                         column_title_side = "top",
                         cell_size = cell_size,
                         na.value = na.value,
                         use_viridis = use_viridis,
                         viridis_color_map = viridis_color_map,
                         viridis_direction = viridis_direction,
                         legend.position = legend.position,
                         legend.length = heatmap.legend.length,
                         range.data = range.data,
                         legend.width = heatmap.legend.width,
                         legend.framecolor = heatmap.legend.framecolor,
                         data_range = "both",
                         symmetrical_scale = enforce_symmetry)
    list.heatmaps[[variant]] <- out[["heatmap"]]
    list.legends[[variant]] <- out[["legend"]]
  }


  # Compute joint heatmap.
  grDevices::pdf(NULL)
  ht_list <- NULL
  # Append heatmaps vertically.
  suppressWarnings({
    for (heatmap in list.heatmaps){
      if (isTRUE(flip)){
        ht_list <- ht_list + heatmap
      } else if (isFALSE(flip)){
        ht_list <- ht_list %v% heatmap
      }
    }
  })

  # Control gap between legends.
  ComplexHeatmap::ht_opt(message = FALSE)

  # Draw final heatmap.
  h <- ComplexHeatmap::draw(ht_list,
                            heatmap_legend_list = list.legends[[1]],
                            heatmap_legend_side = if (legend.position %in% c("top", "bottom")){"bottom"} else {"right"},
                            padding = ggplot2::unit(c(5, 5, 5, 5), "mm"),
                            ht_gap = ggplot2::unit(heatmap_gap, "cm"))
  grDevices::dev.off()

  out.list <- list()
  out.list[["heatmap"]] <- h

  if (isTRUE(plot_FeaturePlots)){
    list.FeaturePlots <- list()
    for (sig in names(input_gene_list)){
      p <- do_FeaturePlot(sample = sample,
                          features = sig,
                          assay = if (is.null(assay)){Seurat::DefaultAssay(sample)} else {assay},
                          reduction = if(is.null(reduction)){Seurat::Reductions(sample)[length(Seurat::Reductions(sample))]} else {reduction},
                          slot = if (is.null(slot)){"data"} else {slot},
                          pt.size = pt.size,
                          order = FALSE,
                          border.size = border.size,
                          enforce_symmetry = enforce_symmetry,
                          plot_cell_borders = plot_cell_borders,
                          font.size = font.size,
                          font.type = font.type,
                          legend.position = legend.position,
                          legend.type = legend.type,
                          legend.framecolor = legend.framecolor,
                          legend.tickcolor = legend.tickcolor,
                          legend.framewidth = legend.framewidth,
                          legend.tickwidth = legend.tickwidth,
                          legend.length = legend.length,
                          legend.width = legend.width,
                          legend.title = paste0(sig, " enrichment"),
                          viridis_color_map = viridis_color_map,
                          viridis_direction = viridis_direction,
                          min.cutoff = if (is.null(min.cutoff)) {NA} else {min.cutoff},
                          max.cutoff = if (is.null(max.cutoff)) {NA} else {max.cutoff})
      list.FeaturePlots[[sig]] <- p
    }
    out.list[["FeaturePlots"]] <- list.FeaturePlots
  }

  for (var in group.by){
    list.group.by = list()
    if (isTRUE(plot_GeyserPlots)){
      list.GeyserPlots <- list()
      for (sig in names(input_gene_list)){
        p <- do_GeyserPlot(sample = sample,
                           assay = if (is.null(assay)){Seurat::DefaultAssay(sample)} else {assay},
                           slot = if (is.null(slot)){"data"} else {slot},
                           features = sig,
                           group.by = var,
                           pt.size = pt.size,
                           border.size = border.size,
                           enforce_symmetry = enforce_symmetry,
                           scale_type = geyser_scale_type,
                           order_by_mean = geyser_order_by_mean,
                           plot_cell_borders = plot_cell_borders,
                           font.size = font.size,
                           font.type = font.type,
                           legend.position = legend.position,
                           legend.type = legend.type,
                           legend.framecolor = legend.framecolor,
                           legend.tickcolor = legend.tickcolor,
                           legend.framewidth = legend.framewidth,
                           legend.tickwidth = legend.tickwidth,
                           legend.length = legend.length,
                           legend.width = legend.width,
                           xlab = if (is.null(group.by)) {"Clusters"} else {group.by},
                           ylab = paste0(sig, " enrichment"),
                           legend.title = paste0(sig, " enrichment"),
                           rotate_x_axis_labels = rotate_x_axis_labels,
                           viridis_color_map = viridis_color_map,
                           viridis_direction = viridis_direction,
                           min.cutoff = min.cutoff,
                           max.cutoff = max.cutoff)
        list.GeyserPlots[[sig]] <- p
      }
      list.group.by[[var]] <- list.GeyserPlots
    }
    out.list[["GeyserPlots"]] <- list.group.by
  }

  if (isFALSE(plot_GeyserPlots) & isFALSE(plot_FeaturePlots)){
    return_object <- h
  } else {
    return_object <- out.list
  }

  return(return_object)
}
