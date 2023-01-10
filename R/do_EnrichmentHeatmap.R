#' Create enrichment scores heatmaps.
#'
#' This function computes the enrichment scores for the cells using \link[Seurat]{AddModuleScore} and then aggregates the scores by the metadata variables provided by the user and displays it as a heatmap, computed by \link[ComplexHeatmap]{Heatmap}.
#'
#' @inheritParams doc_function
#' @param enforce_symmetry \strong{\code{\link[base]{logical}}} | Whether the geyser and feature plot has a symmetrical color scale.
#' @param flavor \strong{\code{\link[base]{character}}} | One of: Seurat, UCell. Compute the enrichment scores using \link[Seurat]{AddModuleScore} or \link[UCell]{AddModuleScore_UCell}.
#' @param ncores \strong{\code{\link[base]{numeric}}} | Number of cores used to run UCell scoring.
#' @param storeRanks \strong{\code{\link[base]{logical}}} | Whether to store the ranks for faster UCell scoring computations. Might require large amounts of RAM.
#' @param geyser_order_by_mean,boxplot_order_by_mean \strong{\code{\link[base]{logical}}} | Whether to order the X axis by the mean of the values.
#' @param geyser_scale_type \strong{\code{\link[base]{character}}} | Type of scale to use. Either "continuous" or "categorical.
#' @param violin_plot_boxplot \strong{\code{\link[base]{logical}}} | Whether to plot the boxplots inside the violin plots.
#' @param violin_boxplot_width \strong{\code{\link[base]{numeric}}} | Width of the boxplots in the violin plots.
#' @param plot_FeaturePlots,plot_GeyserPlots,plot_BeeSwarmPlots,plot_BoxPlots,plot_ViolinPlots \strong{\code{\link[base]{logical}}} | Compute extra visualizations for each of the gene lists.
#' @param return_object \strong{\code{\link[base]{logical}}} | Return the Seurat object with the enrichment scores stored.
#' @param return_matrix \strong{\code{\link[base]{logical}}} | Return the enrichment matrix used for the heatmaps for each value in group.by.
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
                                 row_names_rot = 0,
                                 column_names_rot = 45,
                                 cell_size = 8,
                                 na.value = "grey75",
                                 legend.position = "bottom",
                                 use_viridis = TRUE,
                                 viridis_color_map = "G",
                                 viridis_direction = 1,
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
                                 column_title = NULL,
                                 row_title = NULL,
                                 nbin = 24,
                                 ctrl = 100,
                                 flavor = "Seurat",
                                 legend.title = if (flavor != "AUCell") {"Enrichment"} else {"AUC"},
                                 ncores = 1,
                                 storeRanks = TRUE,
                                 min.cutoff = NULL,
                                 max.cutoff = NULL,
                                 plot_FeaturePlots = FALSE,
                                 plot_GeyserPlots = FALSE,
                                 plot_BeeSwarmPlots = FALSE,
                                 plot_BoxPlots = FALSE,
                                 plot_ViolinPlots = FALSE,
                                 pt.size = 1,
                                 plot_cell_borders = TRUE,
                                 border.size = 2,
                                 geyser_order_by_mean = TRUE,
                                 geyser_scale_type = "continuous",
                                 boxplot_order_by_mean = TRUE,
                                 violin_plot_boxplot = TRUE,
                                 violin_boxplot_width = 0.2,
                                 return_object = FALSE,
                                 return_matrix = FALSE){
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
                       "geyser_order_by_mean" = geyser_order_by_mean,
                       "plot_BoxPlots" = plot_BoxPlots,
                       "plot_BeeSwarmPlots" = plot_BeeSwarmPlots,
                       "plot_ViolinPlots" = plot_ViolinPlots,
                       "boxplot_order_by_mean" = boxplot_order_by_mean,
                       "violin_plot_boxplot" = violin_plot_boxplot)
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
                       "max.cutoff" = max.cutoff,
                       "violin_boxplot_width" = violin_boxplot_width)
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
    stop("If setting flavor = UCell, do not use assay parameter. Instead, make sure the assay you want is set as default.", call. = FALSE)
  }

  if (!(is.null(slot)) & flavor == "Seurat"){
    stop("If setting flavor = Seurat, do not use slot parameter. The slot is determined by default in Seurat.", call. = FALSE)
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
                                      input_gene_list = input_list,
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
  list.matrices <- list()
  if (is.null(group.by)){
    assertthat::assert_that(!("Groups" %in% colnames(sample@meta.data)),
                            msg = "Please make sure you provide a value for group.by or do not have a metadata column named `Groups`.")

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

  # Fix for automatic row and column titles.
  if (is.null(column_title)){
    if (length(group.by) == 1){
      column_title <- ifelse(isTRUE(flip), "Groups", "List of genes")
    } else {
      if (isTRUE(flip)){
        column_title <- rep("Groups", length(group.by))
      } else {
        column_title <- c("List of genes", rep("", length(group.by) - 1))
      }
    }
  } else {
    assertthat::assert_that(length(column_title) == length(group.by),
                            msg = "Please provide as many different column titles as unique values in group.by.")
  }

  if (is.null(row_title)){
    if (length(group.by) == 1){
      row_title <- ifelse(isFALSE(flip), "Groups", "List of genes")
    } else {
      if (isFALSE(flip)){
        row_title <- rep("Groups", length(group.by))
      } else {
        row_title <- c("List of genes", rep("", length(group.by) - 1))
      }
    }
  } else {
    assertthat::assert_that(length(row_title) == length(group.by),
                            msg = "Please provide as many different row titles as unique values in group.by.")
  }


  counter <- 0
  for (variant in group.by){
    counter <- counter + 1
    data <- sample@meta.data %>%
            dplyr::select(dplyr::all_of(c(variant, names(input_list)))) %>%
            dplyr::group_by(.data[[variant]]) %>%
            dplyr::summarize(dplyr::across(.cols = names(input_list), mean)) %>%
            tibble::column_to_rownames(var = variant) %>%
            as.matrix()

    row_title_use <- row_title[counter]
    column_title_use <- column_title[counter]


    out <- heatmap_inner(if (isTRUE(flip)){t(data)} else {data},
                         legend.title = legend.title,
                         column_title = column_title_use,
                         row_title = row_title_use,
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
                         symmetrical_scale = enforce_symmetry,
                         zeros_are_white = if (flavor %in% c("AUCell", "UCell") & viridis_direction == -1) {TRUE} else {FALSE})
    list.heatmaps[[variant]] <- out[["heatmap"]]
    list.legends[[variant]] <- out[["legend"]]
    list.matrices[[variant]] <- data
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

  # Plot extra FeaturePlots.
  if (isTRUE(plot_FeaturePlots)){
    list.FeaturePlots <- list()
    for (sig in names(input_list)){
      p <- SCpubr::do_FeaturePlot(sample = sample,
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
                                  legend.title = if (flavor != "AUCell") {paste0(sig, " enrichment")} else {paste0(sig, " AUC")},
                                  viridis_color_map = viridis_color_map,
                                  viridis_direction = viridis_direction,
                                  min.cutoff = if (is.null(min.cutoff)) {NA} else {min.cutoff},
                                  max.cutoff = if (is.null(max.cutoff)) {NA} else {max.cutoff})
      list.FeaturePlots[[sig]] <- p
    }
    out.list[["FeaturePlots"]] <- list.FeaturePlots
  }

  # Plot extra GeyserPlots.
  if (isTRUE(plot_GeyserPlots)){
    list.group.by = list()
    for (var in group.by){
      list.GeyserPlots <- list()

      for (sig in names(input_list)){
        p <- SCpubr::do_GeyserPlot(sample = sample,
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
                                   xlab = if (is.null(group.by)) {"Clusters"} else {var},
                                   ylab = if (flavor != "AUCell") {paste0(sig, " enrichment")} else {paste0(sig, " AUC")},
                                   legend.title = if (flavor != "AUCell") {paste0(sig, " enrichment")} else {paste0(sig, " AUC")},
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

  # Plot extra BeeSwarmPlots.
  if (isTRUE(plot_BeeSwarmPlots)){
    list.group.by = list()
    for(var in group.by){
      list.BeeSwarmPlots <- list()
      for (sig in names(input_list)){
        p <- SCpubr::do_BeeSwarmPlot(sample = sample,
                                     feature_to_rank = sig,
                                     assay = if (is.null(assay)){Seurat::DefaultAssay(sample)} else {assay},
                                     slot = if (is.null(slot)){"data"} else {slot},
                                     group.by = var,
                                     pt.size = pt.size,
                                     border.size = border.size,
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
                                     continuous_feature = TRUE,
                                     min.cutoff = min.cutoff,
                                     max.cutoff = max.cutoff,
                                     ylab = if (is.null(group.by)) {"Clusters"} else {var},
                                     xlab = if (flavor != "AUCell") {paste0("Ranking of ", sig, " enrichment")} else {paste0("Ranking of ", sig, " AUC")},
                                     viridis_color_map = viridis_color_map,
                                     viridis_direction = viridis_direction)
        list.BeeSwarmPlots[[sig]] <- p
      }
      list.group.by[[var]] <- list.BeeSwarmPlots
    }
    out.list[["BeeSwarmPlots"]] <- list.group.by
  }

  # Plot extra BoxPlots.
  if (isTRUE(plot_BoxPlots)){
    list.group.by = list()
    for(var in group.by){
      list.BoxPlots <- list()
      for (sig in names(input_list)){
        p <- SCpubr::do_BoxPlot(sample = sample,
                                feature = sig,
                                assay = if (is.null(assay)){Seurat::DefaultAssay(sample)} else {assay},
                                slot = if (is.null(slot)){"data"} else {slot},
                                group.by = var,
                                legend.position = legend.position,
                                font.size = font.size,
                                font.type = font.type,
                                xlab = if (is.null(group.by)) {"Clusters"} else {var},
                                ylab = if (flavor != "AUCell") {paste0(sig, " enrichment")} else {paste0(sig, " AUC")},
                                order = boxplot_order_by_mean)
        list.BoxPlots[[sig]] <- p
      }
      list.group.by[[var]] <- list.BoxPlots
    }
    out.list[["BoxPlots"]] <- list.group.by
  }

  # Plot extra ViolinPlots.
  if (isTRUE(plot_ViolinPlots)){
    list.group.by = list()
    for (var in group.by){
      list.ViolinPlots <- list()
      for (sig in names(input_list)){
        p <- SCpubr::do_ViolinPlot(sample = sample,
                                   features = sig,
                                   plot_boxplot = violin_plot_boxplot,
                                   boxplot_width = violin_boxplot_width,
                                   assay = if (is.null(assay)){Seurat::DefaultAssay(sample)} else {assay},
                                   slot = if (is.null(slot)){"data"} else {slot},
                                   group.by = var,
                                   legend.position = legend.position,
                                   font.size = font.size,
                                   font.type = font.type,
                                   xlab = if (is.null(group.by)) {"Clusters"} else {var},
                                   ylab = if (flavor != "AUCell") {paste0(sig, " enrichment")} else {paste0(sig, " AUC")})
        list.ViolinPlots[[sig]] <- p
      }
      list.group.by[[var]] <- list.ViolinPlots
    }
    out.list[["ViolinPlots"]] <- list.group.by
  }

  # Return the scored Seurat object.
  if (isTRUE(return_object)){
    if (isTRUE("Groups" %in% colnames(sample@meta.data))){
      sample$Groups <- NULL
    }
    out.list[["object"]] <- sample
  }

  # Return the scoring matrix.
  if (isTRUE(return_matrix)){
    out.list[["matrices"]] <- list.matrices
  }


  if (isFALSE(plot_GeyserPlots) & isFALSE(plot_FeaturePlots) & isFALSE(plot_BeeSwarmPlots) & isFALSE(plot_BoxPlots) & isFALSE(return_object) & isFALSE(return_matrix)){
    return_me <- h
  } else {
    return_me <- out.list
  }

  return(return_me)
}
