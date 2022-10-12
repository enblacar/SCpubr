#' Create enrichment scores heatmaps.
#'
#' This function computes the enrichment scores for the cells using \link[Seurat]{AddModuleScore} and then aggregates the scores by the metadata variables provided by the user and displays it as a heatmap, computed by \link[ComplexHeatmap]{Heatmap}.
#'
#' @inheritParams doc_function
#' @param symmetrical_scale \strong{\code{\link[base]{logical}}} | Whether to make the color scale symmetrical. Works best when use_viridis = FALSE.
#'
#' @return A ComplexHeatmap object.
#' @export
#'
#' @example /man/examples/examples_do_EnrichmentHeatmap.R
do_EnrichmentHeatmap <- function(sample,
                                 input_gene_list,
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
                                 symmetrical_scale = FALSE,
                                 heatmap_gap = 0.5,
                                 row_names_side = "right",
                                 row_title_side = "left",
                                 row_title_rot = 90,
                                 column_title = if (isFALSE(flip)){"List of Genes"} else {"Groups"},
                                 row_title = if (isFALSE(flip)){"Groups"} else {"List of Genes"},
                                 nbin = 24,
                                 ctrl = 100){
  # Checks for packages.
  check_suggests(function_name = "do_EnrichmentHeatmap")
  # Check if the sample provided is a Seurat object.
  check_Seurat(sample = sample)

  # Check logical parameters.
  logical_list <- list("cluster_cols" = cluster_cols,
                       "cluster_rows" = cluster_rows,
                       "use_viridis" = use_viridis,
                       "symmetrical_scale" = symmetrical_scale)
  check_type(parameters = logical_list, required_type = "logical", test_function = is.logical)
  # Check numeric parameters.
  numeric_list <- list("row_names_rot" = row_names_rot,
                       "column_names_rot" = column_names_rot,
                       "cell_size" = cell_size,
                       "viridis_direction" = viridis_direction,
                       "row_title_rot" = row_title_rot,
                       "nbin" = nbin,
                       "ctrl" = ctrl)
  check_type(parameters = numeric_list, required_type = "numeric", test_function = is.numeric)
  # Check character parameters.
  character_list <- list("input_gene_list" = input_gene_list,
                         "column_title" = column_title,
                         "row_title" = row_title,
                         "legend.title" = legend.title,
                         "group.by" = group.by,
                         "na.value" = na.value,
                         "legend.position" = legend.position,
                         "viridis_color_map" = viridis_color_map,
                         "row_names_side" = row_names_side,
                         "row_title_side" = row_title_side)
  check_type(parameters = character_list, required_type = "character", test_function = is.character)
  check_colors(na.value)
  check_colors(heatmap.legend.framecolor, parameter_name = "heatmap.legend.framecolor")

  check_parameters(parameter = viridis_direction, parameter_name = "viridis_direction")
  check_parameters(parameter = legend.position, parameter_name = "legend.position")

  `%v%` <- ComplexHeatmap::`%v%`
  `%>%` <- magrittr::`%>%`

  if (is.character(input_gene_list)){
    # If input_gene_list is a character of genes.
    input_list <- list("Input" = input_gene_list)
  } else if (is.list(input_gene_list)){
    input_list <- input_gene_list
    assertthat::assert_that(!is.null(names(input_list)),
                            msg = "Please provide a named list. This is, each gene list has to come with a name.")
  }

  # Compute the enrichment scores.
  sample <- compute_enrichment_scores(sample = sample, input_gene_list = input_gene_list, verbose = verbose, nbin = nbin, ctrl = ctrl)

  list.heatmaps <- list()
  list.legends <- list()
  if (is.null(group.by)){
    aggr_entities <- levels(sample)
    sample@meta.data[, "dummy"] <- sample@active.ident
    group.by <- "dummy"
  }

  assertthat::assert_that(isTRUE(group.by %in% colnames(sample@meta.data)),
                          msg = "Please provide a name to group.by that is a categorical column in sample@meta.data.")

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
                         symmetrical_scale = symmetrical_scale)
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


  return(h)
}
