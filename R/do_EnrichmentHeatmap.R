#' Create enrichment scores heatmaps.
#'
#' This function computes the enrichment scores for the cells using \link[Seurat]{AddModuleScore} and then aggregates the scores by the metadata variables provided by the user and displays it as a heatmap, computed by \link[ComplexHeatmap]{Heatmap}.
#'
#' @param sample Your Seurat object.
#' @param list_genes Named list of genes to query to the Seurat object.
#' @param group.by Metadata variable to group the values by.
#' @param column_title,row_title Title for the columns and rows. Only works with single heatmaps.
#' @param verbose Extended output.
#' @param transpose Logical. Transpose the resulting heatmap.
#' @param cluster_cols,cluster_rows Logical. Cluster the columns or rows.
#' @param column_names_rot,row_names_rot Numeric. Degree in which to
#' @param legend_name Text for the legend title.
#' @param colors.use Vector of 2 colors to use to generate the color scale.
#' @param cell_size Numeric. Size of each cell in the heatmap.
#' @param na.value Value for NAs.
#' @param legend.position Character. Location of the legend.
#' @param scale_direction Numeric. Direction of the viridis scales. Either -1 or 1.
#' @param viridis_color_map Character. Viridis color map for the heatmap. One of A, B, C, D, E, F, G, H.
#' @param use_viridis Logical. Whether to use viridis color scales or a bicolor (blue-white-red) color scale.
#' @param heatmap.legend.length,heatmap.legend.width Numeric. Width and length of the legend in the heatmap.
#' @param heatmap.legend.framecolor Character. Color of the edges and ticks of the legend in the heatmap.
#' @param symmetrical_scale Logical. Whether to make the scale symmetrical. Works best when use_viridis = FALSE.
#' @param heatmap_gap Numeric. Gap in cm between legends or heatmaps.
#' @param row_names_side Character. Either left or right.
#' @param row_title_side Character. Either left or right.
#' @param row_title_rotation Numeric. Degree of rotation of the row titles.
#'
#' @return A ComplexHeatmap object.
#' @export
#'
#' @example /man/examples/examples_do_EnrichmentHeatmap.R
do_EnrichmentHeatmap <- function(sample,
                                 list_genes,
                                 group.by = NULL,
                                 column_title = "List of Genes",
                                 row_title = "Groups",
                                 verbose = FALSE,
                                 transpose = FALSE,
                                 cluster_cols = FALSE,
                                 cluster_rows = FALSE,
                                 legend_name = "Enrichment",
                                 row_names_rot = 0,
                                 column_names_rot = 90,
                                 colors.use = NULL,
                                 cell_size = 5,
                                 na.value = "grey75",
                                 legend.position = "bottom",
                                 use_viridis = TRUE,
                                 viridis_color_map = "G",
                                 scale_direction = -1,
                                 heatmap.legend.length = 75,
                                 heatmap.legend.width = 5,
                                 heatmap.legend.framecolor = "black",
                                 symmetrical_scale = FALSE,
                                 heatmap_gap = 0.5,
                                 row_names_side = "right",
                                 row_title_side = "left",
                                 row_title_rotation = 90){
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
                       "scale_direction" = scale_direction,
                       "row_title_rotation" = row_title_rotation)
  check_type(parameters = numeric_list, required_type = "numeric", test_function = is.numeric)
  # Check character parameters.
  character_list <- list("list_genes" = list_genes,
                         "column_title" = column_title,
                         "row_title" = row_title,
                         "legend_name" = legend_name,
                         "colors.use" = colors.use,
                         "group.by" = group.by,
                         "na.value" = na.value,
                         "legend.position" = legend.position,
                         "viridis_color_map" = viridis_color_map,
                         "row_names_side" = row_names_side,
                         "row_title_side" = row_title_side)
  check_type(parameters = character_list, required_type = "character", test_function = is.character)
  check_colors(na.value)
  check_colors(heatmap.legend.framecolor, parameter_name = "heatmap.legend.framecolor")

  `%v%` <- ComplexHeatmap::`%v%`
  `%>%` <- purrr::`%>%`

  if (is.character(list_genes)){
    # If list_genes is a character of genes.
    input_list <- list("Input" = list_genes)
  } else if (is.list(list_genes)){
    input_list <- list_genes
    if (is.null(names(input_list))){
      stop("Please provide a named list. This is, each gene list has to come with a name.", call. = F)
    }
  }

  # Compute the enrichment scores.
  sample <- compute_enrichment_scores(sample = sample, list_genes = list_genes, verbose = verbose)

  list.heatmaps <- list()
  list.legends <- list()
  if (is.null(group.by)){
    aggr_entities <- levels(sample)
    sample@meta.data[, "dummy"] <- sample@active.ident
    group.by <- "dummy"
  }



  max_value_list <- c()
  min_value_list <- c()
  # Compute the max values for all heatmaps.
  for (variant in group.by){
    max_value <- sample@meta.data %>%
                 dplyr::select(dplyr::all_of(c(variant, names(list_genes)))) %>%
                 dplyr::group_by(.data[[variant]]) %>%
                 dplyr::summarize(dplyr::across(.cols = names(list_genes), mean)) %>%
                 tibble::column_to_rownames(var = variant) %>%
                 as.matrix() %>%
                 max()

    min_value <- sample@meta.data %>%
                 dplyr::select(dplyr::all_of(c(variant, names(list_genes)))) %>%
                 dplyr::group_by(.data[[variant]]) %>%
                 dplyr::summarize(dplyr::across(.cols = names(list_genes), mean)) %>%
                 tibble::column_to_rownames(var = variant) %>%
                 as.matrix() %>%
                 min()

    max_value_list <- c(max_value_list, max_value)
    min_value_list <- c(min_value_list, min_value)
  }
  range.data <- c(min(min_value_list), max(max_value_list))


  for (variant in group.by){
    data <- sample@meta.data %>%
      dplyr::select(dplyr::all_of(c(variant, names(list_genes)))) %>%
      dplyr::group_by(.data[[variant]]) %>%
      dplyr::summarize(dplyr::across(.cols = names(list_genes), mean)) %>%
      tibble::column_to_rownames(var = variant) %>%
      as.matrix()

    if (isTRUE(transpose)){
      data <- t(data)
      row_title_use <- column_title
      column_title_use <- row_title
    } else {
      row_title_use <- row_title
      column_title_use <- column_title
    }


    out <- heatmap_inner(data,
                         legend_name = legend_name,
                         column_title = column_title,
                         row_title = row_title,
                         cluster_columns = cluster_cols,
                         cluster_rows = cluster_rows,
                         column_names_rot = column_names_rot,
                         row_names_rot = row_names_rot,
                         row_names_side = row_names_side,
                         row_title_side = row_title_side,
                         row_title_rotation = row_title_rotation,
                         column_title_side = "top",
                         colors.use = colors.use,
                         cell_size = cell_size,
                         na.value = na.value,
                         use_viridis = use_viridis,
                         viridis_color_map = viridis_color_map,
                         viridis_direction = scale_direction,
                         legend.position = legend.position,
                         legend.length = heatmap.legend.length,
                         range.data = range.data,
                         legend.width = heatmap.legend.width,
                         legend.framecolor = heatmap.legend.framecolor,
                         data_range = "both",
                         symmetrical_scale = symmetrical_scale)
    list.heatmaps[[variant]] <- out$heatmap
    list.legends[[variant]] <- out$legend
  }


  # Compute joint heatmap.
  grDevices::pdf(NULL)
  ht_list <- NULL
  # Append heatmaps vertically.
  suppressWarnings({
    for (heatmap in list.heatmaps){
      ht_list = ht_list %v% heatmap
    }
  })

  # Control gap between legends.
  ComplexHeatmap::ht_opt(message = FALSE)

  # Draw final heatmap.
  h <- ComplexHeatmap::draw(ht_list,
                            heatmap_legend_list = list.legends[1],
                            heatmap_legend_side = if (legend.position %in% c("top", "bottom")){"bottom"} else {"right"},
                            padding = ggplot2::unit(c(5, 5, 5, 5), "mm"),
                            ht_gap = ggplot2::unit(heatmap_gap, "cm"))
  grDevices::dev.off()


  return(h)
}
