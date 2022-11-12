#' Create heatmaps of averaged expression by groups.
#'
#' This function generates a heatmap with averaged expression values by the unique groups of the metadata variables provided by the user.
#'
#' @inheritParams doc_function
#'
#'
#' @return A ComplexHeatmap object.
#' @export
#'
#' @examples NULL
do_ExpressionHeatmap <- function(sample,
                                 features,
                                 group.by = NULL,
                                 assay = NULL,
                                 slot = "data",
                                 flip = FALSE,
                                 column_title = if (isFALSE(flip)){"Genes"} else {"Groups"},
                                 row_title = if (isFALSE(flip)){"Groups"} else {"Genes"},
                                 cluster_cols = FALSE,
                                 cluster_rows = FALSE,
                                 legend.title = "Avg. Expression",
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
                                 rotate_x_axis_labels = 45,
                                 enforce_symmetry = FALSE,
                                 heatmap_gap = 0.5,
                                 row_names_side = "right",
                                 row_title_side = "left",
                                 row_title_rot = 90,
                                 min.cutoff = NULL,
                                 max.cutoff = NULL){


  check_suggests(function_name = "do_EnrichmentHeatmap")
  # Check if the sample provided is a Seurat object.
  check_Seurat(sample = sample)

  # Check logical parameters.
  logical_list <- list("flip" = flip,
                       "cluster_cols" = cluster_cols,
                       "cluster_rows" = cluster_rows,
                       "use_viridis" = use_viridis,
                       "enforce_symmetry" = enforce_symmetry)
  check_type(parameters = logical_list, required_type = "logical", test_function = is.logical)
  # Check numeric parameters.
  numeric_list <- list("row_names_rot" = row_names_rot,
                       "column_names_rot" = column_names_rot,
                       "cell_size" = cell_size,
                       "viridis_direction" = viridis_direction,
                       "row_title_rot" = row_title_rot,
                       "rotate_x_axis_labels" = rotate_x_axis_labels,
                       "min.cutoff" = min.cutoff,
                       "max.cutoff" = max.cutoff,
                       "heatmap.legend.length" = heatmap.legend.length,
                       "heatmap.legend.width" = heatmap.legend.width,
                       "heatmap_gap" = heatmap_gap)
  check_type(parameters = numeric_list, required_type = "numeric", test_function = is.numeric)

  # Check character parameters.
  character_list <- list("features" = features,
                         "column_title" = column_title,
                         "row_title" = row_title,
                         "legend.title" = legend.title,
                         "legend.position" = legend.position,
                         "heatmap.legend.framecolor" = heatmap.legend.framecolor,
                         "group.by" = group.by,
                         "na.value" = na.value,
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

  assay <- if (is.null(assay)){Seurat::DefaultAssay(sample)} else {assay}

  if (is.null(group.by)){
    sample@meta.data[, "Groups"] <- sample@active.ident
    group.by <- "Groups"
  }

  features <- check_feature(sample = sample, features = features, permissive = TRUE)

  list.data <- list()
  list.heatmaps <- list()
  list.legends <- list()

  # Check the different values of group.by.
  for (variant in group.by){
    assertthat::assert_that(variant %in% colnames(sample@meta.data),
                            msg = paste0("Value ", variant, " in group.by is not in the sample metadata."))

    assertthat::assert_that(class(sample@meta.data[, variant]) %in% c("character", "factor"),
                            msg = paste0("Value ", variant, " in group.by is not a character or factor column in the metadata."))
  }


  # Generate the heatmap data.
  for (variant in group.by){
    data <- sample@meta.data %>%
      tibble::rownames_to_column(var = "cell") %>%
      dplyr::select(dplyr::all_of(c("cell", variant))) %>%
      tibble::tibble() %>%
      dplyr::left_join(y = {Seurat::GetAssayData(object = sample,
                                                 assay = assay,
                                                 slot = slot)[features, ] %>%
          as.data.frame() %>%
          as.matrix() %>%
          t() %>%
          as.data.frame() %>%
          tibble::rownames_to_column(var = "cell") %>%
          tibble::tibble()},
          by = "cell") %>%
      dplyr::select(-dplyr::all_of(c("cell"))) %>%
      dplyr::group_by(.data[[variant]]) %>%
      dplyr::summarise_at(.vars = features,
                          .funs = mean) %>%
      as.data.frame() %>%
      tibble::column_to_rownames(var = variant) %>%
      as.matrix()
    list.data[[variant]] <- data
  }

  # Apply cutoffs if necessary.

  max_value_list <- c()
  min_value_list <- c()
  # Compute the max values for all heatmaps.
  for (variant in group.by){
    data <- list.data[[variant]]

    max_value_list <- c(max_value_list, max(data))
    min_value_list <- c(min_value_list, min(data))
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

  # Plot the heatmaps.
  for (variant in group.by){
    data <- list.data[[variant]]
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

  return(h)
}
