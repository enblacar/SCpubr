#' Create correlation matrix heatmaps.
#'
#' @inheritParams doc_function
#' @param mode \strong{\code{\link[base]{character}}} | Different types of correlation matrices can be computed. Right now, the only possible value is "hvg", standing for Highly Variable Genes. The sample is subset for the HVG and the data is re-scaled. Scale data is used for the correlation.
#'
#' @return A ComplexHeatmap object.
#' @export
#'
#' @example /man/examples/examples_do_CorrelationPlot.R
do_CorrelationPlot <- function(sample,
                               mode = "hvg",
                               assay = NULL,
                               group.by = NULL,
                               column_title = "",
                               row_title = "",
                               cluster_cols = TRUE,
                               cluster_rows = TRUE,
                               legend.title = "Pearson coef.",
                               row_names_rot = 0,
                               column_names_rot = 0,
                               use_viridis = FALSE,
                               symmetrical_scale = TRUE,
                               viridis_color_map = "G",
                               viridis_direction = 1,
                               cell_size = 8,
                               na.value = "grey75",
                               legend.position = "bottom",
                               legend.length = 75,
                               legend.width = 5,
                               legend.framecolor = "black",
                               min.cutoff = NA,
                               max.cutoff = NA){

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
                       "min.cutoff" = min.cutoff,
                       "max.cutoff" = max.cutoff)
  check_type(parameters = numeric_list, required_type = "numeric", test_function = is.numeric)
  # Check character parameters.
  character_list <- list("mode" = mode,
                         "assay" = assay,
                         "column_title" = column_title,
                         "row_title" = row_title,
                         "legend.title" = legend.title,
                         "viridis_color_map" = viridis_color_map,
                         "group.by" = group.by,
                         "na.value" = na.value)
  check_type(parameters = character_list, required_type = "character", test_function = is.character)
  check_colors(na.value)

  check_parameters(parameter = legend.position, parameter_name = "legend.position")

  `%v%` <- ComplexHeatmap::`%v%`
  `%>%` <- magrittr::`%>%`

  out <- check_and_set_assay(sample = sample, assay = assay)
  sample <- out[["sample"]]
  assay <- out[["assay"]]

  if (mode == "hvg"){
    if (is.null(group.by)){
      assertthat::assert_that(!("Groups" %in% colnames(sample@meta.data)),
                              msg = "Please make sure you provide a value for group.by or do not have a metadata column named `Groups`.")

      sample@meta.data[, "Groups"] <- sample@active.ident
      group.by <- "Groups"
    }

    # Generate a correlation matrix of the HVG.
    variable_genes <- Seurat::VariableFeatures(sample)

    # Subset sample according to the variable genes.
    sample <- sample[variable_genes, ]
    # Scale the data
    sample <- Seurat::ScaleData(sample, verbose = FALSE)

    expr_mat <- data.frame("rownames" = rownames(sample))

    # Retrieve correlation matrix.
    out <- sample@meta.data %>%
           dplyr::select(dplyr::all_of(c(group.by))) %>%
           tibble::rownames_to_column(var = "cell") %>%
           dplyr::left_join(y = {Seurat::GetAssayData(object = sample,
                                                      assay = assay,
                                                      slot = "scale.data") %>%
                                 as.matrix() %>%
                                 t() %>%
                                 as.data.frame() %>%
                                 tibble::rownames_to_column(var = "cell") %>%
                                 tidyr::pivot_longer(-"cell",
                                                     names_to = "gene",
                                                     values_to = "expression")},
                            by = "cell") %>%
           dplyr::select(-"cell") %>%
           dplyr::group_by(.data[[group.by]], .data[["gene"]]) %>%
           dplyr::summarise(mean_expression = mean(.data[["expression"]])) %>%
           tidyr::pivot_wider(names_from = dplyr::all_of(c(group.by)),
                              values_from = "mean_expression") %>%
           as.data.frame() %>%
           tibble::column_to_rownames(var = "gene") %>%
           as.matrix() %>%
           stats::cor() %>%
           round(digits = 2)

    range.data <- c(-1, 1)
    if (!is.na(min.cutoff) & !is.na(max.cutoff)){
      assertthat::assert_that(min.cutoff < max.cutoff,
                              msg = paste0("The value provided for min.cutoff (", min.cutoff, ") has to be lower than the value provided to max.cutoff (", max.cutoff, "). Please select another value."))
    }
    outlier.data <- FALSE
    if (!is.na(min.cutoff)){
      assertthat::assert_that(min.cutoff >= range.data[1],
                              msg = paste0("The value provided for min.cutoff (", min.cutoff, ") is lower than the minimum value in the enrichment matrix (", range.data[1], "). Please select another value."))

      assertthat::assert_that(min.cutoff <= range.data[2],
                              msg = paste0("The value provided for min.cutoff (", min.cutoff, ") is higher than the highest value in the enrichment matrix (", range.data[2], "). Please select another value."))

      range.data <- c(min.cutoff, range.data[2])
      outlier.data <- TRUE
    }

    if (!is.na(max.cutoff)){
      assertthat::assert_that(max.cutoff <= range.data[2],
                              msg = paste0("The value provided for max.cutoff (", max.cutoff, ") is higher than the maximum value in the enrichment matrix (", range.data[2], "). Please select another value."))

      assertthat::assert_that(max.cutoff >= range.data[1],
                              msg = paste0("The value provided for max.cutoff (", max.cutoff, ") is lower than the minimum value in the enrichment matrix (", range.data[1], "). Please select another value."))
      range.data <- c(range.data[1], max.cutoff)
      outlier.data <- TRUE
    }

    if (isTRUE(symmetrical_scale)){
      range.data <- c(-max(abs(range.data)), max(abs(range.data)))
    }

    h <- out %>%
           heatmap_inner(legend.title = legend.title,
                         column_title = column_title,
                         row_title = row_title,
                         row_title_side = "right",
                         row_title_rotation = 0,
                         column_names_rot = column_names_rot,
                         row_names_rot = row_names_rot,
                         cluster_columns = cluster_cols,
                         cluster_rows = cluster_rows,
                         use_viridis = use_viridis,
                         viridis_color_map = viridis_color_map,
                         viridis_direction = viridis_direction,
                         cell_size = cell_size,
                         na.value = na.value,
                         range.data = range.data,
                         legend.position = legend.position,
                         legend.length = legend.length,
                         legend.width = legend.width,
                         legend.framecolor = legend.framecolor,
                         outlier.data = outlier.data,
                         outlier.data.up = if(!is.null(max.cutoff)){TRUE} else {FALSE},
                         outlier.data.down = if(!is.null(min.cutoff)) {TRUE} else {FALSE})

    ComplexHeatmap::ht_opt("HEATMAP_LEGEND_PADDING" = ggplot2::unit(8, "mm"))
    suppressWarnings({
      grDevices::pdf(NULL)
      h <- ComplexHeatmap::draw(h[["heatmap"]],
                                heatmap_legend_list = h[["legend"]],
                                heatmap_legend_side = if (legend.position %in% c("top", "bottom")){"bottom"} else {"right"},
                                padding = ggplot2::unit(c(5, 5, 5, 5), "mm"))
      grDevices::dev.off()
    })

  }
  return(h)
}
