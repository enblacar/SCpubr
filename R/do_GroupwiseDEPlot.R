#' Compute a heatmap with the results of a group-wise DE analysis.
#'
#' @inheritParams doc_function
#' @param de_genes \strong{\code{\link[tibble]{tibble}}} | DE genes matrix resulting of running `Seurat::FindAllMarkers()`.
#' @param top_genes \strong{\code{\link[base]{numeric}}} | Top N differentially expressed (DE) genes by group to retrieve.
#' @param viridis_map_pvalues,viridis_map_logfc,viridis_map_expression \strong{\code{\link[base]{character}}} | Viridis color map for the heatmap of p-values, logFC or expression. One of A, B, C, D, E, F, G, H.
#' @param row_title_p_values \strong{\code{\link[base]{character}}} | Row title for the p-value heatmap. Blank by default.
#' @param row_title_logfc \strong{\code{\link[base]{character}}} | Row title for the logfc heatmap. Clusters by default.
#' @param row_title_expression \strong{\code{\link[base]{character}}} | Vector of titles of equal length as group.by.
#'
#' @return A heatmap composed of 3 main panels: -log10(adjusted p-value), log2(FC) and mean expression by cluster.
#' @export
#'
#' @example /man/examples/examples_do_GroupwiseDEPlot.R
do_GroupwiseDEPlot <- function(sample,
                               de_genes,
                               group.by = NULL,
                               viridis_map_pvalues = "B",
                               viridis_map_logfc = "D",
                               viridis_map_expression = "G",
                               heatmap.legend.length = 75,
                               heatmap.legend.width = 5,
                               heatmap.legend.framecolor = "black",
                               top_genes = 5,
                               viridis_direction = -1,
                               row_title_p_values = "",
                               row_title_logfc = "Clusters",
                               row_title_expression = if (is.null(group.by)){""} else {rep("", length(group.by))},
                               column_title = "DE genes",
                               heatmap_gap = 0.5,
                               legend_gap = 1,
                               assay = NULL,
                               slot = "data",
                               legend.position = "bottom",
                               row_names_side = "right",
                               row_title_side = "left",
                               row_title_rot = 90,
                               column_names_rot = 45,
                               cell_size = 6,
                               min.cutoff = NULL,
                               max.cutoff = NULL){
  check_suggests(function_name = "do_GroupwiseDEPlot")
  # Check if the sample provided is a Seurat object.
  check_Seurat(sample = sample)

  # Check numeric parameters.
  numeric_list <- list("heatmap.legend.length" = heatmap.legend.length,
                       "heatmap.legend.width" = heatmap.legend.width,
                       "cell_size" = cell_size,
                       "viridis_direction" = viridis_direction,
                       "top_genes" = top_genes,
                       "heatmap_gap" = heatmap_gap,
                       "legend_gap" = legend_gap,
                       "row_title_rot" = row_title_rot,
                       "column_names_rot" = column_names_rot,
                       "min.cutoff" = min.cutoff,
                       "max.cutoff" = max.cutoff)
  check_type(parameters = numeric_list, required_type = "numeric", test_function = is.numeric)
  # Check character parameters.
  character_list <- list("group.by" = group.by,
                         "viridis_map_pvalues" = viridis_map_pvalues,
                         "viridis_map_logfc" = viridis_map_logfc,
                         "viridis_map_expression" = viridis_map_expression,
                         "heatmap.legend.framecolor" = heatmap.legend.framecolor,
                         "row_title_p_values" = row_title_p_values,
                         "row_title_logfc" = row_title_logfc,
                         "row_title_expression" = row_title_expression,
                         "column_title" = column_title,
                         "assay" = assay,
                         "slot" = slot,
                         "legend.position" = legend.position,
                         "row_names_side" = row_names_side,
                         "row_title_side" = row_title_side)
  check_type(parameters = character_list, required_type = "character", test_function = is.character)

  check_colors(heatmap.legend.framecolor, parameter_name = "heatmap.legend.framecolor")

  `%v%` <- ComplexHeatmap::`%v%`
  `%>%` <- magrittr::`%>%`
  `.` <- plyr::.()

  check_parameters(parameter = legend.position, parameter_name = "legend.position")
  check_parameters(parameter = viridis_direction, parameter_name = "viridis_direction")
  check_parameters(parameter = viridis_map_pvalues, parameter_name = "viridis_color_map")
  check_parameters(parameter = viridis_map_expression, parameter_name = "viridis_color_map")
  check_parameters(parameter = viridis_map_logfc, parameter_name = "viridis_color_map")

  # Check the assay.
  out <- check_and_set_assay(sample = sample, assay = assay)
  sample <- out[["sample"]]
  assay <- out[["assay"]]

  if (!is.null(group.by)){
    for (value in group.by){
      assertthat::assert_that(isTRUE(value %in% colnames(sample@meta.data)),
                              msg = "Please provide a value for group.by that corresponds to a metadata column.")

      assertthat::assert_that(isTRUE(class(sample@meta.data[, value]) %in% c("factor", "character")),
                              msg = "Please provide a value for group.by that corresponds to a metadata column and this is either a factor or a character column.")
    }
  } else if (is.null(group.by)) {
    sample@meta.data[, "Groups"] <- Seurat::Idents(sample)
    group.by = "Groups"
  }
  assertthat::assert_that(length(group.by) == length(row_title_expression),
                          msg = "Please provide the same number of row titles as the number of items in group.by.")



  magnitude <- ifelse(slot == "data", "avg_log2FC", "avg_diff")
  assertthat::assert_that(min(de_genes[, magnitude]) >= 0,
                          msg = "Please provide a de_genes object in which avg_log2FC/avg.diff has only positive values (including 0).")
  # Compute the top N genes per cluster.
  genes.use <- de_genes %>%
    dplyr::arrange(.data$p_val_adj, dplyr::desc(.data[[magnitude]])) %>%
    dplyr::group_by(.data$cluster) %>%
    dplyr::slice_head(n = top_genes) %>%
    dplyr::pull("gene") %>%
    unique()

  # Compute heatmap of log2FC.
  logfc_out <- de_genes %>%
    dplyr::arrange(.data$p_val_adj, dplyr::desc(.data[[magnitude]])) %>%
    dplyr::group_by(.data$cluster) %>%
    dplyr::slice_head(n = top_genes) %>%
    dplyr::select(dplyr::all_of(c("gene", "cluster", magnitude))) %>%
    tidyr::pivot_wider(names_from = "gene",
                       values_from = dplyr::all_of(c(magnitude))) %>%
    as.data.frame() %>%
    tibble::column_to_rownames(var = "cluster") %>%
    as.matrix() %>%
    replace(is.na(.), 0) %>%
    heatmap_inner(cluster_columns = FALSE,
                  cluster_rows = FALSE,
                  legend.title = ifelse(slot == "data", "Avg. log2(FC)", "Avg. Diff."),
                  data_range = "only_pos",
                  row_title = row_title_logfc,
                  row_names_side = row_names_side,
                  legend.position = legend.position,
                  legend.length = heatmap.legend.length,
                  legend.width = heatmap.legend.width,
                  legend.framecolor = heatmap.legend.framecolor,
                  use_viridis = TRUE,
                  viridis_color_map = viridis_map_logfc,
                  viridis_direction = viridis_direction,
                  zeros_are_white = TRUE,
                  row_title_rotation = row_title_rot,
                  row_title_side = row_title_side,
                  column_names_rot = column_names_rot,
                  cell_size = cell_size)

  # Compute heatmap of -log10FC.
  pvalue_out <- de_genes %>%
    dplyr::arrange(.data$p_val_adj, dplyr::desc(.data[[magnitude]])) %>%
    dplyr::group_by(.data$cluster) %>%
    dplyr::slice_head(n = top_genes) %>%
    dplyr::select(dplyr::all_of(c("gene", "cluster", "p_val_adj"))) %>%
    dplyr::mutate("p_val_adj" = replace(.data$p_val_adj, .data$p_val_adj == 0, .Machine$double.xmin)) %>%
    dplyr::mutate("log10pval" = -log10(.data$p_val_adj)) %>%
    dplyr::select(-"p_val_adj") %>%
    tidyr::pivot_wider(names_from = "gene",
                       values_from = "log10pval") %>%
    as.data.frame() %>%
    tibble::column_to_rownames(var = "cluster") %>%
    as.matrix() %>%
    replace(is.na(.), 0)

  pvalue_out <- heatmap_inner(data = pvalue_out,
                              cluster_columns = FALSE,
                              cluster_rows = FALSE,
                              row_names_side = row_names_side,
                              row_title = row_title_p_values,
                              legend.title = "-log10(Adjusted P-value)",
                              data_range = "only_pos",
                              legend.position = legend.position,
                              legend.length = heatmap.legend.length,
                              legend.width = heatmap.legend.width,
                              legend.framecolor = heatmap.legend.framecolor,
                              na.value = "grey75",
                              column_title = column_title,
                              outlier.data = FALSE,
                              range.data = NULL,
                              outlier.up.label = "Inf",
                              use_viridis = TRUE,
                              viridis_color_map = viridis_map_pvalues,
                              viridis_direction = viridis_direction,
                              zeros_are_white = TRUE,
                              row_title_rotation = row_title_rot,
                              row_title_side = row_title_side,
                              column_names_rot = column_names_rot,
                              cell_size = cell_size)

  # Compute heatmap of expression.
  list.expression.heatmaps <- list()
  list.expression.legends <- list()

  max_value_list <- c()
  min_value_list <- c()
  # Compute the max values for all heatmaps.
  for (variable in group.by){
    max_value <- sample@meta.data %>%
      dplyr::select(dplyr::all_of(c(variable))) %>%
      tibble::rownames_to_column(var = "cell") %>%
      dplyr::left_join(y = {Seurat::GetAssayData(object = sample,
                                                 slot = slot,
                                                 assay = assay)[genes.use, ] %>%
                                                 as.matrix() %>%
                                                 t() %>%
                                                 as.data.frame() %>%
                                                 tibble::rownames_to_column(var = "cell")},
                                                 by = "cell") %>%
      dplyr::select(-"cell") %>%
      tidyr::pivot_longer(cols = -dplyr::all_of(c(variable)),
                          names_to = "gene",
                          values_to = "expression") %>%
      dplyr::group_by(.data[[variable]], .data$gene) %>%
      dplyr::summarise(mean_expression = mean(.data$expression)) %>%
      tidyr::pivot_wider(names_from = "gene",
                         values_from = "mean_expression") %>%
      as.data.frame() %>%
      tibble::column_to_rownames(var = variable) %>%
      dplyr::select(dplyr::all_of(genes.use)) %>%
      as.matrix() %>%
      max()

    min_value <- sample@meta.data %>%
      dplyr::select(dplyr::all_of(c(variable))) %>%
      tibble::rownames_to_column(var = "cell") %>%
      dplyr::left_join(y = {Seurat::GetAssayData(object = sample,
                                                 slot = slot,
                                                 assay = assay)[genes.use, ] %>%
                                                 as.matrix() %>%
                                                 t() %>%
                                                 as.data.frame() %>%
                                                 tibble::rownames_to_column(var = "cell")},
                                                 by = "cell") %>%
      dplyr::select(-"cell") %>%
      tidyr::pivot_longer(cols = -dplyr::all_of(c(variable)),
                          names_to = "gene",
                          values_to = "expression") %>%
      dplyr::group_by(.data[[variable]], .data$gene) %>%
      dplyr::summarise(mean_expression = mean(.data$expression)) %>%
      tidyr::pivot_wider(names_from = "gene",
                         values_from = "mean_expression") %>%
      as.data.frame() %>%
      tibble::column_to_rownames(var = variable) %>%
      dplyr::select(dplyr::all_of(genes.use)) %>%
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
  counter <- 0
  for (variable in group.by){
    counter <- counter + 1

    expression_out <- sample@meta.data %>%
      dplyr::select(dplyr::all_of(c(variable))) %>%
      tibble::rownames_to_column(var = "cell") %>%
      dplyr::left_join(y = {Seurat::GetAssayData(object = sample,
                                                 slot = slot,
                                                 assay = assay)[genes.use, ] %>%
          as.matrix() %>%
          t() %>%
          as.data.frame() %>%
          tibble::rownames_to_column(var = "cell")},
          by = "cell") %>%
      dplyr::select(-"cell") %>%
      tidyr::pivot_longer(cols = -dplyr::all_of(c(variable)),
                          names_to = "gene",
                          values_to = "expression") %>%
      dplyr::group_by(.data[[variable]], .data$gene) %>%
      dplyr::summarise(mean_expression = mean(.data$expression)) %>%
      tidyr::pivot_wider(names_from = "gene",
                         values_from = "mean_expression") %>%
      as.data.frame() %>%
      tibble::column_to_rownames(var = variable) %>%
      dplyr::select(dplyr::all_of(genes.use)) %>%
      as.matrix() %>%
      heatmap_inner(cluster_columns = FALSE,
                    cluster_rows = FALSE,
                    row_names_side = row_names_side,
                    legend.title = "Mean expression",
                    data_range = "both",
                    row_title = row_title_expression[counter],
                    legend.position = legend.position,
                    legend.length = heatmap.legend.length,
                    legend.width = heatmap.legend.width,
                    legend.framecolor = heatmap.legend.framecolor,
                    use_viridis = TRUE,
                    viridis_color_map = viridis_map_expression,
                    viridis_direction = viridis_direction,
                    zeros_are_white = TRUE,
                    range.data = range.data,
                    row_title_rotation = row_title_rot,
                    row_title_side = row_title_side,
                    column_names_rot = column_names_rot,
                    cell_size = cell_size)
    list.expression.heatmaps[[variable]] <- expression_out$heatmap
    list.expression.legends[[variable]] <- expression_out$legend
  }



  # Compute joint heatmap.
  grDevices::pdf(NULL)
  ht_list <- NULL
  list_heatmaps <- c(pvalue_out$heatmap,
                     logfc_out$heatmap,
                     unlist(list.expression.heatmaps))
  list_legends <- c(pvalue_out$legend,
                    logfc_out$legend,
                    list.expression.legends[1])

  # Append heatmaps vertically.
  suppressWarnings({
    for (heatmap in list_heatmaps){
      ht_list <- ht_list %v% heatmap
    }
  })

  # Control gap between legends.
  ComplexHeatmap::ht_opt(legend_gap = ggplot2::unit(rep(legend_gap, 2), "cm"),
                         message = FALSE)

  # Draw final heatmap.
  h <- ComplexHeatmap::draw(ht_list,
                            heatmap_legend_list = list_legends,
                            heatmap_legend_side = if (legend.position %in% c("top", "bottom")){"bottom"} else {"right"},
                            padding = ggplot2::unit(c(5, 5, 5, 5), "mm"),
                            ht_gap = ggplot2::unit(heatmap_gap, "cm"))
  grDevices::dev.off()

  # Return the final heatmap.
  return(h)
}
