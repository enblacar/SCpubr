#' Compute a heatmap with the results of a group-wise DE analysis.
#'
#' @param sample Seurat objects.
#' @param de_genes DE genes matrix resulting of running `Seurat::FindAllMarkers()`.
#' @param assay Character. Assay used to compute the DE genes.
#' @param slot Character. Slot from the assay used to compute the DE genes. Normally, this is "data".
#' @param group.by Character. Metadata variable (s) that was used with `Seurat::FindAllMarkers()`. Generally, it is "seurat_clusters", but it is normally the one used as current identities of the Seurat object. One can add further variables in a vector, and they will be plotted on top of each other.
#' @param top_genes Numeric. Top N DE genes by group to retrieve.
#' @param scale_direction Numeric. Direction of the viridis scales. Either -1 or 1.
#' @param viridis_map_pvalues,viridis_map_logfc,viridis_map_expression Character. Viridis color map for the heatmap of p-values, logFC or expression. One of A, B, C, D, E, F, G, H.
#' @param heatmap.legend.length,heatmap.legend.width Numeric. Width and length of the legend in the heatmap.
#' @param heatmap.legend.framecolor Character. Color of the edges and ticks of the legend in the heatmap.
#' @param column_title Character. Title for the column.
#' @param row_title_p_values Character. Row title for the p-value heatmap. Blank by default.
#' @param row_title_logfc Character. Row title for the logfc heatmap. Clusters by default.
#' @param row_title_expression Character. Vector of titles of equal length as group.by.
#' @param heatmap_gap,legend_gap Numeric. Gap in cm between legends or heatmaps.
#' @param legend.position Character. Where to place the legends. Either "top", "botom", "right", "left".
#' @param row_names_side Character. Either left or right.
#' @param row_title_side Character. Either left or right.
#' @param row_title_rotation Numeric. Degree of rotation of the row titles.
#' @param cell_size Numeric. Size of the cells in the heatmap.
#'
#' @return A heatmap composed of 3 main panels: -log10(adjusted p-value), log2(FC) and mean expression by cluster.
#' @export
#'
#' @example /man/examples/examples_do_GroupwiseDEPlot.R
do_GroupwiseDEPlot <- function(sample,
                               de_genes,
                               group.by = "seurat_clusters",
                               viridis_map_pvalues = "B",
                               viridis_map_logfc = "D",
                               viridis_map_expression = "G",
                               heatmap.legend.length = 75,
                               heatmap.legend.width = 5,
                               heatmap.legend.framecolor = "black",
                               top_genes = 5,
                               scale_direction = -1,
                               row_title_p_values = "",
                               row_title_logfc = "Clusters",
                               row_title_expression = c(""),
                               column_title = "DE genes",
                               heatmap_gap = 0.5,
                               legend_gap = 1,
                               assay = "SCT",
                               slot = "data",
                               legend.position = "bottom",
                               row_names_side = "right",
                               row_title_side = "left",
                               row_title_rotation = 90,
                               cell_size = 5){
  # Checks for packages.
  check_suggests(function_name = "do_GroupwiseDEPlot")
  # Check if the sample provided is a Seurat object.
  check_Seurat(sample = sample)

  # Check logical parameters.
  # logical_list <- list()
  # check_type(parameters = logical_list, required_type = "logical", test_function = is.logical)
  # Check numeric parameters.
  numeric_list <- list("heatmap.legend.length" = heatmap.legend.length,
                       "heatmap.legend.width" = heatmap.legend.width,
                       "cell_size" = cell_size,
                       "scale_direction" = scale_direction,
                       "top_genes" = top_genes,
                       "heatmap_gap" = heatmap_gap,
                       "legend_gap" = legend_gap,
                       "row_title_rotation" = row_title_rotation)
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

   if (length(group.by) != length(row_title_expression)){
     stop("Please provide the same number of row titles as the number of items in group.by.", call. = FALSE)
   }

  if (!(scale_direction %in% c(1, -1))){
    stop("Please provide either -1 or 1 to scale_direction.", call. = FALSE)
  }


  magnitude <- ifelse(slot == "data", "avg_log2FC", "avg_diff")
  # Compute the top N genes per cluster.
  genes.use <- de_genes %>%
               dplyr::arrange(.data$p_val_adj, dplyr::desc(.data[[magnitude]])) %>%
               dplyr::group_by(.data$cluster) %>%
               dplyr::slice_head(n = top_genes) %>%
               dplyr::pull(.data$gene) %>%
               unique()

  # Compute heatmap of log2FC.
  logfc_out <- de_genes %>%
               dplyr::arrange(.data$p_val_adj, dplyr::desc(.data[[magnitude]])) %>%
               dplyr::group_by(.data$cluster) %>%
               dplyr::slice_head(n = top_genes) %>%
               dplyr::select(.data$gene, .data$cluster, .data[[magnitude]]) %>%
               tidyr::pivot_wider(names_from = .data$gene,
                                  values_from = .data[[magnitude]]) %>%
               as.data.frame() %>%
               tibble::column_to_rownames(var = "cluster") %>%
               as.matrix() %>%
               replace(is.na(.), 0) %>%
               heatmap_inner(cluster_columns = F,
                             cluster_rows = F,
                             legend_name = ifelse(slot == "data", "Avg. log2(FC)", "Avg. Diff."),
                             data_range = "only_pos",
                             row_title = row_title_logfc,
                             row_names_side = row_names_side,
                             legend.position = legend.position,
                             legend.length = heatmap.legend.length,
                             legend.width = heatmap.legend.width,
                             legend.framecolor = heatmap.legend.framecolor,
                             use_viridis = TRUE,
                             viridis_color_map = viridis_map_logfc,
                             viridis_direction = scale_direction,
                             zeros_are_white = TRUE,
                             row_title_rotation = row_title_rotation,
                             row_title_side = row_title_side,
                             cell_size = cell_size)

  # Compute heatmap of -log10FC.
  pvalue_out <- de_genes %>%
                dplyr::arrange(.data$p_val_adj, dplyr::desc(.data[[magnitude]])) %>%
                dplyr::group_by(.data$cluster) %>%
                dplyr::slice_head(n = top_genes) %>%
                dplyr::select(.data$gene, .data$cluster, .data$p_val_adj) %>%
                dplyr::mutate("p_val_adj" = replace(.data$p_val_adj, .data$p_val_adj == 0, .Machine$double.xmin)) %>%
                dplyr::mutate("log10pval" = -log10(.data$p_val_adj)) %>%
                dplyr::select(-.data$p_val_adj) %>%
                tidyr::pivot_wider(names_from = .data$gene,
                                   values_from = .data$log10pval) %>%
                as.data.frame() %>%
                tibble::column_to_rownames(var = "cluster") %>%
                as.matrix() %>%
                replace(is.na(.), 0)

  pvalue_out <- heatmap_inner(data = pvalue_out,
                              cluster_columns = F,
                              cluster_rows = F,
                              row_names_side = row_names_side,
                              row_title = row_title_p_values,
                              legend_name = "-log10(Adjusted P-value)",
                              data_range = "only_pos",
                              legend.position = legend.position,
                              legend.length = heatmap.legend.length,
                              legend.width = heatmap.legend.width,
                              legend.framecolor = heatmap.legend.framecolor,
                              na.value = "grey75",
                              column_title = column_title,
                              outlier.data = if(sum(pvalue_out == -log10(.Machine$double.xmin)) > 0) {TRUE} else {FALSE},
                              range.data = if(sum(pvalue_out == -log10(.Machine$double.xmin)) > 0) {max(pvalue_out[pvalue_out < 307])} else {NULL},
                              outlier.up.label = "Inf",
                              use_viridis = TRUE,
                              viridis_color_map = viridis_map_pvalues,
                              viridis_direction = scale_direction,
                              zeros_are_white = TRUE,
                              row_title_rotation = row_title_rotation,
                              row_title_side = row_title_side,
                              cell_size = cell_size)

  # Compute heatmap of expression.
  list.expression.heatmaps <- list()
  list.expression.legends <- list()

  max_value_list <- c()
  # Compute the max values for all heatmaps.
  for (variable in group.by){
    max_value <- sample@meta.data %>%
                 dplyr::select(.data[[variable]]) %>%
                 tibble::rownames_to_column(var = "cell") %>%
                 dplyr::left_join(y = {Seurat::GetAssayData(object = sample,
                                                            slot = slot,
                                                            assay = assay)[genes.use, ] %>%
                                       as.matrix() %>%
                                       t() %>%
                                       as.data.frame() %>%
                                       tibble::rownames_to_column(var = "cell")},
                                       by = "cell") %>%
                 dplyr::select(-.data$cell) %>%
                 tidyr::pivot_longer(cols = -.data[[variable]],
                                     names_to = "gene",
                                     values_to = "expression") %>%
                 dplyr::group_by(.data[[variable]], .data$gene) %>%
                 dplyr::summarise(mean_expression = mean(.data$expression)) %>%
                 tidyr::pivot_wider(names_from = .data$gene,
                                    values_from = .data$mean_expression) %>%
                 as.data.frame() %>%
                 tibble::column_to_rownames(var = variable) %>%
                 dplyr::select(dplyr::all_of(genes.use)) %>%
                 as.matrix() %>%
                 max()

    max_value_list <- c(max_value_list, max_value)
  }

  counter <- 0
  for (variable in group.by){
    counter <- counter + 1
    data_range <- if(slot == "data") {"only_pos"} else if (slot == "scale.data"){"both"}
    expression_out <- sample@meta.data %>%
                      dplyr::select(.data[[variable]]) %>%
                      tibble::rownames_to_column(var = "cell") %>%
                      dplyr::left_join(y = {Seurat::GetAssayData(object = sample,
                                                                 slot = slot,
                                                                 assay = assay)[genes.use, ] %>%
                                            as.matrix() %>%
                                            t() %>%
                                            as.data.frame() %>%
                                            tibble::rownames_to_column(var = "cell")},
                                            by = "cell") %>%
                      dplyr::select(-.data$cell) %>%
                      tidyr::pivot_longer(cols = -.data[[variable]],
                                          names_to = "gene",
                                          values_to = "expression") %>%
                      dplyr::group_by(.data[[variable]], .data$gene) %>%
                      dplyr::summarise(mean_expression = mean(.data$expression)) %>%
                      tidyr::pivot_wider(names_from = .data$gene,
                                         values_from = .data$mean_expression) %>%
                      as.data.frame() %>%
                      tibble::column_to_rownames(var = variable) %>%
                      dplyr::select(dplyr::all_of(genes.use)) %>%
                      as.matrix() %>%
                      heatmap_inner(cluster_columns = F,
                                    cluster_rows = F,
                                    row_names_side = row_names_side,
                                    legend_name = "Mean expression",
                                    data_range = data_range,
                                    row_title = row_title_expression[counter],
                                    legend.position = legend.position,
                                    legend.length = heatmap.legend.length,
                                    legend.width = heatmap.legend.width,
                                    legend.framecolor = heatmap.legend.framecolor,
                                    use_viridis = TRUE,
                                    viridis_color_map = viridis_map_expression,
                                    viridis_direction = scale_direction,
                                    zeros_are_white = TRUE,
                                    range.data = max(max_value_list),
                                    row_title_rotation = row_title_rotation,
                                    row_title_side = row_title_side,
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
      ht_list = ht_list %v% heatmap
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
