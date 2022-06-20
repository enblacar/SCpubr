#' Create correlation matrix heatmaps.
#'
#'
#' @param sample Your Seurat object.
#' @param assay Assay to retrieve the data from.
#' @param mode Different types of correlation matrices can be computed. Right now, the only possible value is "hvg", standing for Highly Variable Genes. The sample is subset for the HVG and the data is re-scaled. Scale data is used for the correlation.
#' @param group.by Metadata variable to group the values by.
#' @param column_title,row_title Title for the columns and rows. Only works with single heatmaps.
#' @param cluster_cols,cluster_rows Logical. Cluster the columns or rows.
#' @param column_names_rot,row_names_rot Numeric. Degree in which to
#' @param legend_name Text for the legend title.
#' @param colors.use Vector of 2 colors to use to generate the color scale.
#' @param cell_size Numeric. Size of each cell in the heatmap.
#' @param na.value Value for NAs.
#' @return A ComplexHeatmap object.
#' @export
#'
#' @example /man/examples/examples_do_CorrelationPlot.R
do_CorrelationPlot <- function(sample,
                               mode = "hvg",
                               assay = NULL,
                               group.by = NULL,
                               column_title = NULL,
                               row_title = NULL,
                               cluster_cols = TRUE,
                               cluster_rows = TRUE,
                               legend_name = "Pearson coef.",
                               row_names_rot = 0,
                               column_names_rot = 90,
                               colors.use = NULL,
                               cell_size = 5,
                               na.value = "grey75"){

  # Checks for packages.
  check_suggests(function_name = "do_CorrelationPlot")
  # Check if the sample provided is a Seurat object.
  check_Seurat(sample = sample)
  # Check logical parameters.
  logical_list <- list("cluster_cols" = cluster_cols,
                       "cluster_rows" = cluster_rows)
  check_type(parameters = logical_list, required_type = "logical", test_function = is.logical)
  # Check numeric parameters.
  numeric_list <- list("row_names_rot" = row_names_rot,
                       "column_names_rot" = column_names_rot,
                       "cell_size" = cell_size)
  check_type(parameters = numeric_list, required_type = "numeric", test_function = is.numeric)
  # Check character parameters.
  character_list <- list("mode" = mode,
                         "assay" = assay,
                         "column_title" = column_title,
                         "row_title" = row_title,
                         "legend_name" = legend_name,
                         "colors.use" = colors.use,
                         "group.by" = group.by,
                         "na.value" = na.value)
  check_type(parameters = character_list, required_type = "character", test_function = is.character)
  check_colors(na.value)

  `%v%` <- ComplexHeatmap::`%v%`
  `%>%` <- purrr::`%>%`

  out <- check_and_set_assay(sample = sample, assay = assay)
  sample <- out[["sample"]]
  assay <- out[["assay"]]

  if (mode == "hvg"){
    if (is.null(group.by)){
      aggr_entities <- levels(sample)
      sample@meta.data[, "dummy"] <- sample@active.ident
      group.by <- "dummy"
    } else {
      if (is.factor(sample@meta.data[, group.by])){
        aggr_entities <- levels(sample@meta.data[, group.by])
      } else {
        aggr_entities <- sort(unique(sample@meta.data[, group.by]))
      }
    }
    # Generate a correlation matrix of the HVG.
    variable_genes <- Seurat::VariableFeatures(sample)

    # Subset sample according to the variable genes.
    sample.variable <- sample[variable_genes, ]
    # Scale the data
    sample.variable <- Seurat::ScaleData(sample.variable, verbose = F)

    expr_mat <- data.frame("rownames" = rownames(sample.variable))


    # Iterate over each marker gene list.
    for (celltype in aggr_entities){

      # Subset only the cells for the cluster.
      subset <- subset(sample.variable, !!rlang::sym(group.by) == celltype)
      # Retrieve which cells are assigned to the cluster.
      expr_scores <- rowMeans(as.matrix(subset@assays[[assay]]@scale.data))

      # Append the scores.
      expr_mat[[celltype]] <- expr_scores
    }
    subset <- NULL
    rownames(expr_mat) <- expr_mat$rownames
    expr_mat$rownames <- NULL

    cor_mat <- round(stats::cor(expr_mat), digits = 2)

    range <- max(abs(cor_mat))

    row_title <- {
      if (!(is.null(row_title))){
        row_title
      } else {
        ""
      }}
    column_title <- {
      if (!(is.null(column_title))){
        column_title
      } else {
        ""
      }}
    out <- heatmap_inner(cor_mat,
                         legend_name = legend_name,
                         column_title = column_title,
                         row_title = row_title,
                         row_title_side = "right",
                         row_title_rotation = 0,
                         column_names_rot = column_names_rot,
                         row_names_rot = row_names_rot,
                         range.data = range,
                         cluster_columns = cluster_cols,
                         cluster_rows = cluster_rows,
                         colors.use = colors.use,
                         cell_size = cell_size,
                         na.value = na.value)
    h <- out[["heatmap"]]
    h_legend <- out[["legend"]]
    ComplexHeatmap::ht_opt("HEATMAP_LEGEND_PADDING" = ggplot2::unit(8, "mm"))
    suppressWarnings({
      grDevices::pdf(NULL)
      h <- ComplexHeatmap::draw(h,
                                heatmap_legend_list = h_legend,
                                padding = ggplot2::unit(c(5, 5, 5, 5), "mm"))
      grDevices::dev.off()
    })

  }
  return(h)
}
