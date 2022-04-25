do_EnrichmentHeatmap <- function(sample,
                                 list_genes,
                                 group.by,
                                 column_title = NULL,
                                 row_title = NULL,
                                 verbose = FALSE){



  if (is.character(list_genes)){
    # If list_genes is a character of genes.
    input_list <- list_genes
    names(input_list) <- list_genes
  } else if (is.list(list_genes)){
    input_list <- list_genes
    if (is.null(names(input_list))){
      stop("Please provide a named list. This is, each gene list has to come with a name.", call. = F)
    }
  }


  # Start the process.
  if (is.factor(sample@meta.data[, group.by])){
    aggr_entities <- levels(sample@meta.data[, group.by])
  } else {
    aggr_entities <- sort(unique(sample@meta.data[, group.by]))
  }
  scoring <- data.frame("rownames" = aggr_entities)


  # Iterate over each marker gene list.
  for (celltype in names(input_list)){
    list_markers <- list(input_list[[celltype]])

    # Compute Seurat AddModuleScore as well.

    test_out <- NA
    control_genes <- 100
    number_bin <- 24

    if (verbose){
      sample <- Seurat::AddModuleScore(sample,
                                       list_markers,
                                       name = celltype,
                                       search = TRUE,
                                       ctrl = control_genes,
                                       nbin = number_bin,
                                       verbose = F)
    } else {
      sample <- suppressMessages(suppressWarnings(Seurat::AddModuleScore(sample,
                                                                         list_markers,
                                                                         name = celltype,
                                                                         search = TRUE,
                                                                         ctrl = control_genes,
                                                                         nbin = number_bin,
                                                                         verbose = F)))
    }


    # Retrieve the scores.
    col_name <- paste0(celltype, "1")

    # Modify the name that Seurat::AddModuleScore gives by default.
    sample@meta.data[, celltype] <- sample@meta.data[, col_name]
    # Remove old metadata.
    sample@meta.data[, col_name] <- NULL

    # Generate empty vectors for the aggregated scores.
    list_score_seurat <- c()


    # Iterate over each cluster.
    for (cluster_name in aggr_entities){
      # Retrieve which cells are assigned to the cluster.
      scores_seurat <- sample@meta.data[sample@meta.data[, group.by] == cluster_name, celltype]

      # Append to the vector the mean for each cell type.
      list_score_seurat <- append(list_score_seurat, mean(scores_seurat))
    }
    # Get the name of the column together with the number of genes used for the enrichment scoring.
    scoring[celltype] <- list_score_seurat
  }
  test_out <- NULL

  # Remove the rownames column in the object, since you set them to be the rownames of the dataframe.
  rownames(scoring) <- scoring$rownames
  scoring$rownames <- NULL

  # Transform the data frames into a Matrix object.
  scoring <- as.matrix(scoring)

  range <- max(abs(scoring))

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
  out <- SCpubr:::heatmap_inner(scoring,
                                legend_name = "Enrichment",
                                column_title = column_title,
                                row_title = row_title
                                )
  h <- out[["heatmap"]]
  h_legend <- out[["legend"]]

  ComplexHeatmap::ht_opt("HEATMAP_LEGEND_PADDING" = ggplot2::unit(8, "mm"))

  grDevices::pdf(NULL)
  h <- ComplexHeatmap::draw(h,
                            heatmap_legend_list = h_legend,
                            padding = ggplot2::unit(c(5, 5, 5, 5), "mm"))
  dev.off()
  return(h)
}
