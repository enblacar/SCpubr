do_CorrelationPlot <- function(sample,
                               mode = NULL){



  if (mode == "hvg"){
    # Generate a correlation matrix of the HVG.
    variable_genes <- Seurat::VariableFeatures(sample)

    # Subset sample according to the variable genes.
    sample.variable <- sample[variable_genes, ]


    expr_mat <- data.frame("rownames" = rownames(sample.variable))


    # Iterate over each marker gene list.
    for (celltype in levels(sample.variable)){

      # Subset only the cells for the cluster.
      subset <- sample.variable[, sample.variable$New_NMF_labelling == celltype]
      # Retrieve which cells are assigned to the cluster.
      expr_scores <- rowMeans(as.matrix(subset@assays$SCT@data))

      # Append the scores.
      expr_mat[[sprintf("%s", celltype)]] <- expr_scores
    }
    subset <- NULL
    rownames(expr_mat) <- expr_mat$rownames
    expr_mat$rownames <- NULL

    cor_mat <- round(cor(expr_mat), digits = 2)

    range <- max(abs(cor_mat))

    out <- SCpubr:::heatmap_inner(cor_mat)
    h <- out[["heatmap"]]
    h_legend <- out[["legend"]]

    ComplexHeatmap::ht_opt("HEATMAP_LEGEND_PADDING" = ggplot2::unit(8, "mm"))
    h <- ComplexHeatmap::draw(h,
                              heatmap_legend_list = h_legend,
                              padding = ggplot2::unit(c(5, 5, 5, 5), "mm"))
  }
  range <- max(abs(data.use))
  h <- pheatmap::pheatmap(data.use,
                          color = viridis::cividis(100),
                          breaks = seq(-range, range, length.out = 100),
                          cellwidth = 15,
                          cellheight = 15)


  # Generate a correlation matrix of the HVG.
  variable_genes <- Seurat::VariableFeatures(sample)
}
