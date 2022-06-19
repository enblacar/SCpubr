\dontrun{
  # Define list of genes.
  genes <- list("Naive CD4+ T" = c("IL7R", "CCR7"),
                "CD14+ Mono" = c("CD14", "LYZ"),
                "Memory CD4+" = c("S100A4"),
                "B" = c("MS4A1"),
                "CD8+ T" = c("CD8A"),
                "FCGR3A+ Mono" = c("FCGR3A", "MS4A7"),
                "NK" = c("GNLY", "NKG7"),
                "DC" = c("FCER1A", "CST3"),
                "Platelet" = c("PPBP"))

  # Default parameters.
  p <- SCpubr::do_EnrichmentHeatmap(sample = sample,
                                    list_genes = genes)
  p

  # Custom aggregated values.
  p <- SCpubr::do_EnrichmentHeatmap(sample = sample,
                                    list_genes = genes,
                                    group.by = "orig.ident")
  p

  # Transposing the matrix.
  p <- SCpubr::do_EnrichmentHeatmap(sample = sample,
                                    list_genes = genes,
                                    transpose = TRUE)
  p

  # Rotating the labels.
  p <- SCpubr::do_EnrichmentHeatmap(sample = sample,
                                    list_genes = genes,
                                    transpose = TRUE,
                                    column_names_rot = 0)
  p


  # Splitting by a second metadata variable.
  clusters <- c("1", "3", "5", "7", "9")
  sample$custom_group <- ifelse(sample$seurat_clusters %in% clusters,
                                "Group A",
                                "Group B")
  p <- SCpubr::do_EnrichmentHeatmap(sample = sample,
                                    list_genes = genes,
                                    transpose = TRUE,
                                    column_names_rot = 0,
                                    cluster_cols = F,
                                    cluster_rows = T,
                                    split.by = "custom_group")
  p

  # Splitting by a second metadata variable and joining vertically.
  p <- SCpubr::do_EnrichmentHeatmap(sample = sample,
                                    list_genes = genes,
                                    transpose = TRUE,
                                    column_names_rot = 0,
                                    cluster_cols = F,
                                    cluster_rows = T,
                                    split.by = "custom_group",
                                    split.horizontal = F)
  p

  # Modifying the tile size.
  p <- SCpubr::do_EnrichmentHeatmap(sample = sample,
                                    list_genes = genes,
                                    transpose = TRUE,
                                    column_names_rot = 0,
                                    cluster_cols = F,
                                    cluster_rows = T,
                                    split.by = "custom_group",
                                    cell_size = 7)
  p

  # Modifying the color scale.
  p <- SCpubr::do_EnrichmentHeatmap(sample = sample,
                                    list_genes = genes,
                                    transpose = TRUE,
                                    column_names_rot = 0,
                                    cluster_cols = F,
                                    cluster_rows = T,
                                    split.by = "custom_group",
                                    colors.use = colortools::opposite("steelblue", plot = F))
  p
}
