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
                                    input_gene_list = genes)
  p

  # Custom aggregated values.
  p <- SCpubr::do_EnrichmentHeatmap(sample = sample,
                                    input_gene_list = genes,
                                    group.by = "orig.ident")
  p

  # Transposing the matrix.
  p <- SCpubr::do_EnrichmentHeatmap(sample = sample,
                                    input_gene_list = genes,
                                    transpose = TRUE)
  p

  # Rotating the labels.
  p <- SCpubr::do_EnrichmentHeatmap(sample = sample,
                                    input_gene_list = genes,
                                    transpose = TRUE,
                                    column_names_rot = 0)
  p

  # Modifying the tile size.
  p <- SCpubr::do_EnrichmentHeatmap(sample = sample,
                                    input_gene_list = genes,
                                    transpose = TRUE,
                                    column_names_rot = 0,
                                    cluster_cols = F,
                                    cluster_rows = T,
                                    split.by = "custom_group",
                                    cell_size = 7)
  p


  # Symmetrical scale viriis.
  p <- SCpubr::do_EnrichmentHeatmap(sample = sample,
                                    input_gene_list = genes,
                                    transpose = TRUE,
                                    column_names_rot = 0,
                                    cluster_cols = F,
                                    cluster_rows = T,
                                    cell_size = 7,
                                    symmetrical_scale = TRUE)
  p


  # Modifying the symmetrical scale non viridis.
  p <- SCpubr::do_EnrichmentHeatmap(sample = sample,
                                    input_gene_list = genes,
                                    transpose = TRUE,
                                    column_names_rot = 0,
                                    cluster_cols = F,
                                    cluster_rows = T,
                                    cell_size = 7,
                                    symmetrical_scale = TRUE,
                                    use_viridis = FALSE)
  p
}
