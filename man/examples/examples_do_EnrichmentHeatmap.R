\donttest{
  # Check Suggests.
  value <- SCpubr:::check_suggests(function_name = "do_EnrichmentHeatmap", passive = TRUE)

  if (isTRUE(value)){
    # Define your Seurat object.
    sample <- readRDS(system.file("extdata/seurat_dataset_example.rds", package = "SCpubr"))

    # Genes have to be unique.
    genes <- list("Naive CD4+ T" = rownames(sample)[1:2],
                  "EPC1+ Mono" = rownames(sample)[3:4],
                  "Memory CD4+" = rownames(sample)[5],
                  "B" = rownames(sample)[6],
                  "CD8+ T" = rownames(sample)[7],
                  "FCGR3A+ Mono" = rownames(sample)[8:9],
                  "NK" = rownames(sample)[10:11],
                  "DC" = rownames(sample)[12:13],
                  "Platelet" = rownames(sample)[14])

    # Default parameters.
    p <- SCpubr::do_EnrichmentHeatmap(sample = sample,
                                      input_gene_list = genes,
                                      nbin = 1,
                                      ctrl = 10)
    p

    # Custom aggregated values.
    p <- SCpubr::do_EnrichmentHeatmap(sample = sample,
                                      input_gene_list = genes,
                                      group.by = "orig.ident",
                                      nbin = 1,
                                      ctrl = 10)
    p

    # Transposing the matrix.
    p <- SCpubr::do_EnrichmentHeatmap(sample = sample,
                                      input_gene_list = genes,
                                      flip = TRUE,
                                      nbin = 1,
                                      ctrl = 10)
    p

    # Rotating the labels.
    p <- SCpubr::do_EnrichmentHeatmap(sample = sample,
                                      input_gene_list = genes,
                                      flip = TRUE,
                                      column_names_rot = 0,
                                      nbin = 1,
                                      ctrl = 10)
    p

    # Modifying the tile size.
    p <- SCpubr::do_EnrichmentHeatmap(sample = sample,
                                      input_gene_list = genes,
                                      flip = TRUE,
                                      column_names_rot = 0,
                                      cluster_cols = FALSE,
                                      cluster_rows = TRUE,
                                      cell_size = 7,
                                      nbin = 1,
                                      ctrl = 10)
    p


    # Symmetrical scale viriis.
    p <- SCpubr::do_EnrichmentHeatmap(sample = sample,
                                      input_gene_list = genes,
                                      flip = TRUE,
                                      column_names_rot = 0,
                                      cluster_cols = FALSE,
                                      cluster_rows = TRUE,
                                      cell_size = 7,
                                      symmetrical_scale = TRUE,
                                      nbin = 1,
                                      ctrl = 10)
    p


    # Modifying the symmetrical scale non viridis.
    p <- SCpubr::do_EnrichmentHeatmap(sample = sample,
                                      input_gene_list = genes,
                                      flip = TRUE,
                                      column_names_rot = 0,
                                      cluster_cols = FALSE,
                                      cluster_rows = TRUE,
                                      cell_size = 7,
                                      symmetrical_scale = TRUE,
                                      use_viridis = FALSE,
                                      nbin = 1,
                                      ctrl = 10)
    p
  } else if (isFALSE(value)){
    message("This function can not be used without its suggested packages.")
    message("Check out which ones are needed using `SCpubr::state_dependencies()`.")
  }
}
