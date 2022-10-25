\donttest{
  # Check Suggests.
  value <- SCpubr:::check_suggests(function_name = "do_EnrichmentHeatmap", passive = TRUE)

  if (isTRUE(value)){
    # Consult the full documentation in https://enblacar.github.io/SCpubr-book/

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

  } else if (isFALSE(value)){
    message("This function can not be used without its suggested packages.")
    message("Check out which ones are needed using `SCpubr::state_dependencies()`.")
  }
}
