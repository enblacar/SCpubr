\donttest{
  # Check Suggests.
  value <- SCpubr:::check_suggests(function_name = "do_CopyNumberVariantPlot", passive = TRUE)

  if (isTRUE(value)){
    # Consult the full documentation in https://enblacar.github.io/SCpubr-book/

    # This function expects that you have run inferCNV on your
    # own and you have access to the output object.

    # Define your Seurat object.
    sample <- readRDS(system.file("extdata/seurat_dataset_example.rds",
                                  package = "SCpubr"))

    # Define your inferCNV object.
    infercnv_object <- readRDS(system.file("extdata/infercnv_object_example.rds",
                                           package = "SCpubr"))

    # Get human chromosome locations.
    chromosome_locations = SCpubr::human_chr_locations

    # Compute for a single chromosome.
    out <- SCpubr::do_CopyNumberVariantPlot(sample = sample,
                                            infercnv_object = infercnv_object,
                                            using_metacells = FALSE,
                                            chromosome_locations = chromosome_locations,
                                            chromosome_focus = "1")
    # Retrieve the UMAP for 1p region.
    out$`1p_umap`
    # Retrieve the dot plot for 1p region.
    out$`1p_dotplot`

  } else if (isFALSE(value)){
    message("This function can not be used without its suggested packages.")
    message("Check out which ones are needed using `SCpubr::state_dependencies()`.")
  }
}
