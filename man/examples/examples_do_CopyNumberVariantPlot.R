\donttest{
  # Check Suggests.
  value <- SCpubr:::check_suggests(function_name = "do_CopyNumberVariantPlot", passive = TRUE)

  if (isTRUE(value)){
    # This function expects that you have run inferCNV on your
    # own and you have access to the output object.

    # Define your Seurat object.
    sample <- readRDS(system.file("extdata/seurat_dataset_example.rds",
                                  package = "SCpubr"))

    # Define your inferCNV object.
    infercnv_object <- readRDS(system.file("extdata/infercnv_object_example.rds",
                                           package = "SCpubr"))
    file_name <- "extdata/infercnv_object_metacells_example.rds"
    infercnv_object_metacells <- readRDS(system.file(file_name,
                                                     package = "SCpubr"))

    # Get human chromosome locations.
    chromosome_locations = SCpubr::human_chr_locations

    # Get metacell mapping if used metacells with inferCNV.
    metacell_mapping <- readRDS(system.file("extdata/metacell_mapping_example.rds",
                                            package = "SCpubr"))

    # Using a normal run - compute for all chromosome arms.
    out <- SCpubr::do_CopyNumberVariantPlot(sample = sample,
                                            infercnv_object = infercnv_object,
                                            using_metacells = FALSE,
                                            chromosome_locations = chromosome_locations)
    # Retrieve the UMAP for 1p region.
    out$`1p_umap`
    # Retrieve the dot plot for 1p region.
    out$`1p_dotplot`

    # Compute for a single chromosoome.
    out <- SCpubr::do_CopyNumberVariantPlot(sample = sample,
                                            infercnv_object = infercnv_object,
                                            using_metacells = FALSE,
                                            chromosome_focus = "2",
                                            chromosome_locations = chromosome_locations)

    # Compute for a run using metacells.
    ## metacell_mapping has to be a vector that, for each cell,
    ## tells which is the corresponding metacell.
    out <- SCpubr::do_CopyNumberVariantPlot(sample = sample,
                                            infercnv_object = infercnv_object_metacells,
                                            using_metacells = TRUE,
                                            metacell_mapping = metacell_mapping,
                                            chromosome_locations = chromosome_locations)
  } else if (isFALSE(value)){
    message("This function can not be used without its suggested packages.")
    message("Check out which ones are needed using `SCpubr::state_dependencies()`.")
  }
}
