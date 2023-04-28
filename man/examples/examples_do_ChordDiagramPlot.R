\donttest{
  # Check Suggests.
  value <- SCpubr:::check_suggests(function_name = "do_ChordDiagramPlot", passive = TRUE)

  if (isTRUE(value)){
    # Consult the full documentation in https://enblacar.github.io/SCpubr-book/

    # Define your Seurat object.
    sample <- readRDS(system.file("extdata/seurat_dataset_example.rds", package = "SCpubr"))

    # Basic chord diagram.
    sample$assignment <- ifelse(sample$seurat_clusters %in% c("0", "4", "7"), "A", "B")
    sample$assignment[sample$seurat_clusters %in% c("1", "2")] <- "C"
    sample$assignment[sample$seurat_clusters %in% c("10", "5")] <- "D"
    sample$assignment[sample$seurat_clusters %in% c("8", "9")] <- "E"

    p <- SCpubr::do_ChordDiagramPlot(sample = sample,
                                     from = "seurat_clusters",
                                     to = "assignment")

    p

  } else if (base::isFALSE(value)){
    message("This function can not be used without its suggested packages.")
    message("Check out which ones are needed using `SCpubr::state_dependencies()`.")
  }
}

