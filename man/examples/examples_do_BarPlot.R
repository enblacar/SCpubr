\donttest{
  # Check Suggests.
  value <- SCpubr:::check_suggests(function_name = "do_BarPlot", passive = TRUE)

  if (isTRUE(value)){
    # Consult the full documentation in https://enblacar.github.io/SCpubr-book/

    # Define your Seurat object.
    sample <- readRDS(system.file("extdata/seurat_dataset_example.rds", package = "SCpubr"))

    # Basic bar plot, horizontal.
    p1 <- SCpubr::do_BarPlot(sample = sample,
                             group.by = "seurat_clusters",
                             legend.position = "none",
                             plot.title = "Number of cells per cluster")

    # Split by a second variable.
    sample$modified_orig.ident <- sample(x = c("Sample_A", "Sample_B", "Sample_C"),
                                         size = ncol(sample),
                                         replace = TRUE,
                                         prob = c(0.2, 0.7, 0.1))

    p <- SCpubr::do_BarPlot(sample,
                            group.by = "seurat_clusters",
                            split.by = "modified_orig.ident",
                            plot.title = "Number of cells per cluster in each sample",
                            position = "stack")


  } else if (base::isFALSE(value)){
    message("This function can not be used without its suggested packages.")
    message("Check out which ones are needed using `SCpubr::state_dependencies()`.")
  }
}
