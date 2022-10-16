\donttest{
  # Check Suggests.
  value <- SCpubr:::check_suggests(function_name = "do_BarPlot", passive = TRUE)

  if (isTRUE(value)){
    # Define your Seurat object.
    sample <- readRDS(system.file("extdata/seurat_dataset_example.rds", package = "SCpubr"))

    # Basic bar plot, horizontal.
    p1 <- SCpubr::do_BarPlot(sample = sample,
                             group.by = "seurat_clusters",
                             legend.position = "none",
                             plot.title = "Number of cells per cluster")

    # Basic bar plot, vertical.
    p2 <- SCpubr::do_BarPlot(sample = sample,
                             group.by = "seurat_clusters",
                             legend.position = "none",
                             plot.title = "Number of cells per cluster",
                             flip = TRUE)
    p <- p1 | p2
    p

    # Split by a second variable.
    sample$modified_orig.ident <- sample(x = c("Sample_A", "Sample_B", "Sample_C"),
                                         size = ncol(sample),
                                         replace = TRUE,
                                         prob = c(0.2, 0.7, 0.1))

    p1 <- SCpubr::do_BarPlot(sample,
                             group.by = "seurat_clusters",
                             split.by = "modified_orig.ident",
                             plot.title = "Number of cells per cluster in each sample",
                             position = "stack")

    p2 <- SCpubr::do_BarPlot(sample,
                             group.by = "modified_orig.ident",
                             split.by = "seurat_clusters",
                             plot.title = "Number of cells per sample in each cluster",
                             position = "stack")
    p <- p1 | p2
    p

    # Position stack and fill with and without split.by.
    p1 <- SCpubr::do_BarPlot(sample,
                             group.by = "seurat_clusters",
                             plot.title = "Without split.by - position = stack",
                             position = "stack",
                             flip = FALSE)

    p2 <- SCpubr::do_BarPlot(sample,
                             group.by = "seurat_clusters",
                             plot.title = "Without split.by - position = fill",
                             position = "fill",
                             flip = FALSE)

    p3 <- SCpubr::do_BarPlot(sample,
                             group.by = "seurat_clusters",
                             split.by = "modified_orig.ident",
                             plot.title = "With split.by - position = stack",
                             position = "stack",
                             flip = FALSE)

    p4 <- SCpubr::do_BarPlot(sample,
                             group.by = "seurat_clusters",
                             split.by = "modified_orig.ident",
                             plot.title = "With split.by - position = fill",
                             position = "fill",
                             flip = FALSE)
    p <- (p1 | p2) / (p3 | p4)
    p
  } else if (isFALSE(value)){
    message("This function can not be used without its suggested packages.")
    message("Check out which ones are needed using `SCpubr::state_dependencies()`.")
  }
}
