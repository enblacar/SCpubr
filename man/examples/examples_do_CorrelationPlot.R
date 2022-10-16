\donttest{
  # Check Suggests.
  value <- SCpubr:::check_suggests(function_name = "do_CorrelationPlot", passive = TRUE)

  if (isTRUE(value)){
    # Define your Seurat object.
    sample <- readRDS(system.file("extdata/seurat_dataset_example.rds", package = "SCpubr"))

    # Default values.
    p <- SCpubr::do_CorrelationPlot(sample = sample)
    p

    # Custom grouping.
    clusters <- c("1", "3", "5", "7", "9")
    sample$custom_group <- ifelse(sample$seurat_clusters %in% clusters,
                                  "Group A",
                                  "Group B")
    p <- SCpubr::do_CorrelationPlot(sample = sample, group.by = "custom_group")
    p

    # Rotated axis labels.
    p <- SCpubr::do_CorrelationPlot(sample = sample,
                                    column_names_rot = 0)
    p

    # Increase cell size and rotate axis labels.
    p <- SCpubr::do_CorrelationPlot(sample = sample,
                                    column_names_rot = 0,
                                    cell_size = 7)
    p
  } else if (isFALSE(value)){
    message("This function can not be used without its suggested packages.")
    message("Check out which ones are needed using `SCpubr::state_dependencies()`.")
  }
}

