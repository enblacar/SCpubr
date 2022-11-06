\donttest{
  # Check Suggests.
  value <- SCpubr:::check_suggests(function_name = "do_SankeyPlot", passive = TRUE)

  if (isTRUE(value)){
    library(dplyr)
    # Define your Seurat object.
    sample <- readRDS(system.file("extdata/seurat_dataset_example.rds", package = "SCpubr"))

    # Compute basic sankey plot.
    p1 <- SCpubr::do_SankeyPlot(sample = sample,
                                first_group = "seurat_clusters",
                                last_group = "orig.ident",
                                type = "sankey")

    # Compute basic alluvial plot.
    p2 <- SCpubr::do_SankeyPlot(sample = sample,
                                first_group = "seurat_clusters",
                                last_group = "orig.ident",
                                type = "alluvial")

    p <- p1 / p2
    p

    sample$assignment <- ifelse(sample$seurat_clusters %in% c("0", "2", "4"), "A", "B")

    # Add more groups.
    p1 <- SCpubr::do_SankeyPlot(sample = sample,
                                first_group = "seurat_clusters",
                                middle_groups = c("assignment"),
                                last_group = "orig.ident",
                                type = "sankey")

    p2 <- SCpubr::do_SankeyPlot(sample = sample,
                                first_group = "seurat_clusters",
                                middle_groups = c("assignment"),
                                last_group = "orig.ident",
                                type = "alluvial")

    p <- p1 / p2
    p
  } else if (isFALSE(value)){
    message("This function can not be used without its suggested packages.")
    message("Check out which ones are needed using `SCpubr::state_dependencies()`.")
  }
}
