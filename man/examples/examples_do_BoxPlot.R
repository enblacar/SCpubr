\donttest{
  # Check Suggests.
  value <- SCpubr:::check_suggests(function_name = "do_BoxPlot", passive = TRUE)

  if (isTRUE(value)){
    # Consult the full documentation in https://enblacar.github.io/SCpubr-book/

    # Define your Seurat object.
    sample <- readRDS(system.file("extdata/seurat_dataset_example.rds", package = "SCpubr"))

    # Basic box plot.
    p <- SCpubr::do_BoxPlot(sample = sample,
                            feature = "nCount_RNA")
    p

    # Use silhouette style.
    p <- SCpubr::do_BoxPlot(sample = sample,
                            feature = "nCount_RNA",
                            use_silhouette = TRUE)
    p

    # Order by mean values.
    p <- SCpubr::do_BoxPlot(sample = sample,
                            feature = "nCount_RNA",
                            order = TRUE)
    p

    # Apply second grouping.
    sample$orig.ident <- ifelse(sample$seurat_clusters %in% c("0", "1", "2", "3"), "A", "B")
    p <- SCpubr::do_BoxPlot(sample = sample,
                            feature = "nCount_RNA",
                            split.by = "orig.ident")
    p

    # Apply statistical tests.
    p <- SCpubr::do_BoxPlot(sample = sample,
                            feature = "nCount_RNA",
                            group.by = "orig.ident",
                            use_test = TRUE,
                            comparisons = list(c("A", "B")))
    p

  } else if (base::isFALSE(value)){
    message("This function can not be used without its suggested packages.")
    message("Check out which ones are needed using `SCpubr::state_dependencies()`.")
  }
}
