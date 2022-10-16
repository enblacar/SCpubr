\donttest{
  # Check Suggests.
  value <- SCpubr:::check_suggests(function_name = "do_BoxPlot", passive = TRUE)

  if (isTRUE(value)){
    # Define your Seurat object.
    sample <- readRDS(system.file("extdata/seurat_dataset_example.rds", package = "SCpubr"))

    # Basic box plot.
    p <- SCpubr::do_BoxPlot(sample = sample,
                            feature = "nCount_RNA")
    p

    # Generate a custom group.
    sample$custom_group = ifelse(colnames(sample) %in% sample(colnames(sample), 50), "A", "B")

    # Use custom grouping.
    p <- SCpubr::do_BoxPlot(sample = sample,
                            feature = "nCount_RNA",
                            group.by = "custom_group")
    p

    # Flip the box plot.
    p <- SCpubr::do_BoxPlot(sample = sample,
                            feature = "nCount_RNA",
                            flip = TRUE)
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
    p <- SCpubr::do_BoxPlot(sample = sample,
                            feature = "nCount_RNA",
                            split.by = "custom_group")
    p

    # Apply statistical tests.
    sample$orig.ident <- ifelse(sample$seurat_clusters %in% c("0", "1", "2", "3"), "A", "B")
    p <- SCpubr::do_BoxPlot(sample = sample,
                            feature = "nCount_RNA",
                            group.by = "orig.ident",
                            use_test = TRUE,
                            comparisons = list(c("A", "B")))
    p

    # Apply statistical tests and show the p-value.
    p <- SCpubr::do_BoxPlot(sample = sample,
                            feature = "nCount_RNA",
                            group.by = "orig.ident",
                            use_test = TRUE,
                            comparisons = list(c("A", "B")),
                            map_signif_level = FALSE)
    p
  } else if (isFALSE(value)){
    message("This function can not be used without its suggested packages.")
    message("Check out which ones are needed using `SCpubr::state_dependencies()`.")
  }
}
