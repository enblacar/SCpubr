\donttest{
  # Check Suggests.
  value <- SCpubr:::check_suggests(function_name = "do_DimPlot", passive = TRUE)

  if (isTRUE(value)){
    # Consult the full documentation in https://enblacar.github.io/SCpubr-book/

    # Define your Seurat object.
    sample <- readRDS(system.file("extdata/seurat_dataset_example.rds", package = "SCpubr"))

    # Basic DimPlot.
    p <- SCpubr::do_DimPlot(sample = sample)


    # Restrict the amount of identities displayed.
    p <- SCpubr::do_DimPlot(sample = sample,
                            idents.keep = c("1", "3", "5"))

    # Group by another variable rather than `Seurat::Idents(sample)`
    p <- SCpubr::do_DimPlot(sample = sample,
                            group.by = "seurat_clusters")

    # Split the output in as many plots as unique identities.
    p <- SCpubr::do_DimPlot(sample = sample,
                            split.by = "seurat_clusters")

    # Highlight given identities
    p <- SCpubr::do_DimPlot(sample,
                            idents.highlight = c("1", "3"))

  } else if (base::isFALSE(value)){
    message("This function can not be used without its suggested packages.")
    message("Check out which ones are needed using `SCpubr::state_dependencies()`.")
  }
}



