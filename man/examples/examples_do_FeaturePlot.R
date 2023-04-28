\donttest{
  # Check Suggests.
  value <- SCpubr:::check_suggests(function_name = "do_FeaturePlot", passive = TRUE)

  if (isTRUE(value)){
    # Consult the full documentation in https://enblacar.github.io/SCpubr-book/

    # Define your Seurat object.
    sample <- readRDS(system.file("extdata/seurat_dataset_example.rds", package = "SCpubr"))

    # Regular FeaturePlot.
    p <- SCpubr::do_FeaturePlot(sample = sample,
                                features = "nCount_RNA")

    # FeaturePlot with a subset of identities
    # (in Seurat::Idents(sample)) maintaining the original UMAP shape.
    idents.use <- levels(sample)[!(levels(sample) %in% c("2", "5", "8"))]
    p <- SCpubr::do_FeaturePlot(sample = sample,
                                idents.highlight = idents.use,
                                features = c("EPC1"))

    # Splitting the FeaturePlot by a variable and
    # maintaining the color scale and the UMAP shape.
    p <- SCpubr::do_FeaturePlot(sample = sample,
                                features = "EPC1",
                                split.by = "seurat_clusters")

  } else if (base::isFALSE(value)){
    message("This function can not be used without its suggested packages.")
    message("Check out which ones are needed using `SCpubr::state_dependencies()`.")
  }
}
