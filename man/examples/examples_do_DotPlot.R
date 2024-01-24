\donttest{
  # Check Suggests.
  value <- SCpubr:::check_suggests(function_name = "do_DotPlot", passive = TRUE)

  if (isTRUE(value)){
    # Consult the full documentation in https://enblacar.github.io/SCpubr-book/
    
    # Define your Seurat object.
    sample <- readRDS(system.file("extdata/seurat_dataset_example.rds", package = "SCpubr"))

    # Basic Dot plot.
    genes <- rownames(sample)[1:14]
    p <- SCpubr::do_DotPlot(sample = sample,
                            features = genes)

  } else if (base::isFALSE(value)){
    message("This function can not be used without its suggested packages.")
    message("Check out which ones are needed using `SCpubr::state_dependencies()`.")
  }
}
