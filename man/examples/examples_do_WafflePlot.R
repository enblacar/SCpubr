\donttest{
  # Check Suggests.
  value <- SCpubr:::check_suggests(function_name = "do_WafflePlot", passive = TRUE)
  
  if (isTRUE(value)){
    # Consult the full documentation in https://enblacar.github.io/SCpubr-book/
    
    # Define your Seurat object.
    sample <- readRDS(system.file("extdata/seurat_dataset_example.rds", package = "SCpubr"))
    
    # Basic Waffle plot.
    p <- SCpubr::do_WafflePlot(sample = sample,
                               group.by = "seurat_clusters")
    
    
  } else if (base::isFALSE(value)){
    message("This function can not be used without its suggested packages.")
    message("Check out which ones are needed using `SCpubr::state_dependencies()`.")
  }
}
