\donttest{
  # Check Suggests.
  value <- SCpubr:::check_suggests(function_name = "do_ExpressionHeatmap", passive = TRUE)

  if (isTRUE(value)){
    # Consult the full documentation in https://enblacar.github.io/SCpubr-book/

    # Define your Seurat object.
    sample <- readRDS(system.file("extdata/seurat_dataset_example.rds", package = "SCpubr"))

    # Define list of genes.
    genes <- rownames(sample)[1:10]

    # Default parameters.
    p <- SCpubr::do_ExpressionHeatmap(sample = sample,
                                      features = genes,
                                      viridis.direction = -1)
    
    p
  } else if (base::isFALSE(value)){
    message("This function can not be used without its suggested packages.")
    message("Check out which ones are needed using `SCpubr::state_dependencies()`.")
  }
}
