\donttest{
  # Check Suggests.
  value <- SCpubr:::check_suggests(function_name = "do_RankedExpressionPlot", passive = TRUE)
  
  if (isTRUE(value)){
    # Consult the full documentation in https://enblacar.github.io/SCpubr-book/
    
    # Define your Seurat object.
    sample <- readRDS(system.file("extdata/seurat_dataset_example.rds", package = "SCpubr"))
    
    # Genes have to be unique.
    genes <- list("A" = rownames(sample)[1:5],
                  "B" = rownames(sample)[6:10],
                  "C" = rownames(sample)[11:15])
    
    
    # This will query, for the provided components, the expression of the genes 
    # for all cells and plot them in the context of the cells reordered by 
    # the position alongside each dimensional reduction component. 
    p <- SCpubr::do_RankedExpressionPlot(sample = sample,
                                         features = genes,
                                         nbin = 1,
                                         ctrl = 5,
                                         subsample = NA,
                                         dims = 1:2,
                                         verbose = FALSE)
    
    p
    
  } else if (base::isFALSE(value)){
    message("This function can not be used without its suggested packages.")
    message("Check out which ones are needed using `SCpubr::state_dependencies()`.")
  }
}
