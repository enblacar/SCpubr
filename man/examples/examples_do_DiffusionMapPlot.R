\donttest{
  # Check Suggests.
  value <- SCpubr:::check_suggests(function_name = "do_DiffusionMapPlot", passive = TRUE)
  
  if (isTRUE(value)){
    # Consult the full documentation in https://enblacar.github.io/SCpubr-book/
    
    # Define your Seurat object.
    sample <- readRDS(system.file("extdata/seurat_dataset_example.rds", package = "SCpubr"))
    
    # Genes have to be unique.
    genes <- list("A" = rownames(sample)[1:5],
                  "B" = rownames(sample)[6:10],
                  "C" = rownames(sample)[11:15])
    
    # Requisite is that you have a diffusion map reduction stored in the Seurat 
    # object under the name "diffusion".
    
    # This will query, for the provided components, the enrichment of the gene 
    # sets for all cells and plot them in the context of the cells reordered by 
    # the position alonside each DC. 
    p <- SCpubr::do_DiffusionMapPlot(sample = sample,
                                     input_gene_list = genes,
                                     nbin = 1,
                                     ctrl = 5,
                                     flavor = "Seurat",
                                     subsample = NA,
                                     dims = 1:2,
                                     verbose = FALSE)
    
    p
    
  } else if (isFALSE(value)){
    message("This function can not be used without its suggested packages.")
    message("Check out which ones are needed using `SCpubr::state_dependencies()`.")
  }
}
