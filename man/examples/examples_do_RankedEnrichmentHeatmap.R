\donttest{
  # Check Suggests.
  value <- SCpubr:::check_suggests(function_name = "do_RankedEnrichmentHeatmap", passive = TRUE)
  
  if (isTRUE(value)){
    # Consult the full documentation in https://enblacar.github.io/SCpubr-book/
    
    # Define your Seurat object.
    sample <- readRDS(system.file("extdata/seurat_dataset_example.rds", package = "SCpubr"))
    
    # Genes have to be unique.
    genes <- list("Gene set A" = rownames(sample)[1:5],
                  "Gene set B" = rownames(sample)[6:10],
                  "Gene set C" = rownames(sample)[11:15])
    
    
    # This will query, for the provided components, the enrichment of the gene 
    # sets for all cells and plot them in the context of the cells reordered by 
    # the position alongside each dimensional reduction component. 
    p <- SCpubr::do_RankedEnrichmentHeatmap(sample = sample,
                                            input_gene_list = genes,
                                            nbin = 1,
                                            ctrl = 5,
                                            flavor = "Seurat",
                                            subsample = NA,
                                            dims = 1:2,
                                            verbose = FALSE)
    
    p
    
  } else if (base::isFALSE(value)){
    message("This function can not be used without its suggested packages.")
    message("Check out which ones are needed using `SCpubr::state_dependencies()`.")
  }
}
