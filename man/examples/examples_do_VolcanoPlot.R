\donttest{
  # Check Suggests.
  value <- SCpubr:::check_suggests(function_name = "do_VolcanoPlot", passive = TRUE)

  if (isTRUE(value)){
    # Consult the full documentation in https://enblacar.github.io/SCpubr-book/

    # Define your Seurat object.
    sample <- readRDS(system.file("extdata/seurat_dataset_example.rds", package = "SCpubr"))

    # Retrieve DE genes.
    de_genes <- readRDS(system.file("extdata/de_genes_example.rds", package = "SCpubr"))

    # Generate a volcano plot.
    p <- SCpubr::do_VolcanoPlot(sample = sample,
                                de_genes = de_genes)
    p

  } else if (isFALSE(value)){
    message("This function can not be used without its suggested packages.")
    message("Check out which ones are needed using `SCpubr::state_dependencies()`.")
  }
}
