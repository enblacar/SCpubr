\donttest{
  # Check Suggests.
  value <- SCpubr:::check_suggests(function_name = "do_GeyserPlot", passive = TRUE)

  if (isTRUE(value)){
    # Consult the full documentation in https://enblacar.github.io/SCpubr-book/

    # Define your Seurat object.
    sample <- readRDS(system.file("extdata/seurat_dataset_example.rds", package = "SCpubr"))

    # Geyser plot with categorical color scale.
    p <- SCpubr::do_GeyserPlot(sample = sample,
                              features = "nCount_RNA",
                              scale_type = "categorical")
    p

    # Geyser plot with continuous color scale.
    p <- SCpubr::do_GeyserPlot(sample = sample,
                              features = "nCount_RNA",
                              scale_type = "continuous")


    p

  } else if (base::isFALSE(value)){
    message("This function can not be used without its suggested packages.")
    message("Check out which ones are needed using `SCpubr::state_dependencies()`.")
  }
}
