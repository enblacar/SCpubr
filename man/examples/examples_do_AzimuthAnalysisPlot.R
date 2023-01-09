\donttest{
  # Check Suggests.
  value <- SCpubr:::check_suggests(function_name = "do_AzimuthAnalysisPlot", passive = TRUE)

  if (isTRUE(value)){
    # Consult the full documentation in https://enblacar.github.io/SCpubr-book/

    # Define your Seurat object that has undergone Azimuth mapping.
    sample <- readRDS(system.file("extdata/seurat_dataset_example.rds", package = "SCpubr"))

    # Generate an Azimuth report.
    out <- SCpubr::do_AzimuthAnalysisPlot(sample = sample,
                                          annotation.labels = "annotation",
                                          annotation.scoring = "annotation.score",
                                          font.size = 18)

  } else if (isFALSE(value)){
    message("This function can not be used without its suggested packages.")
    message("Check out which ones are needed using `SCpubr::state_dependencies()`.")
  }
}
