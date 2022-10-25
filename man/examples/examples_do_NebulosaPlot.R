\donttest{
  # Check Suggests.
  value <- SCpubr:::check_suggests(function_name = "do_NebulosaPlot", passive = TRUE)

  if (isTRUE(value)){
    # Consult the full documentation in https://enblacar.github.io/SCpubr-book/

    # Define your Seurat object.
    sample <- readRDS(system.file("extdata/seurat_dataset_example.rds", package = "SCpubr"))

    # Basic Nebulosa plot.
    p <- SCpubr::do_NebulosaPlot(sample = sample,
                                 features = "EPC1")

    # Compute joint density.
    p <- SCpubr::do_NebulosaPlot(sample = sample,
                                 features = c("EPC1", "TOX2"),
                                 joint = TRUE)

  } else if (isFALSE(value)){
    message("This function can not be used without its suggested packages.")
    message("Check out which ones are needed using `SCpubr::state_dependencies()`.")
  }
}
