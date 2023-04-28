\donttest{
  # Check Suggests.
  value <- SCpubr:::check_suggests(function_name = "do_PathwayActivityPlot", passive = TRUE)

  if (isTRUE(value)){
    # Consult the full documentation in https://enblacar.github.io/SCpubr-book/

    # Define your Seurat object.
    sample <- readRDS(system.file("extdata/seurat_dataset_example.rds",
                                  package = "SCpubr"))

    # Define your activities object.
    progeny_activities <- readRDS(system.file("extdata/progeny_activities_example.rds",
                                              package = "SCpubr"))

    # General heatmap.
    out <- SCpubr::do_PathwayActivityPlot(sample = sample,
                                          activities = progeny_activities)
    p <- out$heatmaps$average_scores
    p

  } else if (base::isFALSE(value)){
    message("This function can not be used without its suggested packages.")
    message("Check out which ones are needed using `SCpubr::state_dependencies()`.")
  }
}
