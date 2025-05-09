\dontrun{
  # Check Suggests.
  value <- SCpubr:::check_suggests(function_name = "do_SavePlot", passive = TRUE)

  if (isTRUE(value)){
    # Define your Seurat object.
    sample <- readRDS(system.file("extdata/seurat_dataset_example.rds", package = "SCpubr"))

    # Generate a plot.
    p <- SCpubr::do_DimPlot(sample = sample)

    # Default parameters.
    SCpubr::do_SavePlot(plot = p)

    # Specifying the name and folder.
    SCpubr::do_SavePlot(plot = p,
                        figure_path = "/path/to/my/figures/",
                        file_name = "my_figure")

    # Specify to also create a new folder.
    SCpubr::do_SavePlot(plot = p,
                        figure_path = "/path/to/my/figures/",
                        file_name = "my_figure",
                        create_path = TRUE)

    # Set dimensions for the figure.
    SCpubr::do_SavePlot(plot = p,
                        figure_path = "/path/to/my/figures/",
                        file_name = "my_figure",
                        create_path = TRUE,
                        width = 8,
                        height = 8)

    # Set quality (dpi).
    SCpubr::do_SavePlot(plot = p,
                        figure_path = "/path/to/my/figures/",
                        file_name = "my_figure",
                        create_path = TRUE,
                        width = 8,
                        height = 8,
                        dpi = 300)
  } else if (base::isFALSE(value)){
    message("This function can not be used without its suggested packages.")
    message("Check out which ones are needed using `SCpubr::state_dependencies()`.")
  }
}

