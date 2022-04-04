\dontrun{
  # Generate a plot.
  p <- SCpubr::do_DimPlot(sample = sample)

  # Default parameters.
  SCpubr::savePlot(plot = p)

  # Specifying the name and folder.
  SCpubr::savePlot(plot = p,
                   figure_path = "/path/to/my/figures/",
                   file_name = "my_figure")

  # Specify to also create a new folder.
  SCpubr::savePlot(plot = p,
                   figure_path = "/path/to/my/figures/",
                   file_name = "my_figure",
                   create_path = TRUE)

  # Set dimensions for the figure.
  SCpubr::savePlot(plot = p,
                   figure_path = "/path/to/my/figures/",
                   file_name = "my_figure",
                   create_path = TRUE,
                   width = 8,
                   height = 8)

  # Set quality (dpi).
  SCpubr::savePlot(plot = p,
                   figure_path = "/path/to/my/figures/",
                   file_name = "my_figure",
                   create_path = TRUE,
                   width = 8,
                   height = 8,
                   dpi = 300)
}
