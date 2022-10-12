\dontrun{
  # Generate a plot.
  p <- SCpubr::do_DimPlot(sample = sample)

  # Default parameters.
  SCpubr::save_Plot(plot = p)

  # Specifying the name and folder.
  SCpubr::save_Plot(plot = p,
                    figure_path = "/path/to/my/figures/",
                    file_name = "my_figure")

  # Specify to also create a new folder.
  SCpubr::save_Plot(plot = p,
                    figure_path = "/path/to/my/figures/",
                    file_name = "my_figure",
                    create_path = TRUE)

  # Set dimensions for the figure.
  SCpubr::save_Plot(plot = p,
                    figure_path = "/path/to/my/figures/",
                    file_name = "my_figure",
                    create_path = TRUE,
                    width = 8,
                    height = 8)

  # Set quality (dpi).
  SCpubr::save_Plot(plot = p,
                    figure_path = "/path/to/my/figures/",
                    file_name = "my_figure",
                    create_path = TRUE,
                    width = 8,
                    height = 8,
                    dpi = 300)
}