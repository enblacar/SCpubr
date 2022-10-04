\dontrun{
  # Basic violin plot.
  p <- SCpubr::do_ViolinPlot(sample = sample,
                             feature = "nCount_RNA")
  p

  # Remove the box plots.
  p <- SCpubr::do_ViolinPlot(sample = sample,
                             feature = "nCount_RNA",
                             plot_boxplot = FALSE)
  p

  # Rotate x axis labels.
  p <- SCpubr::do_ViolinPlot(sample = sample,
                             feature = c("nCount_RNA"),
                             rotate_x_axis_labels = FALSE)
  p

  # Add horizontal lines.
  p <- SCpubr::do_ViolinPlot(sample = sample,
                             feature = "nCount_RNA",
                             y_cut = 25000)
  p

  # Increase line width.
  p1 <- SCpubr::do_ViolinPlot(sample = sample,
                              feature = "nCount_RNA")

  p2 <- SCpubr::do_ViolinPlot(sample = sample,
                              feature = "nCount_RNA",
                              line_width = 2)

  p <- p1 / p2
  p

  # Decrease boxplot width.
  p1 <- SCpubr::do_ViolinPlot(sample = sample,
                              feature = "nCount_RNA")

  p2 <- SCpubr::do_ViolinPlot(sample = sample,
                              feature = "nCount_RNA",
                              boxplot_width = 0.1)

  p <- p1 / p2
  p
}