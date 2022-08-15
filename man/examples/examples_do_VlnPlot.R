\donttest{
  # Basic violin plots.
  p <- SCpubr::do_VlnPlot(sample = sample,
                          features = "nCount_RNA")

  # Remove the inner boxplot.
  p <- SCpubr::do_VlnPlot(sample = sample,
                          features = "nCount_RNA",
                          plot_boxplot = FALSE)
  # Rotate x axis labels.
  p <- SCpubr::do_VlnPlot(sample = sample,
                          features = c("nCount_RNA"),
                          rotate_x_labels = TRUE)

  # Rotate x axis labels of only one panel.
  p <- SCpubr::do_VlnPlot(sample = sample,
                          features = c("nCount_RNA", "nFeature_RNA"),
                          rotate_x_labels = c(FALSE, TRUE),
                          ncol = 1)

  # Add a horizontal line if we want to show a given cutoff.
  p <- SCpubr::do_VlnPlot(sample = sample,
                          features = "nCount_RNA",
                          y_cut = 30000)

  # Add a horizontal line to one of the panels.
  p <- SCpubr::do_VlnPlot(sample = sample,
                          features = c("nCount_RNA", "nFeature_RNA"),
                          y_cut = c(NA, 5000),
                          individual.titles = c("UMIs", NA),
                          ncol = 1)

  # Modifying default colors.
  colors <- c("0" = "#001219",
              "1" = "#005f73",
              "2" = "#0a9396",
              "3" = "#94d2bd",
              "4" = "#e9d8a6",
              "5" = "#ee9b00",
              "6" = "#ca6702",
              "7" = "#bb3e03",
              "8" = "#ae2012",
              "9" = "#9b2226")

  p <- SCpubr::do_VlnPlot(sample = sample,
                          features = "nCount_RNA",
                          colors.use = colors)
}
