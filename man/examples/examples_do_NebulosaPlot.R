\donttest{
  # Basic Nebulosa plot.
  p <- SCpubr::do_NebulosaPlot(sample = sample,
                              features = "CD14")

  # Querying multiple features.
  p <- SCpubr::do_NebulosaPlot(sample,
                               features = c("CD14", "CD8A"))

  # Compute joint density.
  p <- SCpubr::do_NebulosaPlot(sample = sample,
                               features = c("CD14", "CD8A"),
                               joint = TRUE)

  # Return only the joint density panel and modify the plot.title.
  p <- SCpubr::do_NebulosaPlot(sample = sample,
                               features = c("CD14", "CD8A"),
                               joint = TRUE)

  # Modify individual titles and also add a general one.
  # Use NA when you don't want to modify the panel title.
  p <- SCpubr::do_NebulosaPlot(sample = sample,
                               features = c("CD14", "CD8A"),
                               joint = TRUE,
                               individual.titles = c("Plot A",
                                                     NA,
                                                     "Combined density"),
                               plot.title = "Density analysis")

  # Modify viridis color maps.
  p <- SCpubr::do_NebulosaPlot(sample = sample,
                               features = "CD8A",
                               viridis_color_map = "A",
                               plot.title = "Magma")

}
