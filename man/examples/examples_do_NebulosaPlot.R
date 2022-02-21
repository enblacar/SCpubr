\dontrun{
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
  p1 <- SCpubr::do_NebulosaPlot(sample = sample, features = "CD8A", viridis_color_map = "A", plot.title = "Magma")
  p2 <- SCpubr::do_NebulosaPlot(sample = sample, features = "CD8A", viridis_color_map = "B", plot.title = "Inferno")
  p3 <- SCpubr::do_NebulosaPlot(sample = sample, features = "CD8A", viridis_color_map = "C", plot.title = "Plasma")
  p4 <- SCpubr::do_NebulosaPlot(sample = sample, features = "CD8A", viridis_color_map = "D", plot.title = "Viridis")
  p5 <- SCpubr::do_NebulosaPlot(sample = sample, features = "CD8A", viridis_color_map = "E", plot.title = "Cividis")
  p6 <- SCpubr::do_NebulosaPlot(sample = sample, features = "CD8A", viridis_color_map = "F", plot.title = "Rocket")
  p7 <- SCpubr::do_NebulosaPlot(sample = sample, features = "CD8A", viridis_color_map = "G", plot.title = "Mako")
  p8 <- SCpubr::do_NebulosaPlot(sample = sample, features = "CD8A", viridis_color_map = "H", plot.title = "Turbo")
}
