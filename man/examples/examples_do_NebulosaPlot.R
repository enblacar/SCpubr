\donttest{
  # Define your Seurat object.
  sample <- readRDS(system.file("extdata/seurat_dataset_example.rds", package = "SCpubr"))

  # Basic Nebulosa plot.
  p <- SCpubr::do_NebulosaPlot(sample = sample,
                               features = "EPC1")

  # Compute joint density.
  p <- SCpubr::do_NebulosaPlot(sample = sample,
                               features = c("EPC1", "TOX2"),
                               joint = TRUE)

  # Return only the joint density panel and modify the plot.title.
  p <- SCpubr::do_NebulosaPlot(sample = sample,
                               features = c("EPC1", "TOX2"),
                               joint = TRUE)


  # Modify viridis color maps.
  p <- SCpubr::do_NebulosaPlot(sample = sample,
                               features = "TOX2",
                               viridis_color_map = "A",
                               plot.title = "Magma")

}
