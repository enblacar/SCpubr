\donttest{
  # Define your Seurat object.
  sample <- readRDS(system.file("extdata/seurat_dataset_example.rds", package = "SCpubr"))

  # Regular FeaturePlot.
  p <- SCpubr::do_FeaturePlot(sample = sample,
                              features = "nCount_RNA")

  # Add a title, subtitle and caption.
  p <- SCpubr::do_FeaturePlot(sample = sample,
                              features = "nCount_RNA",
                              plot.title = "Number of UMIs",
                              plot.subtitle = "Number of unique mRNAs per cell.",
                              plot.caption = "Plot generated using SCpubr.")

  # Add a subtitle.
  p <- SCpubr::do_FeaturePlot(sample = sample,
                              features = "nCount_RNA",
                              plot.title = "Number of UMIs",
                              plot.subtitle = "Number of unique mRNAs per cell.")

  # Add a caption
  p <- SCpubr::do_FeaturePlot(sample = sample,
                              features = "nCount_RNA",
                              plot.title = "Number of UMIs",
                              plot.subtitle = "Number of unique mRNAs per cell.",
                              plot.caption = "Plot generated using SCpubr.")

  # Plot multiple features and control the output columns.
  p <- SCpubr::do_FeaturePlot(sample, features = c("nCount_RNA",
                                                   "nFeature_RNA",
                                                   "EPC1"),
                              plot.title = "My very important features",
                              ncol = 3)

  # FeaturePlot with a subset of cells maintaining the original UMAP shape.
  cells.plot <- colnames(sample[, !(sample$seurat_clusters %in% c("2", "5", "8"))])
  p <- SCpubr::do_FeaturePlot(sample = sample,
                              cells.highlight = cells.plot,
                              features = c("EPC1"))

  # FeaturePlot with a subset of identities
  # (in Seurat::Idents(sample)) maintaining the original UMAP shape.
  idents.use <- levels(sample)[!(levels(sample) %in% c("2", "5", "8"))]
  p <- SCpubr::do_FeaturePlot(sample = sample,
                              idents.highlight = idents.use,
                              features = c("EPC1"))

  # Splitting the FeaturePlot by a variable and
  # maintaining the color scale and the UMAP shape.
  p <- SCpubr::do_FeaturePlot(sample = sample,
                              features = "EPC1",
                              split.by = "seurat_clusters")

  # Splitting the FeaturePlot by a variable
  # and subset only some of the unique groups.
  p <- SCpubr::do_FeaturePlot(sample,
                              features = c("EPC1", "nCount_RNA", "PC_1"),
                              split.by = "seurat_clusters",
                              split.by.idents = c("2", "5"))

  # Modify the viridis color maps.
  p <- SCpubr::do_FeaturePlot(sample = sample,
                              features = "EPC1",
                              viridis_color_map = "A",
                              plot.title = "Magma")
}
