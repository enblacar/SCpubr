\dontrun{
  # Regular FeaturePlot.
  p <- SCpubr::do_FeaturePlot(sample = sample,
                              features = "nCount_RNA",
                              plot.title = "Number of UMIs")

  # Add a title.
  p <- SCpubr::do_FeaturePlot(sample = sample,
                              features = "nCount_RNA",
                              plot.title = "Number of UMIs")

  # Plot multiple features and control the output columns.
  p <- SCpubr::do_FeaturePlot(sample, features = c("nCount_RNA",
                                                   "nFeature_RNA",
                                                   "percent.mt",
                                                   "CD14"),
                              plot.title = "My very important features",
                              ncol = 2)

  # Plot multiple features and control each panel plot title.
  p <- SCpubr::do_FeaturePlot(sample = sample,
                              features = c("nCount_RNA",
                                           "nFeature_RNA",
                                           "percent.mt",
                                           "CD14"),
                              plot.title = "My very important features",
                              individual.titles = c("Plot A",
                                                    "Plot_B",
                                                    NA,
                                                    "Plot_D"),
                              ncol = 2)

  # FeaturePlot with a subset of cells maintaining the original UMAP shape.
  cells.plot <- colnames(sample[, !(sample$seurat_clusters %in% c("2", "5", "8"))])
  p <- SCpubr::do_FeaturePlot(sample = sample,
                              cells.highlight = cells.plot,
                              features = c("CD14"))

  # FeaturePlot with a subset of identities
  # (in Seurat::Idents(sample)) maintaining the original UMAP shape.
  idents.use <- levels(sample)[!(levels(sample) %in% c("2", "5", "8"))]
  p <- SCpubr::do_FeaturePlot(sample = sample,
                              idents.highlight = idents.use,
                              features = c("CD14"))

  # Splitting the FeaturePlot by a variable and
  # maintaining the color scale and the UMAP shape.
  p <- SCpubr::do_FeaturePlot(sample = sample,
                              features = "LYN",
                              split.by = "new_clusters")

  # Splitting the FeaturePlot by a variable
  # and subset only some of the unique groups.
  p <- SCpubr::do_FeaturePlot(sample,
                              features = c("LYN", "nCount_RNA", "PC_1"),
                              split.by = "seurat_clusters",
                              split.by.idents = c("2", "5"))

  # Modify the viridis color maps.
  p <- SCpubr::do_FeaturePlot(sample = sample,
                              features = "CD14",
                              viridis_color_map = "A",
                              plot.title = "Magma")



}
