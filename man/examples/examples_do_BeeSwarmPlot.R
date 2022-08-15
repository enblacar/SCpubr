\donttest{
  # Basic Bee Swarm plot - categorical coloring.
  # This will color based on the unique values of seurat_clusters.
  p <- SCpubr::do_BeeSwarmPlot(sample = sample,
                              feature_to_rank = "PC_1",
                              group.by = "seurat_clusters",
                              continuous_feature = F)

  # Basic Bee Swarm plot - continuous coloring.
  # This will color based on the PC_1 values.
  p <- SCpubr::do_BeeSwarmPlot(sample = sample,
                               feature_to_rank = "PC_1",
                               group.by = "seurat_clusters",
                               continuous_feature = T)

  # Change default colors - categorical variables.
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

  p <- SCpubr::do_BeeSwarmPlot(sample = sample,
                               feature_to_rank = "Monocyte_signature",
                               group.by = "seurat_clusters",
                               colors.use = colors)

  # Change default colors - continuous variables.
  p <- SCpubr::do_BeeSwarmPlot(sample = sample,
                               feature_to_rank = "Monocyte_signature",
                               group.by = "seurat_clusters",
                               continuous_feature = TRUE,
                               fontsize = 10,
                               viridis_color_map = "A",
                               plot.title = "Magma")


}
