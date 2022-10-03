\dontrun{
  # Basic Bee Swarm plot - categorical coloring.
  # This will color based on the unique values of seurat_clusters.
  p <- SCpubr::do_BeeSwarmPlot(sample = sample,
                               feature_to_rank = "PC_1",
                               group.by = "seurat_clusters",
                               continuous_feature = FFALSE)

  # Basic Bee Swarm plot - continuous coloring.
  # This will color based on the PC_1 values.
  p <- SCpubr::do_BeeSwarmPlot(sample = sample,
                               feature_to_rank = "PC_1",
                               group.by = "seurat_clusters",
                               continuous_feature = TRUE)
}
