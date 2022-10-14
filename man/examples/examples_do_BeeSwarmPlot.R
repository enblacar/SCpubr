\donttest{
  # Define your Seurat object.
  sample <- readRDS(system.file("extdata/seurat_dataset_example.rds", package = "SCpubr"))

  # Basic Bee Swarm plot - categorical coloring.
  # This will color based on the unique values of seurat_clusters.
  p <- SCpubr::do_BeeSwarmPlot(sample = sample,
                               feature_to_rank = "PC_1",
                               group.by = "seurat_clusters",
                               continuous_feature = FALSE)

  # Basic Bee Swarm plot - continuous coloring.
  # This will color based on the PC_1 values.
  p <- SCpubr::do_BeeSwarmPlot(sample = sample,
                               feature_to_rank = "PC_1",
                               group.by = "seurat_clusters",
                               continuous_feature = TRUE)

}
