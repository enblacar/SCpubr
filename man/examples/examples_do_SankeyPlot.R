\dontrun{
  # Generate a more fine-grained clustering.
  sample <- Seurat::FindSubCluster(sample, cluster = c("0", "5"), graph.name = "SCT_snn")
  sample$sub.cluster <- paste0("sub_", sample$sub.cluster)

  # Compute basic sankey plot.
  p1 <- SCpubr::do_SankeyPlot(sample = sample,
                              first_group = "sub.cluster",
                              last_group = "seurat_clusters",
                              type = "sankey")

  # Compute basic alluvial plot.
  p2 <- SCpubr::do_SankeyPlot(sample = sample,
                              first_group = "sub.cluster",
                              last_group = "seurat_clusters",
                              type = "alluvial")

  p <- p1 / p2
  p

  sample$assignment <- ifelse(sample$seurat_clusters %in% c("0", "2", "4"), "A", "B")

  # Add more groups.
  p1 <- SCpubr::do_SankeyPlot(sample = sample,
                              first_group = "sub.cluster",
                              middle_groups = c("seurat_clusters", "assignment"),
                              last_group = "orig.ident",
                              type = "sankey")

  p2 <- SCpubr::do_SankeyPlot(sample = sample,
                              first_group = "sub.cluster",
                              middle_groups = c("seurat_clusters", "assignment"),
                              last_group = "orig.ident",
                              type = "alluvial")

  p <- p1 / p2
  p

  # Control the color and fill of the nodes.
  p1 <- SCpubr::do_SankeyPlot(sample = sample,
                              first_group = "sub.cluster",
                              middle_groups = c("seurat_clusters", "assignment"),
                              last_group = "orig.ident",
                              type = "sankey",
                              node.fill = "grey95",
                              node.color = "black")

  p2 <- SCpubr::do_SankeyPlot(sample = sample,
                              first_group = "sub.cluster",
                              middle_groups = c("seurat_clusters", "assignment"),
                              last_group = "orig.ident",
                              type = "alluvial",
                              node.fill = "grey95",
                              node.color = "black")

  p <- p1 / p2
  p

  # Control the width of the nodes.
  p1 <- SCpubr::do_SankeyPlot(sample = sample,
                              first_group = "sub.cluster",
                              middle_groups = c("seurat_clusters", "assignment"),
                              last_group = "orig.ident",
                              type = "sankey",
                              node.fill = "grey95",
                              node.color = "black")

  p2 <- SCpubr::do_SankeyPlot(sample = sample,
                              first_group = "sub.cluster",
                              middle_groups = c("seurat_clusters", "assignment"),
                              last_group = "orig.ident",
                              type = "sankey",
                              node.fill = "grey95",
                              node.color = "black",
                              width = 0.5)

  p <- p1 / p2
  p

  # Control the alignment of the labels.
  p1 <- SCpubr::do_SankeyPlot(sample = sample,
                              first_group = "sub.cluster",
                              middle_groups = c("seurat_clusters", "assignment"),
                              last_group = "orig.ident",
                              type = "sankey",
                              node.fill = "grey95",
                              node.color = "black",
                              width = 0.5)

  p2 <- SCpubr::do_SankeyPlot(sample = sample,
                              first_group = "sub.cluster",
                              middle_groups = c("seurat_clusters", "assignment"),
                              last_group = "orig.ident",
                              type = "sankey",
                              node.fill = "grey95",
                              node.color = "black",
                              width = 0.5,
                              hjust = 0.5)

  p <- p1 / p2
  p

  # Use text or labels for the nodes.
  p1 <- SCpubr::do_SankeyPlot(sample = sample,
                              first_group = "sub.cluster",
                              middle_groups = c("seurat_clusters", "assignment"),
                              last_group = "orig.ident",
                              type = "sankey")

  p2 <- SCpubr::do_SankeyPlot(sample = sample,
                              first_group = "sub.cluster",
                              middle_groups = c("seurat_clusters", "assignment"),
                              last_group = "orig.ident",
                              type = "sankey",
                              use_labels = TRUE,
                              text_color = "white")

  p <- p1 / p2
  p

  # Modify the space between nodes.
  p1 <- SCpubr::do_SankeyPlot(sample = sample,
                              first_group = "sub.cluster",
                              middle_groups = c("seurat_clusters", "assignment"),
                              last_group = "orig.ident",
                              type = "sankey",
                              space =  1000)

  p2 <- SCpubr::do_SankeyPlot(sample = sample,
                              first_group = "sub.cluster",
                              middle_groups = c("seurat_clusters", "assignment"),
                              last_group = "orig.ident",
                              type = "sankey",
                              space = 5000)

  p <- p1 / p2
  p

  # Modify default colors.
  colors.first <- SCpubr::do_ColorPalette(colors.use = "steelblue",
                                          n = length(unique(sample$sub.cluster)))
  names(colors.first) <- unique(sample$sub.cluster)

  p <- SCpubr::do_SankeyPlot(sample = sample,
                             first_group = "sub.cluster",
                             middle_groups = c("seurat_clusters", "assignment"),
                             last_group = "orig.ident",
                             type = "sankey",
                             colors.first = )
  p
}
