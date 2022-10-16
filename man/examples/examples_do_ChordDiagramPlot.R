\donttest{
  # Check Suggests.
  value <- SCpubr:::check_suggests(function_name = "do_ChordDiagramPlot", passive = TRUE)

  if (isTRUE(value)){
    # Define your Seurat object.
    sample <- readRDS(system.file("extdata/seurat_dataset_example.rds", package = "SCpubr"))

    # Basic chord diagram.
    sample$assignment <- ifelse(sample$seurat_clusters %in% c("0", "4", "7"), "A", "B")
    sample$assignment[sample$seurat_clusters %in% c("1", "2")] <- "C"
    sample$assignment[sample$seurat_clusters %in% c("10", "5")] <- "D"
    sample$assignment[sample$seurat_clusters %in% c("8", "9")] <- "E"

    p <- SCpubr::do_ChordDiagramPlot(sample = sample,
                                     from = "seurat_clusters",
                                     to = "assignment")

    p

    # Increase gap between from and to groups.
    p <- SCpubr::do_ChordDiagramPlot(sample = sample,
                                     from = "seurat_clusters",
                                     to = "assignment",
                                     big.gap = 40)

    p

    # Increase gap width groups in from and to.
    p <- SCpubr::do_ChordDiagramPlot(sample = sample,
                                     from = "seurat_clusters",
                                     to = "assignment",
                                     small.gap = 5)

    p

    # Control the alignment of the diagram.
    p <- SCpubr::do_ChordDiagramPlot(sample = sample,
                                     from = "seurat_clusters",
                                     to = "assignment",
                                     alignment = "horizontal")

    p

    p <- SCpubr::do_ChordDiagramPlot(sample = sample,
                                     from = "seurat_clusters",
                                     to = "assignment",
                                     alignment = "vertical")

    p

    # Control direction of links.
    p <- SCpubr::do_ChordDiagramPlot(sample = sample,
                                     from = "seurat_clusters",
                                     to = "assignment",
                                     directional = 0,
                                     direction.type = "diffHeight")

    p

    p <- SCpubr::do_ChordDiagramPlot(sample = sample,
                                     from = "seurat_clusters",
                                     to = "assignment",
                                     directional = 1)

    p

    p <- SCpubr::do_ChordDiagramPlot(sample = sample,
                                     from = "seurat_clusters",
                                     to = "assignment",
                                     directional = -1)

    p

    p <- SCpubr::do_ChordDiagramPlot(sample = sample,
                                     from = "seurat_clusters",
                                     to = "assignment",
                                     directional = 2,
                                     direction.type = "diffHeight")

    p

    # Add more padding to the labels.
    p <- SCpubr::do_ChordDiagramPlot(sample = sample,
                                     from = "seurat_clusters",
                                     to = "assignment",
                                     padding_labels = 8)

    p

    # Scale the size of the nodes.
    p <- SCpubr::do_ChordDiagramPlot(sample = sample,
                                     from = "seurat_clusters",
                                     to = "assignment",
                                     scale = TRUE,
                                     padding_labels = 8)

    p

    # Prevent self linking.
    sample$seurat_clusters2 <- sample$seurat_clusters
    p <- SCpubr::do_ChordDiagramPlot(sample = sample,
                                     from = "seurat_clusters",
                                     to = "seurat_clusters2",
                                     self.link = 1,
                                     scale = TRUE)

    p

    # Allow self linking.
    p <- SCpubr::do_ChordDiagramPlot(sample = sample,
                                     from = "seurat_clusters",
                                     to = "seurat_clusters2",
                                     self.link = 2,
                                     scale = TRUE)

    p

    # Set triangle arrows.
    p <- SCpubr::do_ChordDiagramPlot(sample = sample,
                                     from = "seurat_clusters",
                                     to = "assignment",
                                     link.arr.type = "triangle")

    p

    # Set big arrows.
    p <- SCpubr::do_ChordDiagramPlot(sample = sample,
                                     from = "seurat_clusters",
                                     to = "assignment",
                                     link.arr.type = "big.arrow")

    p
  } else if (isFALSE(value)){
    message("This function can not be used without its suggested packages.")
    message("Check out which ones are needed using `SCpubr::state_dependencies()`.")
  }
}

