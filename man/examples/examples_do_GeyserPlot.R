\donttest{
  # Check Suggests.
  value <- SCpubr:::check_suggests(function_name = "do_GeyserPlot", passive = TRUE)

  if (isTRUE(value)){
    # Define your Seurat object.
    sample <- readRDS(system.file("extdata/seurat_dataset_example.rds", package = "SCpubr"))

    # Geyser plot with categorical color scale.
    p1 <- SCpubr::do_GeyserPlot(sample = sample,
                                features = "nCount_RNA",
                                scale_type = "categorical")

    # Geyser plot with continuous color scale.
    p2 <- SCpubr::do_GeyserPlot(sample = sample,
                                features = "nCount_RNA",
                                scale_type = "continuous")


    p <- p1 / p2
    p

    # Geyser plot with categorical color scale without ordering by mean.
    p1 <- SCpubr::do_GeyserPlot(sample = sample,
                                features = "nCount_RNA",
                                scale_type = "categorical",
                                order_by_mean = FALSE)

    # Geyser plot with continuous color scale without ordering by mean.
    p2 <- SCpubr::do_GeyserPlot(sample = sample,
                                features = "nCount_RNA",
                                scale_type = "continuous",
                                order_by_mean = FALSE)


    p <- p1 / p2
    p

    # Geyser plot with continuous color scale.
    p1 <- SCpubr::do_GeyserPlot(sample = sample,
                                features = "UMAP_2",
                                scale_type = "continuous",
                                enforce_symmetry = FALSE)

    # Geyser plot with continuous and symmetrical color scale.
    p2 <- SCpubr::do_GeyserPlot(sample = sample,
                                features = "UMAP_2",
                                scale_type = "continuous",
                                enforce_symmetry = TRUE)


    p <- p1 / p2
    p

    # Geyser plot with categorical color scale default X axis grouping.
    p1 <- SCpubr::do_GeyserPlot(sample = sample,
                                features = "UMAP_2",
                                scale_type = "categorical",
                                group.by = NULL,
                                xlab = "Seurat clusters")

    # Geyser plot with categorical color scale and custom grouping.
    p2 <- SCpubr::do_GeyserPlot(sample = sample,
                                features = "UMAP_2",
                                scale_type = "categorical",
                                group.by = "orig.ident",
                                xlab = "Individual sample")


    p <- p1 / p2
    p

    # We only have one value in orig.ident.
    # Let's modify it so that it resembles a multi-sample Seurat object.
    sample$modified_orig.ident <- sample(x = c("Sample_A", "Sample_B", "Sample_C"),
                                         size = ncol(sample),
                                         replace = TRUE,
                                         prob = c(0.2, 0.7, 0.1))

    # Geyser plot with categorical color scale split by seurat clusters.
    p1 <- SCpubr::do_GeyserPlot(sample = sample,
                                features = "UMAP_2",
                                scale_type = "categorical",
                                group.by = "modified_orig.ident",
                                split.by = "seurat_clusters")

    # Geyser plot with continuous color scale split by seurat clusters.
    p2 <- SCpubr::do_GeyserPlot(sample = sample,
                                features = "UMAP_2",
                                scale_type = "continuous",
                                group.by = "modified_orig.ident",
                                split.by = "seurat_clusters")


    p <- p1 / p2
    p

    # Geyser plot with categorical color scale split by seurat clusters rotating labels.
    p1 <- SCpubr::do_GeyserPlot(sample = sample,
                                features = "UMAP_2",
                                scale_type = "categorical",
                                group.by = "modified_orig.ident",
                                split.by = "seurat_clusters",
                                rotate_x_axis_labels = TRUE)

    # Geyser plot with continuous color scale split by seurat clusters rotating labels.
    p2 <- SCpubr::do_GeyserPlot(sample = sample,
                                features = "UMAP_2",
                                scale_type = "continuous",
                                enforce_symmetry = TRUE,
                                group.by = "modified_orig.ident",
                                split.by = "seurat_clusters",
                                rotate_x_axis_labels = TRUE)


    p <- p1 / p2
    p

    # Geyser plot with categorical color scale using color.by.
    p1 <- SCpubr::do_GeyserPlot(sample = sample,
                                features = "UMAP_2",
                                scale_type = "categorical",
                                color.by = "modified_orig.ident")

    # Geyser plot with continuous color scale using color.by.
    p2 <- SCpubr::do_GeyserPlot(sample = sample,
                                features = "UMAP_2",
                                scale_type = "continuous",
                                enforce_symmetry = TRUE,
                                color.by = "nCount_RNA")


    p <- p1 / p2
    p

    # Geyser plot with categorical color scale using color.by and split.by.
    p1 <- SCpubr::do_GeyserPlot(sample = sample,
                                features = "UMAP_2",
                                scale_type = "categorical",
                                group.by = "orig.ident",
                                split.by = "seurat_clusters",
                                color.by = "modified_orig.ident")

    # Geyser plot with continuous color scale using color.by and split.by.
    p2 <- SCpubr::do_GeyserPlot(sample = sample,
                                features = "UMAP_2",
                                scale_type = "continuous",
                                enforce_symmetry = TRUE,
                                color.by = "nCount_RNA",
                                group.by = "orig.ident",
                                split.by = "seurat_clusters")


    p <- p1 / p2
    p

    # Geyser plot with different jitter.
    p0 <- SCpubr::do_GeyserPlot(sample = sample,
                                features = "UMAP_2",
                                scale_type = "continuous",
                                enforce_symmetry = TRUE,
                                jitter = 0.01)

    p1 <- SCpubr::do_GeyserPlot(sample = sample,
                                features = "UMAP_2",
                                scale_type = "continuous",
                                enforce_symmetry = TRUE,
                                jitter = 0.1)

    p2 <- SCpubr::do_GeyserPlot(sample = sample,
                                features = "UMAP_2",
                                scale_type = "continuous",
                                enforce_symmetry = TRUE,
                                jitter = 0.2)

    p3 <- SCpubr::do_GeyserPlot(sample = sample,
                                features = "UMAP_2",
                                scale_type = "continuous",
                                enforce_symmetry = TRUE,
                                jitter = 0.3)

    p4 <- SCpubr::do_GeyserPlot(sample = sample,
                                features = "UMAP_2",
                                scale_type = "continuous",
                                enforce_symmetry = TRUE,
                                jitter = 0.4)

    p5 <- SCpubr::do_GeyserPlot(sample = sample,
                                features = "UMAP_2",
                                scale_type = "continuous",
                                enforce_symmetry = TRUE,
                                jitter = 0.49)


    p <- p0 / p1 / p2 / p3 / p4 / p5
    p
  } else if (isFALSE(value)){
    message("This function can not be used without its suggested packages.")
    message("Check out which ones are needed using `SCpubr::state_dependencies()`.")
  }
}
