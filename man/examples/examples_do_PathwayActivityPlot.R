\donttest{
  # Check Suggests.
  value <- SCpubr:::check_suggests(function_name = "do_PathwayActivityPlot", passive = TRUE)

  if (isTRUE(value)){
    # Define your Seurat object.
    sample <- readRDS(system.file("extdata/seurat_dataset_example.rds",
                                  package = "SCpubr"))

    # Define your activities object.
    progeny_activities <- readRDS(system.file("extdata/progeny_activities_example.rds",
                                              package = "SCpubr"))

    # General heatmap.
    out <- SCpubr::do_PathwayActivityPlot(sample = sample,
                                          activities = progeny_activities)
    p <- out$heatmaps$average_scores
    p

    # Retrieve feature plots.
    out <- SCpubr::do_PathwayActivityPlot(sample = sample,
                                          activities = progeny_activities,
                                          plot_FeaturePlots = TRUE)
    p1 <- SCpubr::do_DimPlot(sample)
    p2 <- out$feature_plots$EGFR
    p <- p1 | p2
    p

    # Retrieve Geyser plots.
    out <- SCpubr::do_PathwayActivityPlot(sample = sample,
                                          activities = progeny_activities,
                                          plot_GeyserPlots = TRUE)
    p1 <- SCpubr::do_DimPlot(sample)
    p2 <- out$geyser_plots$EGFR
    p <- p1 | p2
    p

    # Use non-symmetrical color scale.
    out <- SCpubr::do_PathwayActivityPlot(sample = sample,
                                          activities = progeny_activities,
                                          plot_GeyserPlots = TRUE,
                                          plot_FeaturePlots = TRUE,
                                          enforce_symmetry = FALSE)
    p1 <- out$feature_plots$EGFR
    p2 <- out$geyser_plots$EGFR

    out <- SCpubr::do_PathwayActivityPlot(sample = sample,
                                          activities = progeny_activities,
                                          plot_GeyserPlots = TRUE,
                                          plot_FeaturePlots = TRUE,
                                          enforce_symmetry = TRUE)
    p3 <- out$feature_plots$EGFR
    p4 <- out$geyser_plots$EGFR

    p <- (p1 | p2) / (p3 | p4)
    p

    # Not order Geyser plot by mean values.
    out <- SCpubr::do_PathwayActivityPlot(sample = sample,
                                          activities = progeny_activities,
                                          plot_GeyserPlots = TRUE,
                                          enforce_symmetry = TRUE,
                                          geyser_order_by_mean = FALSE)
    p1 <- out$geyser_plots$EGFR

    out <- SCpubr::do_PathwayActivityPlot(sample = sample,
                                          activities = progeny_activities,
                                          plot_GeyserPlots = TRUE,
                                          enforce_symmetry = TRUE,
                                          geyser_order_by_mean = TRUE)
    p2 <- out$geyser_plots$EGFR

    p <- p1 | p2
    p

    # Plot a third variable in Geyser plots.
    out <- SCpubr::do_PathwayActivityPlot(sample = sample,
                                          activities = progeny_activities,
                                          plot_GeyserPlots = TRUE,
                                          geyser_color.by = "seurat_clusters",
                                          geyser_scale_type = "categorical")
    p1 <- out$geyser_plots$EGFR

    out <- SCpubr::do_PathwayActivityPlot(sample = sample,
                                          activities = progeny_activities,
                                          plot_GeyserPlots = TRUE,
                                          geyser_color.by = "nCount_RNA",
                                          geyser_scale_type = "continuous")
    p2 <- out$geyser_plots$EGFR

    p <- p1 | p2
    p

    # Split the heatmap by another variable.
    sample$split.me <- ifelse(sample$seurat_clusters %in% c("0", "3", "7"), "Group A","Group B")

    out <- SCpubr::do_PathwayActivityPlot(sample = sample,
                                          activities = progeny_activities,
                                          split.by = "split.me")
    p <- out$heatmaps$average_scores
    p
  } else if (isFALSE(value)){
    message("This function can not be used without its suggested packages.")
    message("Check out which ones are needed using `SCpubr::state_dependencies()`.")
  }
}
