\dontrun{
  # Define your sample and assay.
  # sample <- your_seurat_object
  # assay <- "your_normalized_data_assay"

  # Retrieve prior knowledge network.
  # network <- decoupleR::get_progeny(organism = "human")
  #
  # # Run weighted means algorithm.
  # activities <- decoupleR::run_wmean(mat = as.matrix(sample@assays[[assay]]@data),
  #                                    network = network,
  #                                    .source = "source",
  #                                    .targe = "target",
  #                                    .mor = "weight",
  #                                    times = 100,
  #                                    minsize = 5)

  # General heatmap.
  out <- SCpubr::do_PathwayActivityPlot(sample = sample,
                                        activities = activities)
  p <- out$heatmaps$average_scores
  p

  # Retrieve feature plots.
  out <- SCpubr::do_PathwayActivityPlot(sample = sample,
                                        activities = activities,
                                        plot_FeaturePlots = TRUE)
  p1 <- SCpubr::do_DimPlot(sample)
  p2 <- out$feature_plots$EGFR
  p <- p1 | p2
  p

  # Retrieve Geyser plots.
  out <- SCpubr::do_PathwayActivityPlot(sample = sample,
                                        activities = activities,
                                        plot_GeyserPlots = TRUE)
  p1 <- SCpubr::do_DimPlot(sample)
  p2 <- out$geyser_plots$EGFR
  p <- p1 | p2
  p

  # Use non-symmetrical color scale.
  out <- SCpubr::do_PathwayActivityPlot(sample = sample,
                                        activities = activities,
                                        plot_GeyserPlots = TRUE,
                                        plot_FeaturePlots = TRUE,
                                        symmetrical_scale = FALSE)
  p1 <- out$feature_plots$EGFR
  p2 <- out$geyser_plots$EGFR

  out <- SCpubr::do_PathwayActivityPlot(sample = sample,
                                        activities = activities,
                                        plot_GeyserPlots = TRUE,
                                        plot_FeaturePlots = TRUE,
                                        symmetrical_scale = TRUE)
  p3 <- out$feature_plots$EGFR
  p4 <- out$geyser_plots$EGFR

  p <- (p1 | p2) / (p3 | p4)
  p

  # Not order Geyser plot by mean values.
  out <- SCpubr::do_PathwayActivityPlot(sample = sample,
                                        activities = activities,
                                        plot_GeyserPlots = TRUE,
                                        symmetrical_scale = TRUE,
                                        geyser_order_by_mean = FALSE)
  p1 <- out$geyser_plots$EGFR

  out <- SCpubr::do_PathwayActivityPlot(sample = sample,
                                        activities = activities,
                                        plot_GeyserPlots = TRUE,
                                        symmetrical_scale = TRUE,
                                        geyser_order_by_mean = TRUE)
  p2 <- out$geyser_plots$EGFR

  p <- p1 | p2
  p

  # Plot a third variable in Geyser plots.
  out <- SCpubr::do_PathwayActivityPlot(sample = sample,
                                        activities = activities,
                                        plot_GeyserPlots = TRUE,
                                        geyser_color.by = "seurat_clusters",
                                        geyser_scale_type = "categorical")
  p1 <- out$geyser_plots$EGFR

  out <- SCpubr::do_PathwayActivityPlot(sample = sample,
                                        activities = activities,
                                        plot_GeyserPlots = TRUE,
                                        geyser_color.by = "nCount_RNA",
                                        geyser_scale_type = "continuous")
  p2 <- out$geyser_plots$EGFR

  p <- p1 | p2
  p

  # Split the heatmap by another variable.
  sample$split.me <- ifelse(sample$seurat_clusters %in% c("0", "3", "7"), "Group A","Group B")

  out <- SCpubr::do_PathwayActivityPlot(sample = sample,
                                        activities = activities,
                                        split.by = "split.me")
  p <- out$heatmaps$average_scores
  p
}