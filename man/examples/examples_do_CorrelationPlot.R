\dontrun{
  # Default values.
  p <- SCpubr::do_CorrelationPlot(sample = sample)
  p

  # Custom grouping.
  sample$custom_group <- ifelse(sample$seurat_clusters %in% c("1", "3", "5", "7", "9"), "Group A", "Group B")
  p <- SCpubr::do_CorrelationPlot(sample = sample, group.by = "custom_group")
  p

  # Rotated axis labels.
  p <- SCpubr::do_CorrelationPlot(sample = sample,
                                  column_names_rot = 0)
  p

  # Increase cell size.
  p <- SCpubr::do_CorrelationPlot(sample = sample,
                                  column_names_rot = 0,
                                  cell_size = 7)
  p

  # Modified color scale.
  p <- SCpubr::do_CorrelationPlot(sample = sample,
                                  column_names_rot = 0,
                                  colors.use = colortools::opposite("steelblue", plot = F))
  p
}
