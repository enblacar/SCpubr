\dontrun{
  # Basic bar plot showing the number of items in a variable and a title.
  # Legend T/F toggles removes or keeps the legend.
  p <- SCpubr::do_BarPlot(sample = sample,
                          features = "seurat_clusters",
                          legend = F,
                          plot.title = "Number of cells per cluster")

  # Same but horizontal bars.
  p <- SCpubr::do_BarPlot(sample = sample,
                          features = "seurat_clusters",
                          legend = F,
                          horizontal = T)

  # Grouping by a second variable.
  # Position = stack will compute as separate bars.
  p <- SCpubr::do_BarPlot(sample,
                          features = "orig.ident",
                          group.by = "seurat_clusters",
                          position = "stack",
                          legend = T,
                          horizontal = F)
  # Position = fill will compute the proportions of each group in the feature.
  p <- SCpubr::do_BarPlot(sample,
                          features = "orig.ident",
                          group.by = "seurat_clusters",
                          position = "fill",
                          legend = T,
                          horizontal = F)

  # Grouping by a second variable and reordering based on an item in the second variable.
  # Reordering with position = fill will order the bars by descending proportion of order.by.
  p <- SCpubr::do_BarPlot(sample,
                          features = "modified_orig.ident",
                          group.by = "modified_seurat_clusters",
                          plot.title = "Number of cells per cluster and sample",
                          order.by = "1",
                          position = "fill",
                          legend = T,
                          horizontal = F)
  # Reordering with position = stack will order the bars by descending number of order.by.
  p <- SCpubr::do_BarPlot(sample,
                          features = "modified_orig.ident",
                          group.by = "modified_seurat_clusters",
                          plot.title = "Number of cells per sample",
                          order.by = "1",
                          position = "stack",
                          legend = T,
                          horizontal = F)

  # Adding custom colors.
  # Colors must be a named vector of HEX codes of the same lengths as bars/groups to color.
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

  p <- SCpubr::do_BarPlot(sample = sample,
                          features = "seurat_clusters",
                          legend = F,
                          plot.title = "Number of cells per cluster",
                          horizontal = T,
                          colors.use = colors)
}
