\dontrun{
  # Basic box plot.
  p <- SCpubr::do_BoxPlot(sample = sample,
                          feature = "nCount_RNA")
  p

  # Generate a custom group.
  sample$custom_group = ifelse(colnames(sample) %in% sample(colnames(sample), 4000), "A", "B")

  # Use custom grouping.
  p <- SCpubr::do_BoxPlot(sample = sample,
                          feature = "nCount_RNA",
                          group.by = "custom_group")
  p

  # Flip the box plot.
  p <- SCpubr::do_BoxPlot(sample = sample,
                          feature = "nCount_RNA",
                          flip = TRUE)
  p

  # Use silhouette style.
  p <- SCpubr::do_BoxPlot(sample = sample,
                          feature = "nCount_RNA",
                          use_silhouette = TRUE)
  p

  # Order by mean values.
  p <- SCpubr::do_BoxPlot(sample = sample,
                          feature = "nCount_RNA",
                          order = TRUE)
  p

  # Apply second grouping.
  p <- SCpubr::do_BoxPlot(sample = sample,
                          feature = "nCount_RNA",
                          split.by = "custom_group")
  p

  # Apply statistical tests.
  p <- SCpubr::do_BoxPlot(sample = sample,
                          feature = "nCount_RNA",
                          use_test = TRUE,
                          comparisons = list(c("0", "1"),
                                             c("3", "4"),
                                             c("5", "9")))
  p

  # Apply statistical tests and show the p-value.
  p <- SCpubr::do_BoxPlot(sample = sample,
                          feature = "nCount_RNA",
                          use_test = TRUE,
                          comparisons = list(c("0", "1"),
                                             c("3", "4"),
                                             c("5", "9")),
                          map_signif_level = FALSE)
  p
}
