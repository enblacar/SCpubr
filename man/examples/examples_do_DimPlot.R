\donttest{
  # Define your Seurat object.
  sample <- readRDS(system.file("extdata/seurat_dataset_example.rds", package = "SCpubr"))

  # Basic DimPlot.
  p <- SCpubr::do_DimPlot(sample = sample)

  # Control dimensions.
  p <- SCpubr::do_DimPlot(sample = sample,
                          plot.title = "My awesome SC data set",
                          dims = c(2, 1))

  # Include a plot title.
  p <- SCpubr::do_DimPlot(sample = sample,
                          plot.title = "My awesome SC data set")

  # Include a plot subtitle.
  p <- SCpubr::do_DimPlot(sample = sample,
                          plot.subtitle = "My awesome SC data set")

  # Include a plot caption
  p <- SCpubr::do_DimPlot(sample = sample,
                          plot.caption = "My awesome SC data set")

  # Control legend position and number of columns.
  p <- SCpubr::do_DimPlot(sample = sample,
                          plot.title = "My awesome SC data set",
                          legend.position = "left",
                          legend.ncol = 2)

  # Control legend position and number of rows
  p <- SCpubr::do_DimPlot(sample = sample,
                          plot.title = "My awesome SC data set",
                          legend.position = "left",
                          legend.nrow = 2)

  # Use labels instead of legend.
  p <- SCpubr::do_DimPlot(sample = sample,
                          label = TRUE,
                          legend.position = "none")

  # Changing the order of plotting.
  # Using order with one identity value.
  p <- SCpubr::do_DimPlot(sample = sample,
                          shuffle = FALSE,
                          order = "5")

  # Using order with all identity values.
  p <- SCpubr::do_DimPlot(sample = sample,
                          shuffle = FALSE,
                          order = c("5", "8", "4",
                                    "9", "3", "1",
                                    "6", "0", "7", "2"))

  # Restrict the amount of identities displayed.
  p <- SCpubr::do_DimPlot(sample = sample,
                          idents.keep = c("1", "3", "5"))


  # Group by another variable rather than `Seurat::Idents(sample)`
  p <- SCpubr::do_DimPlot(sample = sample,
                          group.by = "seurat_clusters")

  # Split the output in as many plots as unique identities.
  p <- SCpubr::do_DimPlot(sample = sample,
                          split.by = "seurat_clusters")

  # Restrict the amount of identities displayed by split.by.
  p <- SCpubr::do_DimPlot(sample = sample,
                          split.by = "seurat_clusters",
                          idents.keep = c("1", "3", "5"))

  # Modify default colors.
  # Colors need to a named vector of the same length as the identities.
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
  p <- SCpubr::do_DimPlot(sample = sample,
                          colors.use = colors)

  # Highlight cells.
  # Select 1000 random cells out of clusters 1, 5 and 7.
  cells.use <- sample(colnames(sample[, sample$seurat_clusters %in% c("1", "5", "7")]), 50)
  p <- SCpubr::do_DimPlot(sample = sample,
                          cells.highlight = cells.use)

  # Highlight cells and use custom color.
  # Select 1000 random cells out of clusters 1, 5 and 7.
  cells.use <- sample(colnames(sample[, sample$seurat_clusters %in% c("1", "5", "7")]), 50)
  p <- SCpubr::do_DimPlot(sample = sample,
                          cells.highlight = cells.use,
                          colors.use = "black")

  # Change the size of the highlighted cells.
  # Select 1000 random cells out of clusters 1, 5 and 7.
  cells.use <- sample(colnames(sample[, sample$seurat_clusters %in% c("1", "5", "7")]), 50)
  p <- SCpubr::do_DimPlot(sample = sample,
                          cells.highlight = cells.use,
                          colors.use = "black",
                          sizes.highlight = 1)

  # Highlight given identities
  p <- SCpubr::do_DimPlot(sample,
                          idents.highlight = c("1", "3"))

  # Highlight given cells and given identities.
  cells.use <- sample(colnames(sample[, sample$seurat_clusters %in% c("1", "5", "7")]), 50)
  p <- SCpubr::do_DimPlot(sample,
                          cells.highlight = cells.use,
                          idents.highlight = c("2", "4"))
}


