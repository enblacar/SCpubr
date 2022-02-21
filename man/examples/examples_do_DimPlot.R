\dontrun{
# Your Seurat object.
sample <- my_seurat_object

# Basic DimPlot.
p <- SCpubr::do_DimPlot(sample = sample)

# Include a plot title.
p <- SCpubr::do_DimPlot(sample = sample,
                        plot.title = "My awesome SC data set")

# Control legend position and number of columns.
p <- SCpubr::do_DimPlot(sample = sample,
                        plot.title = "My awesome SC data set",
                        legend.position = "left",
                        legend.ncol = 2)

# Use labels instead of legend.
p <- SCpubr::do_DimPlot(sample = sample,
                        label = TRUE,
                        legend = FALSE)

# Group by another variable rather than `Seurat::Idents(sample)`
p <- SCpubr::do_DimPlot(sample = sample,
                        group.by = "orig.ident")

# Split the output in as many plots as unique identities.
p <- SCpubr::do_DimPlot(sample = sample,
                        split.by = "orig.ident")

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
cells.use <- sample(colnames(sample[, sample$seurat_clusters %in% c("1", "5", "7")]), 1000)
p <- Seurat::DimPlot(sample,
                     cells.highlight = cells.use)

# Highlight cells and use custom color.
# Select 1000 random cells out of clusters 1, 5 and 7.
cells.use <- sample(colnames(sample[, sample$seurat_clusters %in% c("1", "5", "7")]), 1000)
p <- Seurat::DimPlot(sample,
                     cells.highlight = cells.use,
                     colors.use = "black")

# Change the size of the highlighted cells.
# Select 1000 random cells out of clusters 1, 5 and 7.
cells.use <- sample(colnames(sample[, sample$seurat_clusters %in% c("1", "5", "7")]), 1000)
p <- Seurat::DimPlot(sample,
                     cells.highlight = cells.use,
                     colors.use = "black",
                     sizes.highlight = 1)
}
