# Define your Seurat object.
sample <- readRDS(system.file("extdata/seurat_dataset_example.rds", package = "SCpubr"))

# Compute basic sankey plot.
p1 <- SCpubr::do_SankeyPlot(sample = sample,
                            first_group = "seurat_clusters",
                            last_group = "orig.ident",
                            type = "sankey")

# Compute basic alluvial plot.
p2 <- SCpubr::do_SankeyPlot(sample = sample,
                            first_group = "seurat_clusters",
                            last_group = "orig.ident",
                            type = "alluvial")

p <- p1 / p2
p

sample$assignment <- ifelse(sample$seurat_clusters %in% c("0", "2", "4"), "A", "B")

# Add more groups.
p1 <- SCpubr::do_SankeyPlot(sample = sample,
                            first_group = "seurat_clusters",
                            middle_groups = c("assignment"),
                            last_group = "orig.ident",
                            type = "sankey")

p2 <- SCpubr::do_SankeyPlot(sample = sample,
                            first_group = "seurat_clusters",
                            middle_groups = c("assignment"),
                            last_group = "orig.ident",
                            type = "alluvial")

p <- p1 / p2
p

# Control the color and fill of the nodes.
p1 <- SCpubr::do_SankeyPlot(sample = sample,
                            first_group = "seurat_clusters",
                            middle_groups = c("assignment"),
                            last_group = "orig.ident",
                            type = "sankey",
                            node.fill = "grey95",
                            node.color = "black")

p2 <- SCpubr::do_SankeyPlot(sample = sample,
                            first_group = "seurat_clusters",
                            middle_groups = c("assignment"),
                            last_group = "orig.ident",
                            type = "alluvial",
                            node.fill = "grey95",
                            node.color = "black")

p <- p1 / p2
p

# Control the width of the nodes.
p1 <- SCpubr::do_SankeyPlot(sample = sample,
                            first_group = "seurat_clusters",
                            middle_groups = c("assignment"),
                            last_group = "orig.ident",
                            type = "sankey",
                            node.fill = "grey95",
                            node.color = "black")

p2 <- SCpubr::do_SankeyPlot(sample = sample,
                            first_group = "seurat_clusters",
                            middle_groups = c("assignment"),
                            last_group = "orig.ident",
                            type = "sankey",
                            node.fill = "grey95",
                            node.color = "black",
                            width = 0.5)

p <- p1 / p2
p

# Control the alignment of the labels.
p1 <- SCpubr::do_SankeyPlot(sample = sample,
                            first_group = "seurat_clusters",
                            middle_groups = c("assignment"),
                            last_group = "orig.ident",
                            type = "sankey",
                            node.fill = "grey95",
                            node.color = "black",
                            width = 0.5)

p2 <- SCpubr::do_SankeyPlot(sample = sample,
                            first_group = "seurat_clusters",
                            middle_groups = c("assignment"),
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
                            first_group = "seurat_clusters",
                            middle_groups = c("assignment"),
                            last_group = "orig.ident",
                            type = "sankey")

p2 <- SCpubr::do_SankeyPlot(sample = sample,
                            first_group = "seurat_clusters",
                            middle_groups = c("assignment"),
                            last_group = "orig.ident",
                            type = "sankey",
                            use_labels = TRUE,
                            text_color = "white")

p <- p1 / p2
p

# Modify the space between nodes.
p1 <- SCpubr::do_SankeyPlot(sample = sample,
                            first_group = "seurat_clusters",
                            middle_groups = c("assignment"),
                            last_group = "orig.ident",
                            type = "sankey",
                            space =  20)

p2 <- SCpubr::do_SankeyPlot(sample = sample,
                            first_group = "seurat_clusters",
                            middle_groups = c("assignment"),
                            last_group = "orig.ident",
                            type = "sankey",
                            space = 40)

p <- p1 / p2
p
