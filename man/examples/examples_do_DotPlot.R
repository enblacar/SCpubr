\donttest{
  # Define your Seurat object.
  sample <- readRDS(system.file("extdata/seurat_dataset_example.rds", package = "SCpubr"))

  # Basic Dot plot.
  p <- SCpubr::do_DotPlot(sample = sample,
                          features = "EPC1")

  # Querying multiple features.
  genes <- rownames(sample)[1:14]
  p <- SCpubr::do_DotPlot(sample = sample,
                          features = genes)

  # Inverting the axes.
  p <- SCpubr::do_DotPlot(sample = sample,
                          features = genes,
                          cluster.idents = TRUE,
                          plot.title = "Clustered",
                          flip = TRUE)

  # Modifying default colors.
  # Two colors to generate a gradient.
  p <- SCpubr::do_DotPlot(sample = sample,
                          features = genes,
                          colors.use = c("#001219", "#e9d8a6"))

  # Querying multiple features as a named list - splitting by each item in list.
  # Genes have to be unique.
  genes <- list("Naive CD4+ T" = rownames(sample)[1:2],
                "EPC1+ Mono" = rownames(sample)[3:4],
                "Memory CD4+" = rownames(sample)[5],
                "B" = rownames(sample)[6],
                "CD8+ T" = rownames(sample)[7],
                "FCGR3A+ Mono" = rownames(sample)[8:9],
                "NK" = rownames(sample)[10:11],
                "DC" = rownames(sample)[12:13],
                "Platelet" = rownames(sample)[14])

  p <- SCpubr::do_DotPlot(sample = sample,
                          features = genes)

  # Clustering the identities.
  p <- SCpubr::do_DotPlot(sample = sample,
                          features = genes,
                          cluster.idents = TRUE,
                          plot.title = "Clustered")
}
