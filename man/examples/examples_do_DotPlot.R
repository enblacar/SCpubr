\dontrun{
  # Basic Dot plot.
  p <- SCpubr::do_DotPlot(sample = sample,
                          features = "CD14")

  # Querying multiple features.
  genes <- c("IL7R", "CCR7", "CD14", "LYZ",
             "S100A4", "MS4A1", "CD8A", "FCGR3A",
             "MS4A7", "GNLY", "NKG7", "FCER1A",
             "CST3", "PPBP")
  p <- SCpubr::do_DotPlot(sample = sample,
                          features = genes)

  # Querying multiple features as a named list - splitting by each item in list.
  # Genes have to be unique.
  genes <- list("Naive CD4+ T" = c("IL7R", "CCR7"),
                "CD14+ Mono" = c("CD14", "LYZ"),
                "Memory CD4+" = c("S100A4"),
                "B" = c("MS4A1"),
                "CD8+ T" = c("CD8A"),
                "FCGR3A+ Mono" = c("FCGR3A", "MS4A7"),
                "NK" = c("GNLY", "NKG7"),
                "DC" = c("FCER1A", "CST3"),
                "Platelet" = c("PPBP"))
  p <- SCpubr::do_DotPlot(sample = sample,
                          features = genes)

  # Clustering the identities.
  genes <- list("Naive CD4+ T" = c("IL7R", "CCR7"),
                "CD14+ Mono" = c("CD14", "LYZ"),
                "Memory CD4+" = c("S100A4"),
                "B" = c("MS4A1"),
                "CD8+ T" = c("CD8A"),
                "FCGR3A+ Mono" = c("FCGR3A", "MS4A7"),
                "NK" = c("GNLY", "NKG7"),
                "DC" = c("FCER1A", "CST3"),
                "Platelet" = c("PPBP"))
  p <- SCpubr::do_DotPlot(sample = sample,
                          features = genes,
                          cluster.idents = TRUE,
                          plot.title = "Clustered")

  # Inverting the axes.
  genes <- c("IL7R", "CCR7", "CD14", "LYZ",
             "S100A4", "MS4A1", "CD8A", "FCGR3A",
             "MS4A7", "GNLY", "NKG7", "FCER1A",
             "CST3", "PPBP")
  p <- SCpubr::do_DotPlot(sample = sample,
                          features = genes,
                          cluster.idents = TRUE,
                          plot.title = "Clustered",
                          flip = T)

  # Modifying default colors.
  # Two colors to generate a gradient.
  genes <- c("IL7R", "CCR7", "CD14", "LYZ",
             "S100A4", "MS4A1", "CD8A", "FCGR3A",
             "MS4A7", "GNLY", "NKG7", "FCER1A",
             "CST3", "PPBP")
  p <- SCpubr::do_DotPlot(sample = sample,
                          features = genes,
                          colors.use = c("#001219", "#e9d8a6"))
}
