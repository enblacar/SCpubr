\dontrun{
  # Define some gene sets to query. It has to be a named list.
  gene_set <- list("A" = Seurat::VariableFeatures(sample)[1:10],
                   "B" = Seurat::VariableFeatures(sample)[11:20],
                   "C" = Seurat::VariableFeatures(sample)[21:30],
                   "D" = Seurat::VariableFeatures(sample)[31:40])

  # Using two variables: A scatter plot X vs Y.
  p <- SCpubr::do_CellularStatesPlot(sample = sample,
                                     gene_list = gene_set,
                                     x1 = "A",
                                     y1 = "B")
  p

  # Using two variables: Enforce symmetry on the plot.
  p <- SCpubr::do_CellularStatesPlot(sample = sample,
                                     gene_list = gene_set,
                                     x1 = "A",
                                     y1 = "B",
                                     enforce_symmetry = TRUE)
  p

  # Using three variables. Figure from: https://www.nature.com/articles/nature20123.
  p <- SCpubr::do_CellularStatesPlot(sample = sample,
                                     gene_list = gene_set,
                                     x1 = "A",
                                     y1 = "B",
                                     x2 = "C")
  p

  # Using three variables. Enforcing symmetry will align X axis with 0 in the center.
  p <- SCpubr::do_CellularStatesPlot(sample = sample,
                                     gene_list = gene_set,
                                     x1 = "A",
                                     y1 = "B",
                                     x2 = "C",
                                     enforce_symmetry = TRUE)
  p

  # Using four variables. Figure from: https://pubmed.ncbi.nlm.nih.gov/31327527/
  p <- SCpubr::do_CellularStatesPlot(sample = sample,
                                     gene_list = gene_set,
                                     x1 = "A",
                                     y1 = "C",
                                     x2 = "B",
                                     y2 = "D")
  p

  # Using four variables. Enforcing symmetry will make the scatter plot squared.
  p <- SCpubr::do_CellularStatesPlot(sample = sample,
                                     gene_list = gene_set,
                                     x1 = "A",
                                     y1 = "C",
                                     x2 = "B",
                                     y2 = "D",
                                     enforce_simmetry = TRUE)
  p
}
