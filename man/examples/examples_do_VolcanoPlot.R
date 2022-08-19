\dontrun{
  # Generate a volcano plot.
  p <- SCpubr::do_VolcanoPlot(sample = sample,
                              de_genes = de_genes)
  p

  # Modify cutoffs.
  p <- SCpubr::do_VolcanoPlot(sample = sample,
                              de_genes = de_genes,
                              pval_cutoff = 1e-15,
                              FC_cutoff = 0.2)
  p

  # Modify number of gene tags.
  p <- SCpubr::do_VolcanoPlot(sample = sample,
                              de_genes = de_genes,
                              n_genes = 15)
  p
}
