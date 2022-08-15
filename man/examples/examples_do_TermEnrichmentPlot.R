\donttest{
  genes <- c("ABCB1", "ABCG2", "AHR", "AKT1", "AR")
  # General use case.
  p <- SCpubr::do_TermEnrichmentPlot(genes = genes,
                                     dbs_use = "C",
                                     ncol = 2)
  p

  # Modify total number of terms displayed.
  p <- SCpubr::do_TermEnrichmentPlot(genes = genes,
                                     dbs_use = "GO_Biological_Process_2021",
                                     nterms = 15)
  p

  # Modify the length of the terms displayed.
  p <- SCpubr::do_TermEnrichmentPlot(genes = genes,
                                     dbs_use = "GO_Biological_Process_2021",
                                     nterms = 15,
                                     nchar_wrap = 30)
  p
}
