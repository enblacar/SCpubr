\dontrun{
  # This function expects that you have run inferCNV on your
  # own and you have access to the output object.

  # Using a normal run - compute for all chromosome arms.
  out <- SCpubr::do_CopyNumberVariantPlot(sample = sample,
                                          infercnv_object = infercnv_object,
                                          using_metacells = FALSE)
  # Retrieve the UMAP for 1p region.
  out$`1p_umap`
  # Retrieve the dot plot for 1p region.
  out$`1p_dotplot`

  # Compute for a single chromosoome.
  out <- SCpubr::do_CopyNumberVariantPlot(sample = sample,
                                          infercnv_object = infercnv_object,
                                          using_metacells = FALSE,
                                          chromosome_focus = "2")

  # Compute for a run using metacells.
  ## metacell_mapping has to be a vector that, for each cell,
  ## tells which is the corresponding metacell.
  out <- SCpubr::do_CopyNumberVariantPlot(sample = sample,
                                          infercnv_object = infercnv_object_metacells,
                                          using_metacells = TRUE,
                                          metacell_mapping = metacell_mapping)

}
