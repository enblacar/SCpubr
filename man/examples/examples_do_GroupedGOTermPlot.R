\donttest{
  # Check Suggests.
  value <- SCpubr:::check_suggests(function_name = "do_GroupedGOTermPlot", passive = TRUE)

  if (isTRUE(value)){
    # Consult the full documentation in https://enblacar.github.io/SCpubr-book/

    # Need to load this library or equivalent.
    suppressMessages(library("org.Hs.eg.db"))

    # Define list of genes to query.
    genes.use <- c("CCR7", "CD14", "LYZ",
                   "S100A4", "MS4A1",
                   "MS4A7", "GNLY", "NKG7", "FCER1A",
                   "CST3", "PPBP")

    # Compute the grouped GO terms.
    out <- SCpubr::do_GroupedGOTermPlot(genes = genes.use,
                                        org.db = org.Hs.eg.db)
  } else if (isFALSE(value)){
    message("This function can not be used without its suggested packages.")
    message("Check out which ones are needed using `SCpubr::state_dependencies()`.")
  }
}
