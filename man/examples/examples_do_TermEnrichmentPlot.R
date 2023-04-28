\donttest{
  # Check Suggests.
  value <- SCpubr:::check_suggests(function_name = "do_TermEnrichmentPlot", passive = TRUE)

  if (isTRUE(value)){
    # Consult the full documentation in https://enblacar.github.io/SCpubr-book/

    # Define your enriched terms.
    enriched_terms <- readRDS(system.file("extdata/enriched_terms_example.rds", package = "SCpubr"))
    enriched_terms$GO_Cellular_Component_2021 <- NULL
    enriched_terms$Azimuth_Cell_Types_2021 <- NULL

    # Default plot.
    p <- SCpubr::do_TermEnrichmentPlot(enriched_terms = enriched_terms)
    p
  } else if (base::isFALSE(value)){
    message("This function can not be used without its suggested packages.")
    message("Check out which ones are needed using `SCpubr::state_dependencies()`.")
  }
}
