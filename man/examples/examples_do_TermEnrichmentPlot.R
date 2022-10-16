\donttest{
  # Check Suggests.
  value <- SCpubr:::check_suggests(function_name = "do_TermEnrichmentPlot", passive = TRUE)

  if (isTRUE(value)){
    # Define your enriched terms.
    enriched_terms <- readRDS(system.file("extdata/enriched_terms_example.rds", package = "SCpubr"))
    enriched_terms$GO_Cellular_Component_2021 <- NULL
    enriched_terms$Azimuth_Cell_Types_2021 <- NULL

    # Default plot.
    p <- SCpubr::do_TermEnrichmentPlot(enriched_terms = enriched_terms)
    p

    # Increased number of terms.
    p <- SCpubr::do_TermEnrichmentPlot(enriched_terms = enriched_terms,
                                       nterms = 15)
    p

    # Control the length of the terms.
    p1 <- SCpubr::do_TermEnrichmentPlot(enriched_terms = enriched_terms,
                                        nterms = 15)
    p2 <- SCpubr::do_TermEnrichmentPlot(enriched_terms = enriched_terms,
                                        nterms = 15,
                                        nchar_wrap = 30)
    p <- p1 / p2
    p

    # Modify font size of the terms.
    p1 <- SCpubr::do_TermEnrichmentPlot(enriched_terms = enriched_terms)
    p2 <- SCpubr::do_TermEnrichmentPlot(enriched_terms = enriched_terms,
                                        text_labels_size = 6)

    p <- p1 / p2
    p
  } else if (isFALSE(value)){
    message("This function can not be used without its suggested packages.")
    message("Check out which ones are needed using `SCpubr::state_dependencies()`.")
  }
}
