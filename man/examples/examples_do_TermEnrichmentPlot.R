\dontrun{
  # Run this to generate a suitable input for the function.

  # Set necessary enrichR global options. This is copied from EnrichR code to avoid having to load the package.
  suppressMessages({
    options(enrichR.base.address = "https://maayanlab.cloud/Enrichr/")
    options(enrichR.live = TRUE)
    options(modEnrichR.use = TRUE)
    options(enrichR.sites.base.address = "https://maayanlab.cloud/")
    options(enrichR.sites = c("Enrichr", "FlyEnrichr", "WormEnrichr", "YeastEnrichr", "FishEnrichr"))

    # Set the search to Human genes.
    enrichR::setEnrichrSite(site = "Enrichr")

    websiteLive <- TRUE
    dbs <- enrichR::listEnrichrDbs()
    # Get all the possible databases to query.
    dbs <- sort(dbs$libraryName)
  })

  # Choose the dataset to query against.
  dbs_use <- c("GO_Biological_Process_2021",
               "GO_Cellular_Component_2021",
               "Azimuth_Cell_Types_2021")

  # List of genes to use as input.
  genes <- c("ABCB1", "ABCG2", "AHR", "AKT1", "AR")

  # Retrieve the enriched terms, this will be used as input for the function.
  enriched_terms <- enrichR::enrichr(genes, dbs_use)

  # Default plot.
  p <- SCpubr::do_TermEnrichmentPlot(enriched_terms = enriched_terms$GO_Biological_Process_2021)
  p

  # Increased number of terms.
  p <- SCpubr::do_TermEnrichmentPlot(enriched_terms = enriched_terms$GO_Biological_Process_2021,
                                     nterms = 15)
  p

  # Control the length of the terms.
  p1 <- SCpubr::do_TermEnrichmentPlot(enriched_terms = enriched_terms$GO_Biological_Process_2021,
                                      nterms = 15)
  p2 <- SCpubr::do_TermEnrichmentPlot(enriched_terms = enriched_terms$GO_Biological_Process_2021,
                                      nterms = 15,
                                      nchar_wrap = 30)
  p <- p1 / p2
  p

  # Modify font size of the terms.
  p1 <- SCpubr::do_TermEnrichmentPlot(enriched_terms = enriched_terms$GO_Biological_Process_2021,)
  p2 <- SCpubr::do_TermEnrichmentPlot(enriched_terms = enriched_terms$GO_Biological_Process_2021,
                                      text_labels_size = 6)

  p <- p1 / p2
  p
}
