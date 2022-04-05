#' Compute the enriched terms for a given list of genes.
#'
#' @param genes Character vector with the genes to query.
#' @param dbs_use Character vector with the desired EnrichR databases to query. This can be consulted with `enrichR::listEnrichrDbs()`. Alternatively, one of the following values:
#' - "A": Will perform a query to 4 databases for cell types (Azimuth, Descartes, PanglaoDB and Descartes) and 4 databases for functional terms (MsigDB, GO-BP, GO-MF and KEGG). This is the default option if this parameter is not provided.
#' - "B": Performs a query for the cell type databases (Azimuth, Descartes, PanglaoDB and Descartes).
#' - "C": Performs a query for the functional terms (MsigDB, GO-BP, GO-MF and KEGG).
#' @param ncol Number of columns to group the output plots. Defaults to a around half of the length of the databases provided, if this is more than 1.
#' @param nchar_wrap Number of characters to use as a limit to wrap the term names. The higher this value, the longer the lines would be for each term in the plots. Defaults to 60.
#' @param nterms Number of terms to report for each database. Terms are arranged by adjusted p-value and selected from lowest to highest. Defaults to 5.
#' @param size Size of the dots. Defaults to 4.
#' @param fontsize Base font size for the plot. Defaults to 14.
#' @param site Site to query the genes against. Can be one of: "Enrichr", "FlyEnrichr", "WormEnrichr", "YeastEnrichr", "FishEnrichr".
#' @param colors.use Character vector of 2 colors (low and high ends of the color scale) to generate the gradient.
#'
#' @return A ggplot2 object with enriched terms.
#' @export
#'
#' @examples
do_EnrichedTermPlot <- function(genes,
                                dbs_use = NULL,
                                ncol = NULL,
                                nchar_wrap = 60,
                                nterms = 5,
                                size = 4,
                                fontsize = 14,
                                site = "Enrichr",
                                colors.use = c("#001FA9", "#00A6A9")){

  sink(tempfile())
  on.exit(sink())
  invisible(force({
    # Checks for packages.
    check_suggests(function_name = "do_EnrichedTermPlot")

    # Define pipe operator internally.
    `%>%` <- purrr::`%>%`


    # Check numeric parameters.
    numeric_list <- list("ncol" = ncol,
                         "nchar_wrap" = nchar_wrap,
                         "nterms" = nterms,
                         "size" = size,
                         "fontsize" = fontsize)
    check_type(parameters = numeric_list, required_type = "numeric", test_function = is.numeric)
    # Check character parameters.
    character_list <- list("genes" = genes,
                           "dbs_use" = dbs_use,
                           "site" = site)
    check_type(parameters = character_list, required_type = "character", test_function = is.character)

    # Define fontsize parameters.
    plot.title.fontsize <- fontsize + 2
    axis.text.fontsize <- fontsize
    axis.title.fontsize <- fontsize + 1
    legend.text.fontsize <- fontsize - 2
    legend.title.fontsize <- fontsize - 2

    # Check colors.
    check_colors(colors.use)

    # Set necessary enrichR global options. This is copied from EnrichR code to avoid having to load the package.
    suppressMessages({
      options(enrichR.base.address = "https://maayanlab.cloud/Enrichr/")
      options(enrichR.live = TRUE)
      options(modEnrichR.use = TRUE)
      options(enrichR.sites.base.address = "https://maayanlab.cloud/")
      options(enrichR.sites = c("Enrichr", "FlyEnrichr", "WormEnrichr", "YeastEnrichr", "FishEnrichr"))
      if (getOption("modEnrichR.use")) {
        enrichR::listEnrichrSites()
      } else {
        enrichR::getEnrichr(url=paste0(getOption("enrichR.base.address"), "datasetStatistics"))
        packageStartupMessage("Enrichr ... ", appendLF = FALSE)
        if (getOption("enrichR.live")) packageStartupMessage("Connection is Live!")
      }
      # Set the search to Human genes.
      enrichR::setEnrichrSite(site = site)

      websiteLive <- TRUE
      dbs <- enrichR::listEnrichrDbs()

      if (is.null(dbs_use)){
        dbs_use <- c("Azimuth_Cell_Types_2021",
                     "CellMarker_Augmented_2021",
                     "PanglaoDB_Augmented_2021",
                     "Descartes_Cell_Types_and_Tissue_2021",
                     "MSigDB_Hallmark_2020",
                     "GO_Biological_Process_2021",
                     "GO_Molecular_Function_2021",
                     "KEGG_2021_Human")
      } else {
        if (dbs_use == "A"){
          dbs_use <- c("Azimuth_Cell_Types_2021",
                       "CellMarker_Augmented_2021",
                       "PanglaoDB_Augmented_2021",
                       "Descartes_Cell_Types_and_Tissue_2021",
                       "MSigDB_Hallmark_2020",
                       "GO_Biological_Process_2021",
                       "GO_Molecular_Function_2021",
                       "KEGG_2021_Human")
        } else if (dbs_use == "B"){
          dbs_use <- c("Azimuth_Cell_Types_2021",
                       "CellMarker_Augmented_2021",
                       "PanglaoDB_Augmented_2021",
                       "Descartes_Cell_Types_and_Tissue_2021")
        } else if (dbs_use == "C"){
          dbs_use <- c("MSigDB_Hallmark_2020",
                       "GO_Biological_Process_2021",
                       "GO_Molecular_Function_2021",
                       "KEGG_2021_Human")
        } else {
          dbs_use <- dbs_use
        }
      }



      enriched <- enrichR::enrichr(genes, dbs_use)
    })



    list_enrichr <- list()
    for (database in names(enriched)){
      # Retrieve the data.
      data <- enriched[[database]] %>%
        dplyr::rowwise() %>%
        dplyr::mutate(Count = {length(unlist(stringr::str_split(.data$Genes, ";")))}) %>%
        dplyr::ungroup() %>%
        dplyr::arrange(Adjusted.P.value) %>%
        dplyr::select(c("Adjusted.P.value", "Term", "Count")) %>%
        dplyr::slice_head(n = nterms) %>%
        dplyr::mutate(Term = ifelse(nchar(Term) >= nchar_wrap, modify_string(Term), Term)) %>%
        dplyr::mutate(Term = factor(Term, levels = rev(Term)))

      # Generate the plot.
      p <- ggplot2::ggplot(data, mapping = ggplot2::aes(x = .data$Count, y = .data$Term, color = .data$Adjusted.P.value))  +
           ggplot2::scale_color_gradient(low = colors.use[1], colors.use[2], trans = 'reverse') +
           ggpubr::theme_pubr(legend = "right") +
           ggplot2::geom_point(size = size) +
           ggplot2::geom_point(shape = 1, size = size, colour = "black") +
           ggplot2::xlab("Genes supporting the Term") +
           ggplot2::ylab("Enriched Term") +
           ggplot2::ggtitle(database) +
           ggplot2::theme(axis.text = ggplot2::element_text(size = axis.text.fontsize, face = "bold"),
                          axis.title = ggplot2::element_text(size = axis.title.fontsize, face = "bold"),
                          plot.title = ggplot2::element_text(size = plot.title.fontsize, face = "bold", hjust = 0.5),
                          legend.title =  ggplot2::element_text(size = legend.title.fontsize, face = "bold")) +
           ggplot2::guides(size = ggplot2::guide_legend(title = "Adj. P-value"))

      # Store in the list.
      list_enrichr[[database]] <- p
    }

    # Check ncol.
    if (is.null(ncol)){
      if (length(dbs_use) > 1){
        ncol <- round(length(dbs_use) / 2, 0)
      } else {
        ncol <- 1
      }
    }

    # Put plots together.
    p <- patchwork::wrap_plots(list_enrichr, ncol = ncol, byrow = TRUE)  &
      ggplot2::theme(legend.text = ggplot2::element_text(size = legend.text.fontsize, face = "bold"))
    return(p)
  }))
}

