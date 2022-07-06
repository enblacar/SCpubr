#' Compute the enriched terms for a given list of genes.
#'
#' @param genes Character vector with the genes to query.
#' @param dbs_use Character vector with the desired EnrichR databases to query. This can be consulted with `enrichR::listEnrichrDbs()`. Alternatively, one of the following values:
#' - "A": Will perform a query to 4 databases for cell types (Azimuth, Descartes, PanglaoDB and Descartes) and 4 databases for functional terms (MsigDB, GO-BP, GO-MF and KEGG). This is the default option if this parameter is not provided.
#' - "B": Performs a query for the cell type databases (Azimuth, Descartes, PanglaoDB and Descartes).
#' - "C": Performs a query for the functional terms (MsigDB, GO-BP, GO-MF and KEGG).
#' @param nchar_wrap Number of characters to use as a limit to wrap the term names. The higher this value, the longer the lines would be for each term in the plots. Defaults to 60.
#' @param nterms Number of terms to report for each database. Terms are arranged by adjusted p-value and selected from lowest to highest. Defaults to 5.
#' @param font.size Base font size for the plot. Defaults to 14.
#' @param font.type Character. Base font for the plot. One of mono, serif or sans.
#' @param site Site to query the genes against. Can be one of: "Enrichr", "FlyEnrichr", "WormEnrichr", "YeastEnrichr", "FishEnrichr".
#' @param colors.use Character vector of 2 colors (low and high ends of the color scale) to generate the gradient.
#' @param legend.length,legend.width Length and width of the legend. Will adjust automatically depending on legend side.
#' @param legend.framewidth,legend.tickwidth Width of the lines of the box in the legend.
#' @param legend.framecolor,legend.tickcolor Color of the lines of the box in the legend.
#' @param legend.type Character. Type of legend to display. One of: normal, colorbar, colorsteps. Colorsteps will default to colorbar if only one value is present for the scale.
#' @param legend.position Position of the legend in the plot.
#' @param text_labels_size Numeric. Controls how big or small labels are in the plot.
#' @return A ggplot2 object with enriched terms.
#' @export
#'
#' @example man/examples/examples_do_TermEnrichmentPlot.R
do_TermEnrichmentPlot <- function(genes,
                                  dbs_use,
                                  nchar_wrap = 20,
                                  nterms = 10,
                                  font.size = 14,
                                  font.type = "sans",
                                  site = "Enrichr",
                                  legend.position = "bottom",
                                  legend.type = "colorsteps",
                                  colors.use = NULL,
                                  text_labels_size = 4,
                                  legend.length = 30,
                                  legend.width = 1,
                                  legend.framewidth = 1.5,
                                  legend.tickwidth = 1.5,
                                  legend.framecolor = "grey50",
                                  legend.tickcolor = "white"){

  sink(tempfile())
  on.exit(sink())
  invisible(force({
    # Checks for packages.
    check_suggests(function_name = "do_TermEnrichmentPlot")


    # Define pipe operator internally.
    `%>%` <- purrr::`%>%`


    # Check numeric parameters.
    numeric_list <- list("nchar_wrap" = nchar_wrap,
                         "nterms" = nterms,
                         "font.size" = font.size,
                         "legend.framewidth" = legend.framewidth,
                         "legend.tickwidth" = legend.tickwidth,
                         "legend.length" = legend.length,
                         "legend.width" = legend.width,
                         "text_labels_size" = text_labels_size)
    check_type(parameters = numeric_list, required_type = "numeric", test_function = is.numeric)
    # Check logical parameters.
    #logical_list <- list("joint_plot" = joint_plot)
    #check_type(parameters = logical_list, required_type = "logical", test_function = is.logical)
    # Check character parameters.
    character_list <- list("genes" = genes,
                           "dbs_use" = dbs_use,
                           "site" = site,
                           "colors.use" = colors.use,
                           "legend.position" = legend.position,
                           "legend.framecolor" = legend.framecolor,
                           "legend.tickcolor" = legend.tickcolor,
                           "legend.type" = legend.type,
                           "font.type" = font.type)
    check_type(parameters = character_list, required_type = "character", test_function = is.character)

    # Check colors.
    if (!is.null(colors.use)){
      check_colors(colors.use)
      if (length(colors.use) != 2){
        stop("Provide only 2 colors.", call. = F)
      }
    } else {
      colors.use <- c("#bdc3c7", "#2c3e50")
    }

    # Check font.type.
    if (!(font.type %in% c("sans", "serif", "mono"))){
      stop("Please select one of the following for font.type: sans, serif, mono.", call. = F)
    }

    # Check the legend.type.
    if (!(legend.type %in% c("normal", "colorbar", "colorsteps"))){
      stop("Please select one of the following for legend.type: normal, colorbar, colorsteps.", call. = FALSE)
    }

    # Check the legend.position.
    if (!(legend.position %in% c("top", "bottom", "left", "right"))){
      stop("Please select one of the following for legend.position: top, bottom, left, right.", call. = FALSE)
    }

    # Define legend parameters. Width and height values will change depending on the legend orientation.
    if (legend.position %in% c("top", "bottom")){
      legend.barwidth <- legend.length
      legend.barheight <- legend.width
    } else if (legend.position %in% c("left", "right")){
      legend.barwidth <- legend.width
      legend.barheight <- legend.length
    }

    # Check the colors provided to legend.framecolor and legend.tickcolor.
    check_colors(legend.framecolor, parameter_name = "legend.framecolor")
    check_colors(legend.tickcolor, parameter_name = "legend.tickcolor")

    # Set necessary enrichR global options. This is copied from EnrichR code to avoid having to load the package.
    suppressMessages({
      options(enrichR.base.address = "https://maayanlab.cloud/Enrichr/")
      options(enrichR.live = TRUE)
      options(modEnrichR.use = TRUE)
      options(enrichR.sites.base.address = "https://maayanlab.cloud/")
      options(enrichR.sites = c("Enrichr", "FlyEnrichr", "WormEnrichr", "YeastEnrichr", "FishEnrichr"))

      # Set the search to Human genes.
      enrichR::setEnrichrSite(site = site)

      websiteLive <- TRUE
      dbs <- enrichR::listEnrichrDbs()
      dbs <- sort(dbs$libraryName)

      if (length(dbs_use) > 1){
        stop("Please provide only a database or one of the combinations (A, B, C).", call. = F)
      }
       else {
         if (dbs_use %in% c("A", "B", "C")){
           dbs_use <- dbs_use
         } else if (!(dbs_use %in% dbs)){
           stop("Please provide a database that is in enrichR. Please run sort(enrichR::listEnrichrDbs()[, 'libraryName']) to retrieve the full list of options.", call. = F)
         } else {
           dbs_use <- dbs_use
         }
      }
    })

    if (!(dbs_use %in% c("A", "B", "C"))){
      enriched <- enrichR::enrichr(genes, dbs_use)
      # Retrieve the data.
      data <- enriched[[dbs_use]] %>%
        dplyr::rowwise() %>%
        dplyr::mutate(Count = {length(unlist(stringr::str_split(.data$Genes, ";")))}) %>%
        dplyr::ungroup() %>%
        dplyr::arrange(.data$Adjusted.P.value) %>%
        dplyr::select(c("Adjusted.P.value", "Term", "Count")) %>%
        dplyr::slice_head(n = nterms) %>%
        dplyr::rowwise() %>% # Apply changes row-wise.
        dplyr::distinct(.data$Term, .keep_all = TRUE) %>% # Remove duplicated entries.
        dplyr::mutate(Term = factor(.data$Term, levels = .data$Term))


      max_value <- max(data[, "Count"]) + 10
      min_value <- -10

      limits <- c(min_value, max_value)

      # This chunk was extracted and adapted from: https://r-graph-gallery.com/296-add-labels-to-circular-barplot.html
      # on 05-07-2022
      # Get the name and the y position of each label
      label_data <- data

      # calculate the ANGLE of the labels
      number_of_bar <- nrow(label_data)
      label_data$id <- seq(1, number_of_bar)
      angle <-  90 - 360 * (label_data$id -0.5) / number_of_bar
      # calculate the alignment of labels: right or left
      # If I am on the left part of the plot, my labels have currently an angle < -90
      label_data$hjust <- ifelse(angle < -90, 1, 0)

      # flip angle BY to make them readable
      label_data$angle <- ifelse(angle < -90, angle + 180, angle)
      # Generate the plot. Based and adapted from: https://r-graph-gallery.com/web-circular-barplot-with-R-and-ggplot2.html
      # on 05-07-2022
      p <- ggplot2::ggplot(data) +
        # This generates lines for each value in Count. Once it's radial, they become circles.
        ggplot2::geom_hline(mapping = ggplot2::aes(yintercept = .data$Count),
                            color = "grey75",
                            linetype = "dashed") +
        # Geom col.
        ggplot2::geom_col(mapping = ggplot2::aes(x = .data$Term,
                                                 y = .data$Count,
                                                 fill = .data$Adjusted.P.value),
                          position = "dodge2",
                          alpha = 1) +
        # Add radial lines that will span further in the Y axis pointing to the term labels.
        ggplot2::geom_segment(mapping = ggplot2::aes(x = .data$Term,
                                                     y = max(.data$Count),
                                                     xend = .data$Term,
                                                     yend = max_value),
                              linetype = "dashed",
                              color = "grey75") +
        # Add "ticks" to each bar in the circular plot.
        ggplot2::geom_segment(mapping = ggplot2::aes(x = .data$Term,
                                                     y = 0,
                                                     xend = .data$Term,
                                                     yend = -0.5),
                              color = "black") +
        # Turn bar plot into circular plot.
        ggplot2::coord_polar(clip = "off") +
        # Set up the y axis limits so that we have empty space in the middle.
        ggplot2::ylim(limits) +
        ggplot2::ggtitle(stringr::str_replace_all(dbs_use, "_", " ")) +
        ggplot2::theme(plot.title = ggplot2::element_text(face = "bold", hjust = 0.5)) +
        # Add the enriched terms as labels with white background at the outside of the plot.
        ggplot2::geom_label(data = label_data,
                            mapping = ggplot2::aes(x = .data$id,
                                                   y = max_value,
                                                   label = stringr::str_wrap(.data$Term, nchar_wrap),
                                                   hjust = .data$hjust),
                            color = "black",
                            fill = "white",
                            label.size = NA,
                            label.padding = ggplot2::unit(0.50, "lines"),
                            fontface = "bold",
                            alpha = 1,
                            angle = 0,
                            size = text_labels_size,
                            inherit.aes = F) +
        # Add the number of genes in each term below y = 0.
        ggplot2::geom_text(data = label_data,
                           mapping = ggplot2::aes(x = .data$id,
                                                  y = -1.5,
                                                  label = .data$Count,
                                                  hjust = 0.5,
                                                  vjust = 0.5),
                           color = "black",
                           fontface = "bold",
                           size = text_labels_size,
                           alpha = 1,
                           angle = 0,
                           inherit.aes = F) +
        # Add black line at y = 0.
        ggplot2::geom_hline(yintercept = 0,
                            color = "black") +
        # Add X axis title in the center of the plot.
        ggplot2::annotate(geom = "text",
                          x = data$Term[1],
                          y = limits[1],
                          angle = 0,
                          hjust = 0.5,
                          vjust = 0.5,
                          size = text_labels_size,
                          label = stringr::str_wrap("Genes in each term", 9),
                          fontface = "bold") +
        ggplot2::scale_fill_gradient("Adj. P-value",
                                     low = colors.use[1],
                                     high = colors.use[2])
      # Add fill scale.
      if (length(unique(data$Adjusted.P.value)) == 1) {
        if (legend.type == "normal"){
          p <- p +
            ggplot2::guides(fill = ggplot2::guide_colorbar(title.position = "top",
                                                           title.hjust = 0.5))
        } else if (legend.type == "colorbar" | legend.type == "colorsteps"){
          p <- p +
            ggplot2::guides(fill = ggplot2::guide_colorbar(title.position = "top",
                                                           barwidth = legend.barwidth,
                                                           barheight = legend.barheight,
                                                           title.hjust = 0.5,
                                                           ticks.linewidth = legend.tickwidth,
                                                           frame.linewidth = legend.framewidth,
                                                           frame.colour = legend.framecolor,
                                                           ticks.colour = legend.tickcolor))
        }
      } else {
        if (legend.type == "normal"){
          p <- p +
            ggplot2::guides(fill = ggplot2::guide_colorbar(title.position = "top",
                                                           title.hjust = 0.5))
        } else if (legend.type == "colorbar"){
          p <- p +
            ggplot2::guides(fill = ggplot2::guide_colorbar(title.position = "top",
                                                           barwidth = legend.barwidth,
                                                           barheight = legend.barheight,
                                                           title.hjust = 0.5,
                                                           ticks.linewidth = legend.tickwidth,
                                                           frame.linewidth = legend.framewidth,
                                                           frame.colour = legend.framecolor,
                                                           ticks.colour = legend.tickcolor))
        } else if (legend.type == "colorsteps"){
          p <- p +
            ggplot2::guides(fill = ggplot2::guide_colorsteps(title.position = "top",
                                                             barwidth = legend.barwidth,
                                                             barheight = legend.barheight,
                                                             title.hjust = 0.5,
                                                             ticks.linewidth = legend.tickwidth,
                                                             frame.linewidth = legend.framewidth,
                                                             frame.colour = legend.framecolor,
                                                             ticks.colour = legend.tickcolor))
        }
      }
      p <- p +
        ggplot2::theme_minimal(base_size = font.size) +
        ggplot2::theme(axis.title = ggplot2::element_blank(),
                       axis.ticks = ggplot2::element_blank(),
                       axis.text = ggplot2::element_blank(),
                       panel.grid.major = ggplot2::element_blank(),
                       plot.title.position = "plot",
                       plot.title = ggplot2::element_text(face = "bold", hjust = 0.5),
                       plot.subtitle = ggplot2::element_text(hjust = 0),
                       plot.caption = ggplot2::element_text(hjust = 1),
                       panel.grid = ggplot2::element_blank(),
                       text = ggplot2::element_text(family = font.type),
                       plot.caption.position = "plot",
                       legend.text = ggplot2::element_text(face = "bold"),
                       legend.position = legend.position,
                       legend.title = ggplot2::element_text(face = "bold"),
                       legend.justification = "center",
                       plot.background = ggplot2::element_rect(fill = "white", color = "white"),
                       plot.margin = ggplot2::margin(t = 10, r = 200, b = 10, l = 200))

    } else if (dbs_use == "A"){
      p <- list("Azimuth_Cell_Types_2021" = SCpubr::do_TermEnrichmentPlot(genes = genes,
                                                                          dbs_use = "Azimuth_Cell_Types_2021",
                                                                          nchar_wrap = nchar_wrap,
                                                                          nterms = nterms,
                                                                          font.size = font.size,
                                                                          site = site,
                                                                          legend.position = legend.position,
                                                                          colors.use = colors.use,
                                                                          legend.length = legend.length,
                                                                          legend.width = legend.width,
                                                                          legend.framewidth = legend.framewidth,
                                                                          legend.tickwidth = legend.tickwidth,
                                                                          legend.framecolor = legend.framecolor,
                                                                          legend.tickcolor = legend.tickcolor),
                "CellMarker_Augmented_2021" = SCpubr::do_TermEnrichmentPlot(genes = genes,
                                                                          dbs_use = "CellMarker_Augmented_2021",
                                                                          nchar_wrap = nchar_wrap,
                                                                          nterms = nterms,
                                                                          font.size = font.size,
                                                                          site = site,
                                                                          legend.position = legend.position,
                                                                          colors.use = colors.use,
                                                                          legend.length = legend.length,
                                                                          legend.width = legend.width,
                                                                          legend.framewidth = legend.framewidth,
                                                                          legend.tickwidth = legend.tickwidth,
                                                                          legend.framecolor = legend.framecolor,
                                                                          legend.tickcolor = legend.tickcolor),
                "PanglaoDB_Augmented_2021" = SCpubr::do_TermEnrichmentPlot(genes = genes,
                                                                          dbs_use = "PanglaoDB_Augmented_2021",
                                                                          nchar_wrap = nchar_wrap,
                                                                          nterms = nterms,
                                                                          font.size = font.size,
                                                                          site = site,
                                                                          legend.position = legend.position,
                                                                          colors.use = colors.use,
                                                                          legend.length = legend.length,
                                                                          legend.width = legend.width,
                                                                          legend.framewidth = legend.framewidth,
                                                                          legend.tickwidth = legend.tickwidth,
                                                                          legend.framecolor = legend.framecolor,
                                                                          legend.tickcolor = legend.tickcolor),
                "Descartes_Cell_Types_and_Tissue_2021" = SCpubr::do_TermEnrichmentPlot(genes = genes,
                                                                          dbs_use = "Descartes_Cell_Types_and_Tissue_2021",
                                                                          nchar_wrap = nchar_wrap,
                                                                          nterms = nterms,
                                                                          font.size = font.size,
                                                                          site = site,
                                                                          legend.position = legend.position,
                                                                          colors.use = colors.use,
                                                                          legend.length = legend.length,
                                                                          legend.width = legend.width,
                                                                          legend.framewidth = legend.framewidth,
                                                                          legend.tickwidth = legend.tickwidth,
                                                                          legend.framecolor = legend.framecolor,
                                                                          legend.tickcolor = legend.tickcolor),
                "MSigDB_Hallmark_2020" = SCpubr::do_TermEnrichmentPlot(genes = genes,
                                                                          dbs_use = "MSigDB_Hallmark_2020",
                                                                          nchar_wrap = nchar_wrap,
                                                                          nterms = nterms,
                                                                          font.size = font.size,
                                                                          site = site,
                                                                          legend.position = legend.position,
                                                                          colors.use = colors.use,
                                                                          legend.length = legend.length,
                                                                          legend.width = legend.width,
                                                                          legend.framewidth = legend.framewidth,
                                                                          legend.tickwidth = legend.tickwidth,
                                                                          legend.framecolor = legend.framecolor,
                                                                          legend.tickcolor = legend.tickcolor),
                "GO_Biological_Process_2021" = SCpubr::do_TermEnrichmentPlot(genes = genes,
                                                                          dbs_use = "GO_Biological_Process_2021",
                                                                          nchar_wrap = nchar_wrap,
                                                                          nterms = nterms,
                                                                          font.size = font.size,
                                                                          site = site,
                                                                          legend.position = legend.position,
                                                                          colors.use = colors.use,
                                                                          legend.length = legend.length,
                                                                          legend.width = legend.width,
                                                                          legend.framewidth = legend.framewidth,
                                                                          legend.tickwidth = legend.tickwidth,
                                                                          legend.framecolor = legend.framecolor,
                                                                          legend.tickcolor = legend.tickcolor),
                "GO_Molecular_Function_2021" = SCpubr::do_TermEnrichmentPlot(genes = genes,
                                                                          dbs_use = "GO_Molecular_Function_2021",
                                                                          nchar_wrap = nchar_wrap,
                                                                          nterms = nterms,
                                                                          font.size = font.size,
                                                                          site = site,
                                                                          legend.position = legend.position,
                                                                          colors.use = colors.use,
                                                                          legend.length = legend.length,
                                                                          legend.width = legend.width,
                                                                          legend.framewidth = legend.framewidth,
                                                                          legend.tickwidth = legend.tickwidth,
                                                                          legend.framecolor = legend.framecolor,
                                                                          legend.tickcolor = legend.tickcolor),
                "KEGG_2021_Human" = SCpubr::do_TermEnrichmentPlot(genes = genes,
                                                                          dbs_use = "KEGG_2021_Human",
                                                                          nchar_wrap = nchar_wrap,
                                                                          nterms = nterms,
                                                                          font.size = font.size,
                                                                          site = site,
                                                                          legend.position = legend.position,
                                                                          colors.use = colors.use,
                                                                          legend.length = legend.length,
                                                                          legend.width = legend.width,
                                                                          legend.framewidth = legend.framewidth,
                                                                          legend.tickwidth = legend.tickwidth,
                                                                          legend.framecolor = legend.framecolor,
                                                                          legend.tickcolor = legend.tickcolor))

    } else if (dbs_use == "B"){
      p <- list("Azimuth_Cell_Types_2021" = SCpubr::do_TermEnrichmentPlot(genes = genes,
                                                                          dbs_use = "Azimuth_Cell_Types_2021",
                                                                          nchar_wrap = nchar_wrap,
                                                                          nterms = nterms,
                                                                          font.size = font.size,
                                                                          site = site,
                                                                          legend.position = legend.position,
                                                                          colors.use = colors.use,
                                                                          legend.length = legend.length,
                                                                          legend.width = legend.width,
                                                                          legend.framewidth = legend.framewidth,
                                                                          legend.tickwidth = legend.tickwidth,
                                                                          legend.framecolor = legend.framecolor,
                                                                          legend.tickcolor = legend.tickcolor),
                "CellMarker_Augmented_2021" = SCpubr::do_TermEnrichmentPlot(genes = genes,
                                                                            dbs_use = "CellMarker_Augmented_2021",
                                                                            nchar_wrap = nchar_wrap,
                                                                            nterms = nterms,
                                                                            font.size = font.size,
                                                                            site = site,
                                                                            legend.position = legend.position,
                                                                            colors.use = colors.use,
                                                                            legend.length = legend.length,
                                                                            legend.width = legend.width,
                                                                            legend.framewidth = legend.framewidth,
                                                                            legend.tickwidth = legend.tickwidth,
                                                                            legend.framecolor = legend.framecolor,
                                                                            legend.tickcolor = legend.tickcolor),
                "PanglaoDB_Augmented_2021" = SCpubr::do_TermEnrichmentPlot(genes = genes,
                                                                           dbs_use = "PanglaoDB_Augmented_2021",
                                                                           nchar_wrap = nchar_wrap,
                                                                           nterms = nterms,
                                                                           font.size = font.size,
                                                                           site = site,
                                                                           legend.position = legend.position,
                                                                           colors.use = colors.use,
                                                                           legend.length = legend.length,
                                                                           legend.width = legend.width,
                                                                           legend.framewidth = legend.framewidth,
                                                                           legend.tickwidth = legend.tickwidth,
                                                                           legend.framecolor = legend.framecolor,
                                                                           legend.tickcolor = legend.tickcolor),
                "Descartes_Cell_Types_and_Tissue_2021" = SCpubr::do_TermEnrichmentPlot(genes = genes,
                                                                                       dbs_use = "Descartes_Cell_Types_and_Tissue_2021",
                                                                                       nchar_wrap = nchar_wrap,
                                                                                       nterms = nterms,
                                                                                       font.size = font.size,
                                                                                       site = site,
                                                                                       legend.position = legend.position,
                                                                                       colors.use = colors.use,
                                                                                       legend.length = legend.length,
                                                                                       legend.width = legend.width,
                                                                                       legend.framewidth = legend.framewidth,
                                                                                       legend.tickwidth = legend.tickwidth,
                                                                                       legend.framecolor = legend.framecolor,
                                                                                       legend.tickcolor = legend.tickcolor))
    } else if (dbs_use == "C"){
      p <- list("MSigDB_Hallmark_2020" = SCpubr::do_TermEnrichmentPlot(genes = genes,
                                                                       dbs_use = "MSigDB_Hallmark_2020",
                                                                       nchar_wrap = nchar_wrap,
                                                                       nterms = nterms,
                                                                       font.size = font.size,
                                                                       site = site,
                                                                       legend.position = legend.position,
                                                                       colors.use = colors.use,
                                                                       legend.length = legend.length,
                                                                       legend.width = legend.width,
                                                                       legend.framewidth = legend.framewidth,
                                                                       legend.tickwidth = legend.tickwidth,
                                                                       legend.framecolor = legend.framecolor,
                                                                       legend.tickcolor = legend.tickcolor),
                "GO_Biological_Process_2021" = SCpubr::do_TermEnrichmentPlot(genes = genes,
                                                                             dbs_use = "GO_Biological_Process_2021",
                                                                             nchar_wrap = nchar_wrap,
                                                                             nterms = nterms,
                                                                             font.size = font.size,
                                                                             site = site,
                                                                             legend.position = legend.position,
                                                                             colors.use = colors.use,
                                                                             legend.length = legend.length,
                                                                             legend.width = legend.width,
                                                                             legend.framewidth = legend.framewidth,
                                                                             legend.tickwidth = legend.tickwidth,
                                                                             legend.framecolor = legend.framecolor,
                                                                             legend.tickcolor = legend.tickcolor),
                "GO_Molecular_Function_2021" = SCpubr::do_TermEnrichmentPlot(genes = genes,
                                                                             dbs_use = "GO_Molecular_Function_2021",
                                                                             nchar_wrap = nchar_wrap,
                                                                             nterms = nterms,
                                                                             font.size = font.size,
                                                                             site = site,
                                                                             legend.position = legend.position,
                                                                             colors.use = colors.use,
                                                                             legend.length = legend.length,
                                                                             legend.width = legend.width,
                                                                             legend.framewidth = legend.framewidth,
                                                                             legend.tickwidth = legend.tickwidth,
                                                                             legend.framecolor = legend.framecolor,
                                                                             legend.tickcolor = legend.tickcolor),
                "KEGG_2021_Human" = SCpubr::do_TermEnrichmentPlot(genes = genes,
                                                                  dbs_use = "KEGG_2021_Human",
                                                                  nchar_wrap = nchar_wrap,
                                                                  nterms = nterms,
                                                                  font.size = font.size,
                                                                  site = site,
                                                                  legend.position = legend.position,
                                                                  colors.use = colors.use,
                                                                  legend.length = legend.length,
                                                                  legend.width = legend.width,
                                                                  legend.framewidth = legend.framewidth,
                                                                  legend.tickwidth = legend.tickwidth,
                                                                  legend.framecolor = legend.framecolor,
                                                                  legend.tickcolor = legend.tickcolor))
    }
    return(p)
  }))
}

