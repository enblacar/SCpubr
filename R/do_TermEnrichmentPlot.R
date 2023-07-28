#' Display the enriched terms for a given list of genes.
#'
#' @inheritParams doc_function
#' @param enriched_terms \strong{\code{\link[base]{list}}} | List containing the output(s) of running Enrichr.
#' @param nchar_wrap \strong{\code{\link[base]{numeric}}} | Number of characters to use as a limit to wrap the term names. The higher this value, the longer the lines would be for each term in the plots. Defaults to 60.
#' @param nterms \strong{\code{\link[base]{numeric}}} | Number of terms to report for each database. Terms are arranged by adjusted p-value and selected from lowest to highest. Defaults to 5.
#' \itemize{
#'   \item \emph{\code{Enrichr}}.
#'   \item \emph{\code{FlyEnrichr}}.
#'   \item \emph{\code{WormEnrichr}}.
#'   \item \emph{\code{YeastEnrichr}}.
#'   \item \emph{\code{FishEnrichr}}.
#' }
#' @param colors.use \strong{\code{\link[base]{character}}} | Character vector of 2 colors (low and high ends of the color scale) to generate the gradient.
#' @param text_labels_size \strong{\code{\link[base]{numeric}}} | Controls how big or small labels are in the plot.
#' @return A ggplot2 object with enriched terms.
#' @export
#'
#' @example man/examples/examples_do_TermEnrichmentPlot.R
do_TermEnrichmentPlot <- function(enriched_terms,
                                  nchar_wrap = 20,
                                  nterms = 10,
                                  font.size = 14,
                                  font.type = "sans",
                                  plot.title = NULL,
                                  plot.subtitle = NULL,
                                  plot.caption = NULL,
                                  legend.position = "bottom",
                                  legend.type = "colorbar",
                                  colors.use = NULL,
                                  text_labels_size = 4,
                                  legend.length = 30,
                                  legend.width = 1,
                                  legend.framewidth = 0.5,
                                  legend.tickwidth = 0.5,
                                  legend.framecolor = "grey50",
                                  legend.tickcolor = "white",
                                  plot.title.face = "bold",
                                  plot.subtitle.face = "plain",
                                  plot.caption.face = "italic",
                                  axis.title.face = "bold",
                                  axis.text.face = "plain",
                                  legend.title.face = "bold",
                                  legend.text.face = "plain"){
    # Add lengthy error messages.
    withr::local_options(.new = list("warning.length" = 8170))
  
    check_suggests(function_name = "do_TermEnrichmentPlot")
    # Define pipe operator internally.
    `%>%` <- magrittr::`%>%`


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

    # Check character parameters.
    character_list <- list("colors.use" = colors.use,
                           "legend.position" = legend.position,
                           "legend.framecolor" = legend.framecolor,
                           "legend.tickcolor" = legend.tickcolor,
                           "legend.type" = legend.type,
                           "font.type" = font.type,
                           "plot.title" = plot.title,
                           "plot.subtitle" = plot.subtitle,
                           "plot.caption" = plot.caption,
                           "plot.title.face" = plot.title.face,
                           "plot.subtitle.face" = plot.subtitle.face,
                           "plot.caption.face" = plot.caption.face,
                           "axis.title.face" = axis.title.face,
                           "axis.text.face" = axis.text.face,
                           "legend.title.face" = legend.title.face,
                           "legend.text.face" = legend.text.face)
    check_type(parameters = character_list, required_type = "character", test_function = is.character)

    # Check colors.
    if (!is.null(colors.use)){
      check_colors(colors.use)
      assertthat::assert_that(length(colors.use) == 2,
                              msg = "Provide only 2 colors to colors.use.")
    } else {
      colors.use <- c("#bdc3c7", "#2c3e50")
    }

    check_parameters(parameter = font.type, parameter_name = "font.type")
    check_parameters(parameter = legend.type, parameter_name = "legend.type")
    check_parameters(parameter = legend.position, parameter_name = "legend.position")
    check_parameters(plot.title.face, parameter_name = "plot.title.face")
    check_parameters(plot.subtitle.face, parameter_name = "plot.subtitle.face")
    check_parameters(plot.caption.face, parameter_name = "plot.caption.face")
    check_parameters(axis.title.face, parameter_name = "axis.title.face")
    check_parameters(axis.text.face, parameter_name = "axis.text.face")
    check_parameters(legend.title.face, parameter_name = "legend.title.face")
    check_parameters(legend.text.face, parameter_name = "legend.text.face")

    # Check the colors provided to legend.framecolor and legend.tickcolor.
    check_colors(legend.framecolor, parameter_name = "legend.framecolor")
    check_colors(legend.tickcolor, parameter_name = "legend.tickcolor")

    output_list <- list()
    for (results in names(enriched_terms)){
      # Retrieve the data.
      data <- enriched_terms[[results]] %>%
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
           ggplot2::labs(title = ifelse(is.null(plot.title), stringr::str_replace_all(results, "_", " "), plot.title),
                         subtitle = plot.subtitle,
                         caption = plot.caption) +
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
                               inherit.aes = FALSE) +
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
                              inherit.aes = FALSE) +
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
          p <- modify_continuous_legend(p = p,
                                        legend.aes = "color",
                                        legend.type = legend.type,
                                        legend.position = legend.position,
                                        legend.length = legend.length,
                                        legend.width = legend.width,
                                        legend.framecolor = legend.framecolor,
                                        legend.tickcolor = legend.tickcolor,
                                        legend.framewidth = legend.framewidth,
                                        legend.tickwidth = legend.tickwidth)
          # Force to colorbar if only 1 p-value is present.
        } else if (legend.type == "colorbar" | legend.type == "colorsteps"){
          p <- modify_continuous_legend(p = p,
                                        legend.aes = "color",
                                        legend.type = "colorbar",
                                        legend.position = legend.position,
                                        legend.length = legend.length,
                                        legend.width = legend.width,
                                        legend.framecolor = legend.framecolor,
                                        legend.tickcolor = legend.tickcolor,
                                        legend.framewidth = legend.framewidth,
                                        legend.tickwidth = legend.tickwidth)
        }
      } else {
        p <- modify_continuous_legend(p = p,
                                      legend.aes = "fill",
                                      legend.type = legend.type,
                                      legend.position = legend.position,
                                      legend.length = legend.length,
                                      legend.width = legend.width,
                                      legend.framecolor = legend.framecolor,
                                      legend.tickcolor = legend.tickcolor,
                                      legend.framewidth = legend.framewidth,
                                      legend.tickwidth = legend.tickwidth)
      }
      p <- p +
           ggplot2::theme_minimal(base_size = font.size) +
           ggplot2::theme(axis.title = ggplot2::element_blank(),
                          axis.ticks = ggplot2::element_blank(),
                          axis.text = ggplot2::element_blank(),
                          panel.grid.major = ggplot2::element_blank(),
                          plot.title.position = "plot",
                          plot.title = ggplot2::element_text(face = plot.title.face,
                                                             hjust = ifelse(is.null(plot.title), 0.5, 0)),
                          plot.subtitle = ggplot2::element_text(face = plot.subtitle.face, hjust = 0),
                          plot.caption = ggplot2::element_text(face = plot.caption.face, hjust = 1),
                          legend.text = ggplot2::element_text(face = legend.text.face),
                          legend.title = ggplot2::element_text(face = legend.title.face),
                          panel.grid = ggplot2::element_blank(),
                          text = ggplot2::element_text(family = font.type),
                          plot.caption.position = "plot",
                          legend.position = legend.position,
                          legend.justification = "center",
                          plot.margin = ggplot2::margin(t = 10, r = 200, b = 10, l = 200),
                          plot.background = ggplot2::element_rect(fill = "white", color = "white"),
                          panel.background = ggplot2::element_rect(fill = "white", color = "white"),
                          legend.background = ggplot2::element_rect(fill = "white", color = "white"))

      output_list[[results]] <- p
    }

    return_object <- if (length(output_list) > 1) {output_list} else {p}
    return(return_object)
}

