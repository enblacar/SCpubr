#' Display the enriched terms for a given list of genes.
#'
#' @inheritParams doc_function
#' @param mat \strong{\code{\link[base]{list}}} | Result of over-representation test with clusterProfiler. Accepts only one result, be aware of that if you compute the test for all GO ontologies. Accessed through \strong{\code{mat@result}}.
#' @param n.chars \strong{\code{\link[base]{numeric}}} | Number of characters to use as a limit to wrap the term names. The higher this value, the longer the lines would be for each term in the plots. Defaults to 40.
#' @param n.terms \strong{\code{\link[base]{numeric}}} | Number of terms to display. Defaults to 25.
#' @return A dotplot object with enriched terms.
#' @export
#'
#' @example man/examples/examples_do_TermEnrichmentPlot.R
do_TermEnrichmentPlot <- function(mat,
                                  n.chars = 40,
                                  n.terms = 25,
                                  font.size = 14,
                                  font.type = "sans",
                                  plot.title = NULL,
                                  plot.subtitle = NULL,
                                  plot.caption = NULL,
                                  use_viridis = FALSE,
                                  viridis.palette = "G",
                                  viridis.direction = -1,
                                  sequential.palette = "YlGnBu",
                                  sequential.direction = 1,
                                  dot.scale = 8, 
                                  legend.type = "colorbar",
                                  legend.position = "bottom",
                                  legend.framewidth = 0.5,
                                  legend.tickwidth = 0.5,
                                  legend.length = 20,
                                  legend.width = 1,
                                  legend.framecolor = "grey50",
                                  legend.tickcolor = "white",
                                  number.breaks = 5,
                                  xlab = NULL,
                                  ylab = NULL,
                                  na.value = "grey75",
                                  grid.color = "grey90",
                                  grid.type = "dashed",
                                  plot.title.face = "bold",
                                  plot.subtitle.face = "plain",
                                  plot.caption.face = "italic",
                                  axis.title.face = "bold",
                                  axis.text.face = "plain",
                                  axis.text.x.angle = 45,
                                  legend.title.face = "bold",
                                  legend.text.face = "plain"){
    # Add lengthy error messages.
    withr::local_options(.new = list("warning.length" = 8170))
  
    check_suggests(function_name = "do_TermEnrichmentPlot")
    # Define pipe operator internally.
    `%>%` <- magrittr::`%>%`

    # Check logical parameters
    logical_list <- list("use_viridis" = use_viridis)
    check_type(parameters = logical_list, required_type = "logical", test_function = is.logical)
    # Check numeric parameters.
    numeric_list <- list("n.chars" = n.chars,
                         "n.terms" = n.terms,
                         "font.size" = font.size,
                         "legend.framewidth" = legend.framewidth,
                         "legend.tickwidth" = legend.tickwidth,
                         "legend.length" = legend.length,
                         "legend.width" = legend.width,
                         "viridis.direction" = viridis.direction,
                         "sequential.direction" = sequential.direction,
                         "number.breaks" = number.breaks,
                         "axis.text.x.angle" = axis.text.x.angle)
    check_type(parameters = numeric_list, required_type = "numeric", test_function = is.numeric)

    # Check character parameters.
    character_list <- list("viridis.palette" = viridis.palette,
                           "sequential.palette" = sequential.palette,
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
                           "legend.text.face" = legend.text.face,
                           "xlab" = xlab,
                           "ylab" = ylab,
                           "na.value" = na.value)
    check_type(parameters = character_list, required_type = "character", test_function = is.character)


    check_parameters(parameter = font.type, parameter_name = "font.type")
    check_parameters(parameter = legend.type, parameter_name = "legend.type")
    check_parameters(parameter = legend.position, parameter_name = "legend.position")
    check_parameters(parameter = plot.title.face, parameter_name = "plot.title.face")
    check_parameters(parameter = plot.subtitle.face, parameter_name = "plot.subtitle.face")
    check_parameters(parameter = plot.caption.face, parameter_name = "plot.caption.face")
    check_parameters(parameter = axis.title.face, parameter_name = "axis.title.face")
    check_parameters(parameter = axis.text.face, parameter_name = "axis.text.face")
    check_parameters(parameter = legend.title.face, parameter_name = "legend.title.face")
    check_parameters(parameter = legend.text.face, parameter_name = "legend.text.face")
    check_parameters(parameter = viridis.palette, parameter_name = "viridis.palette")
    check_parameters(parameter = sequential.palette, parameter_name = "sequential.palette")
    check_parameters(parameter = viridis.direction, parameter_name = "viridis.direction")
    check_parameters(parameter = sequential.direction, parameter_name = "sequential.direction")
    check_parameters(parameter = grid.type, parameter_name = "grid.type")
    
    # Check the colors provided to legend.framecolor and legend.tickcolor.
    check_colors(legend.framecolor, parameter_name = "legend.framecolor")
    check_colors(legend.tickcolor, parameter_name = "legend.tickcolor")
    check_colors(na.value, parameter_name = "na.value")
    check_colors(grid.color, parameter_name = "grid.color")
    
    
    # Check correct colnames.
    for (col_name in c("Description", "GeneRatio", "p.adjust", "Count")){
      assertthat::assert_that(col_name %in% colnames(mat),
                              msg = paste0(add_cross(), 
                                           crayon_body("Missing column "),
                                           crayon_key(col_name),
                                           crayon_body(" in "),
                                           crayon_key("mat"),
                                           crayon_body(".")))
    }
    
    
    # Generate color gradient.
    colors.gradient <- compute_continuous_palette(name = ifelse(isTRUE(use_viridis), viridis.palette, sequential.palette),
                                                  use_viridis = use_viridis,
                                                  direction = ifelse(isTRUE(use_viridis), viridis.direction, sequential.direction),
                                                  enforce_symmetry = FALSE)
    
    
    # PLOT
    
    # Start processing the matrix.
    p <- mat %>%
         dplyr::select(dplyr::all_of(c("Description", "GeneRatio", "p.adjust", "Count"))) %>% 
         # Turn character column GeneRatio into actual numeric GeneRatio.
         # -log10 transform p.adjust column.
         dplyr::mutate("GeneRatio" = unname(vapply(X = sapply(X = .data$GeneRatio, 
                                                              FUN = function(x){stringr::str_split(x, "/")}), 
                                                   FUN = function(x){as.numeric(x[1]) / as.numeric(x[2])}, 
                                                   FUN.VALUE = numeric(1))),
                       "p.adjust" = -log10(.data$p.adjust)) %>% 
         tibble::rownames_to_column(var = "Term") %>% 
         # Retrieve most significant ones.
         dplyr::arrange(dplyr::desc(.data$Count), dplyr::desc(.data$p.adjust)) %>% 
         # Turn Description column into a factor to get the values ordered.
         dplyr::mutate("Description" = factor(.data$Description, levels = rev(.data$Description))) %>% 
         tibble::as_tibble() %>% 
         dplyr::slice_head(n = n.terms) %>% 
         # Start plotting.
         ggplot2::ggplot(mapping = ggplot2::aes(x = .data$GeneRatio,
                                                y = .data$Description,
                                                fill = .data$p.adjust,
                                                size = .data$Count)) + 
         # Geom point.
         ggplot2::geom_point(shape = 21, 
                             color = "black") + 
         # Color scale.
         ggplot2::scale_fill_gradientn(colors = colors.gradient,
                                       na.value = na.value,
                                       name = expression(bold(paste("-", log["10"], "(p.adj)"))),
                                       breaks = scales::extended_breaks(n = number.breaks)) +
         # Add wrapping around Y labels.
         ggplot2::scale_y_discrete(labels = stringr::str_wrap(mat$Description, 
                                                              width = n.chars)) + 
         # Add a size scale.
         ggplot2::scale_size_continuous(range = c(3, dot.scale)) + 
         # Add labs.
         ggplot2::labs(title = plot.title,
                       subtitle = plot.subtitle,
                       caption = plot.caption,
                       x = ifelse(is.null(xlab), "Gene Ratio", xlab),
                       y = ifelse(is.null(ylab), "", ylab)) + 
         # Modify the legend aesthetics of size.
         ggplot2::guides(size = ggplot2::guide_legend(title = "Count",
                                                      title.position = "top",
                                                      title.hjust = 0.5,
                                                      ncol = 2,
                                                      override.aes = ggplot2::aes(fill = "black"))) + 
         # Add a base theme.
         ggplot2::theme_minimal(base_size = font.size) + 
         # Add theme customization.
         ggplot2::theme(axis.text.x = ggplot2::element_text(color = "black",
                                                            face = axis.text.face,
                                                            angle = get_axis_parameters(angle = axis.text.x.angle, flip = TRUE)[["angle"]],
                                                            hjust = get_axis_parameters(angle = axis.text.x.angle, flip = TRUE)[["hjust"]],
                                                            vjust = get_axis_parameters(angle = axis.text.x.angle, flip = TRUE)[["vjust"]]),
                        axis.text.y = ggplot2::element_text(face = axis.text.face, color = "black"),
                        axis.ticks = ggplot2::element_line(color = "black"),
                        axis.line.y = ggplot2::element_line(color = "black"),
                        axis.line.x = ggplot2::element_blank(),
                        axis.title = ggplot2::element_text(face = axis.title.face),
                        plot.title = ggplot2::element_text(face = plot.title.face, hjust = 0),
                        plot.subtitle = ggplot2::element_text(face = plot.subtitle.face, hjust = 0),
                        plot.caption = ggplot2::element_text(face = plot.caption.face, hjust = 1),
                        plot.title.position = "plot",
                        panel.grid.major.y = ggplot2::element_line(color = grid.color, linetype = grid.type),
                        panel.grid.major.x = ggplot2::element_blank(),
                        panel.grid = ggplot2::element_blank(),
                        text = ggplot2::element_text(family = font.type),
                        plot.caption.position = "plot",
                        legend.text = ggplot2::element_text(face = legend.text.face),
                        legend.position = legend.position,
                        legend.title = ggplot2::element_text(face = legend.title.face),
                        legend.justification = "center",
                        plot.margin = ggplot2::margin(t = 10, r = 10, b = 10, l = 10),
                        panel.grid.major = ggplot2::element_blank(),
                        plot.background = ggplot2::element_rect(fill = "white", color = "white"),
                        panel.background = ggplot2::element_rect(fill = "white", color = "white"),
                        legend.background = ggplot2::element_rect(fill = "white", color = "white"))
    
    # Modify fill legend to look nice.
    p <- modify_continuous_legend(p = p,
                                  # nocov start
                                  legend.title = expression(bold(paste("-", log["10"], "(p.adj)"))),
                                  # nocov end
                                  legend.aes = "fill",
                                  legend.type = legend.type,
                                  legend.position = legend.position,
                                  legend.length = legend.length,
                                  legend.width = legend.width,
                                  legend.framecolor = legend.framecolor,
                                  legend.tickcolor = legend.tickcolor,
                                  legend.framewidth = legend.framewidth,
                                  legend.tickwidth = legend.tickwidth)
    
    # Return the plot.
    return(p)
}

