#' Visualize Ligand-Receptor analysis output.
#'
#' This function makes use of [liana](https://github.com/saezlab/liana) package to run Ligand-Receptor analysis. Takes the output of liana and generates a dot-plot visualization according to the user's specifications.
#'
#' @param liana_output Object resulting from running \link[liana]{liana_wrap}.
#' @param split.by Character. Whether to further facet the plot on the y axis by common ligand.complex or receptor.complex. Values to provide: NULL, ligand.complex, receptor.complex.
#' @param keep_source,keep_target Character. Identities to keep for the source/target of the interactions. NULL otherwise.
#' @param top_interactions Numeric. Number of unique interactions to retrieve ordered by magnitude and specificity. It does not necessarily mean that the output will contain as many, but rather an approximate value.
#' @param dot_border Logical. Whether to draw a black border in the dots.
#' @param plot_grid Logical. Whether to plot grid lines.
#' @param border.color Character. Color for the border of the dots.
#' @param x_labels_angle Numeric. One of 0 (horizontal), 45 (diagonal), 90 (vertical). Adjusts to 0 if flip = FALSE and 45 if flip = TRUE.
#' @param rotate_strip_text Logical. Whether the text in the strips should be flipped 90 degrees.
#' @param legend.type Character. Type of legend to display. One of: normal, colorbar, colorsteps.
#' @param legend.position Position of the legend in the plot. Will only work if legend is set to TRUE.
#' @param legend.framewidth,legend.tickwidth Width of the lines of the box in the legend.
#' @param legend.framecolor,legend.tickcolor Color of the lines of the box in the legend.
#' @param legend.length,legend.width Length and width of the legend. Will adjust automatically depending on legend side.
#' @param viridis_color_map Character. A capital letter from A to H or the scale name as in \link[viridis]{scale_fill_viridis}.
#' @param flip Logical. Whether to invert the axis.
#' @param font.size Overall font.size of the plot.
#' @param font.type Character. Base font for the plot. One of mono, serif or sans.
#' @param dot.size Numeric. Size aesthetic for the dots.
#' @param grid.color Character. Color of the grid in the panels.
#' @param grid.type Character. One of the possible linetype options: blank, solid, dashed, dotted, dotdash, longdash, twodash.
#' @param significance_threshold Numeric. Value to filter the interactions by significance. Default is 0.05.
#'
#' @return A ggplot2 plot with the results of the Ligand-Receptor analysis.
#' @export
#'
#' @examples
#' \dontrun{
#' TBD
#' }
do_LigandReceptorPlot <- function(liana_output = NULL,
                                  split.by = NULL,
                                  keep_source = NULL,
                                  keep_target = NULL,
                                  top_interactions = 25,
                                  significance_threshold = 0.05,
                                  dot_border = TRUE,
                                  border.color = "black",
                                  x_labels_angle = 0,
                                  rotate_strip_text = FALSE,
                                  legend.position = "bottom",
                                  legend.type = "colorbar",
                                  legend.length = 20,
                                  legend.width = 1,
                                  legend.framecolor = "grey50",
                                  legend.tickcolor = "white",
                                  legend.framewidth = 1.5,
                                  legend.tickwidth = 1.5,
                                  viridis_color_map = "G",
                                  font.size = 14,
                                  dot.size = 1,
                                  font.type = "sans",
                                  flip = FALSE,
                                  plot_grid = FALSE,
                                  grid.color = "grey90",
                                  grid.type = "dotted"){

  # Checks for packages.
  check_suggests(function_name = "do_LigandReceptorPlot")
  `%>%` <- purrr::`%>%`

  # Check logical parameters.
  logical_list <- list("dot_border" = dot_border,
                       "flip" = flip,
                       "rotate_strip_text" = rotate_strip_text,
                       "plot_grid" = plot_grid)
  check_type(parameters = logical_list, required_type = "logical", test_function = is.logical)
  # Check numeric parameters.
  numeric_list <- list("font.size" = font.size,
                       "top_interactions" = top_interactions,
                       "legend.length" = legend.length,
                       "legend.width" = legend.width,
                       "legend.framewidth" = legend.framewidth,
                       "legend.tickwidth" = legend.tickwidth,
                       "dot.size" = dot.size,
                       "x_labels_angle" = x_labels_angle,
                       "significance_threshold" = significance_threshold)
  check_type(parameters = numeric_list, required_type = "numeric", test_function = is.numeric)
  # Check character parameters.
  character_list <- list("split.by" = split.by,
                         "keep_source" = keep_source,
                         "keep_target" = keep_target,
                         "border.color" = border.color,
                         "legend.position" = legend.position,
                         "legend.type" = legend.type,
                         "legend.position" = legend.position,
                         "legend.framecolor" = legend.framecolor,
                         "viridis_color_map" = viridis_color_map,
                         "legend.tickcolor" = legend.tickcolor,
                         "font.type" = font.type,
                         "grid.color" = grid.color,
                         "grid.type" = grid.type)
  check_type(parameters = character_list, required_type = "character", test_function = is.character)

  # Check border color.
  check_colors(border.color, parameter_name = "border.color")

  # Check viridis_color_map.
  check_viridis_color_map(viridis_color_map = viridis_color_map)

  # Check the colors provided to legend.framecolor and legend.tickcolor.
  check_colors(legend.framecolor, parameter_name = "legend.framecolor")
  check_colors(legend.tickcolor, parameter_name = "legend.tickcolor")
  check_colors(grid.color, parameter_name = "grid.color")

  # Check font.type.
  if (!(font.type %in% c("sans", "serif", "mono"))){
    stop("Please select one of the following for font.type: sans, serif, mono.", call. = F)
  }

  # Check font.type.
  if (!(grid.type %in% c("blank", "solid", "dashed", "dotted", "dotdash", "longdash", "twodash"))){
    stop("Please select one of the following for grid.type: blank, solid, dashed, dotted, dotdash, longdash, twodash.", call. = F)
  }

  # Check the legend.type.
  if (!(legend.type %in% c("normal", "colorbar", "colorsteps"))){
    stop("Please select one of the following for legend.type: normal, colorbar, colorsteps.", call. = FALSE)
  }

  # Check the legend.position.
  if (!(legend.position %in% c("top", "bottom", "left", "right"))){
    stop("Please select one of the following for legend.position: top, bottom, left, right.", call. = FALSE)
  }

  if (!is.null(split.by)){
    if (!(split.by %in% c("receptor.complex", "ligand.complex"))){
      stop("Please select one of the following for split.by: ligand.complex, receptor.complex.", call. = F)
    }
  }

  if (!(x_labels_angle %in% c(0, 45, 90))){
    stop("Please provide one of the following for x_labels_angle: 0, 45, 90.")
  }

  # Define legend parameters. Width and height values will change depending on the legend orientation.
  if (legend.position %in% c("top", "bottom")){
    legend.barwidth <- legend.length
    legend.barheight <- legend.width
    size_title <- "Interaction specificity"
    fill.title <- "Expression Magnitude"
  } else if (legend.position %in% c("left", "right")){
    legend.barwidth <- legend.width
    legend.barheight <- legend.length
    size_title <- stringr::str_wrap("Interaction specificity", width = 10)
    fill.title <- stringr::str_wrap("Expression Magnitude", width = 10)
  }

  if (x_labels_angle == 0){
    hjust = 0.5
    vjust = 1
  } else if (x_labels_angle == 45){
    hjust = 1
    vjust = 1
  } else if (x_labels_angle == 90){
    hjust = 1
    vjust = 0.5
  }

  liana_output <- liana_output %>%
                  liana::liana_aggregate(verbose = FALSE) %>%
                  dplyr::mutate(magnitude = .data$sca.LRscore) %>%
                  dplyr::mutate(specificity = .data$natmi.edge_specificity) %>%
                  dplyr::filter(.data$aggregate_rank <= significance_threshold) %>%
                  dplyr::arrange(dplyr::desc(.data$specificity), dplyr::desc(.data$magnitude))


  liana_output <- liana_output %>%
                  # Merge ligand.complex and receptor.complex columns into one, that will be used for the Y axis.
                  tidyr::unite(c("ligand.complex", "receptor.complex"),
                               col = "interaction",
                               sep = "<span style = 'color:grey50;'> | </span>",
                               remove = FALSE) %>%
                  # Merge source and target column into one, for future filtering.
                  tidyr::unite(c("source", "target"),
                               col = "interacting_clusters",
                               remove = FALSE)


  liana_output <- liana_output %>%
                  # Filter based on the top X interactions of ascending sensibilities.
                  # Retrieve the interaction column and return the unique values.
                  dplyr::filter(.data$interaction %in% {liana_output %>%
                                                        dplyr::pull(.data$interaction) %>%
                                                        unique()}[1:top_interactions])
  # If the user wants to trim the matrix and subset interacting entities.
  if (!(is.null(keep_source))){
    liana_output <- liana_output %>%
                    dplyr::filter(.data$source %in% keep_source)
  }

  if (!(is.null(keep_target))){
    liana_output <- liana_output %>%
                    dplyr::filter(.data$target %in% keep_target)
  }

  # Plot.
  if (isTRUE(flip)){
    if (isTRUE(dot_border)){
      p <-  liana_output %>%
            ggplot2::ggplot(mapping = ggplot2::aes(x = .data$interaction,
                                                   y = .data$target,
                                                   fill = .data$magnitude,
                                                   size = .data$specificity,
                                                   group = .data$interacting_clusters)) +
            ggplot2::geom_point(shape = 21)
    } else if (isFALSE(dot_border)) {
      p <-  liana_output %>%
            ggplot2::ggplot(mapping = ggplot2::aes(x = .data$interaction,
                                                   y = .data$target,
                                                   size = .data$specificity,
                                                   group = .data$interacting_clusters)) +
            ggplot2::geom_point(mapping = ggplot2::aes(color = .data$magnitude),
                                shape = 19)
    }
  } else if (isFALSE(flip)){
    if (isTRUE(dot_border)){
      p <-  liana_output %>%
            ggplot2::ggplot(mapping = ggplot2::aes(x = .data$target,
                                                   y = .data$interaction,
                                                   fill = .data$magnitude,
                                                   size = .data$specificity,
                                                   group = .data$interacting_clusters)) +
            ggplot2::geom_point(shape = 21)
    } else if (isFALSE(dot_border)){
      p <-  liana_output %>%
            ggplot2::ggplot(mapping = ggplot2::aes(x = .data$target,
                                                   y = .data$interaction,
                                                   size = .data$specificity,
                                                   group = .data$interacting_clusters)) +
            ggplot2::geom_point(mapping = ggplot2::aes(color = .data$magnitude),
                                shape = 19)
    }
  }
  p <-  p +
        ggplot2::scale_size_continuous(name = size_title,
                                       range = c(1 * dot.size, 5 * dot.size))
  # Settings for bordered dots.
  if (isTRUE(dot_border)){
    # Add color to aesthetics.
    p$layers[[1]]$aes_params$color <- border.color
    p <- p +
         ggplot2::scale_fill_viridis_c(option = viridis_color_map,
                                       name = fill.title)
  } else {
    p <- p +
         ggplot2::scale_color_viridis_c(option = viridis_color_map,
                                        name = fill.title)
  }
  # Continue plotting.
  if (isFALSE(flip)){
    if (isTRUE(rotate_strip_text)){
      strip_text_angle <- 90
    } else {
      strip_text_angle <- 0
    }
    if (is.null(split.by)){
      p <- p +
        ggplot2::facet_grid(. ~ .data$source,
                            space = "free",
                            scales = "free",
                            drop = FALSE)
    } else if (split.by == "ligand.complex"){
      p <- p +
        ggplot2::facet_grid(.data$ligand.complex ~ .data$source,
                            space = "free",
                            scales = "free",
                            drop = FALSE)
    } else if (split.by == "receptor.complex"){
      p <- p +
        ggplot2::facet_grid(.data$receptor.complex ~ .data$source,
                            space = "free",
                            scales = "free",
                            drop = FALSE)
    }
  } else if (isTRUE(flip)) {
    if (isTRUE(rotate_strip_text)){
      strip_text_angle <- 0
    } else {
      strip_text_angle <- 270
    }
    if (is.null(split.by)){
      p <- p +
        ggplot2::facet_grid(.data$source ~ .,
                            space = "free",
                            scales = "free",
                            drop = FALSE)
    } else if (split.by == "ligand.complex"){
      p <- p +
        ggplot2::facet_grid(.data$source ~ .data$ligand.complex,
                            space = "free",
                            scales = "free",
                            drop = FALSE)
    } else if (split.by == "receptor.complex"){
      p <- p +
        ggplot2::facet_grid(.data$source ~ .data$receptor.complex,
                            space = "free",
                            scales = "free",
                            drop = FALSE)
    }
  }



  p <- p +
       ggplot2::labs(title = "Source") +
       ggplot2::xlab(if (isTRUE(flip)){paste("Ligand", "|", "Receptor", sep = " ")} else if (isFALSE(flip)){"Target"}) +
       ggplot2::ylab(if (isFALSE(flip)){paste("Ligand", "|", "Receptor", sep = " ")} else if (isTRUE(flip)){"Target"}) +
       ggplot2::guides(size = ggplot2::guide_legend(title.position = "top",
                                                    title.hjust = 0.5,
                                                    override.aes = ggplot2::aes(fill = "black"))) +
       ggplot2::theme_minimal(base_size = font.size) +
       ggplot2::theme(plot.title = ggplot2::element_text(face = "bold",
                                                         hjust = if (isFALSE(flip)){0.5} else {1},
                                                         vjust = 0,
                                                         size = font.size),
                      plot.subtitle = ggplot2::element_text(hjust = 0),
                      plot.caption = ggplot2::element_text(hjust = 1),
                      plot.title.position = if (isFALSE(flip)){"panel"} else {"plot"},
                      plot.caption.position = "plot",
                      text = ggplot2::element_text(family = font.type),
                      legend.justification = "center",
                      legend.text = ggplot2::element_text(face = "bold"),
                      legend.title = ggplot2::element_text(face = "bold"),
                      legend.position = legend.position,
                      axis.title.x = ggplot2::element_text(face = "bold", hjust = 0.5),
                      axis.title.y = ggplot2::element_text(face = "bold", angle = 90),
                      axis.text.y = ggplot2::element_text(face = "bold"),
                      axis.text = ggplot2::element_text(face = "bold", color = "black"),
                      axis.ticks = ggplot2::element_line(color = "black"),
                      axis.text.x = if (isFALSE(flip)){
                        ggplot2::element_text(angle = x_labels_angle,
                                              hjust = hjust,
                                              vjust = vjust)
                        } else {
                          ggplot2::element_text(angle = x_labels_angle,
                                                   hjust = hjust,
                                                   vjust = vjust)
                        },
                      strip.text.x = if (isFALSE(flip)) {ggplot2::element_text(face = "bold",
                                                                               angle = strip_text_angle)}
                                     else {ggplot2::element_blank()},
                      strip.text.y = if (isFALSE(flip)) {ggplot2::element_blank()}
                                     else {ggplot2::element_text(face = "bold",
                                                                 angle = strip_text_angle)},
                      panel.border = ggplot2::element_rect(color = "black", fill = NA),
                      panel.grid = if (isTRUE(plot_grid)){ggplot2::element_line(color = grid.color, linetype = grid.type)} else {ggplot2::element_blank()},
                      plot.margin = ggplot2::margin(t = 10, r = 10, b = 10, l = 10),
                      plot.background = ggplot2::element_rect(fill = "white", color = "white"),
                      panel.background = ggplot2::element_rect(fill = "white", color = "black", linetype = "solid"),
                      legend.background = ggplot2::element_rect(fill = "white", color = "white"))

  # Adjust for the type of legend and whether it is fill or color.
  p <- modify_continuous_legend(p = p,
                                legend.aes = ifelse(isTRUE(dot_border), "fill", "color"),
                                legend.type = legend.type,
                                legend.position = legend.position,
                                legend.length = legend.length,
                                legend.width = legend.width,
                                legend.framecolor = legend.framecolor,
                                legend.tickcolor = legend.tickcolor,
                                legend.framewidth = legend.framewidth,
                                legend.tickwidth = legend.tickwidth)

  return(p)
}


