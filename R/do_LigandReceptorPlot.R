#' Visualize Ligand-Receptor analysis output.
#'
#' This function makes use of [liana](https://github.com/saezlab/liana) package to run Ligand-Receptor analysis. Takes the output of liana and generates a dot-plot visualization according to the user's specifications.
#'
#' @param sample Seurat object.
#' @param from_output Logical. Instead of providing a Seurat object, the user can provide the output of \link[liana]{liana_wrap}, thus reducing heavily computation time to get the plots.
#' @param liana_output Object resulting from running \link[liana]{liana_wrap}.
#' @param group.by Character. Only necessary when from_output is FALSE. Metadata variable to group the output by. If NULL, falls back to `Seurat::Idents(sample)`.
#' @param split.by Character. Whether to further facet the plot on the y axis by common ligand.complex or receptor.complex. Values to provide: NULL, ligand.complex, receptor.complex.
#' @param keep_source,keep_target Character. Identities to keep for the source/target of the interactions. NULL otherwise.
#' @param top_interactions Numeric. Number of unique interactions to retrieve ordered by magnitude and specificity. It does not necessarily mean that the output will contain as many, but rather an approximate value.
#' @param assay Character. Only necessary when from_output is FALSE. Assay that contains the normalized data.
#' @param dot_border Logical. Whether to draw a black border in the dots.
#' @param border.color Character. Color for the border of the dots.
#' @param x_labels_angle Numeric. One of 0 (horizontal), 45 (diagonal), 90 (vertical). Adjusts to 0 if flip = FALSE and 45 if flip = TRUE.
#' @param legend.type Character. Type of legend to display. One of: normal, colorbar, colorsteps.
#' @param legend.position Position of the legend in the plot. Will only work if legend is set to TRUE.
#' @param legend.framewidth,legend.tickwidth Width of the lines of the box in the legend.
#' @param legend.framecolor,legend.tickcolor Color of the lines of the box in the legend.
#' @param legend.length,legend.width Length and width of the legend. Will adjust automatically depending on legend side.
#' @param viridis_color_map Character. A capital letter from A to H or the scale name as in \link[viridis]{scale_fill_viridis}.
#' @param verbose Logical. Displays progress bars when running \link[liana]{liana_wrap}.
#' @param flip Logical. Whether to invert the axis.
#' @param font.size Overall font.size of the plot.
#' @param font.type Character. Base font for the plot. One of mono, serif or sans.
#' @param dot.size Numeric. Size aesthetic for the dots.
#'
#' @return A ggplot2 plot with the results of the Ligand-Receptor analysis.
#' @export
#'
#' @examples
#' \dontrun{
#' TBD
#' }
do_LigandReceptorPlot <- function(sample = NULL,
                                  from_output = FALSE,
                                  liana_output = NULL,
                                  group.by = NULL,
                                  split.by = NULL,
                                  keep_source = NULL,
                                  keep_target = NULL,
                                  top_interactions = 25,
                                  assay = "SCT",
                                  dot_border = TRUE,
                                  border.color = "black",
                                  x_labels_angle = 0,
                                  legend.position = "bottom",
                                  legend.type = "colorbar",
                                  legend.length = 20,
                                  legend.width = 1,
                                  legend.framecolor = "grey50",
                                  legend.tickcolor = "white",
                                  legend.framewidth = 1.5,
                                  legend.tickwidth = 1.5,
                                  viridis_color_map = "G",
                                  verbose = FALSE,
                                  font.size = 14,
                                  dot.size = 1,
                                  font.type = "sans",
                                  flip = FALSE){

  # Checks for packages.
  check_suggests(function_name = "do_LigandReceptorPlot")

  # Check logical parameters.
  logical_list <- list("from_output" = from_output,
                       "dot_border" = dot_border,
                       "verbose" = verbose,
                       "flip" = flip)
  check_type(parameters = logical_list, required_type = "logical", test_function = is.logical)
  # Check numeric parameters.
  numeric_list <- list("font.size" = font.size,
                       "top_interactions" = top_interactions,
                       "legend.length" = legend.length,
                       "legend.width" = legend.width,
                       "legend.framewidth" = legend.framewidth,
                       "legend.tickwidth" = legend.tickwidth,
                       "dot.size" = dot.size,
                       "x_labels_angle" = x_labels_angle)
  check_type(parameters = numeric_list, required_type = "numeric", test_function = is.numeric)
  # Check character parameters.
  character_list <- list("group.by" = group.by,
                         "split.by" = split.by,
                         "keep_source" = keep_source,
                         "keep_target" = keep_target,
                         "assay" = assay,
                         "border.color" = border.color,
                         "legend.position" = legend.position,
                         "legend.type" = legend.type,
                         "legend.position" = legend.position,
                         "legend.framecolor" = legend.framecolor,
                         "viridis_color_map" = viridis_color_map,
                         "legend.tickcolor" = legend.tickcolor,
                         "font.type" = font.type)
  check_type(parameters = character_list, required_type = "character", test_function = is.character)

  # Check border color.
  check_colors(border.color, parameter_name = "border.color")

  # Check viridis_color_map.
  check_viridis_color_map(viridis_color_map = viridis_color_map, verbose = verbose)

  # Check the colors provided to legend.framecolor and legend.tickcolor.
  check_colors(legend.framecolor, parameter_name = "legend.framecolor")
  check_colors(legend.tickcolor, parameter_name = "legend.tickcolor")

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

  if (!is.null(split.by)){
    if (!(split.by %in% c("receptor", "ligand"))){
      stop("Please select one of the following for split.by: ligand, receptor.", call. = F)
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

  # If the user wants it to be computed as well.
  if (isFALSE(from_output)){
    # Check if the sample provided is a Seurat object.
    check_Seurat(sample = sample)

    if (is.null(assay)){
      stop("Please provide an assay.", call. = F)
    }

    if (is.null(group.by)){
      sample$dummy <- Seurat::Idents(sample)
    } else {
      sample$dummy <- sample@meta.data[, group.by]
    }
    group.by <- "dummy"

    # Run liana.
    liana_output <- liana::liana_wrap(sce = sample,
                                      method = "cellphonedb",
                                      idents_col = group.by,
                                      verbose = verbose,
                                      assay = assay)

    liana_output <- liana_output %>%
                    dplyr::mutate(magnitude = .data$lr.mean) %>%
                    dplyr::mutate(specificity = .data$pvalue) %>%
                    dplyr::arrange(dplyr::desc(.data$magnitude), .data$specificity)

  # If the user provides the output from liana directly.
  } else if (isTRUE(from_output)){
    if ((sum(c("lr.mean", "pvalue") %in% colnames(liana_output))) != 2){
      stop("The column names of the liana object do not match what is expected from an output using cellphonedb as method.", call. = F)
    }
    liana_output <- liana_output %>%
                    dplyr::mutate(magnitude = .data$lr.mean) %>%
                    dplyr::mutate(specificity = .data$pvalue) %>%
                    dplyr::arrange(dplyr::desc(.data$magnitude), .data$specificity)
  }


  liana_output <- liana_output %>%
                  dplyr::filter(.data$specificity <= 0.05) %>%
                  dplyr::mutate(specificity = -log10(.data$specificity + 0.0000000001)) %>% # Log10 transform and add a small quantity for small values.
                  dplyr::arrange(dplyr::desc(.data$magnitude), dplyr::desc(.data$specificity)) %>%
                  # Merge ligand.complex and receptor.complex columns into one, that will be used for the Y axis.
                  tidyr::unite(c("ligand", "receptor"),
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
  p <-  liana_output %>%
        ggplot2::ggplot(mapping = ggplot2::aes(x = if (isFALSE(flip)){.data$target} else {.data$interaction},
                                               y = if (isFALSE(flip)){.data$interaction} else {.data$target},
                                               fill = if(isTRUE(dot_border)){.data$magnitude} else {NULL},
                                               size = .data$specificity,
                                               group = .data$interacting_clusters)) +
        ggplot2::geom_point(mapping = ggplot2::aes(color = if(isTRUE(dot_border)){NULL} else {.data$magnitude}),
                            shape = if(isTRUE(dot_border)){21} else {19}) +
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
    if (is.null(split.by)){
      p <- p +
        ggplot2::facet_grid(. ~ .data$source,
                            space = "free",
                            scales = "free")
    } else if (split.by == "ligand"){
      p <- p +
        ggplot2::facet_grid(.data$ligand ~ .data$source,
                            space = "free",
                            scales = "free")
    } else if (split.by == "receptor"){
      p <- p +
        ggplot2::facet_grid(.data$receptor ~ .data$source,
                            space = "free",
                            scales = "free")
    }
  } else if (isTRUE(flip)) {
    if (is.null(split.by)){
      p <- p +
        ggplot2::facet_grid(.data$source ~ .,
                            space = "free",
                            scales = "free")
    } else if (split.by == "ligand"){
      p <- p +
        ggplot2::facet_grid(.data$source ~ .data$ligand,
                            space = "free",
                            scales = "free")
    } else if (split.by == "receptor"){
      p <- p +
        ggplot2::facet_grid(.data$source ~ .data$receptor,
                            space = "free",
                            scales = "free")
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
                      plot.subtitle = ggtext::element_markdown(hjust = 0),
                      plot.caption = ggtext::element_markdown(hjust = 1),
                      plot.title.position = if (isFALSE(flip)){"panel"} else {"plot"},
                      plot.caption.position = "plot",
                      text = ggplot2::element_text(family = font.type),
                      legend.justification = "center",
                      legend.text = ggplot2::element_text(face = "bold"),
                      legend.title = ggplot2::element_text(face = "bold"),
                      legend.position = legend.position,
                      axis.title.x = ggplot2::element_text(face = "bold", hjust = 0.5),
                      axis.title.y = ggplot2::element_text(face = "bold", angle = 90),
                      axis.text.y = ggtext::element_markdown(face = "bold"),
                      axis.text = ggplot2::element_text(face = "bold", color = "black"),
                      axis.ticks = ggplot2::element_line(color = "black"),
                      axis.text.x = if (isFALSE(flip)){
                        ggplot2::element_text(angle = x_labels_angle,
                                              hjust = hjust,
                                              vjust = vjust)
                        } else {
                          ggtext::element_markdown(angle = x_labels_angle,
                                                   hjust = hjust,
                                                   vjust = vjust)
                        },
                      strip.text.x = if (isFALSE(flip)) {ggplot2::element_text(face = "bold")} else {ggplot2::element_blank()},
                      strip.text.y = if (isFALSE(flip)) {ggplot2::element_blank()} else {ggplot2::element_text(face = "bold")},
                      panel.grid = ggplot2::element_blank(),
                      plot.margin = ggplot2::margin(t = 10, r = 10, b = 10, l = 10),
                      plot.background = ggplot2::element_rect(fill = "white", color = "white"),
                      panel.background = ggplot2::element_rect(fill = "white", color = "black"),
                      legend.background = ggplot2::element_rect(fill = "white", color = "white"))

  # Adjust for the type of legend and whether it is fill or color.
  if (legend.type == "normal"){
    if (isTRUE(dot_border)){
      p <- p +
           ggplot2::guides(fill = ggplot2::guide_colorbar(title.position = "top",
                                                          title.hjust = 0.5))
    } else {
      p <- p +
           ggplot2::guides(color = ggplot2::guide_colorbar(title.position = "top",
                                                           title.hjust = 0.5))
    }

  } else if (legend.type == "colorbar"){
    if (isTRUE(dot_border)){
      p <- p +
           ggplot2::guides(fill = ggplot2::guide_colorbar(title.position = "top",
                                                          barwidth = legend.barwidth,
                                                          barheight = legend.barheight,
                                                          title.hjust = 0.5,
                                                          ticks.linewidth = legend.tickwidth,
                                                          frame.linewidth = legend.framewidth,
                                                          frame.colour = legend.framecolor,
                                                          ticks.colour = legend.tickcolor))
    } else {
      p <- p +
           ggplot2::guides(color = ggplot2::guide_colorbar(title.position = "top",
                                                           barwidth = legend.barwidth,
                                                           barheight = legend.barheight,
                                                           title.hjust = 0.5,
                                                           ticks.linewidth = legend.tickwidth,
                                                           frame.linewidth = legend.framewidth,
                                                           frame.colour = legend.framecolor,
                                                           ticks.colour = legend.tickcolor))
    }
  } else if (legend.type == "colorsteps"){
    if (isTRUE(dot_border)){
      p <- p +
           ggplot2::guides(fill = ggplot2::guide_colorsteps(title.position = "top",
                                                            barwidth = legend.barwidth,
                                                            barheight = legend.barheight,
                                                            title.hjust = 0.5,
                                                            ticks.linewidth = legend.tickwidth,
                                                            frame.linewidth = legend.framewidth,
                                                            frame.colour = legend.framecolor,
                                                            ticks.colour = legend.tickcolor))
    } else {
      p <- p +
           ggplot2::guides(color = ggplot2::guide_colorsteps(title.position = "top",
                                                             barwidth = legend.barwidth,
                                                             barheight = legend.barheight,
                                                             title.hjust = 0.5,
                                                             ticks.linewidth = legend.tickwidth,
                                                             frame.linewidth = legend.framewidth,
                                                             frame.colour = legend.framecolor,
                                                             ticks.colour = legend.tickcolor))
    }
  }
  return(p)
}


