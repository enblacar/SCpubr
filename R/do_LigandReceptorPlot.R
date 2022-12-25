#' Visualize Ligand-Receptor analysis output.
#'
#' This function makes use of [liana](https://github.com/saezlab/liana) package to run Ligand-Receptor analysis. Takes the output of liana and generates a dot-plot visualization according to the user's specifications.
#'
#' @inheritParams doc_function
#' @param liana_output \strong{\code{\link[tibble]{tibble}}} | Object resulting from running \link[liana]{liana_wrap}.
#' @param split.by \strong{\code{\link[base]{character}}} | Whether to further facet the plot on the y axis by common ligand.complex or receptor.complex. Values to provide: NULL, ligand.complex, receptor.complex.
#' @param keep_source,keep_target \strong{\code{\link[base]{character}}} | Identities to keep for the source/target of the interactions. NULL otherwise.
#' @param top_interactions \strong{\code{\link[base]{numeric}}} | Number of unique interactions to retrieve ordered by magnitude and specificity. It does not necessarily mean that the output will contain as many, but rather an approximate value.
#' @param dot_border \strong{\code{\link[base]{logical}}} | Whether to draw a black border in the dots.
#' @param rotate_strip_text \strong{\code{\link[base]{logical}}} | Whether the text in the strips should be flipped 90 degrees.
#' @param dot.size \strong{\code{\link[base]{numeric}}} | Size aesthetic for the dots.
#' @param compute_ChordDiagrams \strong{\code{\link[base]{logical}}} | Whether to also compute Chord Diagrams for both the number of interactions between source and target but also between ligand.complex and receptor.complex.
#' @param add_missing_LR_combinations \strong{\code{\link[base]{logical}}} | Depending on the value provided to \strong{\code{top_interactions}}, there might be some source-target combinations missing. If set to TRUE, those combinations will be brought back to the plot as NA values.
#' @param arrange_interactions_by \strong{\code{\link[base]{character}}} | Select the method in which the output matrix of interactions is arranged for plotting. One of:
#' \itemize{
#'   \item \emph{\code{aggregate_rank}}: Uses the output matrix of \strong{\code{\link[liana]{liana_aggregate}}}. Interactions are ordered based on \strong{\code{aggregate_rank}}.
#'   \item \emph{\code{specificity}}: Uses the \strong{\code{natmi.edge_specificity}} column to arrange the interactions by their specificity.
#'   \item \emph{\code{magnitude}}: Uses the \strong{\code{sca.LRscore}} column to arrange the interactions by their specificity.
#'   \item \emph{\code{both}}: Uses both \strong{\code{sca.LRscore}} and \strong{\code{natmi.edge_specificity}} columns to arrange the interactions by their specificity and magnitude altogether.
#' }
#' @param sort_interactions_alphabetically \strong{\code{\link[base]{logical}}} | Sort the interactions to be plotted alphabetically (\strong{\code{TRUE}}) or keep them in their original order in the matrix (\strong{\code{FALSE}}).

#' @return A ggplot2 plot with the results of the Ligand-Receptor analysis.
#' @export
#'
#' @example /man/examples/examples_do_LigandReceptorPlot.R

do_LigandReceptorPlot <- function(liana_output,
                                  split.by = NULL,
                                  keep_source = NULL,
                                  keep_target = NULL,
                                  top_interactions = 25,
                                  dot_border = TRUE,
                                  border.color = "black",
                                  rotate_x_axis_labels = 45,
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
                                  viridis_direction = 1,
                                  font.size = 14,
                                  dot.size = 1,
                                  font.type = "sans",
                                  flip = FALSE,
                                  plot.grid = TRUE,
                                  grid.color = "grey90",
                                  grid.type = "dotted",
                                  compute_ChordDiagrams = FALSE,
                                  add_missing_LR_combinations = TRUE,
                                  arrange_interactions_by = "both",
                                  sort_interactions_alphabetically = FALSE){

  # Checks for packages.
  check_suggests(function_name = "do_LigandReceptorPlot")
  `%>%` <- magrittr::`%>%`

  # Check logical parameters.
  logical_list <- list("dot_border" = dot_border,
                       "flip" = flip,
                       "rotate_strip_text" = rotate_strip_text,
                       "plot.grid" = plot.grid,
                       "add_missing_LR_combinations" = add_missing_LR_combinations,
                       "sort_interactions_alphabetically" = sort_interactions_alphabetically)
  check_type(parameters = logical_list, required_type = "logical", test_function = is.logical)
  # Check numeric parameters.
  numeric_list <- list("font.size" = font.size,
                       "top_interactions" = top_interactions,
                       "legend.length" = legend.length,
                       "legend.width" = legend.width,
                       "legend.framewidth" = legend.framewidth,
                       "legend.tickwidth" = legend.tickwidth,
                       "dot.size" = dot.size,
                       "rotate_x_axis_labels" = rotate_x_axis_labels,
                       "viridis_direction" = viridis_direction)
  check_type(parameters = numeric_list, required_type = "numeric", test_function = is.numeric)
  # Check character parameters.
  character_list <- list("split.by" = split.by,
                         "keep_source" = keep_source,
                         "keep_target" = keep_target,
                         "border.color" = border.color,
                         "legend.position" = legend.position,
                         "legend.type" = legend.type,
                         "legend.framecolor" = legend.framecolor,
                         "viridis_color_map" = viridis_color_map,
                         "legend.tickcolor" = legend.tickcolor,
                         "font.type" = font.type,
                         "grid.color" = grid.color,
                         "grid.type" = grid.type,
                         "arrange_interactions_by" = arrange_interactions_by)
  check_type(parameters = character_list, required_type = "character", test_function = is.character)

  # Check border color.
  check_colors(border.color, parameter_name = "border.color")

  # Check viridis_color_map.
  check_viridis_color_map(viridis_color_map = viridis_color_map)

  # Check the colors provided to legend.framecolor and legend.tickcolor.
  check_colors(legend.framecolor, parameter_name = "legend.framecolor")
  check_colors(legend.tickcolor, parameter_name = "legend.tickcolor")
  check_colors(grid.color, parameter_name = "grid.color")

  check_parameters(parameter = font.type, parameter_name = "font.type")
  check_parameters(parameter = legend.type, parameter_name = "legend.type")
  check_parameters(parameter = legend.position, parameter_name = "legend.position")
  check_parameters(parameter = viridis_direction, parameter_name = "viridis_direction")
  check_parameters(parameter = viridis_color_map, parameter_name = "viridis_color_map")
  check_parameters(parameter = grid.type, parameter_name = "grid.type")
  check_parameters(parameter = rotate_x_axis_labels, parameter_name = "rotate_x_axis_labels")
  check_parameters(parameter = arrange_interactions_by, parameter_name = "arrange_interactions_by")

  if (!is.null(split.by)){
    assertthat::assert_that(split.by %in% c("receptor.complex", "ligand.complex"),
                            msg = "Please select one of the following for split.by: ligand.complex, receptor.complex.")
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

  liana_output <- liana_output %>%
                  liana::liana_aggregate(verbose = FALSE)

  # This is to later on add missing interacting pairs.
  possible_interacting_clusters <- c()

  # If we are subsetting the final plot.
  if (isTRUE(add_missing_LR_combinations)){
    possible_sources <- if(is.null(keep_source)){sort(unique(liana_output$source))} else {sort(unique(liana_output$source))[sort(unique(liana_output$source)) %in% keep_source]}
    possible_targets <- if(is.null(keep_target)){sort(unique(liana_output$target))} else {sort(unique(liana_output$target))[sort(unique(liana_output$target)) %in% keep_target]}

    for (source in possible_sources){
      for (target in possible_targets){
        name <- paste0(source, "_", target)
        possible_interacting_clusters <- c(possible_interacting_clusters, name)
      }
    }
  }

  liana_output <- liana_output %>%
                  dplyr::mutate(magnitude = .data$sca.LRscore) %>%
                  dplyr::mutate(specificity = .data$natmi.edge_specificity)

  # Differential arrangement of the interactions.
  if (arrange_interactions_by == "aggregate_rank"){
    liana_output <- liana_output %>%
                    dplyr::arrange(.data$aggregate_rank)
  } else if (arrange_interactions_by == "specificity"){
    liana_output <- liana_output %>%
                    dplyr::arrange(dplyr::desc(.data$specificity))
  } else if (arrange_interactions_by == "magnitude"){
    liana_output <- liana_output %>%
                    dplyr::arrange(dplyr::desc(.data$magnitude))
  } else if (arrange_interactions_by == "both"){
    liana_output <- liana_output %>%
                    dplyr::arrange(dplyr::desc(.data$specificity), dplyr::desc(.data$magnitude))
  }

  liana_output <- liana_output %>%
                  # Merge ligand.complex and receptor.complex columns into one, that will be used for the Y axis.
                  tidyr::unite(c("ligand.complex", "receptor.complex"),
                               col = "interaction",
                               sep = " | ",
                               remove = FALSE) %>%
                  # Merge source and target column into one, for future filtering.
                  tidyr::unite(c("source", "target"),
                               col = "interacting_clusters",
                               remove = FALSE)
  # For Chord diagrams.
  output_copy <- liana_output %>% dplyr::filter(.data$aggregate_rank <= 0.05)


  liana_output <- liana_output %>%
                  # Filter based on the top X interactions of ascending sensibilities.
                  dplyr::inner_join(y = {liana_output %>%
                                         dplyr::distinct_at(c("ligand.complex", "receptor.complex")) %>%
                                         dplyr::slice_head(n = top_interactions)},
                                    by = c("ligand.complex", "receptor.complex"))
  if (isTRUE(add_missing_LR_combinations)){
    # Fix to add missing "NA" interactions.
    liana_output <- liana_output %>%
      dplyr::select(dplyr::all_of(c("interacting_clusters",
                                    "source",
                                    "target",
                                    "interaction",
                                    "ligand.complex",
                                    "receptor.complex",
                                    "magnitude",
                                    "specificity")))

    # Iterate over each possible interaction and each interacting pair.
    not_found_interaction_pairs <- possible_interacting_clusters[possible_interacting_clusters %!in% unique(liana_output$interacting_clusters)]
    interactions <- unique(liana_output$interaction)

    # Iterate over each interaction.
    for(interaction in interactions){
      ligand.complex <- stringr::str_split(interaction, pattern = " \\| ")[[1]][1]
      receptor.complex <- stringr::str_split(interaction, pattern = " \\| ")[[1]][2]
      # For each missing interaction pair.
      for (interacting_cluster in not_found_interaction_pairs){
        source <- stringr::str_split(interacting_cluster, pattern = "_")[[1]][1]
        target <- stringr::str_split(interacting_cluster, pattern = "_")[[1]][2]

        # If the interacting pair - interaction is missing, add a mock row with it.
        column <- tibble::tibble("interacting_clusters" = interacting_cluster,
                                 "source" = source,
                                 "target" = target,
                                 "interaction" = interaction,
                                 "ligand.complex" = ligand.complex,
                                 "receptor.complex" = receptor.complex,
                                 "magnitude" = NA,
                                 "specificity" = NA)
        liana_output <- liana_output %>% dplyr::bind_rows(column)
      }
    }
  }

  # If the user wants to trim the matrix and subset interacting entities.
  if (!(is.null(keep_source))){
    liana_output <- liana_output %>%
                    dplyr::filter(.data$source %in% keep_source)
    output_copy <- output_copy %>%
                   dplyr::filter(.data$source %in% keep_source)
  }

  if (!(is.null(keep_target))){
    liana_output <- liana_output %>%
                    dplyr::filter(.data$target %in% keep_target)
    output_copy <- output_copy %>%
                   dplyr::filter(.data$target %in% keep_target)
  }

  # Make source and target factors, so that they do not get dropped by the plot.
  if (isTRUE(sort_interactions_alphabetically)){
    liana_output$source <- factor(liana_output$source, levels = sort(unique(liana_output$source)))
    liana_output$target <- factor(liana_output$target, levels = sort(unique(liana_output$target)))
    liana_output$interaction <- factor(liana_output$interaction, levels = rev(sort(unique(liana_output$interaction))))
  } else if (isFALSE(sort_interactions_alphabetically)){
    liana_output$source <- factor(liana_output$source, levels = sort(unique(liana_output$source)))
    liana_output$target <- factor(liana_output$target, levels = sort(unique(liana_output$target)))
    liana_output$interaction <- factor(liana_output$interaction, levels = rev(unique(liana_output$interaction)))
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
            ggplot2::geom_point(shape = 21,
                                na.rm = TRUE)
    } else if (isFALSE(dot_border)) {
      p <-  liana_output %>%
            ggplot2::ggplot(mapping = ggplot2::aes(x = .data$interaction,
                                                   y = .data$target,
                                                   size = .data$specificity,
                                                   group = .data$interacting_clusters)) +
            ggplot2::geom_point(mapping = ggplot2::aes(color = .data$magnitude),
                                shape = 19,
                                na.rm = TRUE)
    }
  } else if (isFALSE(flip)){
    if (isTRUE(dot_border)){
      p <-  liana_output %>%
            ggplot2::ggplot(mapping = ggplot2::aes(x = .data$target,
                                                   y = .data$interaction,
                                                   fill = .data$magnitude,
                                                   size = .data$specificity,
                                                   group = .data$interacting_clusters)) +
            ggplot2::geom_point(shape = 21,
                                na.rm = TRUE)
    } else if (isFALSE(dot_border)){
      p <-  liana_output %>%
            ggplot2::ggplot(mapping = ggplot2::aes(x = .data$target,
                                                   y = .data$interaction,
                                                   size = .data$specificity,
                                                   group = .data$interacting_clusters)) +
            ggplot2::geom_point(mapping = ggplot2::aes(color = .data$magnitude),
                                shape = 19,
                                na.rm = TRUE)
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
                                       name = fill.title,
                                       direction = viridis_direction,
                                       na.value = NA)
  } else {
    p <- p +
         ggplot2::scale_color_viridis_c(option = viridis_color_map,
                                        name = fill.title,
                                        direction = viridis_direction,
                                        na.value = NA)
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
                      axis.text.x = ggplot2::element_text(color = "black",
                                                          face = "bold",
                                                          angle = get_axis_parameters(angle = rotate_x_axis_labels, flip = flip)[["angle"]],
                                                          hjust = get_axis_parameters(angle = rotate_x_axis_labels, flip = flip)[["hjust"]],
                                                          vjust = get_axis_parameters(angle = rotate_x_axis_labels, flip = flip)[["vjust"]]),
                      strip.text.x = if (isFALSE(flip)) {ggplot2::element_text(face = "bold",
                                                                               angle = strip_text_angle)}
                                     else {ggplot2::element_blank()},
                      strip.text.y = if (isFALSE(flip)) {ggplot2::element_blank()}
                                     else {ggplot2::element_text(face = "bold",
                                                                 angle = strip_text_angle)},
                      panel.border = ggplot2::element_rect(color = "black", fill = NA),
                      panel.grid = if (isTRUE(plot.grid)){ggplot2::element_line(color = grid.color, linetype = grid.type)} else {ggplot2::element_blank()},
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

  if (isTRUE(compute_ChordDiagrams)){
    data <- output_copy %>%
            dplyr::select(dplyr::all_of(c("source", "target"))) %>%
            dplyr::group_by(.data$target, .data$source) %>%
            dplyr::summarise(value = dplyr::n()) %>%
            dplyr::rename("from" = "source",
                          "to" = "target") %>%
            dplyr::select(dplyr::all_of(c("from", "to", "value")))
    p.source_target <- SCpubr::do_ChordDiagramPlot(from_df = TRUE, df = data, link.border.color = "black", z_index = TRUE)

    data <- liana_output %>%
            dplyr::filter(!(is.na(.data$magnitude))) %>%
            dplyr::select(dplyr::all_of(c("ligand.complex", "receptor.complex"))) %>%
            dplyr::group_by(.data$ligand.complex, .data$receptor.complex) %>%
            dplyr::summarise(value = dplyr::n()) %>%
            dplyr::rename("from" = "ligand.complex",
                          "to" = "receptor.complex") %>%
            dplyr::select(dplyr::all_of(c("from", "to", "value")))
    p.ligand_receptor <- SCpubr::do_ChordDiagramPlot(from_df = TRUE, df = data, link.border.color = "black", z_index = TRUE)
    return(list("dotplot" = p,
                "chord_total_interactions" = p.source_target,
                "chord_ligand_receptor" = p.ligand_receptor))
  } else {
    return(p)
  }

}


