#' Visualize Ligand-Receptor analysis output.
#'
#' This function makes use of [liana](https://github.com/saezlab/liana) package to run Ligand-Receptor analysis. Takes the output of liana and generates a dot-plot visualization according to the user's specifications.
#'
#' @inheritParams doc_function
#' @param liana_output \strong{\code{\link[tibble]{tibble}}} | Object resulting from running \link[liana]{liana_wrap} and \link[liana]{liana_aggregate}.
#' @param split.by \strong{\code{\link[base]{character}}} | Whether to further facet the plot on the y axis by common ligand.complex or receptor.complex. Values to provide: NULL, ligand.complex, receptor.complex.
#' @param keep_source,keep_target \strong{\code{\link[base]{character}}} | Identities to keep for the source/target of the interactions. NULL otherwise.
#' @param top_interactions \strong{\code{\link[base]{numeric}}} | Number of unique interactions to retrieve ordered by magnitude and specificity. It does not necessarily mean that the output will contain as many, but rather an approximate value.
#' @param dot_border \strong{\code{\link[base]{logical}}} | Whether to draw a black border in the dots.
#' @param dot.size \strong{\code{\link[base]{numeric}}} | Size aesthetic for the dots.
#' @param sort.by \strong{\code{\link[base]{character}}} | How to arrange the top interactions. Interactions are sorted and then the top N are retrieved and displayed. This takes place after subsetting for \strong{\code{keep_source}} and \strong{\code{keep_target}}  One of:
#' \itemize{
#'   \item \emph{\code{A}}: Sorts by specificity.
#'   \item \emph{\code{B}}: Sorts by magnitude.
#'   \item \emph{\code{C}}: Sorts by specificity, then magnitude (gives extra weight to specificity).
#'   \item \emph{\code{D}}: Sorts by magnitude, then specificity (gives extra weight to magnitude). Might lead to the display of non-significant results.
#'   \item \emph{\code{E}}: Sorts by specificity and magnitude equally.
#' }
#' @param specificity,magnitude \strong{\code{\link[base]{character}}} | Which columns to use for \strong{\code{specificity}} and \strong{\code{magnitude}}.
#' @param invert_specificity,invert_magnitude \strong{\code{\link[base]{logical}}} | Whether to \strong{\code{-log10}} transform \strong{\code{specificity}} and \strong{\code{magnitude}} columns.
#' @param sorting.type.specificity,sorting.type.magnitude \strong{\code{\link[base]{character}}} | Whether the sorting of e \strong{\code{magnitude}} or \strong{\code{specificity}} columns is done in ascending or descending order. This synergises with the value of e \strong{\code{invert_specificity}} and e \strong{\code{invert_magnitude}} parameters.
#' @param compute_ChordDiagrams \strong{\code{\link[base]{logical}}} | Whether to also compute Chord Diagrams for both the number of interactions between source and target but also between ligand.complex and receptor.complex.
#' @param sort_interactions_alphabetically \strong{\code{\link[base]{logical}}} | Sort the interactions to be plotted alphabetically (\strong{\code{TRUE}}) or keep them in their original order in the matrix (\strong{\code{FALSE}}).
#' @param return_interactions \strong{\code{\link[base]{logical}}} | Whether to return the data.frames with the interactions so that they can be plotted as chord plots using other package functions.
#' 
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
                                  magnitude = "sca.LRscore",
                                  specificity = "aggregate_rank",
                                  sort.by = "E",
                                  sorting.type.specificity = "descending",
                                  sorting.type.magnitude = "descending",
                                  border.color = "black",
                                  axis.text.x.angle = 45,
                                  legend.position = "bottom",
                                  legend.type = "colorbar",
                                  legend.length = 20,
                                  legend.width = 1,
                                  legend.framecolor = "grey50",
                                  legend.tickcolor = "white",
                                  legend.framewidth = 0.5,
                                  legend.tickwidth = 0.5,
                                  use_viridis = FALSE,
                                  viridis.palette = "G",
                                  viridis.direction = 1,
                                  sequential.palette = "YlGnBu",
                                  sequential.direction = 1,
                                  font.size = 14,
                                  dot.size = 1,
                                  font.type = "sans",
                                  plot.grid = TRUE,
                                  grid.color = "grey90",
                                  grid.type = "dotted",
                                  compute_ChordDiagrams = FALSE,
                                  sort_interactions_alphabetically = FALSE,
                                  number.breaks = 5,
                                  plot.title.face = "bold",
                                  plot.subtitle.face = "plain",
                                  plot.caption.face = "italic",
                                  axis.title.face = "bold",
                                  axis.text.face = "bold",
                                  legend.title.face = "bold",
                                  legend.text.face = "plain",
                                  return_interactions = FALSE,
                                  invert_specificity = TRUE,
                                  invert_magnitude = FALSE,
                                  verbose = TRUE){
  # Add lengthy error messages.
  withr::local_options(.new = list("warning.length" = 8170))
  
  # Checks for packages.
  check_suggests(function_name = "do_LigandReceptorPlot")
  `%>%` <- magrittr::`%>%`
  `:=` <- rlang::`:=`

  # Check logical parameters.
  logical_list <- list("dot_border" = dot_border,
                       "plot.grid" = plot.grid,
                       "sort_interactions_alphabetically" = sort_interactions_alphabetically,
                       "use_viridis" = use_viridis,
                       "return_interactions" = return_interactions,
                       "invert_specificity" = invert_specificity,
                       "invert_magnitude" = invert_magnitude,
                       "verbose" = verbose)
  check_type(parameters = logical_list, required_type = "logical", test_function = is.logical)
  # Check numeric parameters.
  numeric_list <- list("font.size" = font.size,
                       "top_interactions" = top_interactions,
                       "legend.length" = legend.length,
                       "legend.width" = legend.width,
                       "legend.framewidth" = legend.framewidth,
                       "legend.tickwidth" = legend.tickwidth,
                       "dot.size" = dot.size,
                       "axis.text.x.angle" = axis.text.x.angle,
                       "viridis.direction" = viridis.direction,
                       "number.breaks" = number.breaks,
                       "sequential.direction" = sequential.direction)
  check_type(parameters = numeric_list, required_type = "numeric", test_function = is.numeric)
  # Check character parameters.
  character_list <- list("split.by" = split.by,
                         "keep_source" = keep_source,
                         "keep_target" = keep_target,
                         "border.color" = border.color,
                         "legend.position" = legend.position,
                         "legend.type" = legend.type,
                         "legend.framecolor" = legend.framecolor,
                         "viridis.palette" = viridis.palette,
                         "legend.tickcolor" = legend.tickcolor,
                         "font.type" = font.type,
                         "grid.color" = grid.color,
                         "grid.type" = grid.type,
                         "sequential.palette" = sequential.palette,
                         "plot.title.face" = plot.title.face,
                         "plot.subtitle.face" = plot.subtitle.face,
                         "plot.caption.face" = plot.caption.face,
                         "axis.title.face" = axis.title.face,
                         "axis.text.face" = axis.text.face,
                         "legend.title.face" = legend.title.face,
                         "legend.text.face" = legend.text.face,
                         "sort.by" = sort.by,
                         "sorting.type.specificity" = sorting.type.specificity,
                         "sorting.type.magnitude" = sorting.type.magnitude)
  check_type(parameters = character_list, required_type = "character", test_function = is.character)

  # Check border color.
  check_colors(border.color, parameter_name = "border.color")

  # Check the colors provided to legend.framecolor and legend.tickcolor.
  check_colors(legend.framecolor, parameter_name = "legend.framecolor")
  check_colors(legend.tickcolor, parameter_name = "legend.tickcolor")
  check_colors(grid.color, parameter_name = "grid.color")

  check_parameters(parameter = font.type, parameter_name = "font.type")
  check_parameters(parameter = legend.type, parameter_name = "legend.type")
  check_parameters(parameter = legend.position, parameter_name = "legend.position")
  check_parameters(parameter = viridis.palette, parameter_name = "viridis.palette")
  check_parameters(parameter = grid.type, parameter_name = "grid.type")
  check_parameters(parameter = axis.text.x.angle, parameter_name = "axis.text.x.angle")
  check_parameters(parameter = number.breaks, parameter_name = "number.breaks")
  check_parameters(plot.title.face, parameter_name = "plot.title.face")
  check_parameters(plot.subtitle.face, parameter_name = "plot.subtitle.face")
  check_parameters(plot.caption.face, parameter_name = "plot.caption.face")
  check_parameters(axis.title.face, parameter_name = "axis.title.face")
  check_parameters(axis.text.face, parameter_name = "axis.text.face")
  check_parameters(legend.title.face, parameter_name = "legend.title.face")
  check_parameters(legend.text.face, parameter_name = "legend.text.face")
  check_parameters(viridis.direction, parameter_name = "viridis.direction")
  check_parameters(sequential.direction, parameter_name = "sequential.direction")
  
  
  colors.gradient <- compute_continuous_palette(name = ifelse(isTRUE(use_viridis), viridis.palette, sequential.palette),
                                                use_viridis = use_viridis,
                                                direction = ifelse(isTRUE(use_viridis), viridis.direction, sequential.direction),
                                                enforce_symmetry = FALSE)
  
  if (!is.null(split.by)){
    assertthat::assert_that(split.by %in% c("receptor.complex", "ligand.complex"),
                            msg = paste0(add_cross,
                                         crayon_body("Please select one of the following for "),
                                         crayon_key("split.by"),
                                         crayon_body(": "),
                                         crayon_key("ligand.complex"),
                                         crayon_body(", "),
                                         crayon_key("receptor.complex"),
                                         crayon_body(".")))
  }

  # Define legend parameters. Width and height values will change depending on the legend orientation.
  if (legend.position %in% c("top", "bottom")){
    size_title <- "Interaction specificity"
    fill.title <- "Expression Magnitude"
  } else if (legend.position %in% c("left", "right")){
    size_title <- stringr::str_wrap("Interaction specificity", width = 10)
    fill.title <- stringr::str_wrap("Expression Magnitude", width = 10)
  }
  
  if (isTRUE(verbose)){
    rlang::inform(paste0(add_info(initial_newline = FALSE),
                         crayon_body("Column for specificity: "),
                         crayon_key(specificity)))
    
    rlang::inform(paste0(add_info(initial_newline = FALSE),
                         crayon_body("Column for magnitude: "),
                         crayon_key(magnitude)))
  }
  
  liana_output <- liana_output %>%
                  dplyr::mutate("magnitude" = .data[[magnitude]]) %>%
                  dplyr::mutate("specificity" = .data[[specificity]])
  
  invert_function <- function(x){-log10(x + 1e-10)}
    
  if (isTRUE(invert_specificity)){
    liana_output <- liana_output %>% 
                    dplyr::mutate("specificity" := invert_function(x = .data$specificity))
  }
    
  if (isTRUE(invert_magnitude)){
    liana_output <- liana_output %>% 
                    dplyr::mutate("magnitude" := invert_function(.data$magnitude))
  }
 
  # Differential arrangement of the interactions.
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
  
  # Sort interactions according to user's preference.
  if (sort.by == "A"){
    if (isTRUE(verbose)){
      rlang::inform(paste0(add_info(initial_newline = FALSE),
                           crayon_body("Sorting interactions by: "),
                           crayon_key("specificify")))
    }
    
    if (sorting.type.specificity == "descending"){
      liana_output <- liana_output %>% 
                      dplyr::arrange(dplyr::desc(.data$specificity))
    } else {
      liana_output <- liana_output %>% 
                      dplyr::arrange(.data$specificity)
    }
    
  } else if (sort.by == "B"){
    if (isTRUE(verbose)){
      rlang::inform(paste0(add_info(initial_newline = FALSE),
                           crayon_body("Sorting interactions by: "),
                           crayon_key("magnitude")))
    }
    
    if (sorting.type.magnitude == "descending"){
      liana_output <- liana_output %>% 
                      dplyr::arrange(dplyr::desc(.data$magnitude))
    } else {
      liana_output <- liana_output %>% 
                      dplyr::arrange(.data$magnitude)
    }
  } else if (sort.by == "C"){
    if (isTRUE(verbose)){
      rlang::inform(paste0(add_info(initial_newline = FALSE),
                           crayon_body("Sorting interactions by: "),
                           crayon_key("specificify"),
                           crayon_body(" then "),
                           crayon_key("magnitude"),
                           crayon_body(".")))
    }
    
    if (sorting.type.magnitude == "ascending" & sorting.type.specificity == "ascending"){
      liana_output <- liana_output %>% 
                      dplyr::arrange(.data$specificity, .data$magnitude)
    } else if (sorting.type.magnitude == "descending" & sorting.type.specificity == "ascending"){
      liana_output <- liana_output %>% 
                      dplyr::arrange(.data$specificity, dplyr::desc(.data$magnitude))
    } else if (sorting.type.magnitude == "ascending" & sorting.type.specificity == "descending"){
      liana_output <- liana_output %>% 
                      dplyr::arrange(dplyr::desc(.data$specificity), .data$magnitude)
    } else if (sorting.type.magnitude == "descending" & sorting.type.specificity == "descending"){
      liana_output <- liana_output %>% 
                      dplyr::arrange(dplyr::desc(.data$specificity), dplyr::desc(.data$magnitude))
    }
    
  } else if (sort.by == "D"){
    if (isTRUE(verbose)){
      rlang::inform(paste0(add_info(initial_newline = FALSE),
                           crayon_body("Sorting interactions by: "),
                           crayon_key("magnitude"),
                           crayon_body(" then "),
                           crayon_key("specificity"),
                           crayon_body(".")))
    }
    
    if (sorting.type.magnitude == "ascending" & sorting.type.specificity == "ascending"){
      liana_output <- liana_output %>% 
                      dplyr::arrange(.data$magnitude, .data$specificity)
    } else if (sorting.type.magnitude == "descending" & sorting.type.specificity == "ascending"){
      liana_output <- liana_output %>% 
                      dplyr::arrange(.data$magnitude, dplyr::desc(.data$specificity))
    } else if (sorting.type.magnitude == "ascending" & sorting.type.specificity == "descending"){
      liana_output <- liana_output %>% 
                      dplyr::arrange(dplyr::desc(.data$magnitude), .data$specificity)
    } else if (sorting.type.magnitude == "descending" & sorting.type.specificity == "descending"){
      liana_output <- liana_output %>% 
                      dplyr::arrange(dplyr::desc(.data$magnitude), dplyr::desc(.data$specificity))
    }
  } else if (sort.by == "E"){
    if (isTRUE(verbose)){
      rlang::inform(paste0(add_info(initial_newline = FALSE),
                           crayon_body("Sorting interactions by: "),
                           crayon_key("magnitude"),
                           crayon_body(" and "),
                           crayon_key("specificity"),
                           crayon_body(" with equal weights.")))
    }
    
    if (sorting.type.magnitude == "ascending"){
      liana_output_magnitude <- liana_output %>%
                                dplyr::arrange(.data$magnitude) %>% 
                                tibble::rowid_to_column(var = "magnitude_rank")
    } else {
      liana_output_magnitude <- liana_output %>%
                                dplyr::arrange(dplyr::desc(.data$magnitude)) %>% 
                                tibble::rowid_to_column(var = "magnitude_rank")
    }
    
    if (sorting.type.specificity == "ascending"){
      liana_output_specificity <- liana_output %>%
                                  dplyr::arrange(.data$specificity) %>% 
                                  tibble::rowid_to_column(var = "specificity_rank")
    } else {
      liana_output_specificity <- liana_output %>%
                                  dplyr::arrange(dplyr::desc(.data$specificity)) %>% 
                                  tibble::rowid_to_column(var = "specificity_rank")
    }
    
    liana_output <- liana_output %>% 
                    dplyr::left_join(y = liana_output_specificity %>% dplyr::select(dplyr::all_of(c("interaction", "specificity_rank"))),
                                     by = "interaction",
                                     relationship = "many-to-many") %>% 
                    dplyr::left_join(y = liana_output_magnitude %>% dplyr::select(dplyr::all_of(c("interaction", "magnitude_rank"))),
                                     by = "interaction",
                                     relationship = "many-to-many") %>% 
                    dplyr::mutate("rank" = .data$magnitude_rank + .data$specificity_rank) %>% 
                    dplyr::arrange(.data$rank) %>% 
                    dplyr::select(!dplyr::all_of(c("rank", "magnitude_rank", "specificity_rank")))
    rm(liana_output_magnitude)
    rm(liana_output_specificity)
  }
  
  if (isTRUE(verbose)){
    rlang::inform(paste0(add_info(initial_newline = FALSE),
                         crayon_body("Sorting type specificity: "),
                         crayon_key(sorting.type.specificity)))
    
    rlang::inform(paste0(add_info(initial_newline = FALSE),
                         crayon_body("Sorting type magnitude: "),
                         crayon_key(sorting.type.magnitude)))
    
    rlang::inform(paste0(add_info(initial_newline = FALSE),
                         crayon_body("Plotting the following top interanctions: "),
                         crayon_key(top_interactions)))
  }
  
  liana_output <- liana_output %>%
                  # Filter based on the top X interactions of ascending sensibilities.
                  dplyr::inner_join(y = {liana_output %>%
                                         dplyr::distinct_at(c("ligand.complex", "receptor.complex")) %>%
                                         dplyr::slice_head(n = top_interactions)},
                                    by = c("ligand.complex", "receptor.complex"),
                                    relationship = "many-to-many")
  
  assertthat::assert_that(nrow(liana_output) > 0,
                          msg = paste0(add_cross(), crayon_body("Whith the current presets of "),
                                       crayon_key("keep_source"),
                                       crayon_body(" and "),
                                       crayon_key("keep_target"),
                                       crayon_body(" there are no interactions left.")))

  # Make source and target factors, so that they do not get dropped by the plot.
  if (isTRUE(sort_interactions_alphabetically)){
    liana_output$source <- factor(liana_output$source, levels = sort(unique(liana_output$source)))
    liana_output$target <- factor(liana_output$target, levels = sort(unique(liana_output$target)))
    liana_output$interaction <- factor(liana_output$interaction, levels = rev(sort(unique(liana_output$interaction))))
  } else if (base::isFALSE(sort_interactions_alphabetically)){
    liana_output$source <- factor(liana_output$source, levels = sort(unique(liana_output$source)))
    liana_output$target <- factor(liana_output$target, levels = sort(unique(liana_output$target)))
    liana_output$interaction <- factor(liana_output$interaction, levels = rev(unique(liana_output$interaction)))
  }


  # Plot.
  if (isTRUE(dot_border)){
    p <-  liana_output %>%
          ggplot2::ggplot(mapping = ggplot2::aes(x = .data$target,
                                                 y = .data$interaction,
                                                 fill = .data$magnitude,
                                                 size = .data$specificity,
                                                 group = .data$interacting_clusters)) +
          ggplot2::geom_point(shape = 21,
                              na.rm = TRUE)
  } else if (base::isFALSE(dot_border)){
    p <-  liana_output %>%
          ggplot2::ggplot(mapping = ggplot2::aes(x = .data$target,
                                                 y = .data$interaction,
                                                 size = .data$specificity,
                                                 group = .data$interacting_clusters)) +
          ggplot2::geom_point(mapping = ggplot2::aes(color = .data$magnitude),
                              shape = 19,
                              na.rm = TRUE)
  }
  
  p <- p +
       ggplot2::guides(y.sec = guide_axis_label_trans(~paste0(levels(.data$interaction)))) +
       ggplot2::scale_size_continuous(name = size_title,
                                      range = c(2 * dot.size, 10 * dot.size)) 
  
  # Settings for bordered dots.
  limits <- c(min(liana_output$magnitude, na.rm = TRUE),
              max(liana_output$magnitude, na.rm = TRUE))
  
  scale.setup <- compute_scales(sample = NULL,
                                feature = NULL,
                                assay = NULL,
                                reduction = NULL,
                                slot = NULL,
                                number.breaks = number.breaks,
                                min.cutoff = NA,
                                max.cutoff = NA,
                                flavor = "Seurat",
                                enforce_symmetry = FALSE,
                                from_data = TRUE,
                                limits.use = limits)

  if (isTRUE(dot_border)){
    # Add color to aesthetics.
    p$layers[[1]]$aes_params$color <- border.color
    p <- p + 
         ggplot2::scale_fill_gradientn(colors = colors.gradient,
                                       na.value = NA,
                                       name = fill.title,
                                       breaks = scale.setup$breaks,
                                       labels = scale.setup$labels,
                                       limits = scale.setup$limits)
  } else {
    p <- p + 
         ggplot2::scale_color_gradientn(colors = colors.gradient,
                                       na.value = NA,
                                       name = fill.title,
                                       breaks = scale.setup$breaks,
                                       labels = scale.setup$labels,
                                       limits = scale.setup$limits)
  }
  # Continue plotting.
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




  p <- p +
       ggplot2::labs(title = "Source") +
       ggplot2::xlab("Target") +
       ggplot2::ylab(paste("Ligand", "|", "Receptor", sep = " ")) +
       ggplot2::guides(size = ggplot2::guide_legend(title.position = "top",
                                                    title.hjust = 0.5,
                                                    override.aes = ggplot2::aes(fill = "black"))) +
       ggplot2::theme_minimal(base_size = font.size) +
       ggplot2::theme(plot.title = ggplot2::element_text(face = plot.title.face,
                                                         hjust = 0.5,
                                                         vjust = 0,
                                                         size = font.size),
                      plot.subtitle = ggplot2::element_text(face = plot.subtitle.face, hjust = 0),
                      plot.caption = ggplot2::element_text(face = plot.caption.face, hjust = 1),
                      legend.text = ggplot2::element_text(face = legend.text.face),
                      legend.title = ggplot2::element_text(face = legend.title.face),
                      plot.title.position = "panel",
                      plot.caption.position = "plot",
                      text = ggplot2::element_text(family = font.type),
                      legend.justification = "center",
                      legend.position = legend.position,
                      axis.title.x = ggplot2::element_text(color = "black", face = axis.title.face, hjust = 0.5),
                      axis.title.y.left = ggplot2::element_text(color = "black", face = axis.title.face, angle = 90),
                      axis.title.y.right = ggplot2::element_blank(),
                      axis.text.y.right = ggplot2::element_text(color = "black",
                                                                face = axis.text.face),
                      axis.text.y.left = ggplot2::element_blank(),
                      axis.ticks.x = ggplot2::element_line(color = "black"),
                      axis.ticks.y.left = ggplot2::element_blank(),
                      axis.ticks.y.right = ggplot2::element_line(color = "black"),
                      axis.text.x = ggplot2::element_text(color = "black",
                                                          face = axis.text.face,
                                                          angle = get_axis_parameters(angle = axis.text.x.angle, flip = FALSE)[["angle"]],
                                                          hjust = get_axis_parameters(angle = axis.text.x.angle, flip = FALSE)[["hjust"]],
                                                          vjust = get_axis_parameters(angle = axis.text.x.angle, flip = FALSE)[["vjust"]]),
                      strip.text.x = ggplot2::element_text(face = "bold",
                                                           color = "black"),
                      strip.text.y = ggplot2::element_blank(),
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

  if (isTRUE(return_interactions)){
    data_interactions <- output_copy %>%
                         dplyr::select(dplyr::all_of(c("source", "target"))) %>%
                         dplyr::group_by(.data$target, .data$source) %>%
                         dplyr::summarise(value = dplyr::n()) %>%
                         dplyr::rename("from" = "source",
                                       "to" = "target") %>%
                         dplyr::select(dplyr::all_of(c("from", "to", "value")))
                 

    data_LF <- liana_output %>%
               dplyr::filter(!(is.na(.data$magnitude))) %>%
               dplyr::select(dplyr::all_of(c("ligand.complex", "receptor.complex"))) %>%
               dplyr::group_by(.data$ligand.complex, .data$receptor.complex) %>%
               dplyr::summarise(value = dplyr::n()) %>%
               dplyr::rename("from" = "ligand.complex",
                             "to" = "receptor.complex") %>%
               dplyr::select(dplyr::all_of(c("from", "to", "value")))
    
    return(list("Plot" = p,
                "Group Interactions" = data_interactions,
                "LR Interactions" = data_LF))
  } else {
    return(p)
  }
}


