#' Compute an overview of the GO terms associated with the input list of genes.
#'
#' @inheritParams doc_function
#' @param levels.use \strong{\code{\link[base]{numeric}}} | Vector of numerics corresponding to the GO ontology levels to plot. If NULL will compute all recursively until there are no results.
#' @param reverse.levels \strong{\code{\link[base]{logical}}} | Whether to place the higher levels first when computing the joint heatmap.
#' @param return_matrices \strong{\code{\link[base]{logical}}} | Returns the matrices of grouped GO terms.
#' @param colors.use \strong{\code{\link[base]{character}}} | Named vector with two colors assigned to the names \strong{\code{Present}} and \strong{\code{Absent}}.
#' 
#' @return A list containing all the matrices for the respective GO levels and all the individual and combined heatmaps.
#' @export
#'
#' @example /man/examples/examples_do_GroupedGOTermPlot.R
do_GroupedGOTermPlot <- function(genes,
                                 org.db,
                                 levels.use = NULL,
                                 GO_ontology = "BP",
                                 min.overlap = 3,
                                 flip = TRUE,
                                 colors.use = c("Present" = "#1e3d59", 
                                                 "Absent" = "#bccbcd"),
                                 legend.position = "bottom",
                                 reverse.levels = TRUE,
                                 axis.text.x.angle = 45,
                                 font.size = 10,
                                 font.type = "sans",
                                 plot.title = paste0("GO | ", GO_ontology),
                                 plot.subtitle = NULL,
                                 plot.caption = NULL,
                                 verbose = FALSE,
                                 return_matrices = FALSE,
                                 grid.color = "white",
                                 border.color = "black",
                                 plot.title.face = "bold",
                                 plot.subtitle.face = "plain",
                                 plot.caption.face = "italic",
                                 axis.title.face = "bold",
                                 axis.text.face = "plain",
                                 legend.title.face = "bold",
                                 legend.text.face = "plain"){
  # Add lengthy error messages.
  withr::local_options(.new = list("warning.length" = 8170))
  
  
  
  check_suggests(function_name = "do_GroupedGOTermPlot")
  `%>%` <- magrittr::`%>%`
  # Check logical parameters.
  logical_list <- list("flip" = flip,
                       "return_matrices" = return_matrices)
  check_type(parameters = logical_list, required_type = "logical", test_function = is.logical)
  # Check numeric parameters.
  numeric_list <- list("font.size" = font.size,
                       "axis.text.x.angle" = axis.text.x.angle,
                       "min.overlap" = min.overlap,
                       "levels.use" = levels.use)
  check_type(parameters = numeric_list, required_type = "numeric", test_function = is.numeric)
  # Check character parameters.
  character_list <- list("legend.position" = legend.position,
                         "GO_ontology" = GO_ontology,
                         "colors.use" = colors.use,
                         "genes" = genes,
                         "plot.title" = plot.title,
                         "plot.subtitle" = plot.subtitle,
                         "plot.caption" = plot.caption,
                         "grid.color" = grid.color,
                         "border.color" = border.color,
                         "plot.title.face" = plot.title.face,
                         "plot.subtitle.face" = plot.subtitle.face,
                         "plot.caption.face" = plot.caption.face,
                         "axis.title.face" = axis.title.face,
                         "axis.text.face" = axis.text.face,
                         "legend.title.face" = legend.title.face,
                         "legend.text.face" = legend.text.face)
  check_type(parameters = character_list, required_type = "character", test_function = is.character)
  
  assertthat::assert_that(min.overlap >= 1,
                          msg = paste0(add_cross(), crayon_body("Please provide a "),
                                       crayon_key("value"),
                                       crayon_body(" higher or equal to "),
                                       crayon_key("1"),
                                       crayon_body(" to "),
                                       crayon_key("min.overlap"),
                                       crayon_body(".")))
  
  assertthat::assert_that("OrgDb" %in% class(org.db),
                          msg = paste0(add_cross(), crayon_body("Please provide a valid"),
                                       crayon_key("OrgDb"),
                                       crayon_body(" object to "),
                                       crayon_key("org.db"),
                                       crayon_body(".")))

  
  check_colors(grid.color, parameter_name = "grid.color")
  check_colors(border.color, parameter_name = "border.color")
  
  check_parameters(parameter = legend.position, parameter_name = "legend.position")
  check_parameters(parameter = axis.text.x.angle, parameter_name = "axis.text.x.angle")
  check_parameters(parameter = GO_ontology, parameter_name = "GO_ontology")
  check_parameters(plot.title.face, parameter_name = "plot.title.face")
  check_parameters(plot.subtitle.face, parameter_name = "plot.subtitle.face")
  check_parameters(plot.caption.face, parameter_name = "plot.caption.face")
  check_parameters(axis.title.face, parameter_name = "axis.title.face")
  check_parameters(axis.text.face, parameter_name = "axis.text.face")
  check_parameters(legend.title.face, parameter_name = "legend.title.face")
  check_parameters(legend.text.face, parameter_name = "legend.text.face")
  
  
  # Compute the GO terms associated with each gene.
  matrices <- do_GroupedGO_matrices(genes = genes,
                                    org.db = org.db,
                                    levels.use = levels.use,
                                    ontologies = GO_ontology,
                                    min.overlap = min.overlap,
                                    verbose = verbose)
  
  matrix.out <- data.frame()
  for (name in names(matrices[[GO_ontology]])){
    matrix.use <- matrices[[GO_ontology]][[name]]
    matrix.use <- matrix.use %>% dplyr::mutate('Level' = name)
    matrix.out <- rbind(matrix.out, matrix.use)
  }
  

  list.heatmaps <- list()
  counter <- 0
  for (level in rev(unique(matrix.out$Level))){
    counter <- counter + 1
    matrix.use <- matrix.out %>% dplyr::filter(.data$Level == level)
    
    df <- data.frame(row.names = genes)
    df.order <- data.frame(row.names = genes)
    for (term in matrix.use$Description){
      genes.use <- stringr::str_split(matrix.use %>%
                                        dplyr::filter(.data$Description == term) %>%
                                        dplyr::pull(.data$geneID), pattern = "/")[[1]]
      df[[term]] <- ifelse(rownames(df) %in% genes.use, "Present", "Absent")
      df.order[[term]] <- ifelse(rownames(df.order) %in% genes.use, 1, 0)
    }
    
    # Clustering.
    if (counter == 1){
      if(length(rownames(df.order)) == 1){
        row_order <- rownames(df.order)[1]
      } else {
        row_order <- rownames(df.order)[stats::hclust(stats::dist(df.order, method = "euclidean"), method = "ward.D")$order]
      }
    }
    
    if (length(colnames(df.order)) == 1){
      col_order <- colnames(df.order)[1]
    } else {
      col_order <- colnames(df.order)[stats::hclust(stats::dist(t(df.order), method = "euclidean"), method = "ward.D")$order]
    }
    
    p <- df %>% 
         tibble::rownames_to_column(var = "Gene") %>% 
         tidyr::pivot_longer(cols = -dplyr::all_of("Gene"),
                             names_to = "Description",
                             values_to = "Status") %>% 
         dplyr::mutate("Status" = factor(.data$Status, levels = c("Present", "Absent")),
                       "Gene" = factor(.data$Gene, levels = row_order),
                       "Description" = factor(.data$Description, levels = col_order)) %>% 
         # nocov start
         ggplot2::ggplot(mapping = ggplot2::aes(x = if (isTRUE(flip)){.data$Description} else {.data$Gene},
                                                y = if (isTRUE(flip)){.data$Gene} else {.data$Description},
                                                fill = .data$Status)) +
         # nocov end
         ggplot2::geom_tile(color = grid.color, linewidth = 0.5, na.rm = TRUE) +
         ggplot2::scale_y_discrete(expand = c(0, 0)) +
         ggplot2::scale_x_discrete(expand = c(0, 0),
                                   position = "top") + 
         ggplot2::coord_equal() + 
         ggplot2::scale_fill_manual(values = colors.use)
    if (isTRUE(flip)){
      p <- p + 
        ggplot2::guides(y.sec = guide_axis_label_trans(~paste0(levels(.data$Description))),
                        x.sec = guide_axis_label_trans(~paste0(levels(.data$Gene))))
    } else {
      p <- p + 
        ggplot2::guides(y.sec = guide_axis_label_trans(~paste0(levels(.data$Gene))),
                        x.sec = guide_axis_label_trans(~paste0(levels(.data$Description))))
    }
    
    # Set axis titles.
    if (isTRUE(flip)){
      if (counter == 1){
        xlab <- level
        ylab <- "Gene"
      } else {
        xlab <- level
        ylab <- NULL
      }
    } else {
      if (counter == length(unique(matrix.out$Level))){
        xlab <- "Gene"
        ylab <- level
      } else {
        xlab <- NULL
        ylab <- level
      } 
    }
    
    
    axis.parameters <- handle_axis(flip = flip,
                                   group.by = unique(matrix.out$Level),
                                   group = level,
                                   counter = counter,
                                   axis.text.x.angle = axis.text.x.angle,
                                   plot.title.face = plot.title.face,
                                   plot.subtitle.face = plot.subtitle.face,
                                   plot.caption.face = plot.caption.face,
                                   axis.title.face = axis.title.face,
                                   axis.text.face = axis.text.face,
                                   legend.title.face = legend.title.face,
                                   legend.text.face = legend.text.face)
    
    if (counter == 1){
      legend.position.use <- legend.position
    } else {
      legend.position.use <- "none"
    }
    # Set theme
    p <- p +
         ggplot2::xlab(xlab) +
         ggplot2::ylab(ylab) +
         ggplot2::theme_minimal(base_size = font.size) +
         ggplot2::guides(fill = ggplot2::guide_legend(title.position = "top",
                                                      title.hjust = 0.5)) + 
         ggplot2::theme(axis.ticks.x.bottom = axis.parameters$axis.ticks.x.bottom,
                        axis.ticks.x.top = axis.parameters$axis.ticks.x.top,
                        axis.ticks.y.left = axis.parameters$axis.ticks.y.left,
                        axis.ticks.y.right = axis.parameters$axis.ticks.y.right,
                        axis.text.y.left = axis.parameters$axis.text.y.left,
                        axis.text.y.right = axis.parameters$axis.text.y.right,
                        axis.text.x.top = axis.parameters$axis.text.x.top,
                        axis.text.x.bottom = axis.parameters$axis.text.x.bottom,
                        axis.title.x.bottom = axis.parameters$axis.title.x.bottom,
                        axis.title.x.top = axis.parameters$axis.title.x.top,
                        axis.title.y.right = axis.parameters$axis.title.y.right,
                        axis.title.y.left = axis.parameters$axis.title.y.left,
                        strip.background = axis.parameters$strip.background,
                        strip.clip = axis.parameters$strip.clip,
                        strip.text = axis.parameters$strip.text,
                        legend.position = legend.position.use,
                        axis.line = ggplot2::element_blank(),
                        plot.title = ggplot2::element_text(face = plot.title.face, hjust = 0),
                        plot.subtitle = ggplot2::element_text(face = plot.subtitle.face, hjust = 0),
                        plot.caption = ggplot2::element_text(face = plot.caption.face, hjust = 1),
                        legend.text = ggplot2::element_text(face = legend.text.face),
                        legend.title = ggplot2::element_text(face = legend.title.face),
                        plot.title.position = "plot",
                        panel.grid = ggplot2::element_blank(),
                        panel.grid.minor.y = ggplot2::element_line(color = "white", linewidth = 1),
                        text = ggplot2::element_text(family = font.type),
                        plot.caption.position = "plot",
                        legend.justification = "center",
                        plot.margin = ggplot2::margin(t = 5, 
                                                      r = 0, 
                                                      b = 0, 
                                                      l = 5),
                        panel.border = ggplot2::element_rect(fill = NA, color = border.color, linewidth = 1),
                        panel.grid.major = ggplot2::element_blank(),
                        plot.background = ggplot2::element_rect(fill = "white", color = "white"),
                        panel.background = ggplot2::element_rect(fill = "white", color = "white"),
                        legend.background = ggplot2::element_rect(fill = "white", color = "white"))
    
    list.heatmaps[[level]] <- p
  }
  # Revert the order.
  list.heatmaps <- list.heatmaps[rev(names(list.heatmaps))]
  
  
  # Plot the combined plot
  if (isTRUE(flip)){
    names.use <- rev(unique(matrix.out$Level))
  } else {
    names.use <- unique(matrix.out$Level)
  }
  
  input <- list.heatmaps[names.use]
  
  p <- patchwork::wrap_plots(input,
                             ncol = if (base::isFALSE(flip)){1} else {NULL},
                             nrow = if(isTRUE(flip)) {1} else {NULL},
                             guides = "collect")
  p <- p +
       patchwork::plot_annotation(title = plot.title,
                                  subtitle = plot.subtitle,
                                  caption = plot.caption,
                                  theme = ggplot2::theme(legend.position = legend.position,
                                                         plot.title = ggplot2::element_text(family = font.type,
                                                                                            color = "black",
                                                                                            face = plot.title.face,
                                                                                            hjust = 0),
                                                         plot.subtitle = ggplot2::element_text(family = font.type,
                                                                                               face = plot.subtitle.face,
                                                                                               color = "black",
                                                                                               hjust = 0),
                                                         plot.caption = ggplot2::element_text(family = font.type,
                                                                                              face = plot.caption.face,
                                                                                              color = "black",
                                                                                              hjust = 1),
                                                         plot.caption.position = "plot",
                                                         plot.margin = if (isTRUE(flip)){ggplot2::margin(t = 0, 
                                                                                       b = 0,
                                                                                       r = 1,
                                                                                       l = 2,
                                                                                       unit = "cm")} else {NULL}))
  
  return_me <- list()
  
  if (isTRUE(return_matrices)){
    return_me[["Matrices"]] <- matrix.out
    return_me[["Plot"]] <- return_me
  } else {
    return_me <- p
  }
  
  return(return_me)
}




