#' Display CNV scores from inferCNV as Feature Plots.
#'
#'
#'
#' @param sample Seurat object.
#' @param infercnv_object Output inferCNV object run on the same Seurat object.
#' @param group.by Metadata variable to group cells by. If not provided, defaults to active identitites.
#' @param using_metacells Logical. Whether inferCNV was run using metacells or not.
#' @param metacell_mapping Vector or cell - metacell mapping.
#' @param chromosome_focus Character. Region stating which chromosome to plot. Eg: 1p, 19q. NULL will plot all regions.
#' @param chromosome_locations Tibble. Tibble containing the chromosome regions to use. Can be obtained using utils::data("human_chr_locations", package = "SCpubr")
#' @param font.type Character. Base font for the plot. One of mono, serif or sans.
#' @param legend.type Character. Type of legend to display. One of: normal, colorbar, colorsteps.
#' @param legend.position Position of the legend in the plot. Will only work if legend is set to TRUE.
#' @param legend.framewidth,legend.tickwidth Width of the lines of the box in the legend.
#' @param legend.framecolor,legend.tickcolor Color of the lines of the box in the legend.
#' @param legend.length,legend.width Length and width of the legend. Will adjust automatically depending on legend side.
#' @param pt.size Size of the dots in the dotplot pltos.
#' @param font.size Overall fontsize for the plots.
#' @param border.size Thickness of the border around cells in all plots.
#' @param rotate_x_axis_labels Logical. Rotate x axis labels 90 degrees in the dotplot plots.
#'
#' @return A list containing Feature Plots for different chromosome regions and corresponding dot plots by groups..
#' @export
#'
#' @example man/examples/examples_do_CopyNumberVariantPlot.R
do_CopyNumberVariantPlot <- function(sample,
                                     infercnv_object,
                                     chromosome_locations,
                                     group.by = NULL,
                                     using_metacells = FALSE,
                                     metacell_mapping = NULL,
                                     chromosome_focus = NULL,
                                     legend.type = "colorbar",
                                     legend.position = "bottom",
                                     legend.length = 30,
                                     legend.width = 1,
                                     legend.framewidth = 1.5,
                                     legend.tickwidth = 1.5,
                                     legend.framecolor = "grey50",
                                     legend.tickcolor = "white",
                                     font.size = 14,
                                     pt.size = 1,
                                     font.type = "sans",
                                     border.size = 1.5,
                                     rotate_x_axis_labels = TRUE){


  check_suggests(function_name = "do_CopyNumberVariantPlot")

  # Check logical parameters.
  logical_list <- list("using_metacells" = using_metacells,
                       "rotate_x_axis_labels" = rotate_x_axis_labels)
  check_type(parameters = logical_list, required_type = "logical", test_function = is.logical)
  # Check numeric parameters.
  numeric_list <- list("font.size" = font.size,
                       "legend.length" = legend.length,
                       "legend.width" = legend.width,
                       "legend.framewidth" = legend.framewidth,
                       "legend.tickwidth" = legend.tickwidth,
                       "pt.size" = pt.size,
                       "border.size" = border.size)
  check_type(parameters = numeric_list, required_type = "numeric", test_function = is.numeric)
  # Check character parameters.
  character_list <- list("group.by" = group.by,
                         "legend.type" = legend.type,
                         "legend.position" = legend.position,
                         "legend.framecolor" = legend.framecolor,
                         "legend.tickcolor" = legend.tickcolor,
                         "font.type" = font.type)
  check_type(parameters = character_list, required_type = "character", test_function = is.character)



  `%>%` <- purrr::`%>%`
  check_colors(legend.framecolor, parameter_name = "legend.framecolor")
  check_colors(legend.tickcolor, parameter_name = "legend.tickcolor")

  if (is.null(chromosome_focus)){
    chromosome_list <- c(as.character(seq(1, 22)))
  } else {
    chromosome_list <- chromosome_focus
  }

  # Fix for group.by.
  if (is.null(group.by)){
    group.by <- "dummy"
    sample@meta.data$dummy <- sample@active.ident
  } else {
    sample@meta.data$dummy <- sample@meta.data[, group.by]
    group.by <- "dummy"
  }

  # Retrieve the genes.
  genes <- infercnv_object@gene_order

  # Retrieve chr 1p start and end coordinates.
  chr_locations <- chromosome_locations

  # This list will contain all the outputs.
  return_list <- list()
  events_list <- c()
  for (chromosome in chromosome_list){
    # Retrieve chr locations of the chromosome.
    locations <- chr_locations %>% dplyr::filter(.data$chr == chromosome)

    for (chr_arm in c("p", "q")){
      # Retrieve the start.
      start <- locations %>% dplyr::filter(.data$arm == chr_arm) %>% dplyr::pull(start)
      # Retrieve the end.
      end <- locations %>% dplyr::filter(.data$arm == chr_arm) %>% dplyr::pull(end)

      # Retrieve the genes present in the chromosome arm.
      genes_use <- rownames(genes %>% dplyr::filter(.data$chr == paste0("chr", chromosome),  stop <= end))

      # Retrieve the CNV scores from the inferCNV object.
      CNV_scores <- infercnv_object@expr.data

      # Make it at least 2 genes in the object (this will only be applicable in VERY LOW QUALITY DATASETS)
      if (sum(genes_use %in% rownames(CNV_scores)) > 1){
        # Filter the scores for only the genes in the chromosome arm.
        CNV_scores <- CNV_scores[genes_use[genes_use %in% rownames(CNV_scores)], ]

        CNV_scores_final <- tibble::tibble("score" = colMeans(CNV_scores),
                                           "cells" = colnames(CNV_scores))

        scores_name <- paste0(chromosome, chr_arm)
        # If metacells were used.
        if (isTRUE(using_metacells)){
          # Magic
          sample@meta.data[, "metacell_mapping"] <- metacell_mapping
          sample@meta.data[, scores_name] <- CNV_scores_final$score[match(sample@meta.data$metacell_mapping,CNV_scores_final$cells)]
          # If no metacells were used.
        } else if (isFALSE(using_metacells)){
          sample@meta.data[, scores_name] <- CNV_scores_final$score
        }
        events_list <- append(events_list, scores_name)
      }
    }
  }

  sample$cells <- colnames(sample)
  sample$group <- sample@meta.data$dummy
  for (event in events_list){
    p <- sample@meta.data %>%
            dplyr::select(event, "group", "cells") %>%
            ggplot2::ggplot(mapping = ggplot2::aes(x = .data$group,
                                                   y = !!(rlang::sym(event)),
                                                   color = !!(rlang::sym(event)))) +
            ggplot2::geom_point(position = ggplot2::position_jitter(width = 0.455,
                                                                   seed = 0),
                              size = pt.size * border.size,
                                      color = "black",
                              na.rm = TRUE) +
            ggplot2::geom_point(position = ggplot2::position_jitter(width = 0.45,
                                                                   seed = 0),size = pt.size,
                                na.rm = TRUE) +
            ggdist::stat_pointinterval(position = ggplot2::position_dodge(width = 1),
                                       na.rm = TRUE) +
            ggplot2::scale_color_gradientn(colors = c("#033270", "#4091C9", "#fdf0d5", "#c94040", "#65010C"),
                                           limits = c(0.8, 1.2),
                                           breaks = seq(0.8, 1.2, by = 0.1),
                                           labels = as.character(seq(0.8, 1.2, by = 0.1))) +
            # Add radial lines that will span further in the Y axis pointing to the term labels.
            ggplot2::geom_hline(mapping = ggplot2::aes(yintercept = 1),
                               color = "black",
                               linetype = "dashed") +
            ggplot2::ylim(c(0.8, 1.2)) +
            ggplot2::theme_minimal(base_size = font.size) +
            ggplot2::theme(axis.title = ggplot2::element_blank(),
                           axis.line.x = ggplot2::element_line(color = "black"),
                           axis.text.x = ggplot2::element_text(color = "black",
                                                               face = "bold",
                                                               angle = ifelse(isTRUE(rotate_x_axis_labels), 90, 0),
                                                               hjust = ifelse(isTRUE(rotate_x_axis_labels), 1, 0.5),
                                                               vjust = ifelse(isTRUE(rotate_x_axis_labels), 0.5, 1)),
                           axis.text.y = ggplot2::element_text(color = "black", face = "bold", hjust = 0),
                           axis.ticks = ggplot2::element_line(color = "black"),
                           panel.grid.major = ggplot2::element_blank(),
                           plot.title.position = "plot",
                           plot.title = ggtext::element_markdown(face = "bold", hjust = 0),
                           plot.subtitle = ggtext::element_markdown(hjust = 0),
                           plot.caption = ggtext::element_markdown(hjust = 1),
                           panel.grid = ggplot2::element_blank(),
                           text = ggplot2::element_text(family = font.type),
                           plot.caption.position = "plot",
                           legend.text = ggplot2::element_text(face = "bold"),
                           legend.position = legend.position,
                           legend.title = ggplot2::element_text(face = "bold"),
                           legend.justification = "center",
                           plot.margin = ggplot2::margin(t = 10, r = 40, b = 10, l = 40),
                           plot.background = ggplot2::element_rect(fill = "white", color = "white"),
                           panel.background = ggplot2::element_rect(fill = "white", color = "white"),
                           legend.background = ggplot2::element_rect(fill = "white", color = "white"))

    # Define legend parameters. Width and height values will change depending on the legend orientation.
    if (legend.position %in% c("top", "bottom")){
      legend.barwidth <- legend.length
      legend.barheight <- legend.width
    } else if (legend.position %in% c("left", "right")){
      legend.barwidth <- legend.width
      legend.barheight <- legend.length
    }

    # Modify the way the legend is displayed.
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

    # Plot the scores!
    p.f <- do_FeaturePlot(sample = sample,
                        features = event,
                        plot_cell_borders = T)
    p.f <- add_scale(p = p.f,
                     scale = "color",
                     function_use = ggplot2::scale_color_gradientn(colors = c("#033270", "#4091C9", "#fdf0d5", "#c94040", "#65010C"),
                                                                   limits = c(0.8, 1.2),
                                                                   breaks = seq(0.8, 1.2, by = 0.1),
                                                                   labels = as.character(seq(0.8, 1.2, by = 0.1))))
    return_list[[paste0(event, "_umap")]] <- p.f
    return_list[[paste0(event, "_dotplot")]] <- p
  }
  return(return_list)
}
