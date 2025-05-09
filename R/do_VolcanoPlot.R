#' Compute a Volcano plot out of DE genes.
#'
#' @inheritParams doc_function
#' @param de_genes \strong{\code{\link[tibble]{tibble}}} | Output of `Seurat::FindMarkers()`.
#' @param pval_cutoff \strong{\code{\link[base]{numeric}}} | Cutoff for the p-value.
#' @param FC_cutoff \strong{\code{\link[base]{numeric}}} | Cutoff for the avg_log2FC.
#' @param plot_lines \strong{\code{\link[base]{logical}}} | Whether to plot the division lines.
#' @param line_color \strong{\code{\link[base]{character}}} | Color for the lines.
#' @param line_size \strong{\code{\link[base]{numeric}}} | Size of the lines in the plot.
#' @param add_gene_tags \strong{\code{\link[base]{logical}}} | Whether to plot the top genes.
#' @param add_tag_side \strong{\code{\link[base]{logical}}} | Either "both", "positive" or "negative" to indicate which side of genes to tag
#' @param order_tags_by \strong{\code{\link[base]{character}}} | Either "both", "pvalue" or "logfc".
#' @param tag_size \strong{\code{\link[base]{numeric}}} | Size of the text/label for the tags.
#' @param n_genes \strong{\code{\link[base]{numeric}}} | Number of top genes to plot.
#' @param use_labels \strong{\code{\link[base]{logical}}} | Whether to use labels instead of text for the tags.
#' @param colors.use \strong{\code{\link[base]{character}}} | Color to generate a tetradic color scale with. If NULL, default colors are used.
#'
#' @return A volcano plot as a ggplot2 object.
#' @export
#'
#' @example /man/examples/examples_do_VolcanoPlot.R
do_VolcanoPlot <- function(sample,
                           de_genes,
                           pval_cutoff = 0.05,
                           FC_cutoff = 2,
                           pt.size = 1,
                           border.size = 1.5,
                           border.color = "black",
                           font.size = 14,
                           font.type = "sans",
                           plot.title = NULL,
                           plot.subtitle = NULL,
                           plot.caption = NULL,
                           plot_lines = TRUE,
                           line_color = "grey75",
                           line_size = 0.5,
                           add_gene_tags = TRUE,
                           add_tag_side = "both",
                           order_tags_by = "both",
                           tag_size = 6,
                           n_genes = 5,
                           use_labels = FALSE,
                           colors.use = NULL,
                           plot.title.face = "bold",
                           plot.subtitle.face = "plain",
                           plot.caption.face = "italic",
                           axis.title.face = "bold",
                           axis.text.face = "plain",
                           legend.title.face = "bold",
                           legend.text.face = "plain"){
  # Add lengthy error messages.
  withr::local_options(.new = list("warning.length" = 8170))

  check_suggests(function_name = "do_VolcanoPlot")
  # Check if the sample provided is a Seurat object.
  check_Seurat(sample = sample)

  # Check logical parameters.
  logical_list <- list("add_gene_tags" = add_gene_tags,
                       "plot_lines" = plot_lines,
                       "use_labels" = use_labels)
  check_type(parameters = logical_list, required_type = "logical", test_function = is.logical)
  # Check numeric parameters.
  numeric_list <- list("pval_cutoff" = pval_cutoff,
                       "FC_cutoff" = FC_cutoff,
                       "pt.size" = pt.size,
                       "border.size" = border.size,
                       "font.size" = font.size,
                       "line_size" = line_size,
                       "n_genes" = n_genes,
                       "tag_size" = tag_size)
  check_type(parameters = numeric_list, required_type = "numeric", test_function = is.numeric)
  # Check character parameters.
  character_list <- list("border.color" = border.color,
                         "font.type" = font.type,
                         "line_color" = line_color,
                         "plot.title" = plot.title,
                         "plot.subtitle" = plot.subtitle,
                         "plot.caption" = plot.caption,
                         "add_tag_side" = add_tag_side,
                         "order_tags_by" = order_tags_by,
                         "colors.use" = colors.use,
                         "plot.title.face" = plot.title.face,
                         "plot.subtitle.face" = plot.subtitle.face,
                         "plot.caption.face" = plot.caption.face,
                         "axis.title.face" = axis.title.face,
                         "axis.text.face" = axis.text.face,
                         "legend.title.face" = legend.title.face,
                         "legend.text.face" = legend.text.face)
  check_type(parameters = character_list, required_type = "character", test_function = is.character)

  check_colors(border.color, parameter_name = "border.color")
  check_colors(line_color, parameter_name = "line_color")
  check_colors(colors.use, parameter_name = "colors.use")

  check_parameters(parameter = font.type, parameter_name = "font.type")
  check_parameters(plot.title.face, parameter_name = "plot.title.face")
  check_parameters(plot.subtitle.face, parameter_name = "plot.subtitle.face")
  check_parameters(plot.caption.face, parameter_name = "plot.caption.face")
  check_parameters(axis.title.face, parameter_name = "axis.title.face")
  check_parameters(axis.text.face, parameter_name = "axis.text.face")
  check_parameters(legend.title.face, parameter_name = "legend.title.face")
  check_parameters(legend.text.face, parameter_name = "legend.text.face")


  assertthat::assert_that(order_tags_by %in% c("both", "pvalue", "logfc"),
                          msg = "Please use either both, pvalue or logfc in order_tags_by.")

  assertthat::assert_that(add_tag_side %in% c("both", "positive", "negative"),
                          msg = "Please use either both, positive or negative in add_tag_side")

  `%>%` <- magrittr::`%>%`
  
  if (!is.null(colors.use)){
    colors <- do_ColorPalette(colors.use, tetradic = TRUE)
    names(colors) <- c("A", "C", "B", "D")
  } else {
    colors <- c("A" = "#385f71",
                "B" = "#d7b377",
                "C" = "#f5f0f6",
                "D" = "grey75")
  }
  

  if (!("gene" %in% colnames(de_genes))){
    data <- de_genes %>%
            tibble::rownames_to_column(var = "gene")
  } else {
    data <- de_genes
  }


  data <- data %>%
          tibble::as_tibble() %>%
          dplyr::select(c("p_val_adj", "avg_log2FC", "gene")) %>%
          dplyr::mutate("p_val_adj" = replace(.data$p_val_adj, .data$p_val_adj == 0, .Machine$double.xmin)) %>%
          dplyr::mutate(log_p = -log10(.data$p_val_adj)) %>%
          dplyr::select(-"p_val_adj")

  pval_cutoff <- -log10(pval_cutoff)
  data$color <- NA
  data$color[abs(data$avg_log2FC) >= FC_cutoff & data$log_p >= pval_cutoff] <- "A"
  data$color[abs(data$avg_log2FC) < FC_cutoff & data$log_p >= pval_cutoff] <- "B"
  data$color[abs(data$avg_log2FC) < FC_cutoff & data$log_p < pval_cutoff] <- "C"
  data$color[abs(data$avg_log2FC) >= FC_cutoff & data$log_p < pval_cutoff] <- "D"

  max_value <- max(abs(c(min(data$avg_log2FC), max(data$avg_log2FC))))
  x_lims <- c(-max_value, max_value)

  # Shuffle the data.
  data <- data[sample(rownames(data), nrow(data)), ]
  p <- data %>%
       ggplot2::ggplot(mapping = ggplot2::aes(x = .data$avg_log2FC,
                                              y = .data$log_p)) +
       ggplot2::geom_point(size = pt.size * border.size,
                           color = border.color) +
       ggplot2::geom_point(mapping = ggplot2::aes(color = .data$color),
                           size = pt.size) +
       ggplot2::labs(title = plot.title,
                     subtitle = plot.subtitle,
                     caption = plot.caption) +
       ggplot2::scale_color_manual(values = colors) +
       ggplot2::guides(color = ggplot2::guide_legend(override.aes = list(size = 4),
                                                     title.position = "top",
                                                     title.hjust = 0.5)) +
       ggplot2::xlim(x_lims) +
       ggplot2::xlab(expression(bold(paste("Avg. ", log["2"], "(FC)")))) +
       ggplot2::ylab(expression(bold(paste("-", log["10"], "(p.adj.)"))))

  if (isTRUE(plot_lines)){
    p <- p +
         ggplot2::geom_hline(yintercept = pval_cutoff,
                             color = line_color,
                             linewidth = line_size,
                             linetype = "dashed") +
         ggplot2::geom_vline(xintercept = FC_cutoff,
                             color = line_color,
                             linewidth = line_size,
                             linetype = "dashed") +
         ggplot2::geom_vline(xintercept = -FC_cutoff,
                             color = line_color,
                             linewidth = line_size,
                             linetype = "dashed")
  }

  if (isTRUE(add_gene_tags)){
    if (order_tags_by == "both"){
      data.label <- data %>%
                    dplyr::mutate("abs_avg_log2FC" = abs(.data$avg_log2FC)) %>% 
                    dplyr::arrange(dplyr::desc(.data$log_p),
                                   dplyr::desc(.data$abs_avg_log2FC)) %>%
                    as.data.frame()
      
      data.up <- data.label %>% 
                 dplyr::filter(.data$avg_log2FC > 0) %>% 
                 utils::head(n_genes)
      
      data.down <- data.label %>% 
                   dplyr::filter(.data$avg_log2FC < 0) %>% 
                   utils::head(n_genes)
      
      data.label <- dplyr::bind_rows(data.up, data.down)
                    
    } else if (order_tags_by == "pvalue"){
      data.up <- data %>%
                 dplyr::filter(.data$avg_log2FC > 0) %>%
                 dplyr::arrange(dplyr::desc(.data$log_p),
                                dplyr::desc(.data$avg_log2FC)) %>%
                 as.data.frame() %>%
                 utils::head(n_genes)

      data.down <- data %>%
                   dplyr::filter(.data$avg_log2FC < 0) %>%
                   dplyr::arrange(dplyr::desc(.data$log_p)) %>%
                   as.data.frame() %>%
                   utils::head(n_genes)
      
      data.label <- dplyr::bind_rows(data.up, data.down)
    } else if (order_tags_by == "logfc"){
      data.up <- data %>%
                 dplyr::arrange(dplyr::desc(.data$avg_log2FC)) %>%
                 as.data.frame() %>%
                 utils::head(n_genes)

      data.down <- data %>%
                   dplyr::arrange(.data$avg_log2FC) %>%
                   as.data.frame() %>%
                   utils::head(n_genes)

      data.label <- dplyr::bind_rows(data.up, data.down)
    }

    if (add_tag_side == "positive") {
      data.label <- data.up

    } else if (add_tag_side == "negative") {
      data.label <- data.down

    }

    if (base::isFALSE(use_labels)){
      p <- p +
           ggrepel::geom_text_repel(data = data.label,
                                    mapping = ggplot2::aes(label = .data$gene),
                                    max.overlaps = 1000,
                                    color = "black",
                                    fontface = "bold",
                                    size = tag_size)
    } else if (isTRUE(use_labels)){
      p <- p +
           ggrepel::geom_label_repel(data = data.label,
                                     mapping = ggplot2::aes(label = .data$gene),
                                     max.overlaps = 1000,
                                     color = "black",
                                     fontface = "bold",
                                     size = tag_size)
    }

  }
  p <- p +
       ggplot2::theme_minimal(base_size = font.size) +
       ggplot2::theme(plot.margin = ggplot2::margin(t = 10, r = 10, b = 10, l = 10),
                      plot.title = ggplot2::element_text(face = plot.title.face, hjust = 0),
                      plot.subtitle = ggplot2::element_text(face = plot.subtitle.face, hjust = 0),
                      plot.caption = ggplot2::element_text(face = plot.caption.face, hjust = 1),
                      legend.text = ggplot2::element_text(face = legend.text.face),
                      legend.title = ggplot2::element_text(face = legend.title.face),
                      panel.grid = ggplot2::element_blank(),
                      plot.title.position = "plot",
                      plot.caption.position = "plot",
                      text = ggplot2::element_text(family = font.type),
                      legend.position = "none",
                      legend.justification = "center",
                      axis.title.x = ggplot2::element_text(face = axis.title.face, color = "black"),
                      axis.title.y = ggplot2::element_text(face = axis.title.face, angle = 90, color = "black"),
                      axis.text = ggplot2::element_text(face = axis.text.face, color = "black"),
                      axis.line = ggplot2::element_line(color = "black"),
                      axis.ticks = ggplot2::element_line(color = "black"),
                      plot.background = ggplot2::element_rect(fill = "white", color = "white"),
                      panel.background = ggplot2::element_rect(fill = "white", color = "white"),
                      legend.background = ggplot2::element_rect(fill = "white", color = "white"))
  return(p)
}
