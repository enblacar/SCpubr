#' Create correlation matrix heatmaps.
#'
#' @inheritParams doc_function
#' @param mode \strong{\code{\link[base]{character}}} | Different types of correlation matrices can be computed. Right now, the only possible value is "hvg", standing for Highly Variable Genes. The sample is subset for the HVG and the data is re-scaled. Scale data is used for the correlation.
#' @param cluster \strong{\code{\link[base]{logical}}} | Whether to cluster the elements in the heatmap or not.
#' @param remove.diagonal \strong{\code{\link[base]{logical}}} | Whether to convert diagnoal to NA. Normally this value would be 1, heavily shifting the color scale.
#' @return A ggplot2 object.
#' @export
#'
#' @example /man/examples/examples_do_CorrelationPlot.R
do_CorrelationPlot <- function(sample = NULL,
                               input_gene_list = NULL,
                               cluster = TRUE,
                               remove.diagonal = TRUE,
                               mode = "hvg",
                               assay = NULL,
                               group.by = NULL,
                               legend.title = "Pearson coef.",
                               enforce_symmetry = ifelse(mode == "hvg", TRUE, FALSE),
                               font.size = 14,
                               font.type = "sans",
                               na.value = "grey75",
                               legend.width = 1,
                               legend.length = 20,
                               legend.framewidth = 0.5,
                               legend.tickwidth = 0.5,
                               legend.framecolor = "grey50",
                               legend.tickcolor = "white",
                               legend.type = "colorbar",
                               legend.position = "bottom",
                               min.cutoff = NA,
                               max.cutoff = NA,
                               number.breaks = 5,
                               plot.title = NULL,
                               plot.subtitle = NULL,
                               plot.caption = NULL,
                               diverging.palette = "RdBu",
                               diverging.direction = -1,
                               use_viridis = FALSE,
                               viridis.palette = "G",
                               viridis.direction = -1,
                               sequential.palette = "YlGnBu",
                               sequential.direction = 1,
                               axis.text.x.angle = 45,
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

  `%>%` <- magrittr::`%>%`

  # Check logical parameters.
  logical_list <- list("enforce_symmetry" = enforce_symmetry,
                       "cluster" = cluster,
                       "remove.diagonal" = remove.diagonal)
  check_type(parameters = logical_list, required_type = "logical", test_function = is.logical)
  # Check numeric parameters.
  numeric_list <- list("min.cutoff" = min.cutoff,
                       "max.cutoff" = max.cutoff,
                       "number.breaks" = number.breaks,
                       "legend.width" = legend.width,
                       "legend.length" = legend.length,
                       "legend.tickwidth" = legend.tickwidth,
                       "legend.framewidth" = legend.framewidth,
                       "font.size" = font.size,
                       "axis.text.x.angle" = axis.text.x.angle,
                       "sequential.direction" = sequential.direction,
                       "viridis.direction" = viridis.direction,
                       "diverging.direction" = diverging.direction)
  check_type(parameters = numeric_list, required_type = "numeric", test_function = is.numeric)
  # Check character parameters.
  character_list <- list("mode" = mode,
                         "assay" = assay,
                         "legend.title" = legend.title,
                         "group.by" = group.by,
                         "na.value" = na.value,
                         "legend.framecolor" = legend.framecolor,
                         "legend.tickcolor" = legend.tickcolor,
                         "legend.type" = legend.type,
                         "plot.title" = plot.title,
                         "plot.subtitle" = plot.subtitle,
                         "plot.caption" = plot.caption,
                         "font.type" = font.type,
                         "diverging.palette" = diverging.palette,
                         "sequential.palette" = sequential.palette,
                         "viridis.palette" = viridis.palette,
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

  check_colors(na.value, parameter_name = "na.value")
  check_colors(legend.framecolor, parameter_name = "legend.framecolor")
  check_colors(legend.tickcolor, parameter_name = "legend.tickcolor")
  check_colors(grid.color, parameter_name = "grid.color")
  check_colors(border.color, parameter_name = "border.color")

  check_parameters(parameter = legend.position, parameter_name = "legend.position")
  check_parameters(parameter = font.type, parameter_name = "font.type")
  check_parameters(parameter = legend.type, parameter_name = "legend.type")
  check_parameters(parameter = number.breaks, parameter_name = "number.breaks")
  check_parameters(parameter = diverging.palette, parameter_name = "diverging.palette")
  check_parameters(parameter = sequential.palette, parameter_name = "sequential.palette")
  check_parameters(plot.title.face, parameter_name = "plot.title.face")
  check_parameters(plot.subtitle.face, parameter_name = "plot.subtitle.face")
  check_parameters(plot.caption.face, parameter_name = "plot.caption.face")
  check_parameters(axis.title.face, parameter_name = "axis.title.face")
  check_parameters(axis.text.face, parameter_name = "axis.text.face")
  check_parameters(legend.title.face, parameter_name = "legend.title.face")
  check_parameters(legend.text.face, parameter_name = "legend.text.face")
  check_parameters(viridis.direction, parameter_name = "viridis.direction")
  check_parameters(sequential.direction, parameter_name = "sequential.direction")
  check_parameters(diverging.direction, parameter_name = "diverging.direction")

  if (isTRUE(enforce_symmetry)){
    colors.gradient <- compute_continuous_palette(name = diverging.palette,
                                                  use_viridis = FALSE,
                                                  direction = diverging.direction,
                                                  enforce_symmetry = enforce_symmetry)
  } else {
    colors.gradient <- compute_continuous_palette(name = ifelse(isTRUE(use_viridis), viridis.palette, sequential.palette),
                                                  use_viridis = use_viridis,
                                                  direction = ifelse(isTRUE(use_viridis), viridis.direction, sequential.direction),
                                                  enforce_symmetry = enforce_symmetry)
  }


  if (mode == "hvg"){
    # Check if the sample provided is a Seurat object.
    check_Seurat(sample = sample)
    out <- check_and_set_assay(sample = sample, assay = assay)
    sample <- out[["sample"]]
    assay <- out[["assay"]]

    # Check group.by.
    out <- check_group_by(sample = sample,
                          group.by = group.by,
                          is.heatmap = TRUE)
    sample <- out[["sample"]]
    group.by <- out[["group.by"]]

    # Generate a correlation matrix of the HVG.
    variable_genes <- Seurat::VariableFeatures(sample)

    # Sort them in order (for ATAC experiments).
    suppressWarnings({
    genes <- rownames(SeuratObject::GetAssayData(object = sample,
                                    assay = assay,
                                    slot = "data"))
    genes <- data.frame("Genes" = genes) %>%
             tibble::rowid_to_column(var = "Position") %>%
             tibble::as_tibble() %>%
             dplyr::filter(.data$Genes %in% variable_genes) %>%
             dplyr::arrange(.data$Position) %>%
             dplyr::pull(.data$Genes)
    })

    # Subset sample according to the variable genes.
    sample <- sample[genes, ]
    # Scale the data
    sample <- Seurat::ScaleData(sample, features = genes, verbose = FALSE)

    # Retrieve correlation matrix.
    suppressWarnings({
    out <- sample@meta.data %>%
           dplyr::select(dplyr::all_of(c(group.by))) %>%
           tibble::rownames_to_column(var = "cell") %>%
           dplyr::left_join(y = {SeuratObject::GetAssayData(object = sample,
                                                      assay = assay,
                                                      slot = "scale.data") %>%
                                 as.matrix() %>%
                                 t() %>%
                                 as.data.frame() %>%
                                 tibble::rownames_to_column(var = "cell") %>%
                                 tidyr::pivot_longer(-"cell",
                                                     names_to = "gene",
                                                     values_to = "expression")},
                            by = "cell") %>%
           dplyr::select(-"cell") %>%
           dplyr::group_by(.data[[group.by]], .data[["gene"]]) %>%
           dplyr::summarise(mean_expression = mean(.data[["expression"]])) %>%
           tidyr::pivot_wider(names_from = dplyr::all_of(c(group.by)),
                              values_from = "mean_expression") %>%
           as.data.frame() %>%
           tibble::column_to_rownames(var = "gene") %>%
           as.matrix() %>%
           stats::cor() %>%
           round(digits = 2)
    })
    # Compute hclust.
    if (isTRUE(cluster)){
      order <- rownames(out)[stats::hclust(stats::dist(out, method = "euclidean"), method = "ward.D")$order]
    } else {
      order <- rownames(out)
    }



    out.long <- out %>%
                as.data.frame() %>%
                tibble::rownames_to_column(var = "x") %>%
                tibble::as_tibble() %>%
                tidyr::pivot_longer(cols = -"x",
                                    names_to = "y",
                                    values_to = "score") %>%
                dplyr::mutate("x" = factor(.data$x, levels = order),
                              "y" = factor(.data$y, levels = rev(order))) %>%
                dplyr::mutate("score" = ifelse(as.character(.data$x) == as.character(.data$y), ifelse(isTRUE(remove.diagonal), NA, .data$score), .data$score))

    limits <- c(min(out.long$score, na.rm = TRUE),
                max(out.long$score, na.rm = TRUE))

    # Compute scale limits, breaks, etc.
    scale.setup <- compute_scales(sample = sample,
                                  feature = " ",
                                  assay = "SCT",
                                  reduction = NULL,
                                  slot = "scale.data",
                                  number.breaks = number.breaks,
                                  min.cutoff = min.cutoff,
                                  max.cutoff = max.cutoff,
                                  flavor = "Seurat",
                                  enforce_symmetry = enforce_symmetry,
                                  from_data = TRUE,
                                  limits.use = limits)

    # Modify according to min.cutoff and max.cutoff.
    if (!is.na(min.cutoff)){
      out.long <- out.long %>%
                  dplyr::mutate("score" = ifelse(.data$score < min.cutoff, min.cutoff, .data$score))
    }

    if (!is.na(max.cutoff)){
      out.long <- out.long %>%
                  dplyr::mutate("score" = ifelse(.data$score > max.cutoff, max.cutoff, .data$score))
    }


    p <- ggplot2::ggplot(out.long,
                         mapping = ggplot2::aes(x = .data$x,
                                                y = .data$y,
                                                fill = .data$score)) +
         ggplot2::geom_tile(color = grid.color, linewidth = 0.5) +
         ggplot2::scale_y_discrete(expand = c(0, 0)) +
         ggplot2::scale_x_discrete(expand = c(0, 0),
                                   position = "top") +
         ggplot2::guides(y.sec = guide_axis_label_trans(~paste0(levels(.data$y))),
                         x.sec = guide_axis_label_trans(~paste0(levels(.data$x)))) +
         ggplot2::scale_fill_gradientn(colors = colors.gradient,
                                       na.value = na.value,
                                       name = legend.title,
                                       breaks = scale.setup$breaks,
                                       labels = scale.setup$labels,
                                       limits = scale.setup$limits) +
         ggplot2::coord_equal() +
         ggplot2::xlab(NULL) +
         ggplot2::ylab(NULL) +
         ggplot2::labs(title = plot.title,
                       subtitle = plot.subtitle,
                       caption = plot.caption)

    p <- modify_continuous_legend(p = p,
                                  legend.title = legend.title,
                                  legend.aes = "fill",
                                  legend.type = legend.type,
                                  legend.position = legend.position,
                                  legend.length = legend.length,
                                  legend.width = legend.width,
                                  legend.framecolor = legend.framecolor,
                                  legend.tickcolor = legend.tickcolor,
                                  legend.framewidth = legend.framewidth,
                                  legend.tickwidth = legend.tickwidth)

    p <- p +
         ggplot2::theme_minimal(base_size = font.size) +
         ggplot2::theme(axis.ticks.x.bottom = ggplot2::element_line(color = "black"),
                        axis.ticks.x.top = ggplot2::element_blank(),
                        axis.ticks.y.left = ggplot2::element_blank(),
                        axis.ticks.y.right = ggplot2::element_line(color = "black"),
                        axis.text.y.left = ggplot2::element_blank(),
                        axis.text.y.right = ggplot2::element_text(color = "black",
                                                                  face = axis.text.face),
                        axis.text.x.top = ggplot2::element_blank(),
                        axis.text.x.bottom = ggplot2::element_text(color = "black",
                                                                   face = axis.text.face,
                                                                   angle = get_axis_parameters(angle = axis.text.x.angle, flip = FALSE)[["angle"]],
                                                                   hjust = get_axis_parameters(angle = axis.text.x.angle, flip = FALSE)[["hjust"]],
                                                                   vjust = get_axis_parameters(angle = axis.text.x.angle, flip = FALSE)[["vjust"]]),
                        axis.title.x.bottom = ggplot2::element_blank(),
                        axis.title.x.top = ggplot2::element_text(color = "black",
                                                                 face = axis.title.face),
                        axis.title.y.right = ggplot2::element_blank(),
                        axis.title.y.left = ggplot2::element_text(color = "black",
                                                                  face = axis.title.face),
                        axis.line = ggplot2::element_blank(),
                        plot.title = ggplot2::element_text(face = plot.title.face, hjust = 0),
                        plot.subtitle = ggplot2::element_text(face = plot.subtitle.face, hjust = 0),
                        plot.caption = ggplot2::element_text(face = plot.caption.face, hjust = 1),
                        plot.title.position = "plot",
                        panel.grid = ggplot2::element_blank(),
                        panel.grid.minor.y = ggplot2::element_line(color = "white", linewidth = 1),
                        text = ggplot2::element_text(family = font.type),
                        plot.caption.position = "plot",
                        legend.text = ggplot2::element_text(face = legend.text.face),
                        legend.title = ggplot2::element_text(face = legend.title.face),
                        legend.justification = "center",
                        plot.margin = ggplot2::margin(t = 0, r = 10, b = 0, l = 40),
                        panel.border = ggplot2::element_rect(fill = NA, color = border.color, linewidth = 1),
                        panel.grid.major = ggplot2::element_blank(),
                        legend.position = legend.position,
                        plot.background = ggplot2::element_rect(fill = "white", color = "white"),
                        panel.background = ggplot2::element_rect(fill = "white", color = "white"),
                        legend.background = ggplot2::element_rect(fill = "white", color = "white"))


  } else if (mode == "jaccard"){

    # Compute jaccard indext.
    jaccard <- function(set_1, set_2) {
      # Compute intersection.
      intersection <- length(dplyr::intersect(set_1, set_2))
      # Compute the union.
      union <- length(set_1) + length(set_2) - intersection
      # Jaccard index is just the number of shared genes divided by the number of non-shared genes.
      jaccard_index <- intersection / union
      return(jaccard_index)
    }

    jaccard_scores <- list()
    for(listname_store in names(input_gene_list)){
      vector_scores <- NULL
      for(listname in names(input_gene_list)){
        scores <- jaccard(set_1 = input_gene_list[[listname_store]], set_2 = input_gene_list[[listname]])
        names(scores) <- listname
        vector_scores <- append(vector_scores, round(scores, 2))
      }
      jaccard_scores[[listname_store]] <- vector_scores
    }

    jaccard_matrix <- as.matrix(as.data.frame(jaccard_scores))
    colnames(jaccard_matrix) <- rownames(jaccard_matrix)
    if (isTRUE(cluster)){
      order <- rownames(jaccard_matrix)[stats::hclust(stats::dist(jaccard_matrix, method = "euclidean"), method = "ward.D")$order]
    } else {
      order <- rownames(jaccard_matrix)
    }

    jaccard_matrix <- jaccard_matrix[order, order]
    if (isTRUE(remove.diagonal)){
      jaccard_matrix[jaccard_matrix == 1] <- NA
    }


    data <- jaccard_matrix %>%
            as.data.frame() %>%
            tibble::rownames_to_column(var = "x") %>%
            tidyr::pivot_longer(cols = -dplyr::all_of("x"),
                                names_to = "y",
                                values_to = "score") %>%
            dplyr::mutate("x" = factor(.data$x, levels = order),
                          "y" = factor(.data$y, levels = rev(order)))

    limits <- c(min(data$score, na.rm = TRUE),
                max(data$score, na.rm = TRUE))

    assertthat::assert_that(limits[[1]] != limits[[2]],
                            msg = paste0(add_cross(), crayon_body("The "),
                                         crayon_key(" jaccard similarity matrix "),
                                         crayon_body(" has no different values. Try another gene set.")))

    scale.setup <- compute_scales(sample = NULL,
                                  feature = NULL,
                                  assay = NULL,
                                  reduction = NULL,
                                  slot = NULL,
                                  number.breaks = number.breaks,
                                  min.cutoff = min.cutoff,
                                  max.cutoff = max.cutoff,
                                  flavor = "Seurat",
                                  enforce_symmetry = FALSE,
                                  from_data = TRUE,
                                  limits.use = limits)

    # Modify according to min.cutoff and max.cutoff.
    if (!is.na(min.cutoff)){
      data <- data %>%
              dplyr::mutate("score" = ifelse(.data$score < min.cutoff, min.cutoff, .data$score))
    }

    if (!is.na(max.cutoff)){
      data <- data %>%
              dplyr::mutate("score" = ifelse(.data$score > max.cutoff, max.cutoff, .data$score))
    }

    p <- data %>%
         ggplot2::ggplot(mapping = ggplot2::aes(x = .data$x,
                                                y = .data$y,
                                                fill = .data$score)) +
         ggplot2::geom_tile(color = grid.color, linewidth = 0.5, na.rm = TRUE) +
         ggplot2::scale_y_discrete(expand = c(0, 0)) +
         ggplot2::scale_x_discrete(expand = c(0, 0),
                                   position = "top") +
         ggplot2::coord_equal() +
         ggplot2::guides(y.sec = guide_axis_label_trans(~paste0(levels(.data$y))),
                         x.sec = guide_axis_label_trans(~paste0(levels(.data$x))))

    axis.parameters <- handle_axis(flip = FALSE,
                                   group.by = "A",
                                   group = "A",
                                   counter = 1,
                                   axis.text.x.angle = axis.text.x.angle,
                                   plot.title.face = plot.title.face,
                                   plot.subtitle.face = plot.subtitle.face,
                                   plot.caption.face = plot.caption.face,
                                   axis.title.face = axis.title.face,
                                   axis.text.face = axis.text.face,
                                   legend.title.face = legend.title.face,
                                   legend.text.face = legend.text.face)

    p <- p +
         ggplot2::scale_fill_gradientn(colors = colors.gradient,
                                       na.value = na.value,
                                       name = legend.title,
                                       breaks = scale.setup$breaks,
                                       labels = scale.setup$labels,
                                       limits = scale.setup$limits)

    p <- modify_continuous_legend(p = p,
                                  legend.title = "Jaccard score",
                                  legend.aes = "fill",
                                  legend.type = legend.type,
                                  legend.position = legend.position,
                                  legend.length = legend.length,
                                  legend.width = legend.width,
                                  legend.framecolor = legend.framecolor,
                                  legend.tickcolor = legend.tickcolor,
                                  legend.framewidth = 0.5,
                                  legend.tickwidth = 0.5)

    p <- p +
         ggplot2::xlab("") +
         ggplot2::ylab("") +
         ggplot2::theme_minimal(base_size = font.size) +
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
                        legend.position = "bottom",
                        axis.line = ggplot2::element_blank(),
                        plot.title = ggplot2::element_text(face = plot.title.face, hjust = 0),
                        plot.subtitle = ggplot2::element_text(face = plot.subtitle.face, hjust = 0),
                        plot.caption = ggplot2::element_text(face = plot.caption.face, hjust = 1),
                        plot.title.position = "plot",
                        panel.grid = ggplot2::element_blank(),
                        panel.grid.minor.y = ggplot2::element_line(color = "white", linewidth = 1),
                        text = ggplot2::element_text(family = "sans"),
                        plot.caption.position = "plot",
                        legend.text = ggplot2::element_text(face = legend.text.face),
                        legend.title = ggplot2::element_text(face = legend.title.face),
                        legend.justification = "center",
                        plot.margin = ggplot2::margin(t = 10,
                                                      r = 0,
                                                      b = 0,
                                                      l = 40),
                        panel.border = ggplot2::element_rect(fill = NA, color = border.color, linewidth = 1),
                        panel.grid.major = ggplot2::element_blank(),
                        plot.background = ggplot2::element_rect(fill = "white", color = "white"),
                        panel.background = ggplot2::element_rect(fill = "white", color = "white"),
                        legend.background = ggplot2::element_rect(fill = "white", color = "white"))
  }
  return(p)
}
