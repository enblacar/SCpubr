#' Create correlation matrix heatmaps.
#'
#' @inheritParams doc_function
#' @param mode \strong{\code{\link[base]{character}}} | Different types of correlation matrices can be computed. Right now, the only possible value is "hvg", standing for Highly Variable Genes. The sample is subset for the HVG and the data is re-scaled. Scale data is used for the correlation.
#'
#' @return A ComplexHeatmap object.
#' @export
#'
#' @example /man/examples/examples_do_CorrelationPlot.R
do_CorrelationPlot <- function(sample,
                               mode = "hvg",
                               assay = NULL,
                               group.by = NULL,
                               legend.title = "Pearson coef.",
                               enforce_symmetry = TRUE,
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
                               rotate_x_axis_labels = 45){

  # Check if the sample provided is a Seurat object.
  check_Seurat(sample = sample)
  # Check logical parameters.
  logical_list <- list("enforce_symmetry" = enforce_symmetry)
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
                       "rotate_x_axis_labels" = rotate_x_axis_labels)
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
                         "diverging.palette" = diverging.palette)
  check_type(parameters = character_list, required_type = "character", test_function = is.character)

  check_colors(na.value)
  check_colors(legend.framecolor)
  check_colors(legend.tickcolor)

  check_parameters(parameter = legend.position, parameter_name = "legend.position")
  check_parameters(parameter = font.type, parameter_name = "font.type")
  check_parameters(parameter = legend.type, parameter_name = "legend.type")
  check_parameters(parameter = number.breaks, parameter_name = "number.breaks")
  check_parameters(parameter = diverging.palette, parameter_name = "diverging.palette")

  `%>%` <- magrittr::`%>%`

  out <- check_and_set_assay(sample = sample, assay = assay)
  sample <- out[["sample"]]
  assay <- out[["assay"]]

  if (mode == "hvg"){
    if (is.null(group.by)){
      assertthat::assert_that(!("Groups" %in% colnames(sample@meta.data)),
                              msg = paste0(crayon_body("Please, make sure you provide a value for "),
                                           crayon_key("group.by"),
                                           crayon_body(" and that this is not called "),
                                           crayon_key("Groups")))

      sample@meta.data[, "Groups"] <- sample@active.ident
      group.by <- "Groups"
    }

    # Generate a correlation matrix of the HVG.
    variable_genes <- Seurat::VariableFeatures(sample)

    # Subset sample according to the variable genes.
    sample <- sample[variable_genes, ]
    # Scale the data
    sample <- Seurat::ScaleData(sample, verbose = FALSE)

    expr_mat <- data.frame("rownames" = rownames(sample))

    # Retrieve correlation matrix.
    out <- sample@meta.data %>%
           dplyr::select(dplyr::all_of(c(group.by))) %>%
           tibble::rownames_to_column(var = "cell") %>%
           dplyr::left_join(y = {Seurat::GetAssayData(object = sample,
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
    assertthat::assert_that(sum(dim(out)) > 2,
                            msg = paste0(crayon_body("Please provide a variable to "),
                                         crayon_key("group.by"),
                                         crayon_body(" that does not generate a 1x1 matrix.")))
    # Compute hclust.
    if(length(rownames(out)) == 1){
      order <- rownames(out)[1]
    } else {
      order <- rownames(out)[stats::hclust(stats::dist(out, method = "euclidean"), method = "ward.D")$order]
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
                dplyr::mutate("score" = ifelse(as.character(.data$x) == as.character(.data$y), NA, .data$score))

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

    p <- ggplot2::ggplot(out.long,
                         mapping = ggplot2::aes(x = .data$x,
                                                y = .data$y,
                                                fill = .data$score)) +
         ggplot2::geom_tile(color = "white", linewidth = 0.5) +
         ggplot2::scale_y_discrete(expand = c(0, 0)) +
         ggplot2::scale_x_discrete(expand = c(0, 0),
                                   position = "top") +
         ggplot2::guides(y.sec = guide_axis_label_trans(~paste0(levels(.data$y))),
                         x.sec = guide_axis_label_trans(~paste0(levels(.data$x)))) +
         ggplot2::scale_fill_gradientn(colors = RColorBrewer::brewer.pal(n = 11, name = diverging.palette) %>% rev(),
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
                                                                  face = "bold"),
                        axis.text.x.top = ggplot2::element_blank(),
                        axis.text.x.bottom = ggplot2::element_text(color = "black",
                                                                   face = "bold",
                                                                   angle = get_axis_parameters(angle = rotate_x_axis_labels, flip = FALSE)[["angle"]],
                                                                   hjust = get_axis_parameters(angle = rotate_x_axis_labels, flip = FALSE)[["hjust"]],
                                                                   vjust = get_axis_parameters(angle = rotate_x_axis_labels, flip = FALSE)[["vjust"]]),
                        axis.title.x.bottom = ggplot2::element_blank(),
                        axis.title.x.top = ggplot2::element_text(color = "black",
                                                                 face = "bold"),
                        axis.title.y.right = ggplot2::element_blank(),
                        axis.title.y.left = ggplot2::element_text(color = "black",
                                                                  face = "bold"),
                        axis.line = ggplot2::element_blank(),
                        plot.title = ggplot2::element_text(face = "bold", hjust = 0),
                        plot.subtitle = ggplot2::element_text(hjust = 0),
                        plot.caption = ggplot2::element_text(hjust = 1),
                        plot.title.position = "plot",
                        panel.grid = ggplot2::element_blank(),
                        panel.grid.minor.y = ggplot2::element_line(color = "white", linewidth = 1),
                        text = ggplot2::element_text(family = font.type),
                        plot.caption.position = "plot",
                        legend.text = ggplot2::element_text(face = "bold"),
                        legend.title = ggplot2::element_text(face = "bold"),
                        legend.justification = "center",
                        plot.margin = ggplot2::margin(t = 0, r = 10, b = 0, l = 10),
                        panel.border = ggplot2::element_rect(fill = NA, color = "black", linewidth = 1),
                        panel.grid.major = ggplot2::element_blank(),
                        legend.position = legend.position,
                        plot.background = ggplot2::element_rect(fill = "white", color = "white"),
                        panel.background = ggplot2::element_rect(fill = "white", color = "white"),
                        legend.background = ggplot2::element_rect(fill = "white", color = "white"))


  }
  return(p)
}
