#' Compute a dotplot with the results of a group-wise DE analysis.
#'
#' @inheritParams doc_function
#' @param de_genes \strong{\code{\link[tibble]{tibble}}} | DE genes matrix resulting of running `Seurat::FindAllMarkers()`.
#' @param top_genes \strong{\code{\link[base]{numeric}}} | Top N differentially expressed (DE) genes by group to retrieve.
#' @param p.cutoff \strong{\code{\link[base]{numeric}}} | Cutoff to use for adjusted p.value to filter significant genes.
#' 
#' @return A dotplot composed of 3 main panels: -log10(adjusted p-value), log2(FC) and mean expression by cluster.
#' @export
#'
#' @example /man/examples/examples_do_GroupwiseDEHeatmap.R
do_GroupwiseDEHeatmap <- function(sample,
                               de_genes,
                               group.by = NULL,
                               assay = NULL,
                               slot = "data",
                               number.breaks = 5,
                               dot.scale = 8,
                               top_genes = 5,
                               p.cutoff = 0.05,
                               flip = FALSE,
                               plot.title = NULL,
                               plot.subtitle = NULL,
                               plot.caption = NULL,
                               xlab = NULL,
                               ylab = NULL,
                               use_viridis = FALSE,
                               colors.use = NULL,
                               colorblind = FALSE,
                               viridis.direction = -1,
                               viridis.palette = "G",
                               sequential.direction = 1,
                               sequential.palette = "YlGnBu",
                               diverging.palette = "RdBu",
                               diverging.direction = -1,
                               legend.position = "bottom",
                               legend.title = NULL,
                               legend.width = 1,
                               legend.length = 7.5,
                               legend.framewidth = 0.5,
                               legend.tickwidth = 0.5,
                               legend.framecolor = "grey50",
                               legend.tickcolor = "white",
                               legend.ncol = NULL,
                               legend.nrow = NULL,
                               legend.byrow = FALSE,
                               legend.type = "colorbar",
                               font.size = 14,
                               font.type = "sans",
                               axis.text.x.angle = 45,
                               min.cutoff = NA,
                               max.cutoff = NA,
                               enforce_symmetry = FALSE,
                               na.value = "grey75",
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

  check_suggests(function_name = "do_GroupwiseDEHeatmap")
  # Check if the sample provided is a Seurat object.
  check_Seurat(sample = sample)

  logical_list <- list("use_viridis" = use_viridis,
                       "enforce_symmetry" = enforce_symmetry,
                       "flip" = flip,
                       "legend.byrow" = legend.byrow,
                       "colorblind" = colorblind)
  check_type(parameters = logical_list, required_type = "logical", test_function = is.logical)
  # Check numeric parameters.
  numeric_list <- list("number.breaks" = number.breaks,
                       "top_genes" = top_genes,
                       "viridis.direction" = viridis.direction,
                       "legend.width" = legend.width,
                       "legend.length" = legend.length,
                       "legend.framewidth" = legend.framewidth,
                       "legend.tickwidth" = legend.tickwidth,
                       "font.size" = font.size,
                       "axis.text.x.angle" = axis.text.x.angle,
                       "min.cutoff" = min.cutoff,
                       "max.cutoff" = max.cutoff,
                       "dot.scale" = dot.scale,
                       "legend.nrow" = legend.nrow,
                       "legend.ncol" = legend.ncol)
  check_type(parameters = numeric_list, required_type = "numeric", test_function = is.numeric)
  # Check character parameters.
  character_list <- list("group.by" = group.by,
                         "slot" = slot,
                         "legend.position" = legend.position,
                         "legend.framecolor" = legend.framecolor,
                         "legend.tickcolor" = legend.tickcolor,
                         "legend.type" = legend.type,
                         "font.type" = font.type,
                         "viridis.palette" = viridis.palette,
                         "sequential.palette" = sequential.palette,
                         "na.value" = na.value,
                         "border.color" = border.color,
                         "plot.title.face" = plot.title.face,
                         "plot.subtitle.face" = plot.subtitle.face,
                         "plot.caption.face" = plot.caption.face,
                         "axis.title.face" = axis.title.face,
                         "axis.text.face" = axis.text.face,
                         "legend.title" = legend.title,
                         "legend.title.face" = legend.title.face,
                         "legend.text.face" = legend.text.face,
                         "xlab" = xlab,
                         "ylab" = ylab,
                         "plot.title" = plot.title,
                         "plot.subtitle" = plot.subtitle,
                         "plot.caption" = plot.caption)
  check_type(parameters = character_list, required_type = "character", test_function = is.character)


  `%>%` <- magrittr::`%>%`
  `:=` <- rlang::`:=`

  check_colors(legend.framecolor, parameter_name = "legend.framecolor")
  check_colors(legend.tickcolor, parameter_name = "legend.tickcolor")
  check_colors(na.value, parameter_name = "na.value")
  check_colors(border.color, parameter_name = "border.color")

  check_parameters(parameter = legend.position, parameter_name = "legend.position")
  check_parameters(parameter = viridis.palette, parameter_name = "viridis.palette")
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
  
  # Check the assay.
  out <- check_and_set_assay(sample = sample, assay = assay)
  sample <- out[["sample"]]
  assay <- out[["assay"]]

  # Check group.by.
  out <- check_group_by(sample = sample,
                        group.by = group.by,
                        is.heatmap = FALSE)
  sample <- out[["sample"]]
  group.by <- out[["group.by"]]
  
  values <- unique(sample@meta.data[, group.by])
  data.values <- unique(de_genes$cluster)
  
  # Check that group.by works with input de_genes.
  assertthat::assert_that(sum(values %in% data.values) == length(data.values))


  if (base::isFALSE(enforce_symmetry)){
    colors.gradient.exp <- compute_continuous_palette(name = ifelse(isTRUE(use_viridis), viridis.palette, sequential.palette),
                                                      use_viridis = use_viridis,
                                                      direction = ifelse(isTRUE(use_viridis), viridis.direction, sequential.direction),
                                                      enforce_symmetry = enforce_symmetry)
    
    colors.gradient.fc <- compute_continuous_palette(name = "YlOrRd",
                                                     use_viridis = FALSE,
                                                     direction = ifelse(isTRUE(use_viridis), viridis.direction, sequential.direction),
                                                     enforce_symmetry = enforce_symmetry)
    
    colors.gradient.pval <- compute_continuous_palette(name = "PuBuGn",
                                                       use_viridis = FALSE,
                                                       direction = ifelse(isTRUE(use_viridis), viridis.direction, sequential.direction),
                                                       enforce_symmetry = enforce_symmetry)
  } else {
    colors.gradient.exp <- compute_continuous_palette(name = diverging.palette,
                                                      use_viridis = FALSE,
                                                      direction = diverging.direction,
                                                      enforce_symmetry = enforce_symmetry)
    
    colors.gradient.fc <- compute_continuous_palette(name = "BrBG",
                                                     use_viridis = FALSE,
                                                     direction = diverging.direction,
                                                     enforce_symmetry = enforce_symmetry)
    
    colors.gradient.pval <- compute_continuous_palette(name = "PuOr",
                                                       use_viridis = FALSE,
                                                       direction = diverging.direction,
                                                       enforce_symmetry = enforce_symmetry)
  }
  
  if (is.null(colors.use)){
    colors.use <- generate_color_scale(names_use = if (is.factor(sample@meta.data[, group.by])) {levels(sample@meta.data[, group.by])} else {sort(unique(sample@meta.data[, group.by]))}, colorblind = colorblind)
  } else {
    check_colors(colors.use, parameter_name = "colors.use")
    check_consistency_colors_and_names(sample = sample, colors = colors.use, grouping_variable = group.by)
    colors.use <- colors.use[unique(sample@meta.data[, group.by])]
  }


  magnitude <- ifelse(slot == "data", "avg_log2FC", "avg_diff")
  specificity <- "p_val_adj"
  
  


  data.use <- de_genes %>%
              dplyr::arrange(.data[[specificity]], dplyr::desc(.data[[magnitude]])) %>%
              dplyr::group_by(.data$cluster) %>%
              dplyr::filter(.data[[specificity]] <= p.cutoff) %>% 
              dplyr::slice_head(n = top_genes)

  
  gene.order <- data.use$gene %>% unique()
  
  data.use <- data.use %>% 
              dplyr::rename("{group.by}" := "cluster") %>% 
              dplyr::select(dplyr::all_of(c("gene", group.by, magnitude, specificity))) %>% 
              dplyr::group_by(.data$gene) %>% 
              dplyr::arrange(dplyr::desc(.data[[magnitude]])) %>% 
              dplyr::slice_head(n = 1) %>% 
              dplyr::ungroup() %>% 
              dplyr::mutate("p.adj.log10.minus" := -log10(.data[[specificity]])) %>% 
              dplyr::arrange(dplyr::desc(.data[[magnitude]]), dplyr::desc(.data$p.adj.log10.minus)) %>% 
              dplyr::select(-dplyr::all_of(c(specificity)))
  
  data.use$p.adj.log10.minus[data.use$p.adj.log10.minus == Inf] <- .Machine$double.xmax
  
  
  
  # Workaround parameter deprecation.
  if (base::isTRUE(utils::packageVersion("Seurat") < "4.9.9")){
    data <- Seurat::GetAssayData(object = sample,
                                 assay = assay,
                                 slot = slot)
  } else {
    data <- SeuratObject::LayerData(object = sample,
                                    assay = assay,
                                    layer = slot)
  }
  
  expr.data <- data[data.use$gene, , drop = FALSE] %>% 
               as.data.frame() %>% 
               tibble::rownames_to_column(var = "gene") %>% 
               tidyr::pivot_longer(cols = -"gene",
                                   names_to = "Cell",
                                   values_to = "Expression") %>% 
               dplyr::left_join(y = {sample@meta.data %>% 
                                     dplyr::select(dplyr::all_of(c(group.by))) %>% 
                                     tibble::rownames_to_column(var = "Cell") %>% 
                                     dplyr::mutate("Groups.use" = .data[[group.by]]) %>% 
                                     dplyr::select(-dplyr::all_of(c(group.by)))},
                                by = "Cell")
  
  
  data.use <- expr.data %>% 
              dplyr::mutate("logical" = ifelse(.data$Expression == 0, 0, 1)) %>% 
              dplyr::group_by(.data$gene, .data$Groups.use) %>% 
              dplyr::summarise("Avg.Exp" = mean(.data$Expression, na.rm = TRUE),
                               "N.Exp" = sum(.data$logical),
                               "N" = dplyr::n()) %>% 
              dplyr::mutate("P.Exp" = (.data$N.Exp / .data$N) * 100) %>% 
              dplyr::left_join(y = {data.use %>% 
                                    dplyr::mutate("Groups.use" = .data[[group.by]]) %>% 
                                    dplyr::select(-dplyr::all_of(c(group.by)))},
                               by = c("Groups.use", "gene")) %>% 
              dplyr::select(dplyr::all_of(c("gene", "Groups.use", "Avg.Exp", "P.Exp", magnitude, "p.adj.log10.minus")))
  
  
  if (is.factor(sample@meta.data[, group.by])){
    order.use <- rev(levels(sample@meta.data[, group.by]))
  } else {
    order.use <- rev(sort(unique(sample@meta.data[, group.by])))
  }
  
  data.use <- data.use %>% 
              dplyr::mutate("gene" = factor(.data$gene, levels = gene.order)) %>% 
              dplyr::mutate("Groups.use" = factor(.data$Groups.use, levels = order.use))
  
  # Define cutoffs.
  range.data.exp <- c(min(data.use[, "Avg.Exp"], na.rm = TRUE),
                      max(data.use[, "Avg.Exp"], na.rm = TRUE))
  
  out <- check_cutoffs(min.cutoff = min.cutoff,
                       max.cutoff = max.cutoff,
                       limits = range.data.exp)
  range.data.exp <- out$limits
  
  
  scale.setup.exp <- compute_scales(sample = sample,
                                    feature = NULL,
                                    assay = assay,
                                    reduction = NULL,
                                    slot = slot,
                                    number.breaks = number.breaks,
                                    min.cutoff = min.cutoff,
                                    max.cutoff = max.cutoff,
                                    flavor = "Seurat",
                                    enforce_symmetry = enforce_symmetry,
                                    from_data = TRUE,
                                    limits.use = range.data.exp)
  
  
  range.data.fc <- c(min(data.use[, magnitude], na.rm = TRUE),
                     max(data.use[, magnitude], na.rm = TRUE))
  
  scale.setup.fc <- compute_scales(sample = sample,
                                   feature = NULL,
                                   assay = assay,
                                   reduction = NULL,
                                   slot = slot,
                                   number.breaks = number.breaks,
                                   min.cutoff = NA,
                                   max.cutoff = NA,
                                   flavor = "Seurat",
                                   enforce_symmetry = enforce_symmetry,
                                   from_data = TRUE,
                                   limits.use = range.data.fc)
  
  
  range.data.pval <- c(-log10(p.cutoff),
                       max(data.use[, "p.adj.log10.minus"], na.rm = TRUE))
  
  scale.setup.pval <- compute_scales(sample = sample,
                                     feature = NULL,
                                     assay = assay,
                                     reduction = NULL,
                                     slot = slot,
                                     number.breaks = number.breaks,
                                     min.cutoff = NA,
                                     max.cutoff = NA,
                                     flavor = "Seurat",
                                     enforce_symmetry = enforce_symmetry,
                                     from_data = TRUE,
                                     limits.use = range.data.pval)
  
  # Modify values
  if (!is.na(min.cutoff)){
    data.use$Avg.Exp <- ifelse(data.use$Avg.Exp <= min.cutoff, min.cutoff, data.use$Avg.Exp)
  }
  
  if (!is.na(max.cutoff)){
    data.use$Avg.Exp <- ifelse(data.use$Avg.Exp >= max.cutoff, max.cutoff, data.use$Avg.Exp)
  }
  
  # Plot
  
  
  p1 <- data.use %>%
        ggplot2::ggplot(mapping = ggplot2::aes(x = if (base::isFALSE(flip)){.data$gene} else {.data$Groups.use},
                                               y = if (base::isFALSE(flip)){.data$Groups.use} else {.data$gene},
                                               fill = .data$Avg.Exp,
                                               size = .data$P.Exp)) + 
        ggplot2::geom_point(color = "black", shape = 21) +
        ggplot2::scale_size_continuous(range = c(0, dot.scale)) +
        ggplot2::scale_fill_gradientn(colors = colors.gradient.exp,
                                      na.value = na.value,
                                      name = if (is.null(legend.title)){"Avg. Expression"} else {legend.title},
                                      breaks = scale.setup.exp$breaks,
                                      labels = scale.setup.exp$labels,
                                      limits = scale.setup.exp$limits) + 
        ggplot2::labs(title = plot.title,
                      subtitle = plot.subtitle,
                      caption = plot.caption,
                      x = ifelse(is.null(xlab), "Genes", xlab),
                      y = ifelse(is.null(ylab), "", ylab)) +
        ggplot2::theme_minimal(base_size = font.size) +
        ggplot2::theme(axis.text.x = ggplot2::element_text(color = "black",
                                                           face = axis.text.face,
                                                           angle = get_axis_parameters(angle = axis.text.x.angle, flip = flip)[["angle"]],
                                                           hjust = get_axis_parameters(angle = axis.text.x.angle, flip = flip)[["hjust"]],
                                                           vjust = get_axis_parameters(angle = axis.text.x.angle, flip = flip)[["vjust"]]),
                       axis.text.y = ggplot2::element_text(face = axis.text.face, color = "black"),
                       axis.ticks = ggplot2::element_line(color = "black"),
                       axis.line = ggplot2::element_line(color = "black"),
                       axis.title = ggplot2::element_text(face = axis.title.face),
                       plot.title = ggplot2::element_text(face = plot.title.face, hjust = 0),
                       plot.subtitle = ggplot2::element_text(face = plot.subtitle.face, hjust = 0),
                       plot.caption = ggplot2::element_text(face = plot.caption.face, hjust = 1),
                       plot.title.position = "plot",
                       panel.grid = ggplot2::element_blank(),
                       text = ggplot2::element_text(family = font.type),
                       plot.caption.position = "plot",
                       legend.text = ggplot2::element_text(face = legend.text.face),
                       legend.position = legend.position,
                       legend.title = ggplot2::element_text(face = legend.title.face),
                       legend.justification = "center",
                       plot.margin = ggplot2::margin(t = 20, r = 0, b = 0, l = 0),
                       panel.grid.major = ggplot2::element_blank(),
                       plot.background = ggplot2::element_rect(fill = "white", color = "white"),
                       panel.background = ggplot2::element_rect(fill = "white", color = "white"),
                       legend.background = ggplot2::element_rect(fill = "white", color = "white")) +
        ggplot2::guides(size = ggplot2::guide_legend(title = "Pct. Exp.",
                                                     title.position = "top",
                                                     title.hjust = 0.5,
                                                     ncol = 1,
                                                     nrow = legend.nrow,
                                                     byrow = legend.byrow,
                                                     override.aes = ggplot2::aes(fill = "black")))
  
  # Add leyend modifiers.
  p1 <- modify_continuous_legend(p = p1,
                                 # nocov start
                                 legend.title = if (is.null(legend.title)){"Avg. Exp."} else {legend.title},
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
  
  
  
  
  p2 <- data.use %>%
        dplyr::mutate("X_mask" = "X_label") %>% 
        dplyr::filter(!is.na(.data$p.adj.log10.minus)) %>%
        ggplot2::ggplot(mapping = ggplot2::aes(x = if (base::isFALSE(flip)){.data$gene} else {.data$X_mask},
                                               y = if (base::isFALSE(flip)){.data$X_mask} else {.data$gene},
                                               fill = .data$avg_log2FC)) + 
        ggplot2::geom_point(color = "black", shape = 22, na.rm = TRUE, size = dot.scale) +
        ggplot2::scale_fill_gradientn(colors = colors.gradient.fc,
                                      na.value = na.value,
                                      name = magnitude,
                                      breaks = scale.setup.fc$breaks,
                                      labels = scale.setup.fc$labels,
                                      limits = scale.setup.fc$limits)  +
        ggplot2::labs(title = NULL,
                      subtitle = NULL,
                      caption = NULL,
                      x = NULL,
                      y = NULL) +
        ggplot2::theme_minimal(base_size = font.size) +
        ggplot2::theme(axis.text.x = ggplot2::element_blank(),
                       axis.text.y = ggplot2::element_blank(),
                       axis.ticks = ggplot2::element_blank(),
                       axis.line = ggplot2::element_blank(),
                       axis.title = ggplot2::element_text(face = axis.title.face),
                       plot.title = ggplot2::element_text(face = plot.title.face, hjust = 0),
                       plot.subtitle = ggplot2::element_text(face = plot.subtitle.face, hjust = 0),
                       plot.caption = ggplot2::element_text(face = plot.caption.face, hjust = 1),
                       plot.title.position = "plot",
                       panel.grid = ggplot2::element_blank(),
                       text = ggplot2::element_text(family = font.type),
                       plot.caption.position = "plot",
                       legend.text = ggplot2::element_text(face = legend.text.face),
                       legend.position = legend.position,
                       legend.title = ggplot2::element_text(face = legend.title.face),
                       legend.justification = "center",
                       plot.margin = ggplot2::margin(t = 0, r = 0, b = 0, l = 0, unit = "null"),
                       panel.grid.major = ggplot2::element_blank(),
                       plot.background = ggplot2::element_rect(fill = "white", color = "white"),
                       panel.background = ggplot2::element_rect(fill = "white", color = "white"),
                       legend.background = ggplot2::element_rect(fill = "white", color = "white"))
  
  # Add leyend modifiers.
  p2 <- modify_continuous_legend(p = p2,
                                 # nocov start
                                 legend.title = expression(bold(paste("avg_", log["2"], "(FC)"))),
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
  
  
  p3 <- data.use %>%
        dplyr::mutate("X_mask" = "X_label") %>% 
        dplyr::filter(!is.na(.data$p.adj.log10.minus)) %>%
        ggplot2::ggplot(mapping = ggplot2::aes(x = if (base::isFALSE(flip)){.data$gene} else {.data$X_mask},
                                               y = if (base::isFALSE(flip)){.data$X_mask} else {.data$gene},
                                               fill = .data$p.adj.log10.minus)) + 
        ggplot2::geom_point(color = "black", shape = 22, na.rm = TRUE, size = dot.scale) +
        ggplot2::scale_fill_gradientn(colors = colors.gradient.pval,
                                      na.value = na.value,
                                      name = magnitude,
                                      breaks = scale.setup.pval$breaks,
                                      labels = scale.setup.pval$labels,
                                      limits = scale.setup.pval$limits)  +
        ggplot2::labs(title = NULL,
                      subtitle = NULL,
                      caption = NULL,
                      x = NULL,
                      y = NULL) +
        ggplot2::theme_minimal(base_size = font.size) +
        ggplot2::theme(axis.text.x = ggplot2::element_blank(),
                       axis.text.y = ggplot2::element_blank(),
                       axis.ticks = ggplot2::element_blank(),
                       axis.line = ggplot2::element_blank(),
                       axis.title = ggplot2::element_text(face = axis.title.face),
                       plot.title = ggplot2::element_text(face = plot.title.face, hjust = 0),
                       plot.subtitle = ggplot2::element_text(face = plot.subtitle.face, hjust = 0),
                       plot.caption = ggplot2::element_text(face = plot.caption.face, hjust = 1),
                       plot.title.position = "plot",
                       panel.grid = ggplot2::element_blank(),
                       text = ggplot2::element_text(family = font.type),
                       plot.caption.position = "plot",
                       legend.text = ggplot2::element_text(face = legend.text.face),
                       legend.position = legend.position,
                       legend.title = ggplot2::element_text(face = legend.title.face),
                       legend.justification = "center",
                       plot.margin = ggplot2::margin(t = 0, r = 0, b = 0, l = 0, unit = "null"),
                       panel.grid.major = ggplot2::element_blank(),
                       plot.background = ggplot2::element_rect(fill = "white", color = "white"),
                       panel.background = ggplot2::element_rect(fill = "white", color = "white"),
                       legend.background = ggplot2::element_rect(fill = "white", color = "white"))
  
  # Add leyend modifiers.
  p3 <- modify_continuous_legend(p = p3,
                                 # nocov start
                                 legend.title = expression(bold(paste("-", log["10"], "(p.adj.)"))),
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
  
  p4 <- data.use %>%
        dplyr::mutate("X_mask" = "X_label",
                      "Groups.use" = factor(.data$Groups.use, levels = rev(order.use))) %>%
        dplyr::filter(!is.na(.data$p.adj.log10.minus)) %>% 
        ggplot2::ggplot(mapping = ggplot2::aes(x = if (base::isFALSE(flip)){.data$gene} else {.data$X_mask},
                                               y = if (base::isFALSE(flip)){.data$X_mask} else {.data$gene},
                                               fill = .data$Groups.use)) + 
        ggplot2::geom_point(color = "black", shape = 22, na.rm = TRUE, size = dot.scale) +
        ggplot2::scale_size_continuous(range = c(dot.scale, dot.scale)) +
        ggplot2::scale_fill_manual(values = colors.use,
                                   na.value = na.value,
                                   name = group.by)  +
        ggplot2::labs(title = NULL,
                      subtitle = NULL,
                      caption = NULL,
                      x = NULL,
                      y = NULL) +
        ggplot2::theme_minimal(base_size = font.size) +
        ggplot2::theme(axis.text.x = ggplot2::element_blank(),
                       axis.text.y = ggplot2::element_blank(),
                       axis.ticks = ggplot2::element_blank(),
                       axis.line = ggplot2::element_blank(),
                       axis.title = ggplot2::element_text(face = axis.title.face),
                       plot.title = ggplot2::element_text(face = plot.title.face, hjust = 0),
                       plot.subtitle = ggplot2::element_text(face = plot.subtitle.face, hjust = 0),
                       plot.caption = ggplot2::element_text(face = plot.caption.face, hjust = 1),
                       plot.title.position = "plot",
                       panel.grid = ggplot2::element_blank(),
                       text = ggplot2::element_text(family = font.type),
                       plot.caption.position = "plot",
                       legend.text = ggplot2::element_text(face = legend.text.face),
                       legend.position = legend.position,
                       legend.title = ggplot2::element_text(face = legend.title.face),
                       legend.justification = "center",
                       plot.margin = ggplot2::margin(t = 0, r = 0, b = 0, l = 0, unit = "null"),
                       panel.grid.major = ggplot2::element_blank(),
                       plot.background = ggplot2::element_rect(fill = "white", color = "white"),
                       panel.background = ggplot2::element_rect(fill = "white", color = "white"),
                       legend.background = ggplot2::element_rect(fill = "white", color = "white")) +
         ggplot2::guides(fill = ggplot2::guide_legend(title = group.by,
                                                      title.position = "top",
                                                      title.hjust = 0.5,
                                                      ncol = legend.ncol,
                                                      nrow = legend.nrow,
                                                      byrow = legend.byrow))
  
  if (base::isFALSE(flip)){
    layout <- paste(c(paste(rep("A", 1), collapse = "\n"), 
                      paste(rep("B", 1), collapse = "\n"),
                      paste(rep("C", 1), collapse = "\n"),
                      paste(rep("D", length(order.use)), collapse = "\n"),
                      paste(rep("D", length(order.use)), collapse = "\n"),
                      paste(rep("D", length(order.use)), collapse = "\n"),
                      paste(rep("E", 2), collapse = "\n"),
                      paste(rep("E", 2), collapse = "\n"),
                      paste(rep("E", 2), collapse = "\n")), collapse = "\n")
  } else {
    first <- paste(c(rep("D", length(order.use)), "C", "B", "A"), collapse = "")
    layout <- paste(c(first, first, first, first, first, paste(rep("E", nchar(first)), collapse = "")), collapse = "\n")
  }
  
  
  p <- patchwork::wrap_plots(A = p4,
                             B = p3,
                             C = p2,
                             D = p1,
                             E = patchwork::guide_area(),
                             design = layout, 
                             guides = "collect") &
       patchwork::plot_annotation(theme = ggplot2::theme(legend.position = legend.position))
  
  # Return the final heatmap.
  return(p)
}
