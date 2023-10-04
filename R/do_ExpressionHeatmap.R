#' Create heatmaps of averaged expression by groups.
#'
#' This function generates a heatmap with averaged expression values by the unique groups of the metadata variables provided by the user.
#'
#' @inheritParams doc_function
#'
#'
#' @return A ggplot2 object.
#' @export
#'
#' @example /man/examples/examples_do_ExpressionHeatmap.R
do_ExpressionHeatmap <- function(sample,
                                 features,
                                 group.by = NULL,
                                 assay = NULL,
                                 cluster = TRUE,
                                 features.order = NULL,
                                 groups.order = NULL,
                                 slot = "data",
                                 legend.title = "Avg. Expression",
                                 na.value = "grey75",
                                 legend.position = "bottom",
                                 legend.width = 1,
                                 legend.length = 20,
                                 legend.framewidth = 0.5,
                                 legend.tickwidth = 0.5,
                                 legend.framecolor = "grey50",
                                 legend.tickcolor = "white",
                                 legend.type = "colorbar",
                                 font.size = 14,
                                 font.type = "sans",
                                 axis.text.x.angle = 45,
                                 enforce_symmetry = FALSE,
                                 min.cutoff = NA,
                                 max.cutoff = NA,
                                 diverging.palette = "RdBu",
                                 diverging.direction = -1,
                                 sequential.palette = "YlGnBu",
                                 sequential.direction = 1,
                                 number.breaks = 5,
                                 use_viridis = FALSE,
                                 viridis.palette = "G",
                                 viridis.direction = -1,
                                 flip = FALSE,
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

  check_suggests(function_name = "do_ExpressionHeatmap")
  # Check if the sample provided is a Seurat object.
  check_Seurat(sample = sample)

  # Check logical parameters.
  logical_list <- list("enforce_symmetry" = enforce_symmetry,
                       "use_viridis" = use_viridis,
                       "flip" = flip,
                       "cluster" = cluster)
  check_type(parameters = logical_list, required_type = "logical", test_function = is.logical)
  # Check numeric parameters.
  numeric_list <- list("axis.text.x.angle" = axis.text.x.angle,
                       "min.cutoff" = min.cutoff,
                       "max.cutoff" = max.cutoff,
                       "legend.width" = legend.width,
                       "legend.length" = legend.length,
                       "legend.framewidth" = legend.framewidth,
                       "legend.tickwidth" = legend.tickwidth,
                       "font.size" = font.size,
                       "number.breaks" = number.breaks,
                       "viridis.direction" = viridis.direction,
                       "sequential.direction" = sequential.direction,
                       "diverging.direction" = diverging.direction)
  check_type(parameters = numeric_list, required_type = "numeric", test_function = is.numeric)

  # Check character parameters.
  character_list <- list("features" = features,
                         "legend.type" = legend.type,
                         "font.type" = font.type,
                         "legend.title" = legend.title,
                         "legend.position" = legend.position,
                         "legend.framecolor" = legend.framecolor,
                         "legend.tickcolor" = legend.tickcolor,
                         "group.by" = group.by,
                         "na.value" = na.value,
                         "slot" = slot,
                         "assay" = assay,
                         "group.by" = group.by,
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

  check_colors(na.value)
  check_colors(legend.framecolor)
  check_colors(legend.tickcolor)
  check_colors(grid.color)
  check_colors(border.color)

  check_parameters(parameter = legend.position, parameter_name = "legend.position")
  check_parameters(parameter = font.type, parameter_name = "font.type")
  check_parameters(parameter = legend.type, parameter_name = "legend.type")
  check_parameters(parameter = number.breaks, parameter_name = "number.breaks")
  check_parameters(parameter = diverging.palette, parameter_name = "diverging.palette")
  check_parameters(parameter = sequential.palette, parameter_name = "sequential.palette")
  check_parameters(parameter = viridis.palette, parameter_name = "viridis.palette")
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


  `%>%` <- magrittr::`%>%`

  # Generate the continuous color palette.
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

  assay <- if (is.null(assay)){Seurat::DefaultAssay(sample)} else {assay}

  Seurat::DefaultAssay(sample) <- assay

  # Check group.by.
  out <- check_group_by(sample = sample,
                        group.by = group.by,
                        is.heatmap = TRUE)
  sample <- out[["sample"]]
  group.by <- out[["group.by"]]

  if (is.list(features)){
    warning(paste0(add_warning(), crayon_body("You have provided a "),
                   crayon_key("list"),
                   crayon_body(" to the parameter "),
                   crayon_key("features"),
                   crayon_body("Transforming into a character vector.")), call. = FALSE)
    features <- unname(unlist(features))
  }

  features <- remove_duplicated_features(features)

  # Generate the heatmap data.
  if (sum(!features %in% rownames(sample)) >= 1){
    warning(paste0(add_warning(), crayon_body("The following features are not found in the "),
                   crayon_key("row names"),
                   crayon_body(" of the specified "),
                   crayon_key("assay"),
                   crayon_body(" (default assay if not):\n"),
                   paste(vapply(features[!features %in% rownames(sample)], crayon_key, FUN.VALUE = character(1)), collapse = crayon_body(", ")), "\n"), call. = FALSE)
  }

  features <- features[features %in% rownames(sample)]

  assertthat::assert_that(length(features) >= 1,
                          msg = paste0(add_cross(), crayon_body("None of the provided "),
                                       crayon_key("features"),
                                       crayon_body(" are present in the "),
                                       crayon_key("sample"),
                                       crayon_body(".")))

  matrix.list <- list()
  for (group in group.by){
    # Extract activities from object as a long dataframe
    suppressMessages({
      sample$group.by <- sample@meta.data[, group]
      suppressWarnings({
      df <- SeuratObject::GetAssayData(sample,
                          assay = assay,
                          slot = slot)[features, , drop = FALSE] %>%
            as.matrix() %>%
            t() %>%
            as.data.frame() %>%
            tibble::rownames_to_column(var = "cell") %>%
            dplyr::left_join(y = {sample@meta.data[, "group.by", drop = FALSE] %>%
                                  tibble::rownames_to_column(var = "cell")},
                                  by = "cell") %>%
            dplyr::select(-"cell") %>%
            tidyr::pivot_longer(cols = -"group.by",
                                names_to = "gene",
                                values_to = "expression") %>%
            dplyr::group_by(.data$group.by, .data$gene) %>%
            dplyr::summarise(mean = mean(.data$expression, na.rm = TRUE))
      df.order <- df
      })
    })
    matrix.list[[group]][["df"]] <- df
    matrix.list[[group]][["df.order"]] <- df.order
  }


  counter <- 0
  for (group in group.by){
    counter <- counter + 1

    df <- matrix.list[[group]][["df"]]
    df.order <- matrix.list[[group]][["df.order"]]

    # Transform to wide to retrieve the hclust.
    df.order <- df.order %>%
                tidyr::pivot_wider(id_cols = "group.by",
                                   names_from = 'gene',
                                   values_from = 'mean') %>%
                tibble::column_to_rownames("group.by") %>%
                as.matrix()

    df.order[is.na(df.order)] <- 0
    if(length(rownames(df.order)) == 1){
      row_order <- rownames(df.order)[1]
    } else {
      if (isTRUE(cluster)){
        row_order <- rownames(df.order)[stats::hclust(stats::dist(df.order, method = "euclidean"), method = "ward.D")$order]
      } else {
        row_order <- rownames(df.order)
      }

    }
    if (counter == 1){
      if (length(colnames(df.order)) == 1){
        col_order <- colnames(df.order)[1]
      } else {
        if (isTRUE(cluster)){
          col_order <- colnames(df.order)[stats::hclust(stats::dist(t(df.order), method = "euclidean"), method = "ward.D")$order]
        } else {
          col_order <- colnames(df.order)
        }

      }
    }

    if (!is.null(groups.order) & (group %in% names(groups.order))){
      groups.order.use <- groups.order[[group]]
    } else {
      groups.order.use <- row_order
    }

    data <- df %>%
            dplyr::mutate("gene" = factor(.data$gene, levels = if (is.null(features.order)){rev(col_order)} else {features.order}),
                          "group.by" = factor(.data$group.by, levels = groups.order.use))

    if (!is.na(min.cutoff)){
      data <- data %>%
              dplyr::mutate("mean" = ifelse(.data$mean < min.cutoff, min.cutoff, .data$mean))
    }

    if (!is.na(max.cutoff)){
      data <- data %>%
              dplyr::mutate("mean" = ifelse(.data$mean > max.cutoff, max.cutoff, .data$mean))
    }

    matrix.list[[group]][["data"]] <- data
  }

  # Compute limits.
  min.vector <- NULL
  max.vector <- NULL

  for (group in group.by){
    data <- matrix.list[[group]][["data"]]

    min.vector <- append(min.vector, min(data$mean, na.rm = TRUE))
    max.vector <- append(max.vector, max(data$mean, na.rm = TRUE))
  }

  # Get the absolute limits of the datasets.
  limits <- c(min(min.vector),
              max(max.vector))

  # Compute overarching scales for all heatmaps.
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


  # Plot individual heatmaps.
  counter <- 0
  list.heatmaps <- list()
  for (group in group.by){
    counter <- counter + 1
    data <- matrix.list[[group]][["data"]]

    p <- data %>%
         # nocov start
         ggplot2::ggplot(mapping = ggplot2::aes(x = if (base::isFALSE(flip)){.data$gene} else {.data$group.by},
                                                y = if (base::isFALSE(flip)){.data$group.by} else {.data$gene},
                                                fill = .data$mean)) +
         # nocov end
         ggplot2::geom_tile(color = grid.color, linewidth = 0.5) +
         ggplot2::scale_y_discrete(expand = c(0, 0)) +
         ggplot2::scale_x_discrete(expand = c(0, 0),
                                   position = "top") +
         ggplot2::guides(y.sec = guide_axis_label_trans(~paste0(levels(.data$group.by))),
                         x.sec = guide_axis_label_trans(~paste0(levels(.data$gene)))) +
         ggplot2::coord_equal() +
         ggplot2::scale_fill_gradientn(colors = colors.gradient,
                                       na.value = na.value,
                                       name = legend.title,
                                       breaks = scale.setup$breaks,
                                       labels = scale.setup$labels,
                                       limits = scale.setup$limits)


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

    axis.parameters <- handle_axis(flip = flip,
                                   group.by = group.by,
                                   group = group,
                                   counter = counter,
                                   axis.text.x.angle = axis.text.x.angle,
                                   plot.title.face = plot.title.face,
                                   plot.subtitle.face = plot.subtitle.face,
                                   plot.caption.face = plot.caption.face,
                                   axis.title.face = axis.title.face,
                                   axis.text.face = axis.text.face,
                                   legend.title.face = legend.title.face,
                                   legend.text.face = legend.text.face)
    # nocov start
    # Set axis titles.
    if (base::isFALSE(flip)){
      if (counter == 1){
        if (length(group.by) > 1){
          xlab <- NULL
        } else {
          xlab <- "Genes"
        }

        ylab <- group
      } else {
        if (length(group.by) > 1){
          if (counter == length(group.by)){
            xlab <- "Genes"
          } else {
            xlab <- NULL
          }
        } else {
          xlab <- NULL
        }
        ylab <- group
      }
    } else {
      if (counter == 1){
        ylab <- "Genes"

        xlab <- group
      } else {
        ylab <- NULL
        xlab <- group
      }
    }
    # nocov end

    # Set theme
    p <- p +
         ggplot2::xlab(xlab) +
         ggplot2::ylab(ylab) +
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
                        plot.margin = ggplot2::margin(t = 0, r = 0, b = 0, l = 0),
                        panel.border = ggplot2::element_rect(fill = NA, color = border.color, linewidth = 1),
                        panel.grid.major = ggplot2::element_blank(),
                        legend.position = legend.position,
                        plot.background = ggplot2::element_rect(fill = "white", color = "white"),
                        panel.background = ggplot2::element_rect(fill = "white", color = "white"),
                        legend.background = ggplot2::element_rect(fill = "white", color = "white"))

    list.heatmaps[[group]] <- p
  }

  # Plot the combined plot
  input <- if(base::isFALSE(flip)){list.heatmaps[rev(group.by)]}else{list.heatmaps[group.by]}
  p <- patchwork::wrap_plots(input,
                             ncol = if (base::isFALSE(flip)){1} else {NULL},
                             nrow = if(isTRUE(flip)) {1} else {NULL},
                             guides = "collect")
  p <- p +
       patchwork::plot_annotation(theme = ggplot2::theme(legend.position = legend.position,
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
                                                         plot.caption.position = "plot"))

  return(p)
}
