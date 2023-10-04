#' Plot TF Activities from decoupleR using Dorothea prior knowledge.
#'
#'
#' @inheritParams doc_function
#' @param activities \strong{\code{\link[tibble]{tibble}}} | Result of running decoupleR method with dorothea regulon prior knowledge.
#' @param n_tfs \strong{\code{\link[base]{numeric}}} | Number of top regulons to consider for downstream analysis.
#' @param tfs.use \strong{\code{\link[base]{character}}} | Restrict the analysis to given regulons.
#' @param enforce_symmetry \strong{\code{\link[base]{logical}}} | Whether the geyser and feature plot has a symmetrical color scale.
#'
#' @return A ggplot2 object.
#' @export
#'
#' @example /man/examples/examples_do_TFActivityPlot.R

do_TFActivityPlot <- function(sample,
                              activities,
                              n_tfs = 25,
                              slot = "scale.data",
                              statistic = "norm_wmean",
                              tfs.use = NULL,
                              group.by = NULL,
                              split.by = NULL,
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
                              enforce_symmetry = TRUE,
                              diverging.palette = "RdBu",
                              diverging.direction = -1,
                              use_viridis = FALSE,
                              viridis.palette = "G",
                              viridis.direction = -1,
                              sequential.palette = "YlGnBu",
                              sequential.direction = 1,
                              min.cutoff = NA,
                              max.cutoff = NA,
                              number.breaks = 5,
                              flip = FALSE,
                              return_object = FALSE,
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

  check_suggests(function_name = "do_TFActivityPlot")
  # Check if the sample provided is a Seurat object.
  check_Seurat(sample = sample)

  # Check logical parameters.
  logical_list <- list("enforce_symmetry" = enforce_symmetry,
                       "flip" = flip,
                       "return_object" = return_object,
                       "use_viridis" = use_viridis)
  check_type(parameters = logical_list, required_type = "logical", test_function = is.logical)
  # Check numeric parameters.
  numeric_list <- list("n_tfs" = n_tfs,
                       "font.size" = font.size,
                       "legend.width" = legend.width,
                       "legend.length" = legend.length,
                       "legend.framewidth" = legend.framewidth,
                       "legend.tickwidth" = legend.tickwidth,
                       "axis.text.x.angle" = axis.text.x.angle,
                       "min.cutoff" = min.cutoff,
                       "max.cutoff" = max.cutoff,
                       "number.breaks" = number.breaks,
                       "viridis.direction" = viridis.direction,
                       "sequential.direction" = sequential.direction,
                       "diverging.direction" = diverging.direction)
  check_type(parameters = numeric_list, required_type = "numeric", test_function = is.numeric)
  # Check character parameters.
  character_list <- list("group.by" = group.by,
                         "slot" = slot,
                         "split.by" = split.by,
                         "na.value" = na.value,
                         "legend.position" = legend.position,
                         "legend.framecolor" = legend.framecolor,
                         "font.type" = font.type,
                         "legend.tickcolor" = legend.tickcolor,
                         "legend.type" = legend.type,
                         "tfs.use" = tfs.use,
                         "viridis.palette" = viridis.palette,
                         "sequential.palette" = sequential.palette,
                         "statistic" = statistic,
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

  `%>%` <- magrittr::`%>%`

  check_colors(legend.framecolor, parameter_name = "legend.framecolor")
  check_colors(legend.tickcolor, parameter_name = "legend.tickcolor")
  check_colors(na.value, parameter_name = "na.value")
  check_colors(grid.color, parameter_name = "grid.color")
  check_colors(border.color, parameter_name = "border.color")

  check_parameters(parameter = font.type, parameter_name = "font.type")
  check_parameters(parameter = legend.type, parameter_name = "legend.type")
  check_parameters(parameter = legend.position, parameter_name = "legend.position")
  check_parameters(parameter = axis.text.x.angle, parameter_name = "axis.text.x.angle")
  check_parameters(parameter = number.breaks, parameter_name = "number.breaks")
  check_parameters(parameter = diverging.palette, parameter_name = "diverging.palette")
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


  sample[['dorothea']] <- activities %>%
                          dplyr::filter(.data$statistic == .env$statistic) %>%
                          tidyr::pivot_wider(id_cols = 'source',
                                             names_from = 'condition',
                                             values_from = 'score') %>%
                          tibble::column_to_rownames('source') %>%
                          Seurat::CreateAssayObject()

  Seurat::DefaultAssay(sample) <- "dorothea"
  sample@assays$dorothea@key <- "dorothea_"

  # Scale the data.
  sample <- Seurat::ScaleData(sample, verbose = FALSE)

  if (!is.null(split.by) & !is.null(group.by)){
    assertthat::assert_that(length(group.by) == 1,
                            msg = paste0(add_cross(), crayon_body("When using "),
                                         crayon_key("split.by"),
                                         crayon_body(" make sure you only provide a single value to "),
                                         crayon_key("group.by"),
                                         crayon_body(". Otherwise, the prot will not keep the proportions. This is a design choice. Thanks!")))
  }


  if (is.null(group.by)) {
    sample$Groups <- Seurat::Idents(sample)
    sample$group.by <- sample$Groups
    group.by <- "Groups"
  }


  # Plotting
  list.out <- list()

    matrix.list <- list()
    list.tfs <- list()
    for (group in group.by){
      # Extract activities from object as a long dataframe
      suppressMessages({
        sample$group.by <- sample@meta.data[, group]
        suppressWarnings({
        df <- t(as.matrix(SeuratObject::GetAssayData(object = sample,
                                        assay = "dorothea",
                                        slot = slot))) %>%
              as.data.frame() %>%
              tibble::rownames_to_column(var = "cell") %>%
              dplyr::left_join(y = {sample@meta.data[, "group.by", drop = FALSE] %>%
                                    tibble::rownames_to_column(var = "cell")},
                  by = "cell") %>%
              dplyr::select(-"cell") %>%
              tidyr::pivot_longer(cols = -"group.by",
                                  names_to = "source",
                                  values_to = "score") %>%
              dplyr::group_by(.data$group.by, .data$source) %>%
              dplyr::summarise(mean = mean(.data$score, na.rm = TRUE))
        df.order <- df
        })
        if (!is.null(split.by)){
          sample$split.by <- sample@meta.data[, split.by]
          suppressWarnings({
          df.split <- t(as.matrix(SeuratObject::GetAssayData(object = sample,
                                                assay = "dorothea",
                                                slot = slot))) %>%
                      as.data.frame() %>%
                      tibble::rownames_to_column(var = "cell") %>%
                      dplyr::left_join(y = {sample@meta.data[, c("group.by", "split.by"), drop = FALSE] %>%
                                            tibble::rownames_to_column(var = "cell")},
                                            by = "cell") %>%
                      dplyr::select(-"cell") %>%
                      tidyr::pivot_longer(cols = -c("group.by", "split.by"),
                                          names_to = "source",
                                          values_to = "score") %>%
                      dplyr::group_by(.data$split.by, .data$group.by, .data$source) %>%
                      dplyr::summarise(mean = mean(.data$score, na.rm = TRUE))
          matrix.list[[group]][["df.split"]] <- df.split
          })
        }


      })

      # Get top tfs with more variable means across clusters
      tfs <- df.order %>%
             dplyr::group_by(.data$source) %>%
             dplyr::summarise(std = stats::sd(.data$mean, na.rm = TRUE)) %>%
             dplyr::arrange(-abs(.data$std)) %>%
             dplyr::slice_head(n = n_tfs) %>%
             dplyr::pull(.data$source)
      matrix.list[[group]][["df"]] <- df
      matrix.list[[group]][["df.order"]] <- df.order
      list.tfs[[group]] <- tfs
    }

    shared_tfs <- NULL
    if (is.null(tfs.use)){
      for (group in group.by){
        shared_tfs <- append(shared_tfs, list.tfs[[group]])
      }
      shared_tfs <- unique(shared_tfs)
    } else {
      shared_tfs <- unique(tfs.use[tfs.use %in% rownames(sample)])
    }

    counter <- 0
    for (group in group.by){
      counter <- counter + 1

      df <- matrix.list[[group]][["df"]]
      df.order <- matrix.list[[group]][["df.order"]]

      # Subset long data frame to top tfs and transform to wide matrix
      data <- df %>%
              dplyr::filter(.data$source %in% shared_tfs)

      if (!is.null(split.by)){
        df.split <- matrix.list[[group]][["df.split"]]
        data <- df.split %>%
                dplyr::filter(.data$source %in% shared_tfs)
      }

      # Transform to wide to retrieve the hclust.
      df.order <- df.order %>%
                  dplyr::filter(.data$source %in% shared_tfs) %>%
                  tidyr::pivot_wider(id_cols = "group.by",
                                     names_from = 'source',
                                     values_from = 'mean') %>%
                  tibble::column_to_rownames("group.by") %>%
                  as.matrix()

      df.order[is.na(df.order)] <- 0
      if(length(rownames(df.order)) == 1){
        row_order <- rownames(df.order)[1]
      } else {
        row_order <- rownames(df.order)[stats::hclust(stats::dist(df.order, method = "euclidean"), method = "ward.D")$order]
      }
      if (counter == 1){
        # nocov start
        if (length(colnames(df.order)) == 1){
          col_order <- colnames(df.order)[1]
        # nocov end
        } else {
          col_order <- colnames(df.order)[stats::hclust(stats::dist(t(df.order), method = "euclidean"), method = "ward.D")$order]
        }
      }

      data <- data %>%
              dplyr::mutate("source" = factor(.data$source, levels = rev(col_order)),
                            "group.by" = factor(.data$group.by, levels = row_order))

      matrix.list[[group]][["data.mean"]] <- data

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
      data <- matrix.list[[group]][["data.mean"]]

      min.vector <- append(min.vector, min(data$mean, na.rm = TRUE))
      max.vector <- append(max.vector, max(data$mean, na.rm = TRUE))
    }

    # Get the absolute limits of the datasets.
    limits <- c(min(min.vector),
                max(max.vector))

    # Compute overarching scales for all heatmaps.
    scale.setup <- compute_scales(sample = sample,
                                  feature = " ",
                                  assay = "dorothea",
                                  reduction = NULL,
                                  slot = slot,
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
           ggplot2::ggplot(mapping = ggplot2::aes(x = if(base::isFALSE(flip)){.data$source} else {.data$group.by},
                                                  y = if(base::isFALSE(flip)){.data$group.by} else {.data$source},
                                                  fill = .data$mean)) +
           # nocov end
           ggplot2::geom_tile(color = grid.color, linewidth = 0.5) +
           ggplot2::scale_y_discrete(expand = c(0, 0)) +
           ggplot2::scale_x_discrete(expand = c(0, 0),
                                     position = "top") +
           ggplot2::guides(y.sec = guide_axis_label_trans(~paste0(levels(.data$group.by))),
                           x.sec = guide_axis_label_trans(~paste0(levels(.data$source)))) +
           ggplot2::coord_equal() +
           ggplot2::scale_fill_gradientn(colors = colors.gradient,
                                         na.value = na.value,
                                         name = paste0("Regulon score | ", statistic, ifelse(slot == "scale.data", " | Scaled + Centered", "")),
                                         breaks = scale.setup$breaks,
                                         labels = scale.setup$labels,
                                         limits = scale.setup$limits)


      if (!is.null(split.by)){
        p <- p +
          ggplot2::facet_grid(.data$split.by ~ .,
                              drop = FALSE,
                              switch = "y")
      }

      p <- modify_continuous_legend(p = p,
                                    legend.title = paste0("Regulon score | ", statistic),
                                    legend.aes = "fill",
                                    legend.type = legend.type,
                                    legend.position = legend.position,
                                    legend.length = legend.length,
                                    legend.width = legend.width,
                                    legend.framecolor = legend.framecolor,
                                    legend.tickcolor = legend.tickcolor,
                                    legend.framewidth = legend.framewidth,
                                    legend.tickwidth = legend.tickwidth)

      # nocov start
      # Set axis titles.
      if (base::isFALSE(flip)){
        if (counter == 1){
          if (length(group.by) > 1){
            xlab <- NULL
          } else {
            xlab <- "Regulon"
          }

          ylab <- group
        } else {
          if (length(group.by) > 1){
            if (counter == length(group.by)){
              xlab <- "Regulon"
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
          ylab <- "Regulon"

          xlab <- group
        } else {
          ylab <- NULL
          xlab <- group
        }
      }
      # nocov end


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
                          strip.background = axis.parameters$strip.background,
                          strip.clip = axis.parameters$strip.clip,
                          strip.text = axis.parameters$strip.text,
                          legend.position = if (is.null(split.by)) {legend.position} else {"bottom"},
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
                          plot.margin = ggplot2::margin(t = 0, r = 10, b = 0, l = 10),
                          panel.border = ggplot2::element_rect(fill = NA, color = border.color, linewidth = 1),
                          panel.grid.major = ggplot2::element_blank(),
                          plot.background = ggplot2::element_rect(fill = "white", color = "white"),
                          panel.background = ggplot2::element_rect(fill = "white", color = "white"),
                          legend.background = ggplot2::element_rect(fill = "white", color = "white"),
                          panel.spacing.x = ggplot2::unit(0, "cm"))

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


    list.out[["Heatmap"]] <- p


    if (isTRUE(return_object)){
      list.out[["Object"]] <- sample
      return_me <- list.out
    } else{
      return_me <- list.out[["Heatmap"]]
    }

    return(return_me)
}
