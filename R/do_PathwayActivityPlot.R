#' Plot Pathway Activities from decoupleR using Progeny prior knowledge.
#'
#' @inheritParams doc_function
#' @param activities \strong{\code{\link[tibble]{tibble}}} | Result of running decoupleR method with progeny regulon prior knowledge.
#' @param plot_FeaturePlots \strong{\code{\link[base]{logical}}} | Compute output FeaturePlots for each of the top regulons.
#' @param plot_Heatmaps \strong{\code{\link[base]{logical}}} | Compute output heatmap showcasing the average TF activity per regulon and group.by variable.
#' @param plot_GeyserPlots \strong{\code{\link[base]{logical}}} | Compute output GeyserPlots for each of the top regulons and group.by variable.
#' @param geyser_order_by_mean \strong{\code{\link[base]{logical}}} | Whether to order the X axis by the mean of the values.
#' @param geyser_scale_type \strong{\code{\link[base]{character}}} | Type of scale to use. Either "continuous" or "categorical.
#'
#' @return A list containing several output plots according to the user's specifications.
#' @export
#'
#' @example /man/examples/examples_do_PathwayActivityPlot.R

do_PathwayActivityPlot <- function(sample,
                                   activities,
                                   group.by = NULL,
                                   split.by = NULL,
                                   plot_FeaturePlots = FALSE,
                                   plot_Heatmaps = TRUE,
                                   plot_GeyserPlots = FALSE,
                                   row_title = NULL,
                                   column_title = NULL,
                                   flip = FALSE,
                                   cluster_cols = TRUE,
                                   cluster_rows = TRUE,
                                   row_names_rot = 0,
                                   column_names_rot = 45,
                                   cell_size = 5,
                                   pt.size = 1,
                                   plot_cell_borders = TRUE,
                                   border.size = 2,
                                   na.value = "grey75",
                                   legend.position = "bottom",
                                   heatmap.legend.length = 75,
                                   heatmap.legend.width = 5,
                                   heatmap.legend.framecolor = "black",
                                   legend.width = 1,
                                   legend.length = 20,
                                   legend.framewidth = 1.5,
                                   legend.tickwidth = 1.5,
                                   legend.framecolor = "grey50",
                                   legend.tickcolor = "white",
                                   legend.type = "colorbar",
                                   font.size = 14,
                                   font.type = "sans",
                                   rotate_x_axis_labels = 45,
                                   enforce_symmetry = TRUE,
                                   geyser_order_by_mean = TRUE,
                                   geyser_scale_type = "continuous",
                                   viridis_color_map = "G",
                                   viridis_direction = 1,
                                   min.cutoff = NULL,
                                   max.cutoff = NULL){

  check_suggests(function_name = "do_PathwayActivityPlot")
  # Check if the sample provided is a Seurat object.
  check_Seurat(sample = sample)

  # Check logical parameters.
  logical_list <- list("plot_FeaturePlots" = plot_FeaturePlots,
                       "plot_Heatmaps" = plot_Heatmaps,
                       "plot_GeyserPlots" = plot_GeyserPlots,
                       "flip" = flip,
                       "cluster_cols" = cluster_cols,
                       "cluster_rows" = cluster_rows,
                       "plot_cell_borders" = plot_cell_borders,
                       "geyser_order_by_mean" = geyser_order_by_mean,
                       "enforce_symmetry" = enforce_symmetry)
  check_type(parameters = logical_list, required_type = "logical", test_function = is.logical)
  # Check numeric parameters.
  numeric_list <- list("row_names_rot" = row_names_rot,
                       "column_names_rot" = column_names_rot,
                       "cell_size" = cell_size,
                       "heatmap.legend.length" = heatmap.legend.length,
                       "heatmap.legend.width" = heatmap.legend.width,
                       "pt.size" = pt.size,
                       "border.size" = border.size,
                       "font.size" = font.size,
                       "legend.width" = legend.width,
                       "legend.length" = legend.length,
                       "legend.framewidth" = legend.framewidth,
                       "legend.tickwidth" = legend.tickwidth,
                       "viridis_direction" = viridis_direction,
                       "rotate_x_axis_labels" = rotate_x_axis_labels,
                       "min.cutoff" = min.cutoff,
                       "max.cutoff" = max.cutoff)
  check_type(parameters = numeric_list, required_type = "numeric", test_function = is.numeric)
  # Check character parameters.
  character_list <- list("group.by" = group.by,
                         "split.by" = split.by,
                         "na.value" = na.value,
                         "legend.position" = legend.position,
                         "legend.framecolor" = legend.framecolor,
                         "font.type" = font.type,
                         "legend.tickcolor" = legend.tickcolor,
                         "legend.type" = legend.type,
                         "geyser_scale_type" = geyser_scale_type,
                         "viridis_color_map" = viridis_color_map)
  check_type(parameters = character_list, required_type = "character", test_function = is.character)

  check_colors(legend.framecolor, parameter_name = "legend.framecolor")
  check_colors(legend.tickcolor, parameter_name = "legend.tickcolor")
  check_colors(heatmap.legend.framecolor, parameter_name = "heatmap.legend.framecolor")
  check_colors(na.value, parameter_name = "na.value")

  check_parameters(parameter = font.type, parameter_name = "font.type")
  check_parameters(parameter = legend.type, parameter_name = "legend.type")
  check_parameters(parameter = legend.position, parameter_name = "legend.position")
  check_parameters(parameter = viridis_direction, parameter_name = "viridis_direction")
  check_parameters(parameter = viridis_color_map, parameter_name = "viridis_color_map")
  check_parameters(parameter = rotate_x_axis_labels, parameter_name = "rotate_x_axis_labels")


  `%v%` <- ComplexHeatmap::`%v%`
  `%>%` <- magrittr::`%>%`

  sample[["progeny"]] <- activities %>%
                         dplyr::filter(.data$statistic == "norm_wmean") %>%
                         tidyr::pivot_wider(id_cols = "source",
                                            names_from = "condition",
                                            values_from = "score") %>%
                         tibble::column_to_rownames("source") %>%
                         Seurat::CreateAssayObject()

  Seurat::DefaultAssay(sample) <- "progeny"
  # Scale data.
  sample <- Seurat::ScaleData(object = sample, verbose = FALSE)

  list.out <- list()

  if (isTRUE(plot_FeaturePlots)){
    list.features <- list()
    for (pathway in sort(rownames(sample))){
      p <- do_FeaturePlot(sample = sample,
                          features = pathway,
                          assay = "progeny",
                          reduction = "umap",
                          slot = "scale.data",
                          pt.size = pt.size,
                          order = FALSE,
                          border.size = border.size,
                          enforce_symmetry = enforce_symmetry,
                          plot_cell_borders = plot_cell_borders,
                          font.size = font.size,
                          font.type = font.type,
                          legend.title = paste0(pathway, " activity"),
                          legend.position = legend.position,
                          legend.type = legend.type,
                          legend.framecolor = legend.framecolor,
                          legend.tickcolor = legend.tickcolor,
                          legend.framewidth = legend.framewidth,
                          legend.tickwidth = legend.tickwidth,
                          legend.length = legend.length,
                          legend.width = legend.width,
                          viridis_color_map = viridis_color_map,
                          viridis_direction = viridis_direction,
                          min.cutoff = if (is.null(min.cutoff)) {NA} else {min.cutoff},
                          max.cutoff = if (is.null(max.cutoff)) {NA} else {max.cutoff})

      list.features[[pathway]] <- p
    }
    list.out[["feature_plots"]] <- list.features
  }

  if (isTRUE(plot_GeyserPlots)){
    list.geysers <- list()
    for (pathway in sort(rownames(sample))){
      p <- do_GeyserPlot(sample = sample,
                         assay = "progeny",
                         slot = "scale.data",
                         features = pathway,
                         group.by = group.by,
                         pt.size = pt.size,
                         border.size = border.size,
                         enforce_symmetry = enforce_symmetry,
                         scale_type = geyser_scale_type,
                         order_by_mean = geyser_order_by_mean,
                         plot_cell_borders = plot_cell_borders,
                         font.size = font.size,
                         font.type = font.type,
                         legend.position = legend.position,
                         legend.type = legend.type,
                         legend.framecolor = legend.framecolor,
                         legend.tickcolor = legend.tickcolor,
                         legend.framewidth = legend.framewidth,
                         legend.tickwidth = legend.tickwidth,
                         legend.length = legend.length,
                         legend.width = legend.width,
                         xlab = if (is.null(group.by)) {"Clusters"} else {group.by},
                         ylab = paste0(pathway, " activity"),
                         legend.title = paste0(pathway, " activity"),
                         rotate_x_axis_labels = rotate_x_axis_labels,
                         viridis_color_map = viridis_color_map,
                         viridis_direction = viridis_direction,
                         min.cutoff = min.cutoff,
                         max.cutoff = max.cutoff)
      list.geysers[[pathway]] <- p
    }
    list.out[["geyser_plots"]] <- list.geysers
  }


  if (isTRUE(plot_Heatmaps)){
    # Start the process.
    if (is.null(group.by)){
      sample@meta.data[, "dummy"] <- sample@active.ident
    } else {
      sample@meta.data[, "dummy"] <- sample@meta.data[, group.by]
    }

    group.by <- "dummy"
    if (is.null(split.by)){
      suppressMessages({
        data <- sample@assays$progeny@scale.data %>%
                t() %>%
                as.data.frame() %>%
                tibble::rownames_to_column(var = "cell") %>%
                dplyr::left_join(y = {sample@meta.data[, group.by, drop = FALSE] %>%
                                      tibble::rownames_to_column(var = "cell")},
                                      by = "cell") %>%
                dplyr::select(-dplyr::all_of(c("cell"))) %>%
                tidyr::pivot_longer(cols = -dplyr::all_of(c(group.by)),
                                    names_to = "source",
                                    values_to = "score") %>%
                dplyr::group_by(.data[[group.by]], .data$source) %>%
                dplyr::summarise(mean = mean(.data$score)) %>%
                tidyr::pivot_wider(id_cols = dplyr::all_of(c(group.by)),
                                   names_from = 'source',
                                   values_from = 'mean') %>%
                tibble::column_to_rownames(group.by) %>%
                as.matrix()
      })

      row_title <- {
        if (!(is.null(row_title))){
          row_title
        } else {
          ""
        }}

      column_title <- {
        if (!(is.null(column_title))){
          column_title
        } else {
          ""
        }}

      if (isTRUE(flip)){
        data <- t(data)
      }

      range.data <- c(min(data), max(data))
      if (!is.null(min.cutoff) & !is.null(max.cutoff)){
        assertthat::assert_that(min.cutoff < max.cutoff,
                                msg = paste0("The value provided for min.cutoff (", min.cutoff, ") has to be lower than the value provided to max.cutoff (", max.cutoff, "). Please select another value."))

        assertthat::assert_that(max.cutoff > min.cutoff,
                                msg = paste0("The value provided for max.cutoff (", max.cutoff, ") has to be higher than the value provided to min.cutoff (", min.cutoff, "). Please select another value."))

        assertthat::assert_that(max.cutoff != min.cutoff,
                                msg = paste0("The value provided for max.cutoff (", max.cutoff, ") can not be the same than the value provided to min.cutoff (", min.cutoff, "). Please select another value."))

      }

      if (!is.null(min.cutoff)){
        assertthat::assert_that(min.cutoff >= range.data[1],
                                msg = paste0("The value provided for min.cutoff (", min.cutoff, ") is lower than the minimum value in the enrichment matrix (", range.data[1], "). Please select another value."))
        range.data <- c(min.cutoff, range.data[2])
      }

      if (!is.null(max.cutoff)){
        assertthat::assert_that(max.cutoff <= range.data[2],
                                msg = paste0("The value provided for max.cutoff (", max.cutoff, ") is lower than the maximum value in the enrichment matrix (", range.data[2], "). Please select another value."))
        range.data <- c(range.data[1], max.cutoff)
      }

      out <- heatmap_inner(data,
                           legend.title = "Pathway activity",
                           column_title = column_title,
                           row_title = row_title,
                           cluster_columns = cluster_cols,
                           cluster_rows = cluster_rows,
                           column_names_rot = column_names_rot,
                           row_names_rot = row_names_rot,
                           cell_size = cell_size,
                           na.value = na.value,
                           data_range = "both",
                           range.data = range.data,
                           legend.position = legend.position,
                           legend.length = heatmap.legend.length,
                           legend.width = heatmap.legend.width,
                           legend.framecolor = heatmap.legend.framecolor,
                           symmetrical_scale = enforce_symmetry)
      h <- out[["heatmap"]]
      h_legend <- out[["legend"]]
      ComplexHeatmap::ht_opt("HEATMAP_LEGEND_PADDING" = ggplot2::unit(8, "mm"))
      suppressWarnings({
        grDevices::pdf(NULL)
        h <- ComplexHeatmap::draw(h,
                                  heatmap_legend_list = h_legend,
                                  heatmap_legend_side = if (legend.position %in% c("top", "bottom")){"bottom"} else {"right"},
                                  padding = ggplot2::unit(c(5, 5, 5, 5), "mm"))
        grDevices::dev.off()
      })
    } else {
      split.values <- as.character(sort(unique(sample@meta.data %>% dplyr::pull(!!rlang::sym(split.by)))))
      list.heatmaps <- list()
      # Get the maximum range.
      data <- sample@assays$progeny@scale.data %>%
        t() %>%
        as.data.frame() %>%
        tibble::rownames_to_column(var = "cell") %>%
        dplyr::left_join(y = {sample@meta.data[, c(group.by, split.by)] %>%
            tibble::rownames_to_column(var = "cell")},
            by = "cell") %>%
        dplyr::select(-dplyr::all_of(c("cell", split.by))) %>%
        tidyr::pivot_longer(cols = -dplyr::all_of(c(group.by)),
                            names_to = "source",
                            values_to = "score") %>%
        dplyr::group_by(.data[[group.by]], .data$source) %>%
        dplyr::summarise(mean = mean(.data$score)) %>%
        tidyr::pivot_wider(id_cols = dplyr::all_of(c(group.by)),
                           names_from = 'source',
                           values_from = 'mean') %>%
        tibble::column_to_rownames(group.by) %>%
        as.matrix()
      range.data <- c(min(data), max(data))

      if (!is.null(min.cutoff)){
        assertthat::assert_that(min.cutoff >= range.data[1],
                                msg = paste0("The value provided for min.cutoff (", min.cutoff, ") is lower than the minimum value in the enrichment matrix (", range.data[1], "). Please select another value."))
        range.data <- c(min.cutoff, range.data[2])
      }

      if (!is.null(max.cutoff)){
        assertthat::assert_that(max.cutoff <= range.data[2],
                                msg = paste0("The value provided for max.cutoff (", max.cutoff, ") is lower than the maximum value in the enrichment matrix (", range.data[2], "). Please select another value."))
        range.data <- c(range.data[1], max.cutoff)

      }
      for (split.value in split.values){
        suppressMessages({
          data <- sample@assays$progeny@scale.data %>%
                  t() %>%
                  as.data.frame() %>%
                  tibble::rownames_to_column(var = "cell") %>%
                  dplyr::left_join(y = {sample@meta.data[, c(group.by, split.by)] %>%
                                        tibble::rownames_to_column(var = "cell")},
                                        by = "cell") %>%
                  dplyr::filter(.data[[split.by]] == split.value) %>%  # This is key.
                  dplyr::select(-dplyr::all_of(c("cell", split.by))) %>%
                  tidyr::pivot_longer(cols = -dplyr::all_of(c(group.by)),
                                      names_to = "source",
                                      values_to = "score") %>%
                  dplyr::group_by(.data[[group.by]], .data$source) %>%
                  dplyr::summarise(mean = mean(.data$score)) %>%
                  tidyr::pivot_wider(id_cols = dplyr::all_of(c(group.by)),
                                     names_from = 'source',
                                     values_from = 'mean') %>%
                  tibble::column_to_rownames(group.by) %>%
                  as.matrix()
        })

        row_title <- split.value
        column_title <- ""

        if (isTRUE(flip)){
          data <- t(data)
        }


        out <- heatmap_inner(data,
                             legend.title = "Pathway activity",
                             column_title = column_title,
                             row_title = row_title,
                             cluster_columns = cluster_cols,
                             cluster_rows = cluster_rows,
                             range.data = range.data,
                             column_names_rot = column_names_rot,
                             row_names_rot = row_names_rot,
                             cell_size = cell_size,
                             na.value = na.value,
                             legend.position = legend.position,
                             legend.length = heatmap.legend.length,
                             legend.width = heatmap.legend.width,
                             legend.framecolor = heatmap.legend.framecolor,
                             symmetrical_scale = enforce_symmetry)
        h <- out[["heatmap"]]
        h_legend <- out[["legend"]]
        list.heatmaps[[split.value]] <- h
      }
      ComplexHeatmap::ht_opt("HEATMAP_LEGEND_PADDING" = ggplot2::unit(8, "mm"))
      suppressWarnings({
        grDevices::pdf(NULL)
        ht_list <- NULL
        for (name in names(list.heatmaps)){
          ht_list <- ht_list %v% list.heatmaps[[name]]
        }
        h <- ComplexHeatmap::draw(ht_list,
                                  heatmap_legend_list = h_legend,
                                  heatmap_legend_side = if (legend.position %in% c("top", "bottom")){"bottom"} else {"right"},
                                  padding = ggplot2::unit(c(5, 5, 5, 5), "mm"))
        grDevices::dev.off()
      })
    }
    list.out[["heatmaps"]] <- list("average_scores" = h)
  }
  return(list.out)
}
