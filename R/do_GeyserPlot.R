#' Generate a Geyser plot.
#'
#' A Geyser plot is a custom plot in which we plot continuous values on the Y axis grouped by a categorical value in the X. This is plotted as a dot plot, jittered so that the dots span
#' all the way to the other groups. On top of this, the mean and .66 and .95 of the data is plotted, depicting the overall distribution of the dots. The cells can, then, be colored by
#' a continuous variable (same as Y axis or different) or a categorical one (same as X axis or different).
#'
#' Special thanks to Christina Blume for coming up with the name of the plot.
#'
#' @inheritParams doc_function
#' @param scale_type \strong{\code{\link[base]{character}}} | Type of color scale to use.  One of:
#' \itemize{
#'   \item \emph{\code{categorical}}: Use a categorical color scale based on the values of "group.by".
#'   \item \emph{\code{continuous}}: Use a continuous color scale based on the values of "feature".
#' }
#' @param order_by_mean \strong{\code{\link[base]{logical}}} | Whether to order the groups by the mean of the data (highest to lowest).
#' @param jitter \strong{\code{\link[base]{numeric}}} | Amount of jitter in the plot along the X axis. The lower the value, the more compacted the dots are.
#' @param colors.use \strong{\code{\link[base]{character}}} | Named vector of colors to use. Has to match the unique values of group.by when scale_type is set to categorical.
#'
#' @return Either a plot of a list of plots, depending on the number of features provided.
#' @export
#' @example /man/examples/examples_do_GeyserPlot.R

do_GeyserPlot <- function(sample,
                          features,
                          assay = NULL,
                          slot = "data",
                          group.by = NULL,
                          split.by = NULL,
                          enforce_symmetry = FALSE,
                          scale_type = "continuous",
                          order_by_mean = TRUE,
                          plot_cell_borders = TRUE,
                          jitter = 0.45,
                          pt.size = 1,
                          border.size = 2,
                          border.color = "black",
                          legend.position = "bottom",
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
                          viridis_color_map = "G",
                          viridis_direction = 1,
                          colors.use = NULL,
                          na.value = "grey75",
                          legend.ncol = NULL,
                          legend.nrow = NULL,
                          legend.icon.size = 4,
                          legend.byrow = FALSE,
                          legend.title = NULL,
                          plot.title = NULL,
                          plot.subtitle = NULL,
                          plot.caption = NULL,
                          xlab = "Groups",
                          ylab = feature,
                          flip = FALSE,
                          min.cutoff = NULL,
                          max.cutoff = NULL){

  check_suggests(function_name = "do_GeyserPlot")
  # Check if the sample provided is a Seurat object.
  check_Seurat(sample = sample)

  # Check logical parameters.
  logical_list <- list("enforce_symmetry" = enforce_symmetry,
                       "order_by_mean" = order_by_mean,
                       "plot_cell_borders" = plot_cell_borders,
                       "flip" = flip)
  check_type(parameters = logical_list, required_type = "logical", test_function = is.logical)
  # Check numeric parameters.
  numeric_list <- list("pt.size" = pt.size,
                       "jitter" = jitter,
                       "font.size" = font.size,
                       "legend.framewidth" = legend.framewidth,
                       "legend.tickwidth" = legend.tickwidth,
                       "legend.length" = legend.length,
                       "legend.width" = legend.width,
                       "border.size" = border.size,
                       "legend.ncol" = legend.ncol,
                       "legend.nrow" = legend.nrow,
                       "legend.icon.size" = legend.icon.size,
                       "viridis_direction" = viridis_direction,
                       "rotate_x_axis_labels" = rotate_x_axis_labels)
  check_type(parameters = numeric_list, required_type = "numeric", test_function = is.numeric)
  # Check character parameters.

  character_list <- list("legend.position" = legend.position,
                         "features" = features,
                         "assay" = assay,
                         "group.by" = group.by,
                         "slot" = slot,
                         "split.by" = split.by,
                         "plot.title" = plot.title,
                         "plot.subtitle" = plot.subtitle,
                         "plot.caption" = plot.caption,
                         "scale_type" = scale_type,
                         "viridis_color_map" = viridis_color_map,
                         "legend.framecolor" = legend.framecolor,
                         "legend.tickcolor" = legend.tickcolor,
                         "legend.type" = legend.type,
                         "font.type" = font.type,
                         "border.color" = border.color,
                         "na.value" = na.value)
  # Checks
  check_type(parameters = character_list, required_type = "character", test_function = is.character)

  check_colors(border.color, parameter_name = "border.color")
  check_colors(legend.tickcolor, parameter_name = "legend.tickcolor")
  check_colors(legend.framecolor, parameter_name = "legend.framecolor")
  check_colors(na.value, parameter_name = "na.value")

  check_parameters(parameter = font.type, parameter_name = "font.type")
  check_parameters(parameter = legend.type, parameter_name = "legend.type")
  check_parameters(parameter = legend.position, parameter_name = "legend.position")
  check_parameters(parameter = viridis_direction, parameter_name = "viridis_direction")
  check_parameters(parameter = viridis_color_map, parameter_name = "viridis_color_map")
  check_parameters(parameter = scale_type, parameter_name = "scale_type")
  check_parameters(parameter = rotate_x_axis_labels, parameter_name = "rotate_x_axis_labels")

  `%>%` <- magrittr::`%>%`
  # Check the assay.
  out <- check_and_set_assay(sample = sample, assay = assay)
  sample <- out[["sample"]]
  assay <- out[["assay"]]
  rm(out)

  # Check that split.by is in metadata variables.
  if (!is.null(split.by)){
    assertthat::assert_that(split.by %in% colnames(sample@meta.data),
                            msg = "The variable for split.by has to be on the metadata of the object.")
  }


  # Check that group.by is in metadata variables.
  if (!is.null(group.by)){
    assertthat::assert_that(group.by %in% colnames(sample@meta.data),
                            msg = "The variable for group.by has to be on the metadata of the object.")
  }

  # Check that jitter is in range.
  assertthat::assert_that(jitter > 0 & jitter < 0.5,
                          msg = "Value for jitter has to be betwen 0 and 0.49.")

  # Will contain the output.
  list.out <- list()

  # Assign group.by to a metadata variable.
  if (is.null(group.by)){
    sample@meta.data[, "Groups"] <- sample@active.ident
    group.by <- "Groups"
  }

  # Iterate for each feature.
  for (feature in features){
    # Check the feature.
    check_feature(sample = sample,
                  features = feature)


    # Get a vector of all dimensional reduction compontents.
    dim_colnames <- c()
    for(red in Seurat::Reductions(object = sample)){
      col.names <- colnames(sample@reductions[[red]][[]])
      dim_colnames <- c(dim_colnames, col.names)
      if (feature %in% col.names){
        # Get the reduction in which the feature is, if this is the case.
        reduction <- red
      }
    }


    # Depending on where the feature is, generate a tibble accordingly.
    if (isTRUE(feature %in% colnames(sample@meta.data))){
      data <- sample@meta.data %>%
              dplyr::select(dplyr::all_of(c(group.by, feature))) %>%
              tibble::rownames_to_column(var = "cell") %>%
              tibble::as_tibble()
    } else if (isTRUE(feature %in% rownames(sample))){
      data <- Seurat::GetAssayData(object = sample,
                                   assay = assay,
                                   slot = slot)[feature, , drop = FALSE] %>%
              as.matrix() %>%
              t() %>%
              as.data.frame() %>%
              tibble::rownames_to_column(var = "cell") %>%
              tibble::tibble() %>%
              dplyr::left_join(y = {sample@meta.data %>%
                                    dplyr::select(dplyr::all_of(c(group.by))) %>%
                                    tibble::rownames_to_column(var = "cell")},
                               by = "cell")
    } else if (isTRUE(feature %in% dim_colnames)){
      data <- sample@reductions[[reduction]][[]][, feature, drop = FALSE] %>%
              as.data.frame() %>%
              tibble::rownames_to_column(var = "cell") %>%
              tibble::tibble() %>%
              dplyr::left_join(y = {sample@meta.data %>%
                                    dplyr::select(dplyr::all_of(c(group.by))) %>%
                                    tibble::rownames_to_column(var = "cell")},
                                    by = "cell")
    }

    # If we also want additional split.by.
    if (!(is.null(split.by))){
      data <- data %>%
              dplyr::left_join(y = {sample@meta.data %>%
                                    dplyr::select(dplyr::all_of(c(split.by))) %>%
                                    tibble::rownames_to_column(var = "cell")},
                               by = "cell") %>%
              dplyr::mutate("split.by" = .data[[split.by]]) %>%
              dplyr::select(-dplyr::all_of(c(split.by)))

    }

    # Proceed with the regular plot.
    if (isTRUE(order_by_mean)){
      data <- data %>%
              dplyr::mutate("group.by" = factor(.data[[group.by]], levels = {data %>%
                                                                             dplyr::group_by(.data[[group.by]]) %>%
                                                                             dplyr::summarise("mean" = mean(.data[[feature]])) %>%
                                                                             dplyr::arrange(dplyr::desc(.data$mean)) %>%
                                                                             dplyr::pull(.data[[group.by]]) %>%
                                                                             as.character()}),
                            "values" = .data[[feature]])
    } else if (isFALSE(order_by_mean)){
      data <- data %>%
              dplyr::mutate("group.by" = .data[[group.by]],
                            "values" = .data[[feature]])
    }

    # Get the final column names.
    cols.use <- c("values", "group.by")
    if (!(is.null(split.by))){
      cols.use <- append(cols.use, "split.by")
    }

    data <- data %>%
            dplyr::select(dplyr::all_of(cols.use))

    # Define cutoffs.
    range.data <- c(min(data[, "values"], na.rm = TRUE), max(data[, "values"], na.rm = TRUE))


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


    # Plot.
    p <- ggplot2::ggplot(data = data,
                         mapping = ggplot2::aes(x = .data[["group.by"]],
                                                y = .data[["values"]],
                                                color = if (scale_type == "categorical"){.data[["group.by"]]} else {.data[["values"]]}))

    if (isTRUE(plot_cell_borders)){
      p <- p +
           ggplot2::geom_point(position = ggplot2::position_jitter(width = jitter,
                                                                   seed = 0),
                               size = pt.size * border.size,
                               color = border.color,
                               na.rm = TRUE)
    }

    if (isTRUE(scale_type == "continuous")){
      if (isTRUE(enforce_symmetry)){
        limits <- c(min(data[, "values"], na.rm = TRUE),
                    max(data[, "values"], na.rm = TRUE))
        if (limits[1] != range.data[1]){
          limits <- c(range.data[1], limits[2])
        }

        if (limits[2] != range.data[2]){
          limits <- c(limits[1], range.data[2])
        }
        end_value <- max(abs(limits))
        scale.use <- ggplot2::scale_color_gradientn(colors = c("#033270", "#4091C9", "grey95", "#c94040", "#65010C"),
                                                    limits = c(-end_value, end_value),
                                                    na.value = na.value)

      } else if (isFALSE(enforce_symmetry)){
        scale.use <- ggplot2::scale_color_viridis_c(option = viridis_color_map,
                                                    na.value = na.value,
                                                    direction = viridis_direction,
                                                    limits = range.data)
      }
    } else if (isTRUE(scale_type == "categorical")){
      limits <- c(min(data[, "values"], na.rm = TRUE),
                  max(data[, "values"], na.rm = TRUE))
      if (limits[1] != range.data[1]){
        limits <- c(range.data[1], limits[2])
      }

      if (limits[2] != range.data[2]){
        limits <- c(limits[1], range.data[2])
      }
      end_value <- max(abs(limits))
      if (is.null(colors.use)){
        values <- data %>% dplyr::pull(.data[["values"]])
        names.use <- if (is.factor(values)){levels(values)} else {sort(unique(values))}
        colors.use <- generate_color_scale(names_use = names.use)
      } else {
        check_colors(colors.use)
      }
      scale.use <- ggplot2::scale_color_manual(values = colors.use,
                                               na.value = na.value)
    }

    p <- p +
         ggplot2::geom_point(position = ggplot2::position_jitter(width = jitter,
                                                                 seed = 0),
                             size = pt.size,
                             na.rm = TRUE) +
         ggdist::stat_pointinterval(interval_size_range = c(2, 3),
                                    fatten_point = 1.5,
                                    interval_color = "white",
                                    point_color = "white",
                                    position = ggplot2::position_dodge(width = 1),
                                    na.rm = TRUE,
                                    show.legend = FALSE) +
         ggdist::stat_pointinterval(interval_size_range = c(1, 2),
                                    interval_color = "black",
                                    point_color = "black",
                                    position = ggplot2::position_dodge(width = 1),
                                    na.rm = TRUE,
                                    show.legend = FALSE) +
         scale.use

    if (!(is.null(split.by))){
      p <- p +
           ggplot2::facet_grid(. ~ split.by)
    }
    p <- p +
         ggplot2::scale_y_continuous(labels = scales::label_number(),
                                     limits = if (isTRUE(enforce_symmetry)) {c(-end_value, end_value)} else {range.data}) +
         ggplot2::labs(title = plot.title,
                       subtitle = plot.subtitle,
                       caption = plot.caption) +
         ggplot2::xlab(xlab) +
         ggplot2::ylab(ylab) +
         ggplot2::theme_minimal(base_size = font.size) +
         ggplot2::theme(axis.title = ggplot2::element_text(color = "black",
                                                           face = "bold"),
                        axis.line.x = if (isFALSE(flip)) {ggplot2::element_line(color = "black")} else if (isTRUE(flip)) {ggplot2::element_blank()},
                        axis.line.y = if (isTRUE(flip)) {ggplot2::element_line(color = "black")} else if (isFALSE(flip)) {ggplot2::element_blank()},
                        axis.text.x = ggplot2::element_text(color = "black",
                                                            face = "bold",
                                                            angle = get_axis_parameters(angle = rotate_x_axis_labels, flip = flip)[["angle"]],
                                                            hjust = get_axis_parameters(angle = rotate_x_axis_labels, flip = flip)[["hjust"]],
                                                            vjust = get_axis_parameters(angle = rotate_x_axis_labels, flip = flip)[["vjust"]]),
                        axis.text.y = ggplot2::element_text(color = "black", face = "bold"),
                        axis.ticks = ggplot2::element_line(color = "black"),
                        panel.grid.major = ggplot2::element_blank(),
                        plot.title.position = "plot",
                        plot.title = ggplot2::element_text(face = "bold", hjust = 0),
                        plot.subtitle = ggplot2::element_text(hjust = 0),
                        plot.caption = ggplot2::element_text(hjust = 1),
                        panel.grid = ggplot2::element_blank(),
                        text = ggplot2::element_text(family = font.type),
                        plot.caption.position = "plot",
                        legend.text = ggplot2::element_text(face = "bold"),
                        legend.position = legend.position,
                        legend.title = ggplot2::element_text(face = "bold"),
                        legend.justification = "center",
                        plot.margin = ggplot2::margin(t = 10, r = 10, b = 10, l = 10),
                        plot.background = ggplot2::element_rect(fill = "white", color = "white"),
                        panel.background = ggplot2::element_rect(fill = "white", color = "white"),
                        legend.background = ggplot2::element_rect(fill = "white", color = "white"),
                        strip.text =ggplot2::element_text(color = "black", face = "bold"))

    if (isTRUE(scale_type == "continuous")){
      if (is.null(legend.title)){
        legend.title <- feature
      }
      p <- modify_continuous_legend(p = p,
                                    legend.title = legend.title,
                                    legend.aes = "color",
                                    legend.type = legend.type,
                                    legend.position = legend.position,
                                    legend.length = legend.length,
                                    legend.width = legend.width,
                                    legend.framecolor = legend.framecolor,
                                    legend.tickcolor = legend.tickcolor,
                                    legend.framewidth = legend.framewidth,
                                    legend.tickwidth = legend.tickwidth)
    } else if (isTRUE(scale_type == "categorical")){
      if (is.null(legend.title)){
        legend.title <- ""
      }
      p <- p +
           ggplot2::guides(color = ggplot2::guide_legend(title = legend.title,
                                                         ncol = legend.ncol,
                                                         nrow = legend.nrow,
                                                         byrow = legend.byrow,
                                                         override.aes = list(size = legend.icon.size),
                                                         title.position = "top",
                                                         title.hjust = 0.5))
    }

    list.out[[feature]] <- p
  }

  if (isTRUE(flip)){
    p <- p +
         ggplot2::coord_flip()
  }
  return(if (length(features) > 1) {list.out} else {p})
}
