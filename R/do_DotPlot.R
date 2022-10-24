#' This function is a wrapper for \link[Seurat]{DotPlot}. It provides most of its functionalities while adding extra.
#' You can
#'
#' @inheritParams doc_function
#' @param colors.use \strong{\code{\link[base]{character}}} | Two colors if split.by is not set, which will define a gradient. As many numbers as unique values in split.by, if set, which each own will define its own gradient. Defaults to predefined color scales if not provided.
#' @param dot.scale \strong{\code{\link[base]{numeric}}} | Scale the size of the dots.
#' @param cluster.idents \strong{\code{\link[base]{logical}}} | Whether to cluster the identities based on the expression of the features.
#' @param scale.by \strong{\code{\link[base]{character}}} | How to scale the size of the dots. One of:
#' \itemize{
#'   \item \emph{\code{radius}}: use radius aesthetic.
#'   \item \emph{\code{size}}: use size aesthetic.
#' }
#' @param dot_border \strong{\code{\link[base]{logical}}} | Whether to plot a border around dots.
#'
#' @return A ggplot2 object containing a Dot Plot.
#' @export
#'
#' @example man/examples/examples_do_DotPlot.R
do_DotPlot <- function(sample,
                       features,
                       assay = NULL,
                       group.by = NULL,
                       split.by = NULL,
                       legend.type = "colorbar",
                       legend.position = "bottom",
                       legend.framewidth = 1.5,
                       legend.tickwidth = 1.5,
                       legend.length = 20,
                       legend.width = 1,
                       legend.framecolor = "grey50",
                       legend.tickcolor = "white",
                       dot.scale = 6,
                       colors.use = c("#1BFFFF25", "#2E3192"),
                       plot.title = NULL,
                       plot.subtitle = NULL,
                       plot.caption = NULL,
                       xlab = NULL,
                       ylab = NULL,
                       font.size = 14,
                       font.type = "sans",
                       cluster.idents = FALSE,
                       flip = FALSE,
                       rotate_x_axis_labels = 45,
                       scale.by = "size",
                       use_viridis = FALSE,
                       viridis_color_map = "G",
                       viridis_direction = -1,
                       na.value = "grey75",
                       dot_border = TRUE,
                       plot.grid = TRUE,
                       grid.color = "grey75",
                       grid.type = "dashed"){
    check_suggests(function_name = "do_DotPlot")
    check_Seurat(sample = sample)
    # Check the assay.
    out <- check_and_set_assay(sample, assay = assay)
    sample <- out[["sample"]]
    assay <- out[["assay"]]

    # Check logical parameters.
    logical_list <- list("flip" = flip,
                         "cluster.idents" = cluster.idents,
                         "use_viridis" = use_viridis,
                         "dot_border" = dot_border,
                         "plot.grid" = plot.grid)
    check_type(parameters = logical_list, required_type = "logical", test_function = is.logical)
    # Check numeric parameters.
    numeric_list <- list("dot.scale" = dot.scale,
                         "font.size" = font.size,
                         "legend.framewidth" = legend.framewidth,
                         "legend.tickwidth" = legend.tickwidth,
                         "legend.length" = legend.length,
                         "legend.width" = legend.width,
                         "viridis_direction" = viridis_direction,
                         "rotate_x_axis_labels" = rotate_x_axis_labels)
    check_type(parameters = numeric_list, required_type = "numeric", test_function = is.numeric)
    # Check character parameters.
    character_list <- list("legend.position" = legend.position,
                           "plot.title" = plot.title,
                           "features" = unlist(features),
                           "xlab" = xlab,
                           "ylab" = ylab,
                           "colors.use" = colors.use,
                           "group.by" = group.by,
                           "split.by" = split.by,
                           "scale.by" = scale.by,
                           "legend.framecolor" = legend.framecolor,
                           "legend.tickcolor" = legend.tickcolor,
                           "legend.type" = legend.type,
                           "font.type" = font.type,
                           "viridis_color_map" = viridis_color_map,
                           "grid.color" = grid.color,
                           "grid.type" = grid.type)
    check_type(parameters = character_list, required_type = "character", test_function = is.character)

    # Check the features.
    features <- check_feature(sample = sample, features = features, permissive = TRUE)
    features <- remove_duplicated_features(features = features)

    # Check that flip is not set to TRUE and features is not a named list.
    if (isTRUE(flip)){
      assertthat::assert_that(!is.list(features),
                              msg = "Please provide the genes as a simple character vector or set flip to FALSE.")
    }

    if (is.list(features)){
      assertthat::assert_that(isFALSE(flip),
                              msg = "Please provide the genes as a simple character vector or set flip to FALSE.")
    }

    # Check that split.by is set and the user has not provided a correct vector of colors.
    if (!(is.null(split.by))){
      if (length(colors.use) != length(as.character(unique(Seurat::FetchData(sample, vars = split.by)[, 1])))){
        names.use <- if (is.factor(unique(Seurat::FetchData(sample, vars = split.by)[, 1]))){
          levels(unique(Seurat::FetchData(sample, vars = split.by)[, 1]))
        } else {
          as.character(unique(Seurat::FetchData(sample, vars = split.by)[, 1]))
        }
        colors.use <- generate_color_scale(names.use)
      }
    }
    check_parameters(parameter = font.type, parameter_name = "font.type")
    check_parameters(parameter = legend.type, parameter_name = "legend.type")
    check_parameters(parameter = legend.position, parameter_name = "legend.position")
    check_parameters(parameter = viridis_direction, parameter_name = "viridis_direction")
    check_parameters(parameter = viridis_color_map, parameter_name = "viridis_color_map")
    check_parameters(parameter = grid.type, parameter_name = "grid.type")
    check_parameters(parameter = rotate_x_axis_labels, parameter_name = "rotate_x_axis_labels")

    # Check colors.
    check_colors(colors.use)

    # Check the colors provided to legend.framecolor and legend.tickcolor.
    check_colors(legend.framecolor, parameter_name = "legend.framecolor")
    check_colors(legend.tickcolor, parameter_name = "legend.tickcolor")
    check_colors(na.value, parameter_name = "na.value")
    check_colors(grid.color, parameter_name = "grid.color")

    # Define legend parameters.
    if (legend.position %in% c("top", "bottom")){
      legend.barwidth <- legend.length
      legend.barheight <- legend.width
    } else if (legend.position %in% c("left", "right")){
      legend.barwidth <- legend.width
      legend.barheight <- round(legend.length / 2, 0)
    }


    p <- Seurat::DotPlot(sample,
                         features = features,
                         cols = colors.use,
                         group.by = group.by,
                         split.by = split.by,
                         dot.scale = dot.scale,
                         cluster.idents = cluster.idents,
                         scale.by = scale.by)
    if (isTRUE(dot_border)){
      suppressMessages({
        p <- p +
          ggplot2::geom_point(mapping = ggplot2::aes(x = p[["data"]][["features.plot"]],
                                                     y = p[["data"]][["id"]],
                                                     fill = p[["data"]][["avg.exp.scaled"]],
                                                     size = p[["data"]][["pct.exp"]]),
                              shape = 21) +
          ggplot2::scale_size_continuous(range = c(0, dot.scale))
      })

      p[["layers"]][[1]] <- NULL
    }
    if (isTRUE(use_viridis)){
      if (isFALSE(dot_border)){
        p <- add_scale(p = p,
                       function_use = ggplot2::scale_color_viridis_c(na.value = na.value,
                                                                     option = viridis_color_map,
                                                                     direction = viridis_direction),
                       scale = "color")
      } else if (isTRUE(dot_border)){
        p <- p +
             ggplot2::scale_fill_viridis_c(na.value = na.value,
                                           option = viridis_color_map,
                                           direction = viridis_direction)
      }
    } else if (isFALSE(use_viridis)){
      if (isFALSE(dot_border)){
        p <- add_scale(p = p,
                       function_use = ggplot2::scale_color_gradientn(na.value = na.value,
                                                                     colors = colors.use),
                       scale = "color")
      } else if (isTRUE(dot_border)){
        p <- p +
             ggplot2::scale_fill_gradientn(na.value = na.value,
                                           colors = colors.use)
      }
    }
    p <- p +
         ggplot2::xlab(xlab) +
         ggplot2::ylab(ylab) +
         ggplot2::labs(title = plot.title,
                       subtitle = plot.subtitle,
                       caption = plot.caption) +
         ggplot2::theme_minimal(base_size = font.size) +
         ggplot2::theme(axis.text.x = ggplot2::element_text(color = "black",
                                                            face = "bold",
                                                            angle = get_axis_parameters(angle = rotate_x_axis_labels, flip = flip)[["angle"]],
                                                            hjust = get_axis_parameters(angle = rotate_x_axis_labels, flip = flip)[["hjust"]],
                                                            vjust = get_axis_parameters(angle = rotate_x_axis_labels, flip = flip)[["vjust"]]),
                        axis.text.y = ggplot2::element_text(face = "bold", color = "black"),
                        axis.ticks = ggplot2::element_line(color = "black"),
                        axis.line = ggplot2::element_line(color = "black"),
                        axis.title = ggplot2::element_text(face = "bold"),
                        plot.title = ggplot2::element_text(face = "bold", hjust = 0),
                        plot.subtitle = ggplot2::element_text(hjust = 0),
                        plot.caption = ggplot2::element_text(hjust = 1),
                        plot.title.position = "plot",
                        panel.grid = if (isTRUE(plot.grid)){ggplot2::element_line(color = grid.color, linetype = grid.type)} else {ggplot2::element_blank()},
                        text = ggplot2::element_text(family = font.type),
                        plot.caption.position = "plot",
                        legend.text = ggplot2::element_text(face = "bold"),
                        legend.position = legend.position,
                        legend.title = ggplot2::element_text(face = "bold"),
                        legend.justification = "center",
                        plot.margin = ggplot2::margin(t = 10, r = 10, b = 10, l = 10),
                        panel.grid.major = ggplot2::element_blank(),
                        plot.background = ggplot2::element_rect(fill = "white", color = "white"),
                        panel.background = ggplot2::element_rect(fill = "white", color = "white"),
                        legend.background = ggplot2::element_rect(fill = "white", color = "white"))
    # Add leyend modifiers.
    p <- modify_continuous_legend(p = p,
                                  legend.title = "Avg. Expression",
                                  legend.aes = if (isTRUE(dot_border)) {"fill"} else {"color"},
                                  legend.type = legend.type,
                                  legend.position = legend.position,
                                  legend.length = legend.length,
                                  legend.width = legend.width,
                                  legend.framecolor = legend.framecolor,
                                  legend.tickcolor = legend.tickcolor,
                                  legend.framewidth = legend.framewidth,
                                  legend.tickwidth = legend.tickwidth)

    # Modify size legend.
    if (isTRUE(dot_border)){
      p <- p +
           ggplot2::guides(size = ggplot2::guide_legend(title = "Percent Expressed",
                                                        title.position = "top",
                                                        title.hjust = 0.5,
                                                        override.aes = ggplot2::aes(fill = "black")))
    } else {
      p <- p +
           ggplot2::guides(size = ggplot2::guide_legend(title = "Percent Expressed",
                                                        title.position = "top",
                                                        title.hjust = 0.5))
    }

    if (is.list(features)){
      p <- p + ggplot2::theme(strip.background = ggplot2::element_rect(color = 'black', fill = 'white'),
                                    strip.text = ggplot2::element_text(face = "bold", color = "black"))
    }
    if (flip == TRUE){
        p <- p + ggplot2::coord_flip()
    }
    return(p)

}
