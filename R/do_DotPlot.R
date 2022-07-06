#' Wrapper for \link[Seurat]{DotPlot}.
#'
#'
#' @param sample Seurat object.
#' @param features Features to represent.
#' @param assay Assay to use. Defaults to the current assay.
#' @param group.by Variable you want the cells to be colored for.
#' @param split.by Split into as many plots as unique values in the variable provided.
#' @param colors.use Two colors if split.by is not set, which will define a gradient. As many numbers as unique values in split.by, if set, which each own will define its own gradient. Defaults to predefined color scales if not provided.
#' @param legend Whether to plot the legend or not.
#' @param legend.type Character. Type of legend to display. One of: normal, colorbar, colorsteps.
#' @param legend.position Position of the legend in the plot. Will only work if legend is set to TRUE.
#' @param legend.framewidth,legend.tickwidth Width of the lines of the box in the legend.
#' @param legend.framecolor,legend.tickcolor Color of the lines of the box in the legend.
#' @param legend.length,legend.width Length and width of the legend. Will adjust automatically depending on legend side.
#' @param plot.title,plot.subtitle,plot.caption Title, subtitle or caption to use in the plot.
#' @param plot.title,plot.subtitle,plot.caption Title to use in the plot.
#' @param xlab Title for the X axis.
#' @param ylab Title for the Y axis.
#' @param font.size Base font.size of the plot.
#' @param font.type Character. Base font for the plot. One of mono, serif or sans.
#' @param flip Whether to flip the axis.
#' @param dot.scale Scale the size of the dots.
#' @param cluster.idents Logical. Whether to cluster the identities based on the expression of the features.
#' @param rotate_x_labels Logical. Whether to rotate X axis labels to horizontal or not. If multiple features, a vector of logical values of the same length.
#' @param scale.by Whether to scale the size of the dots by radius or size aesthetic.
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
                       legend = TRUE,
                       legend.type = "colorbar",
                       legend.position = "bottom",
                       legend.framewidth = 1.5,
                       legend.tickwidth = 1.5,
                       legend.length = 20,
                       legend.width = 1,
                       legend.framecolor = "grey50",
                       legend.tickcolor = "white",
                       dot.scale = 6,
                       colors.use = c("#bdc3c7", "#2c3e50"),
                       plot.title = NULL,
                       plot.subtitle = NULL,
                       plot.caption = NULL,
                       xlab = NULL,
                       ylab = NULL,
                       font.size = 14,
                       font.type = "sans",
                       cluster.idents = FALSE,
                       flip = FALSE,
                       rotate_x_labels = NULL,
                       scale.by = "size"){
    # Checks for packages.
    check_suggests(function_name = "do_DotPlot")
    # Check the assay.
    out <- check_and_set_assay(sample, assay = assay)
    sample <- out[["sample"]]
    assay <- out[["assay"]]

    # Check logical parameters.
    logical_list <- list("legend" = legend,
                         "flip" = flip,
                         "cluster.idents" = cluster.idents,
                         "rotate_x_labels" = rotate_x_labels)
    check_type(parameters = logical_list, required_type = "logical", test_function = is.logical)
    # Check numeric parameters.
    numeric_list <- list("dot.scale" = dot.scale,
                         "font.size" = font.size,
                         "legend.framewidth" = legend.framewidth,
                         "legend.tickwidth" = legend.tickwidth,
                         "legend.length" = legend.length,
                         "legend.width" = legend.width)
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
                           "font.type" = font.type)
    check_type(parameters = character_list, required_type = "character", test_function = is.character)

    # Check the features.
    features <- check_feature(sample = sample, features = features, permissive = TRUE)
    features <- remove_duplicated_features(features = features)

    # Check that flip is not set to TRUE and features is not a named list.
    if (isTRUE(flip) & is.list(features)){stop("Please provide the genes as a simple character vector or set flip to FALSE.", call. = F)}
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
    # Check font.type.
    if (!(font.type %in% c("sans", "serif", "mono"))){
      stop("Please select one of the following for font.type: sans, serif, mono.", call. = F)
    }

    # Check colors.
    check_colors(colors.use)

    # Check the colors provided to legend.framecolor and legend.tickcolor.
    check_colors(legend.framecolor, parameter_name = "legend.framecolor")
    check_colors(legend.tickcolor, parameter_name = "legend.tickcolor")

    # Check the legend.type.
    if (!(legend.type %in% c("normal", "colorbar", "colorsteps"))){
      stop("Please select one of the following for legend.type: normal, colorbar, colorsteps.", call. = FALSE)
    }

    # Check the legend.position.
    if (!(legend.position %in% c("top", "bottom", "left", "right"))){
      stop("Please select one of the following for legend.position: top, bottom, left, right.", call. = FALSE)
    }

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
                         scale.by = scale.by) &
         ggplot2::theme_minimal(base_size = font.size) &
         ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1, face = "bold", color = "black"),
                        axis.text.y = ggplot2::element_text(face = "bold", color = "black"),
                        axis.line = ggplot2::element_line(color = "black"),
                        axis.title = ggplot2::element_text(face = "bold"),
                        plot.title = ggtext::element_markdown(face = "bold", hjust = 0),
                        plot.subtitle = ggtext::element_markdown(hjust = 0),
                        plot.caption = ggtext::element_markdown(hjust = 1),
                        plot.title.position = "plot",
                        panel.grid = ggplot2::element_blank(),
                        text = ggplot2::element_text(family = font.type),
                        plot.caption.position = "plot",
                        legend.text = ggplot2::element_text(face = "bold"),
                        legend.position = legend.position,
                        legend.title = ggplot2::element_text(face = "bold"),
                        legend.justification = "center",
                        plot.margin = ggplot2::margin(t = 10, r = 10, b = 10, l = 10),
                        panel.grid.major = ggplot2::element_blank(),
                        plot.background = ggplot2::element_rect(fill = "white", color = "white"),)
    # Add leyend modifiers.
    if (legend.type == "normal"){
      p <- p +
        ggplot2::guides(color = ggplot2::guide_colorbar(title = "Avg. Expression",
                                                        title.position = "top",
                                                        title.hjust = 0.5))
    } else if (legend.type == "colorbar"){
      p <- p +
        ggplot2::guides(color = ggplot2::guide_colorbar(title = "Avg. Expression",
                                                        title.position = "top",
                                                        barwidth = legend.barwidth,
                                                        barheight = legend.barheight,
                                                        title.hjust = 0.5,
                                                        ticks.linewidth = legend.tickwidth,
                                                        frame.linewidth = legend.framewidth,
                                                        frame.colour = legend.framecolor,
                                                        ticks.colour = legend.tickcolor))
    } else if (legend.type == "colorsteps"){
      p <- p +
        ggplot2::guides(color = ggplot2::guide_colorsteps(title = "Avg. Expression",
                                                          title.position = "top",
                                                          barwidth = legend.barwidth,
                                                          barheight = legend.barheight,
                                                          title.hjust = 0.5,
                                                          ticks.linewidth = legend.tickwidth,
                                                          frame.linewidth = legend.framewidth,
                                                          frame.colour = legend.framecolor,
                                                          ticks.colour = legend.tickcolor))
    }

    # Modify size legend.
    p <- p +
         ggplot2::guides(size = ggplot2::guide_legend(title = "Percent Expressed",
                                                      title.position = "top",
                                                      title.hjust = 0.5))

    if (!is.null(xlab)){
      p <- p & ggplot2::xlab(xlab)
    } else {
      p <- p & ggplot2::xlab("")
    }
    if (!is.null(ylab)){
      p <- p & ggplot2::ylab(ylab)
    } else {
      p <- p & ggplot2::ylab("")
    }
    if (!is.null(plot.title)){
      p <- p &
          ggplot2::labs(title = plot.title)
    }


    # Add custom subtitle.
    if (!is.null(plot.subtitle)){
      p <- p +
          ggplot2::labs(subtitle = plot.subtitle)
    }

    # Add custom caption
    if (!is.null(plot.caption)){
      p <- p +
          ggplot2::labs(caption = plot.caption)
      }

    if (!(is.null(rotate_x_labels))){
      if (isTRUE(rotate_x_labels)){
        p <- p & ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 0, hjust = 0.5))
      }
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
