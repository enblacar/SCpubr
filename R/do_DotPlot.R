#' Wrapper for \link[Seurat]{DotPlot}.
#'
#'
#' @param sample Seurat object.
#' @param features Features to represent.
#' @param assay Assay to use. Defaults to the current assay.
#' @param group.by Variable you want the cells to be colored for.
#' @param split.by Split into as many plots as unique values in the variable provided.
#' @param cols Two colors if split.by is not set, which will define a gradient. As many numbers as unique values in split.by, if set, which each own will define its own gradient. Defaults to predefined color scales if not provided.
#' @param legend Whether to plot the legend or not.
#' @param legend.position Position of the legend in the plot. Will only work if legend is set to TRUE.
#' @param plot.title Title to use in the plot.
#' @param xlab Title for the X axis.
#' @param ylab Title for the Y axis.
#' @param axis.text.fonsize Modify the font size for axis texts.
#' @param axis.title.fonsize Modify the font size for axis titles.
#' @param plot.title.fonsize Modify the font size for the plot title.
#' @param legend.text.fontsize Modify the font size for the legend text.
#' @param legend.title.fontsize Modify the font size for the legend title.
#' @param flip Whether to flip the axis.
#' @param dot.scale Scale the size of the dots.
#' @param cluster.idents Logical. Whether to cluster the identities based on the expression of the features.
#'
#' @return A ggplot2 object containing a Dot Plot.
#' @export
#'
#' @examples
#' \dontrun{
#' TBD
#' }
do_DotPlot <- function(sample,
                       features,
                       assay = NULL,
                       group.by = NULL,
                       split.by = NULL,
                       legend = TRUE,
                       dot.scale = 6,
                       cols = c("grey75", "#014f86"),
                       legend.position = "right",
                       plot.title = NULL,
                       xlab = NULL,
                       ylab = NULL,
                       axis.text.fonsize = 16,
                       axis.title.fonsize = 16,
                       plot.title.fonsize = 18,
                       legend.text.fontsize = 12,
                       legend.title.fontsize = 12,
                       cluster.idents = FALSE,
                       flip = FALSE){
    # Checks for packages.
    check_suggests(function_name = "do_DotPlot")
    # Check the assay.
    out <- check_and_set_assay(sample, assay = assay)
    sample <- out[["sample"]]
    assay <- out[["assay"]]

    # Check logical parameters.
    logical_list <- list("legend" = legend,
                         "flip" = flip,
                         "cluster.idents" = cluster.idents)
    check_type(parameters = logical_list, required_type = "logical", test_function = is.logical)
    # Check numeric parameters.
    numeric_list <- list("dot.scale" = dot.scale,
                         "axis.text.fonsize" = axis.text.fonsize,
                         "axis.title.fonsize" = axis.title.fonsize,
                         "plot.title.fonsize" = plot.title.fonsize,
                         "legend.text.fontsize" = legend.text.fontsize,
                         "legend.title.fontsize" = legend.title.fontsize)
    check_type(parameters = numeric_list, required_type = "numeric", test_function = is.numeric)
    # Check character parameters.
    character_list <- list("legend.position" = legend.position,
                           "plot.title" = plot.title,
                           "features" = unlist(features),
                           "xlab" = xlab,
                           "ylab" = ylab,
                           "cols" = cols,
                           "group.by" = group.by,
                           "split.by" = split.by)
    check_type(parameters = character_list, required_type = "character", test_function = is.character)

    # Check the features.
    check_feature(sample = sample, features = features)

    # Check that flip is not set to TRUE and features is not a named list.
    if (isTRUE(flip) & is.list(features)){stop("Please provide the genes as a simple character vector or set flip to FALSE.")}
    # Check that split.by is set and the user has not provided a correct vector of colors.
    if (!(is.null(split.by))){
      if (length(cols) != length(as.character(unique(Seurat::FetchData(sample, vars = split.by)[, 1])))){
        names.use <- if (is.factor(unique(Seurat::FetchData(sample, vars = split.by)[, 1]))){
          levels(unique(Seurat::FetchData(sample, vars = split.by)[, 1]))
        } else {
          as.character(unique(Seurat::FetchData(sample, vars = split.by)[, 1]))
        }
        cols <- generate_color_scale(names.use)
      }
    }

    plot <- Seurat::DotPlot(sample,
                            features = features,
                            cols = cols,
                            group.by = group.by,
                            split.by = split.by,
                            dot.scale = dot.scale,
                            cluster.idents = cluster.idents) +
            ggpubr::theme_pubr(legend = legend.position) +
            ggplot2::theme(axis.text.x = ggplot2::element_text(size = axis.text.fonsize, angle = 90, vjust = 0.5, hjust = 1, face = "bold"),
                           axis.text.y = ggplot2::element_text(size = axis.text.fonsize, face = "bold"),
                           axis.title = ggplot2::element_text(face = "bold", size = axis.title.fonsize),
                           legend.text = ggplot2::element_text(size = legend.text.fontsize, hjust = 1),
                           legend.title = ggplot2::element_text(size = legend.title.fontsize, face = "bold"),
                           plot.title = ggplot2::element_text(size = plot.title.fonsize, face = "bold", hjust = 0.5))

    if (!is.null(xlab)){
      plot <- plot & ggplot2::xlab(xlab)
    } else {
      plot <- plot & ggplot2::xlab("")
    }
    if (!is.null(ylab)){
      plot <- plot & ggplot2::ylab(ylab)
    } else {
      plot <- plot & ggplot2::ylab("")
    }
    if (!is.null(plot.title)){
      plot <- plot + ggplot2::ggtitle(plot.title)
    }

    if (is.list(features)){
      plot <- plot + ggplot2::theme(strip.background = ggplot2::element_rect(color = 'black', fill = 'white'),
                                    strip.text = ggplot2::element_text(face = "bold", size = legend.text.fontsize - 2, color = "black"))
    }
    if (flip == TRUE){
        plot <- plot + ggplot2::coord_flip()
    }
    return(plot)

}
