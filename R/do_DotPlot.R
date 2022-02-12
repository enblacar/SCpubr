#' Wrapper for \link[Seurat]{DotPlot}.
#'
#'
#' @param sample Seurat object.
#' @param features Features to represent.
#' @param assay Assay to use. Defaults to the current assay.
#' @param group.by Variable you want the cells to be colored for.
#' @param split.by Split into as many plots as unique values in the variable provided.
#' @param cols From \link[Seurat]{DotPlot}: Colors to plot: the name of a palette from `RColorBrewer::brewer.pal.info`, a pair of colors defining a gradient, or 3+ colors defining multiple gradients (if split.by is set).
#' @param legend Whether to plot the legend or not.
#' @param legend.position Position of the legend in the plot. Will only work if legend is set to TRUE.
#' @param plot.title Title to use in the plot.
#' @param xlab Title for the X axis.
#' @param ylab Title for the Y axis.
#' @param remove_x_axis Remove X axis labels and ticks from the plot.
#' @param remove_y_axis Remove Y axis labels and ticks from the plot.
#' @param axis.text.fonsize Modify the font size for axis texts.
#' @param axis.title.fonsize Modify the font size for axis titles.
#' @param plot.title.fonsize Modify the font size for the plot title.
#' @param legend.text.fontsize Modify the font size for the legend text.
#' @param legend.title.fontsize Modify the font size for the legend title.
#' @param flip Whether to flip the axis.
#' @param dot.scale Scale the size of the dots.
#' @param cluster.idents Logical. Whether to cluster the identities based on the expression of the features.
#' @param ... Other arguments provided to \link[Seurat]{DotPlot}.
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
                       plot.title = "",
                       xlab = "",
                       ylab = "",
                       axis.text.fonsize = 16,
                       axis.title.fonsize = 16,
                       plot.title.fonsize = 18,
                       legend.text.fontsize = 16,
                       legend.title.fontsize = 16,
                       remove_x_axis = FALSE,
                       remove_y_axis = FALSE,
                       cluster.idents = TRUE,
                       flip = FALSE,
                       ...){
    # Checks for packages.
    check_suggests(function_name = "do_DotPlot")
    # Check the assay.
    out <- check_and_set_assay(sample, assay = assay)
    sample <- out[["sample"]]
    assay <- out[["assay"]]

    # Check logical parameters.
    logical_list <- list("legend" = legend,
                         "remove_x_axis" = remove_x_axis,
                         "remove_y_axis" = remove_y_axis,
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
                           "features" = features,
                           "xlab" = xlab,
                           "ylab" = ylab,
                           "cols" = cols)
    check_type(parameters = character_list, required_type = "character", test_function = is.character)
    if(!(is.null(group.by))){check_type(parameters = list("group.by" = group.by), required_type = "chracter", test_function = is.character)}
    if(!(is.null(split.by))){check_type(parameters = list("split.by" = split.by), required_type = "chracter", test_function = is.character)}

    # Check the features.
    check_feature(sample = sample, features = features)

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
                            ...) +
            ggplot2::ylab(ylab) +
            ggplot2::xlab(xlab) +
            ggpubr::theme_pubr(legend = legend.position) +
            ggplot2::theme(axis.text.x = ggplot2::element_text(size = axis.text.fonsize, angle = 90, vjust = 0.5, hjust = 1, face = "bold"),
                           axis.text.y = ggplot2::element_text(size = axis.text.fonsize, face = "bold"),
                           axis.title = ggplot2::element_text(face = "bold", size = axis.title.fonsize),
                           legend.text = ggplot2::element_text(size = legend.text.fontsize, hjust = 1),
                           legend.title = ggplot2::element_text(size = legend.title.fontsize, face = "bold"),
                           plot.title = ggplot2::element_text(size = plot.title.fonsize, face = "bold", hjust = 0.5))
    if (remove_x_axis == TRUE){
        plot <- plot + ggpubr::rremove("x.text") + ggpubr::rremove("x.ticks")
    }
    if (remove_y_axis == TRUE){
        plot <- plot + ggpubr::rremove("y.text") + ggpubr::rremove("y.ticks")
    }
    if (flip == TRUE){
        plot <- plot + ggplot2::coord_flip()
    }
    return(plot)

}
