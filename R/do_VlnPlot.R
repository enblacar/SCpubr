#' Wrapper for \link[Seurat]{VlnPlot}.
#'
#'
#' @param sample  Seurat object.
#' @param features Features to represent.
#' @param group.by  Variable you want the cells to be colored for.
#' @param split.by  Split into as many plots as unique values in the variable provided.
#' @param cols  From \link[Seurat]{DotPlot}: Colors to plot: the name of a palette from `RColorBrewer::brewer.pal.info`, a pair of colors defining a gradient, or 3+ colors defining multiple gradients (if split.by is set).
#' @param cols.group  Vector of colors matching the unique values in group.by.
#' @param cols.split  Vector of colors matching the unique values in split.by.
#' @param plot_boxplot Logical. Whether to plot a boxplot inside the violin or not.
#' @param legend  Whether to plot the legend or not.
#' @param legend.position Position of the legend in the plot. Will only work if legend is set to TRUE.
#' @param plot.title  Title to use in the plot.
#' @param pt.size  Size of points in the VlnPlot.
#' @param xlab  Title for the X axis.
#' @param ylab  Title for the Y axis.
#' @param remove_x_axis  Remove X axis labels and ticks from the plot.
#' @param remove_y_axis  Remove Y axis labels and ticks from the plot.
#' @param axis.text.fonsize  Modify the fontsize for axis texts.
#' @param axis.title.fonsize  Modify the fontsize for axis titles.
#' @param plot.title.fonsize  Modify the fontsize for the plot title.
#' @param legend.text.fontsize  Modify the fontsize for the legend text.
#' @param legend.title.fontsize  Modify the fontsize for the legend title.
#' @param y_cut  Vector with the values in which the Violins should be cut. Only works for one feature.
#' @param legend.ncol  Number of columns in the legend.


#' @return A ggplot2 object containing a Violin Plot.
#' @export
#'
#' @examples
#' \dontrun{
#' TBD
#' }
do_VlnPlot <- function(sample,
                       features,
                       group.by = NULL,
                       split.by = NULL,
                       legend = TRUE,
                       cols = NULL,
                       cols.group = NULL,
                       cols.split = NULL,
                       pt.size = 0,
                       y_cut = NULL,
                       plot_boxplot = TRUE,
                       legend.position = "bottom",
                       plot.title = "",
                       xlab = "",
                       ylab = "",
                       axis.text.fonsize = 28,
                       axis.title.fonsize = 28,
                       plot.title.fonsize = 30,
                       legend.text.fontsize = 20,
                       legend.title.fontsize = 24,
                       remove_x_axis = TRUE,
                       remove_y_axis = FALSE,
                       legend.ncol = 3){
    # Checks for packages.
    check_suggests(function_name = "do_VlnPlot")
    if (!(is.null(y_cut)) & length(features) > 1){
        stop('You can only provide values for y_cut if only one feature is provided to the function.')
    }

    if (is.null(cols)){
        if (is.null(group.by) & is.null(split.by)){
            colors.use <- "steelblue"
        }

        if (!(is.null(group.by)) & is.null(cols.group)){
            colors.use <- "steelblue"
        } else if (!(is.null(group.by)) & !(is.null(cols.group))){
            colors.use <- cols.group
        }

        if (!(is.null(split.by)) & is.null(cols.split)){
            colors.use <- "steelblue"
        } else if (!(is.null(split.by)) & !(is.null(cols.split))){
            colors.use <- cols.split
        }

    }

    plot <- Seurat::VlnPlot(sample,
                            features = features,
                            cols = colors.use,
                            group.by = group.by,
                            split.by = split.by,
                            pt.size = pt.size) &
        ggplot2::ylab(ylab) &
        ggplot2::xlab(xlab) &
        ggpubr::theme_pubr(legend = legend.position) &
        ggplot2::scale_y_continuous(labels = scales::comma) &
        ggplot2::theme(axis.text.x = ggplot2::element_text(size = axis.text.fonsize, angle = 90, vjust = 0.5, hjust = 1, face = "bold"),
                       axis.text.y = ggplot2::element_text(size = axis.text.fonsize, face = "bold"),
                       axis.title = ggplot2::element_text(face = "bold", size = axis.title.fonsize),
                       legend.text = ggplot2::element_text(size = legend.text.fontsize, hjust = 0, face = "bold"),
                       legend.title = ggplot2::element_text(size = legend.title.fontsize, face = "bold"),
                       plot.title = ggplot2::element_text(size = plot.title.fonsize, face = "bold", hjust = 0.5)) &
        ggplot2::guides(fill = ggplot2::guide_legend(ncol = legend.ncol))

    if (plot_boxplot == TRUE){
      plot <- plot &
        ggplot2::geom_boxplot(width = .1, fill = "white", outlier.colour = NA)
    }

    if (remove_x_axis == TRUE){
        plot <- plot & ggpubr::rremove("x.text") + ggpubr::rremove("x.ticks")
    }
    if (remove_y_axis == TRUE){
        plot <- plot & ggpubr::rremove("y.text") + ggpubr::rremove("y.ticks")
    }

    if (legend == FALSE){
        plot <- plot & Seurat::NoLegend()
    }

    if (!(is.null(y_cut))){
        for (value in y_cut){
            plot <- plot +
                ggplot2::geom_hline(yintercept = value, linetype = "dashed", colour = "grey25", size = 1.5)
        }
    }

    return(plot)

}
