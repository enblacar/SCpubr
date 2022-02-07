#' RankPlot plot.
#'
#'
#' @param sample Seurat object.
#' @param assay Seurat assay to choose.
#' @param feature_to_rank Features for which the cells are going to be ranked.
#' @param continuous_feature Is the feature to rank and color for continuous? I.e: an enrichment score.
#' @param group.by Variable you want the cells to be grouped for.
#' @param colors.use Named vector with the color assignment.
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
#'
#' @return  A ggplot2 object containing a Bee Swarm plot.
#' @export
#'
#' @examples
#' \dontrun{
#' TBD
#' }
do_RankPlot <- function(sample,
                       feature_to_rank,
                       group.by,
                       assay = "SCT",
                       continuous_feature = FALSE,
                       colors.use = NULL,
                       legend.position = "right",
                       plot.title = "",
                       xlab = NULL,
                       ylab = "",
                       axis.text.fonsize = 15,
                       axis.title.fonsize = 15,
                       plot.title.fonsize = 15,
                       legend.text.fontsize = 10,
                       legend.title.fontsize = 12,
                       remove_x_axis = FALSE,
                       remove_y_axis = FALSE,
                       flip = FALSE){

    if (!(feature_to_rank %in% colnames(sample@meta.data)) & (feature_to_rank %in% rownames(sample))){
        message("Feature provided is not in the metadata columns but found as a gene name.")
        sample$rank_me <- sample@assays[[assay]]@data[feature_to_rank, ]
        sample$rank <- rank(sample$rank_me)
    } else if (feature_to_rank %in% colnames(sample@meta.data)) {
        sample$rank_me <- sample@meta.data[, feature_to_rank]
        sample$rank <- rank(sample$rank_me)
    } else {
        stop(paste0("The following feature was not found as either metadata variable nor as a gene name: ", feature_to_rank))
    }
    # Compute the ranking
    sample$ranked_groups <- factor(sample@meta.data[, group.by], levels = sort(unique(sample@meta.data[, group.by])))

    color_by <- ifelse(continuous_feature == T, "rank_me", "ranked_groups")

    plot <- ggplot2::ggplot(sample@meta.data, mapping = ggplot2::aes(x = rank, y = .data$ranked_groups, color = !!rlang::sym(color_by))) +
        ggbeeswarm::geom_quasirandom(groupOnX = FALSE) +
        ggpubr::theme_pubr(legend = legend.position) +
        ggplot2::ggtitle(plot.title) +
        ggpubr::rremove("legend.title") +
        ggplot2::theme(axis.text = ggplot2::element_text(size = axis.text.fonsize,
                                                         face = "bold"),
                       axis.title = ggplot2::element_text(size = axis.title.fonsize,
                                                          face = "bold"),
                       plot.title = ggplot2::element_text(size = plot.title.fonsize,
                                                          face = "bold",
                                                          hjust = 0.5),
                       legend.text = ggplot2::element_text(size = 10, face = "bold", hjust = 1))

    if (continuous_feature == TRUE){
        plot <- plot + viridis::scale_color_viridis()
    } else {
        if (is.null(colors.use)){
            colors.use <- colortools::wheel("#457b9d", length(levels(sample)))
            names(colors.use) <- unique(levels(sample))
        }

        plot <- plot + ggplot2::scale_color_manual(values = colors.use)
    }

    if (continuous_feature == FALSE){
        plot <- plot + ggpubr::rremove("legend")
    }
    if (remove_x_axis == TRUE){
        plot <- plot + ggpubr::rremove("x.text") + ggpubr::rremove("x.ticks")
    }
    if (remove_y_axis == TRUE){
        plot <- plot + ggpubr::rremove("y.text") + ggpubr::rremove("y.ticks")
    }
    if (flip == TRUE){
        plot <- plot + ggplot2::coord_flip() + ggpubr::rremove("y.ticks") +
            ggpubr::rremove("y.text") +
            ggplot2::xlab(ifelse(is.null(ylab), paste0("Ordering of cells across ", feature_to_rank), ylab)) +
            ggplot2::ylab(xlab)

    } else {
        plot <- plot + ggpubr::rremove("x.ticks") +
            ggpubr::rremove("x.text") +
            ggplot2::xlab(ifelse(is.null(xlab), paste0("Ordering of cells across ", feature_to_rank), xlab)) +
            ggplot2::ylab(ylab)

    }
    return(plot)

}
