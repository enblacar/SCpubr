#' BeeSwarm plot.
#'
#'
#' @param sample Seurat object.
#' @param assay Assay to use. Defauls to the current assay.
#' @param reduction Reduction to use. Can be the canonical ones such as "umap", "pca", or any custom ones, such as "diffusion". If you are unsure about which reductions you have, use `Seurat::Reductions(sample)`. Defaults to "umap" if present or to the last computed reduction if the argument is not provided.
#' @param slot Data slot to use. Character. Only one of: counts, data, scale.data. Defaults to "data".
#' @param feature_to_rank Features for which the cells are going to be ranked. Ideal case is that this feature is stored as a metadata column.
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
do_BeeSwarmPlot <- function(sample,
                       feature_to_rank,
                       group.by,
                       assay = NULL,
                       reduction = NULL,
                       slot = NULL,
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
    # Checks for packages.
    check_suggests(function_name = "do_BeeSwarmPlot")
    # Check the assay.
    out <- check_and_set_assay(sample, assay = assay)
    sample <- out[["sample"]]
    assay <- out[["assay"]]
    # Check the reduction.
    reduction <- check_and_set_reduction(sample = sample, reduction = reduction)
    # Check logical parameters.
    logical_list <- list("continuous_feature" = continuous_feature,
                         "remove_x_axis" = remove_x_axis,
                         "remove_y_axis" = remove_y_axis,
                         "flip" = flip)
    check_type(parameters = logical_list, required_type = "logical", test_function = is.logical)
    # Check numeric parameters.
    numeric_list <- list("axis.text.fonsize" = axis.text.fonsize,
                         "axis.title.fonsize" = axis.title.fonsize,
                         "plot.title.fonsize" = plot.title.fonsize,
                         "legend.text.fontsize" = legend.text.fontsize,
                         "legend.title.fontsize" = legend.title.fontsize)
    check_type(parameters = numeric_list, required_type = "numeric", test_function = is.numeric)
    # Check character parameters.
    character_list <- list("legend.position" = legend.position,
                           "plot.title" = plot.title,
                           "feature_to_rank" = feature_to_rank,
                           "group.by" = group.by,
                           "ylab" = ylab,
                           "xlab" = xlab,
                           "slot" = slot)
    check_type(parameters = character_list, required_type = "character", test_function = is.character)
    # Check slot.
    slot <- check_and_set_slot(slot = slot)

    dim_colnames <- check_feature(sample = sample, features = feature_to_rank, dump_reduction_names = TRUE)
    if (feature_to_rank %in% colnames(sample@meta.data)) {
      sample$rank_me <- sample@meta.data[, feature_to_rank]
      sample$rank <- rank(sample$rank_me)
    } else if (feature_to_rank %in% rownames(sample)){
      sample$rank_me <- Seurat::GetAssayData(object = sample, slot = slot)[feature_to_rank, ]
      sample$rank <- rank(sample$rank_me)
    } else if (feature_to_rank %in% dim_colnames){
      for(red in Seurat::Reductions(object = sample)){
        if (feature_to_rank %in% colnames(sample@reductions[[red]][[]])){
          reduction <- red
          sample$rank_me <- sample@reductions[[reduction]][[]][, feature_to_rank]
          sample$rank <- rank(sample$rank_me)
        }
      }
    } else {
        stop(paste0("The following feature was not found: ", feature_to_rank, "\n Please check whether it is a metadata variable, a gene name, or the name of a dimension reduction component (and specified it correctly in the reduction parameter)."))
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
    } else if (continuous_feature == FALSE) {
        if (is.null(colors.use)){colors.use <- generate_color_scale(levels(sample))}
        plot <- plot +
                  ggplot2::scale_color_manual(values = colors.use) +
                  ggpubr::rremove("legend")
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