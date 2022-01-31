#' Butterfly plot.
#'
#' This plot aims to show the relationships between distinct enrichment scores. If 3 variables are provided, the relationship is between the Y axis and the dual X axis.
#' If 4 variables are provided, each corner of the plot represents how enriched the cells are in that given list. How to interpret this? In a 3-variable plot, the Y axis
#' just means one variable. The higher the cells are in the Y axis the more enriched they are in that given variable. The X axis is a dual parameter one. Cells falling
#' into each extreme of the axis are highly enriched for either x1 or x2, while cells falling in between are not enriched for any of the two. In a 4-variable plot, each corner
#' shows the enrichment for one of the 4 given features. Cells will tend to locate in either of the four corners, but there will be cases of cells locating mid-way between two
#' given corners (enriched in both features) or in the middle of the plot (not enriched for any).
#'
#' This plots are based on the following publications:
#' - Neftel, C. et al. An Integrative Model of Cellular States, Plasticity, and Genetics for Glioblastoma. Cell 178, 835-849.e21 (2019). https://doi.org/10.1016/j.cell.2019.06.024
#' - Tirosh, I., Venteicher, A., Hebert, C. et al. Single-cell RNA-seq supports a developmental hierarchy in human oligodendroglioma. Nature 539, 309–313 (2016). https://doi.org/10.1038/nature20123
#'
#' @param sample  Seurat object.
#' @param x1  First feature in the X axis. Will go on the right side if y2 is not provided and top-right quadrant if provided.
#' @param x2  Second feature on the X axis. Will go on the left side if y2 is not provided and top-left quadrant if provided.
#' @param y1  First feature on the Y axis. Will become the Y axis if y2 is not provided and bottom-right quadrant if provided.
#' @param y2  Second feature on the Y axis. Will become the bottom-left quadrant if provided.
#' @param categorical  Do you want to color the cells using a categorical variable?
#' @param categorical_feature  Categorical metadata variable to color the cells for.
#' @param continuous  Do you want to color the cells using a continuous variable?
#' @param continuous_feature  Continuous metadata variable to color the cells for.
#' @param legend.position  Position of the legend in the plot. One of: top, bottom, left, right.
#' @param plot.title  Title to use in the plot.
#' @param xlab  Title for the X axis. Only works if y2 is not set up.
#' @param ylab  Title for the Y axis. Only works if y2 is not set up.
#' @param axis.ticks  Whether to show axis ticks.
#' @param axis.text  Whether to show axis text.
#' @param complex.output  Logical. Returns a patchwork plot with the RankPlot of variables across X and Y axis.
#' @param complex.output.grouping.variable  Variable to group the values in the RankPlot.
#'
#' @return
#' @export
#'
#' @examples
#' \dontrun{
#' TBD
#' }
do_ButterflyPlot <- function(sample,
                        x1,
                        x2,
                        y1,
                        y2 = NULL,
                        categorical = FALSE,
                        categorical_feature = FALSE,
                        colors.use = NULL,
                        continuous = FALSE,
                        continuous_feature = FALSE,
                        legend.position = NULL,
                        plot.title = "",
                        xlab = NULL,
                        ylab = NULL,
                        axis.ticks = TRUE,
                        axis.text = TRUE,
                        complex.output = FALSE,
                        complex.output.grouping.variable = FALSE
                        ){

    # Check.
    if (categorical == TRUE & continuous == TRUE){
        stop("Either select categorical or continuous coloring.")
    }

    # 3-variable plot.
    if (is.null(y2)){
        # Retrieve metadata variables.
        variables_to_retrieve <- if (categorical == TRUE){
            c(x1, x2, y1, categorical_feature)
        } else if (continuous == TRUE){
            c(x1, x2, y1, continuous_feature)
        } else {
            c(x1, x2, y1)
        }

        # And store them as a tibble.
        scores <- sample@meta.data[, variables_to_retrieve]
        scores$cell <- rownames(scores)
        # Shuffle the cells so that we accomplish a random plotting, not sample by sample.
        scores <- scores[sample(scores$cell, nrow(scores)), ]
        scores <- tidyr::tibble(scores)

        # Compute the scores for the X axis.
        x <- unlist(pbapply::pbsapply(1:nrow(scores), function(x) {
            d <- log2(abs(scores[x, x1] - scores[x, x2]) + 1)
            ifelse(scores[x, x1] > scores[x, x2], d, -d)
        }))
        names(x) <- scores$cell

        # Compute the scores for the Y axis.
        y <- log2(scores %>% dplyr::pull(y1) + 1)
        names(y) <- scores$cell


        # Define titles.
        x_lab <- ifelse(is.null(xlab), paste0(x2, "  <---->  ", x1), xlab)
        y_lab <- ifelse(is.null(ylab), y1, ylab)

        # Add the variables as metadata in case they are needed for complex output.
        sample$x_axis <- x
        sample$y_axis <- y

        # Plain plot without coloring.
        if (categorical == FALSE & continuous == FALSE) {
            df <- data.frame("set_x" = x, "set_y" = y)
            plot <- ggplot2::ggplot(df, mapping = ggplot2::aes(x = set_x, y = set_y)) +
                        ggplot2::geom_point() +
                        ggpubr::theme_pubr(legend = "bottom") +
                        ggpubr::rremove("legend.title")

        # Color based on categorical variable.
        } else if (categorical == TRUE){
            df <- data.frame("set_x" = x, "set_y" = y, "color" = scores[, categorical_feature])
            colors.use <- colors.use[names(colors.use) %in% unique(df[, categorical_feature])]
            plot <- ggplot2::ggplot(df, mapping = ggplot2::aes(x = set_x, y = set_y, color = !!(sym(categorical_feature)))) +
                ggplot2::geom_point() +
                ggpubr::theme_pubr(legend = ifelse(is.null(legend.position), "bottom", legend.position)) +
                ggpubr::rremove("legend.title") +
                ggplot2::scale_color_manual(values = colors.use)

        # Color based on a continuous variable.
        } else if (continuous == TRUE){
            df <- data.frame("set_x" = x, "set_y" = y, "color" = scores[, continuous_feature])
            plot <- ggplot2::ggplot(df, mapping = ggplot2::aes(x = set_x, y = set_y, color = !!(sym(continuous_feature)))) +
                ggplot2::geom_point() +
                ggpubr::theme_pubr(legend = ifelse(is.null(legend.position), "right", legend.position)) +
                viridis::scale_color_viridis(name = continuous_feature)
        }

        # Add the labels.
        plot <- plot +
                    ggplot2::xlab(x_lab) +
                    ggplot2::ylab(y_lab) +
                    ggplot2::ggtitle(plot.title)

    # 4-parameter plot.
    } else if (!is.null(y2)){
        # Retrieve metadata variables to plot.
        variables_to_retrieve <- if (categorical == TRUE){
            c(x1, x2, y1, y2, categorical_feature)
        } else if (continuous_feature == TRUE){
            c(x1, x2, y1, y2, continuous_feature)
        } else {
            c(x1, x2, y1, y2)
        }

        # And store them as a tibble.
        scores <- sample@meta.data[, variables_to_retrieve]
        scores$cell <- rownames(scores)
        # Shuffle the cells so that we accomplish a random plotting, not sample by sample.
        scores <- scores[sample(scores$cell, nrow(scores)), ]
        scores <- tidyr::tibble(scores)

        # Compute Y axis values.
        d <- apply(scores, 1, function(x){as.double(max(x[c(x1, x2)])) - as.double(max(x[c(y1, y2)]))})
        names(d) <- scores$cell

        # Compute X axis values.
        x <- unlist(pbapply::pbsapply(1:length(d), function(x) {
            if (d[x] > 0) {
                d <- log2(abs(scores[x, x1] - scores[x, x2]) + 1)
                ifelse(scores[x, x1] > scores[x, x2], d, -d)
            } else {
                d <- log2(abs(scores[x, y1] - scores[x, y2]) + 1)
                ifelse(scores[x, y1] > scores[x, y2], d, -d)
            }
        }))
        names(x) <- scores$cell

        # Define titles for the axis.
        x_lab1 <- paste0(y2, "  <---->  ", y1)
        x_lab2 <- paste0(x2, "  <---->  ", x1)
        y_lab1 <- paste0(y2, "  <---->  ", x2)
        y_lab2 <- paste0(x1, "  <---->  ", y1)

        # Add the axis value as metadata for complex output if needed.
        sample$x_axis <- x
        sample$y_axis <- d

        # Plain plot without coloring.
        if (categorical == FALSE & continuous == FALSE) {
            df <- data.frame("set_x" = x, "set_y" = d)
            plot <- ggplot2::ggplot(df, mapping = ggplot2::aes(x = set_x, y = set_y)) +
                ggplot2::geom_point() +
                ggpubr::theme_pubr(legend = "bottom") +
                ggpubr::rremove("legend.title")

        # Color based on a categorical variable.
        } else if (categorical == TRUE){
            df <- data.frame("set_x" = x, "set_y" = d, "color" = scores[, categorical_feature])
            colors.use <- colors.use[names(colors.use) %in% unique(df[, categorical_feature])]
            plot <- ggplot2::ggplot(df, mapping = ggplot2::aes(x = set_x, y = set_y, color = !!(sym(categorical_feature)))) +
                ggplot2::geom_point() +
                ggpubr::theme_pubr(legend = ifelse(is.null(legend.position), "bottom", legend.position)) +
                ggpubr::rremove("legend.title") +
                ggplot2::scale_color_manual(values = colors.use)

        # Color based on a continuous variable.
        } else if (continuous == TRUE){
            df <- data.frame("set_x" = x, "set_y" = d, "color" = scores[, continuous_feature])
            plot <- ggplot2::ggplot(df, mapping = ggplot2::aes(x = set_x, y = set_y, color = !!(sym(continuous_feature)))) +
                ggplot2::geom_point() +
                ggpubr::theme_pubr(legend = ifelse(is.null(legend.position), "right", legend.position)) +
                viridis::scale_color_viridis(name = continuous_feature)
        }

        # Add the extra axis and axis titles.
        plot <- plot +
            ggplot2::scale_y_continuous(sec.axis = sec_axis(~.*1, name = y_lab2)) +
            ggplot2::scale_x_continuous(sec.axis = sec_axis(~.*1, name = x_lab2)) +
            ggplot2::xlab(x_lab1) +
            ggplot2::ylab(y_lab1) +
            ggplot2::ggtitle(plot.title)
    }

    # Overall formatting for the plot.
    plot <- plot +
        ggplot2::theme(axis.text = ggplot2::element_text(face = "bold"),
                       axis.title = ggplot2::element_text(face = "bold"),
                       plot.title = ggplot2::element_text(face = "bold",
                                                          hjust = 0.5),
                       legend.text = ggplot2::element_text(size = 10, face = "bold", hjust = 1))

    # Remove axis ticks?
    if (axis.ticks == FALSE){
        plot <- plot + ggpubr::rremove("ticks")
    }

    # Remove axis text?
    if (axis.text == FALSE){
        plot <- plot + ggpubr::rremove("axis.text")
    }

    # Complex output?
    if (complex.output == TRUE){
        if (categorical == FALSE & continuous == FALSE) {
            message("No Complex output can be returned as no categorical or continuous variables are selected.")
        } else if (categorical == TRUE){
            p.x <- do_RankPlot(sample, feature_to_rank = "x_axis", group.by = complex.output.grouping.variable, colors.use = colors.use)
            p.y <- do_RankPlot(sample, feature_to_rank = "y_axis", group.by = complex.output.grouping.variable, colors.use = colors.use)
            plot <- plot | p.x | p.y

        # Color based on a continuous variable.
        } else if (continuous == TRUE){
            p.x <- do_RankPlot(sample, feature_to_rank = "x_axis", group.by = !!(sym(complex.output.grouping.variable)), continuous_feature = TRUE)
            p.y <- do_RankPlot(sample, feature_to_rank = "y_axis", group.by = !!(sym(complex.output.grouping.variable)), continuous_feature = TRUE)
            plot <- plot | p.x | p.y
        }
    }

    return(plot)

}
