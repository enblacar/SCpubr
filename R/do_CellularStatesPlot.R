#' Cellular States plot.
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
#' @param gene_list Named list of lists of marker genes to query for enrichment.
#' @param x1  A name of a list from gene_list. First feature in the X axis. Will go on the right side of the X axis if y2 is not provided and top-right quadrant if provided.
#' @param x2  A name of a list from gene_list. Second feature on the X axis. Will go on the left side of the X axis if y2 is not provided and top-left quadrant if provided.
#' @param y1  A name of a list from gene_list. First feature on the Y axis. Will become the Y axis if y2 is not provided and bottom-right quadrant if provided.
#' @param y2  A name of a list from gene_list. Second feature on the Y axis. Will become the bottom-left quadrant if provided.
#' @param colors.use Named vector with the names of the unique values in the categorical variable and values the HEX codes for the colors.
#' @param legend.position  Position of the legend in the plot. One of: top, bottom, left, right.
#' @param group.by Metadata variable to color the cells by. Defaults to current identities.
#' @param plot.title,plot.subtitle,plot.caption  Title to use in the plot.
#' @param xlab  Title for the X axis. Only works if y2 is not set up.
#' @param ylab  Title for the Y axis. Only works if y2 is not set up.
#' @param axis.ticks  Whether to show axis ticks.
#' @param axis.text  Whether to show axis text.
#' @param verbose Verbose function?
#'
#' @return  A ggplot2 object containing a butterfly plot.
#' @export
#' @examples
#' \dontrun{
#' TBD
#' }
do_CellularStatesPlot <- function(sample,
                                  gene_list,
                                  x1,
                                  x2 = NULL,
                                  y1,
                                  y2 = NULL,
                                  group.by = NULL,
                                  colors.use = NULL,
                                  legend.position = NULL,
                                  plot.title = "",
                                  plot.subtitle = NULL,
                                  plot.caption = NULL,
                                  xlab = NULL,
                                  ylab = NULL,
                                  axis.ticks = TRUE,
                                  axis.text = TRUE,
                                  verbose = FALSE){
    # Checks for packages.
    check_suggests(function_name = "do_CellularStatesPlot")
    # Check if the sample provided is a Seurat object.
    check_Seurat(sample = sample)
    # Define pipe operator internally.
    `%>%` <- purrr::`%>%`

    # Check the colors provided.
    if (is.null(colors.use)){
      colors.use <- {
        if (is.null(group.by)){
          generate_color_scale(levels(sample))
        } else if (!(is.null(group.by))){
          data.use <- sample[[]][, group.by, drop = F]
          names.use <- if (is.factor(data.use[, 1])){levels(data.use[, 1])} else {sort(unique(data.use[, 1]))}
          generate_color_scale(names.use)
        }
      }
    } else {
      check_colors(colors.use, parameter_name = "colors.use")
      if (is.null(group.by)){
        colors.use <- check_consistency_colors_and_names(sample = sample, colors = colors.use)
      } else {
        colors.use <- check_consistency_colors_and_names(sample = sample, colors = colors.use, grouping_variable = group.by)
      }
    }

    # Fix for group.by.
    if (is.null(group.by)){
      group.by <- "dummy"
      sample@meta.data$dummy <- sample@active.ident
    } else {
      sample@meta.data$dummy <- sample@meta.data[, group.by]
      group.by <- "dummy"
    }

    # Compute the enrichment scores.
    sample <- compute_enrichment_scores(sample = sample, list_genes = gene_list, verbose = verbose)

    # 2-variable plot.
    if (is.null(y2) & is.null(x2)){
      # Check that the names provided are not repeated.
      if (sum(duplicated(c(x1, y1))) > 0){
        stop("The names of the lists to plot can not be the same.", call. = FALSE)
      }
      # Check that the names provided match the marker genes.
      if (!(x1 %in% names(gene_list))){
        stop(paste0(x1, " is not a name of a list of genes provided to gene_list.", call. = FALSE), call. = F)
      }
      if (!(y1 %in% names(gene_list))){
        stop(paste0(y1, " is not a name of a list of genes provided to gene_list.", call. = FALSE), call. = F)
      }
      # Retrieve metadata variables.
      variables_to_retrieve <- c(x1, y1, group.by)
      # And store them as a tibble.
      scores <- sample@meta.data[, variables_to_retrieve]
      scores$cell <- rownames(scores)
      # Shuffle the cells so that we accomplish a random plotting, not sample by sample.
      scores <- scores[sample(scores$cell, nrow(scores)), ]
      scores <- tidyr::tibble(scores)

      # Compute scores for the X axis.
      x <- scores %>% dplyr::pull(x1)
      # Compute scores for the Y axis.
      y <- scores %>% dplyr::pull(y1)

      names(x) <- scores$cell
      names(y) <- scores$cell

      # Define titles.
      x_lab <- ifelse(is.null(xlab), x1, xlab)
      y_lab <- ifelse(is.null(ylab), y1, ylab)

      # Define limits of polots.
      lim1 <- min(min(x), min(y))
      lim2 <- max(max(x), max(y))
      lim_x <- c(lim1, lim2)
      lim_y <- c(lim1, lim2)

      # Plot
      df <- data.frame("set_x" = x, "set_y" = y, "group.by" = scores$dummy)
      plot <- ggplot2::ggplot(df, mapping = ggplot2::aes(x = .data$set_x, y = .data$set_y, color = .data$group.by)) +
              ggplot2::geom_point() +
              ggpubr::theme_pubr(legend = "bottom") +
              ggpubr::rremove("legend.title") +
              ggplot2::scale_color_manual(values = colors.use) +
              ggplot2::xlab(x_lab) +
              ggplot2::ylab(y_lab) +
              ggplot2::ggtitle(plot.title) +
              ggplot2::coord_fixed(xlim = lim_x, ylim = lim_y)

    # 3-variable plot.
    } else if (is.null(y2) & !(is.null(x2))){
        # Check that the names provided are not repeated.
        if (sum(duplicated(c(x1, y1, x2))) > 0){
          stop("The names of the lists to plot can not be the same.", call. = FALSE)
        }
        # Check that the names provided match the marker genes.
        if (!(x1 %in% names(gene_list))){
          stop(paste0(x1, " is not a name of a list of genes provided to gene_list.", call. = FALSE))
        }
        if (!(x2 %in% names(gene_list))){
          stop(paste0(x2, " is not a name of a list of genes provided to gene_list.", call. = FALSE))
        }
        if (!(y1 %in% names(gene_list))){
          stop(paste0(y1, " is not a name of a list of genes provided to gene_list.", call. = FALSE))
        }
        # Retrieve metadata variables.
        variables_to_retrieve <- c(x1, x2, y1, group.by)
        # And store them as a tibble.
        scores <- sample@meta.data[, variables_to_retrieve]
        scores$cell <- rownames(scores)
        # Shuffle the cells so that we accomplish a random plotting, not sample by sample.
        scores <- tidyr::tibble(scores)

        # Compute the scores for the X axis.
        x <- unlist(sapply(1:nrow(scores), function(x) {
          score_1 <- scores[x, x1] + stats::runif(1, min=0, max=0.15)
          score_2 <- scores[x, x2] + stats::runif(1, min=0, max=0.15)
          d <- max(score_1, score_2)
          ifelse(score_1 > score_2, d, -d)
        }))

        # Compute the scores for the Y axis.
        y <- unlist(sapply(1:nrow(scores), function(x) {
          score_1 <- scores[x, x1] + stats::runif(1, min=0, max=0.15)
          score_2 <- scores[x, x2] + stats::runif(1, min=0, max=0.15)
          d <- max(score_1, score_2)
          y <- scores[x, y1] - d
          y
        }))

        names(x) <- scores$cell
        names(y) <- scores$cell

        # Define titles.
        x_lab <- ifelse(is.null(xlab), paste0(x2, "  <---->  ", x1), xlab)
        y_lab <- ifelse(is.null(ylab), y1, ylab)

        # Define limits of polots.
        lim <- max(abs(x))
        lim_x <- c(-lim, lim)
        lim_y <- NULL

        # Plot.
        df <- data.frame("set_x" = x, "set_y" = y, "group.by" = scores$dummy)
        plot <- ggplot2::ggplot(df, mapping = ggplot2::aes(x = .data$set_x, y = .data$set_y, color = .data$group.by)) +
                ggplot2::geom_point() +
                ggpubr::theme_pubr(legend = "bottom") +
                ggpubr::rremove("legend.title")  +
                ggplot2::scale_color_manual(values = colors.use) +
                ggplot2::xlab(x_lab) +
                ggplot2::ylab(y_lab) +
                ggplot2::ggtitle(plot.title) +
                ggplot2::xlim(lim_x)

    # 4-parameter plot.
    } else if (!is.null(y2) & !(is.null(x2))){
        # Check that the names provided are not repeated.
        if (sum(duplicated(c(x1, y1, x2, y2))) > 0){
          stop("The names of the lists to plot can not be the same.", call. = FALSE)
        }
        # Check that the names provided match the marker genes.
        if (!(x1 %in% names(gene_list))){
          stop(paste0(x1, " is not a name of a list of genes provided to gene_list.", call. = FALSE))
        }
        if (!(x2 %in% names(gene_list))){
          stop(paste0(x2, " is not a name of a list of genes provided to gene_list.", call. = FALSE))
        }
        if (!(y1 %in% names(gene_list))){
          stop(paste0(y1, " is not a name of a list of genes provided to gene_list.", call. = FALSE))
        }
        if (!(y2 %in% names(gene_list))){
          stop(paste0(y2, " is not a name of a list of genes provided to gene_list.", call. = FALSE))
        }
        # Retrieve metadata variables to plot.
        variables_to_retrieve <- c(x1, x2, y1, y2)
        # And store them as a tibble.
        scores <- sample@meta.data[, variables_to_retrieve]
        # Shuffle the cells so that we accomplish a random plotting, not sample by sample.

        # Compute Y axis values.
        d <- apply(scores, 1, function(x){max(x[c(x1, x2)]) - max(x[c(y1, y2)])})

        # Compute X axis values.
        x <- sapply(1:length(d), function(x) {
          if (d[x] > 0) {
            d <- log2(abs(scores[x, x1] - scores[x, x2]) + 1)
            ifelse(scores[x, x1] < scores[x, x2], d, -d)
          } else {
            d <- log2(abs(scores[x, y1] - scores[x, y2]) + 1)
            ifelse(scores[x, y1] < scores[x, y2], d, -d)
          }
        })

        names(x) <- rownames(scores)

        # Define titles for the axis.
        x_lab1 <- paste0(y1, "  <---->  ", y2)
        x_lab2 <- paste0(x1, "  <---->  ", x2)
        y_lab1 <- paste0(y1, "  <---->  ", x1)
        y_lab2 <- paste0(x2, "  <---->  ", y2)


        # Plot.
        df <- data.frame(row.names = rownames(scores))
        df$set_x <- x
        df$set_y <- d
        df$group.by <- sample@meta.data[, group.by]
        plot <- ggplot2::ggplot(df, mapping = ggplot2::aes(x = .data$set_x, y = .data$set_y, color = .data$group.by)) +
                ggplot2::geom_point() +
                ggpubr::theme_pubr(legend = "bottom") +
                ggpubr::rremove("legend.title")  +
                ggplot2::scale_color_manual(values = colors.use)  +
                ggplot2::xlab(x_lab1) +
                ggplot2::ylab(y_lab1) +
                ggplot2::ggtitle(plot.title)
        suppressMessages({
          plot <- plot +
                  ggplot2::scale_y_continuous(sec.axis = ggplot2::sec_axis(~., name = y_lab2)) +
                  ggplot2::scale_x_continuous(sec.axis = ggplot2::sec_axis(~., name = x_lab2))
        })


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

    return(plot)

}
