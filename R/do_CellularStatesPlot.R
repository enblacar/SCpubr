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
#' - Neftel, C. \emph{et al}. An Integrative Model of Cellular States, Plasticity, and Genetics for Glioblastoma. Cell 178, 835-849.e21 (2019). \doi{10.1016/j.cell.2019.06.024}
#' - Tirosh, I., Venteicher, A., Hebert, C. \emph{et al}. Single-cell RNA-seq supports a developmental hierarchy in human oligodendroglioma. Nature 539, 309–313 (2016). \doi{10.1038/nature20123}
#' @inheritParams doc_function
#' @param x1  \strong{\code{\link[base]{character}}} | A name of a list from input_gene_list. First feature in the X axis. Will go on the right side of the X axis if y2 is not provided and top-right quadrant if provided.
#' @param x2  \strong{\code{\link[base]{character}}} | A name of a list from input_gene_list. Second feature on the X axis. Will go on the left side of the X axis if y2 is not provided and top-left quadrant if provided.
#' @param y1  \strong{\code{\link[base]{character}}} | A name of a list from input_gene_list. First feature on the Y axis. Will become the Y axis if y2 is not provided and bottom-right quadrant if provided.
#' @param y2  \strong{\code{\link[base]{character}}} | A name of a list from input_gene_list. Second feature on the Y axis. Will become the bottom-left quadrant if provided.
#' @param axis.ticks  \strong{\code{\link[base]{logical}}} | Whether to show axis ticks.
#' @param axis.text  \strong{\code{\link[base]{logical}}} | Whether to show axis text.
#' @param enforce_symmetry \strong{\code{\link[base]{logical}}} | Whether to enforce the plot to follow a symmetry (3 variables, the X axis has 0 as center, 4 variables, all axis have the same range and the plot is squared).
#' @param plot_features \strong{\code{\link[base]{logical}}} | Whether to also report any other feature onto the primary plot.
#' @param features \strong{\code{\link[base]{character}}} | Additional features to plot.
#' @param plot_enrichment_scores \strong{\code{\link[base]{logical}}} | Whether to report enrichment scores for the input lists as plots.
#'
#' @return  A ggplot2 object containing a butterfly plot.
#' @export
#' @example man/examples/examples_do_CellularStatesPlot.R

do_CellularStatesPlot <- function(sample,
                                  input_gene_list,
                                  x1,
                                  y1,
                                  x2 = NULL,
                                  y2 = NULL,
                                  group.by = NULL,
                                  colors.use = NULL,
                                  legend.position = "bottom",
                                  legend.icon.size = 4,
                                  legend.ncol = NULL,
                                  legend.nrow = NULL,
                                  legend.byrow = FALSE,
                                  plot.title = NULL,
                                  plot.subtitle = NULL,
                                  plot.caption = NULL,
                                  font.size = 14,
                                  font.type = "sans",
                                  xlab = NULL,
                                  ylab = NULL,
                                  axis.ticks = TRUE,
                                  axis.text = TRUE,
                                  verbose = FALSE,
                                  enforce_symmetry = FALSE,
                                  plot_marginal_distributions = FALSE,
                                  marginal.type = "density",
                                  marginal.size = 5,
                                  marginal.group = TRUE,
                                  plot_cell_borders = TRUE,
                                  plot_enrichment_scores = FALSE,
                                  border.size = 2,
                                  border.color = "black",
                                  pt.size = 2,
                                  raster = FALSE,
                                  raster.dpi = 1024,
                                  plot_features = FALSE,
                                  features = NULL,
                                  viridis_color_map = "G",
                                  viridis_direction = 1,
                                  nbin = 24,
                                  ctrl = 100){
    check_suggests(function_name = "do_CellularStatesPlot")
    # Check if the sample provided is a Seurat object.
    check_Seurat(sample = sample)

    # Check logical parameters.
    logical_list <- list("axis.ticks" = axis.ticks,
                         "axis.text" = axis.text,
                         "verbose" = verbose,
                         "enforce_symmetry" = enforce_symmetry,
                         "plot_marginal_distributions" = plot_marginal_distributions,
                         "marginal.group" = marginal.group,
                         "legend.byrow" = legend.byrow,
                         "plot_cell_borders" = plot_cell_borders,
                         "raster" = raster,
                         "plot_features" = plot_features,
                         "plot_enrichment_scores" = plot_enrichment_scores)
    check_type(parameters = logical_list, required_type = "logical", test_function = is.logical)
    # Check numeric parameters.
    numeric_list <- list("font.size" = font.size,
                         "marginal.size" = marginal.size,
                         "legend.icon.size" = legend.icon.size,
                         "legend.ncol" = legend.ncol,
                         "legend.nrow" = legend.nrow,
                         "pt.size" = pt.size,
                         "border.size" = border.size,
                         "raster.dpi" = raster.dpi,
                         "viridis_direction" = viridis_direction,
                         "nbin" = nbin,
                         "ctrl" = ctrl)
    check_type(parameters = numeric_list, required_type = "numeric", test_function = is.numeric)
    # Check character parameters.
    character_list <- list("input_gene_list" = input_gene_list,
                           "x1" = x1,
                           "x2" = x2,
                           "y1" = y1,
                           "y2" = y2,
                           "group.by" = group.by,
                           "ylab" = ylab,
                           "xlab" = xlab,
                           "legend.position" = legend.position,
                           "plot.title" = plot.title,
                           "plot.subtitle" = plot.subtitle,
                           "plot.caption" = plot.caption,
                           "font.type" = font.type,
                           "marginal.type" = marginal.type,
                           "border.color" = border.color,
                           "features" = features,
                           "viridis_color_map" = viridis_color_map)
    check_type(parameters = character_list, required_type = "character", test_function = is.character)

    # Define pipe operator internally.
    `%>%` <- magrittr::`%>%`

    # Check the colors provided.
    if (is.null(colors.use)){
      colors.use <- {
        if (is.null(group.by)){
          generate_color_scale(levels(sample))
        } else if (!(is.null(group.by))){
          data.use <- sample[[]][, group.by, drop = FALSE]
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
    # Check border color.
    check_colors(border.color, parameter_name = "border.color")

    # Check group.by
    if (is.null(group.by)){
      assertthat::assert_that(!("Groups" %in% colnames(sample@meta.data)),
                              msg = "Please make sure you provide a value for group.by or do not have a metadata column named `Groups`.")

      sample@meta.data[, "Groups"] <- sample@active.ident
      group.by <- "Groups"
    }

    check_parameters(parameter = font.type, parameter_name = "font.type")
    check_parameters(parameter = legend.position, parameter_name = "legend.position")
    check_parameters(parameter = marginal.type, parameter_name = "marginal.type")
    check_parameters(parameter = viridis_color_map, parameter_name = "viridis_color_map")
    check_parameters(parameter = viridis_direction, parameter_name = "viridis_direction")


    # Compute the enrichment scores.
    sample <- compute_enrichment_scores(sample = sample, input_gene_list = input_gene_list, verbose = verbose, nbin = nbin, ctrl = ctrl)

    # 2-variable plot.
    if (is.null(y2) & is.null(x2)){
      # Check that the names provided are not repeated.
      assertthat::assert_that(sum(duplicated(c(x1, y1))) == 0,
                              msg = "The names of the lists to plot can not be the same.")
      # Check that the names provided match the marker genes.
      assertthat::assert_that(x1 %in% names(input_gene_list),
                              msg = paste0(x1, " is not a name of a list of genes provided to input_gene_list."))

      assertthat::assert_that(y1 %in% names(input_gene_list),
                              msg = paste0(y1, " is not a name of a list of genes provided to input_gene_list."))

      # Retrieve metadata variables.
      variables_to_retrieve <- c(x1, y1, group.by)
      # And store them as a tibble.
      scores <- sample@meta.data[, variables_to_retrieve]
      scores[["cell"]] <- rownames(scores)
      # Shuffle the cells so that we accomplish a random plotting, not sample by sample.
      scores <- scores[sample(scores[["cell"]], nrow(scores)), ]
      scores <- tidyr::tibble(scores)

      # Compute scores for the X axis.
      x <- scores %>% dplyr::pull(x1)
      # Compute scores for the Y axis.
      y <- scores %>% dplyr::pull(y1)

      names(x) <- scores[["cell"]]
      names(y) <- scores[["cell"]]

      # Define titles.
      x_lab <- ifelse(is.null(xlab), x1, xlab)
      y_lab <- ifelse(is.null(ylab), y1, ylab)


      # Plot
      df <- data.frame("set_x" = x, "set_y" = y, "group.by" = scores[[group.by]])
      p <- ggplot2::ggplot(data = df,
                           mapping = ggplot2::aes(x = .data[["set_x"]],
                                                  y = .data[["set_y"]],
                                                  color = .data[["group.by"]]))

      if (isFALSE(raster)){
        p <- p +
             ggplot2::geom_point(size = pt.size)
      } else if (isTRUE(raster)){
        p <- p +
             scattermore::geom_scattermore(size = pt.size,
                                           pointsize = pt.size,
                                           pixels = c(raster.dpi, raster.dpi))
      }
      p <- p +
           ggplot2::scale_color_manual(values = colors.use) +
           ggplot2::guides(color = ggplot2::guide_legend(title = "")) +
           ggplot2::xlab(x_lab) +
           ggplot2::ylab(y_lab) +
           ggplot2::labs(title = plot.title,
                         subtitle = plot.subtitle,
                         caption = plot.caption)

      if (isTRUE(enforce_symmetry)){
        # Define limits of polots.
        lim1 <- min(min(x), min(y))
        lim2 <- max(max(x), max(y))
        lim_x <- c(lim1, lim2)
        lim_y <- c(lim1, lim2)

        p <- p  +
             ggplot2::coord_fixed(xlim = lim_x, ylim = lim_y)
      }


    # 3-variable plot.
    } else if (is.null(y2) & !(is.null(x2))){
        # Check that the names provided are not repeated.
        assertthat::assert_that(sum(duplicated(c(x1, y1, x2))) == 0,
                                msg = "The names of the lists to plot can not be the same.")
        # Check that the names provided match the marker genes.
        assertthat::assert_that(x1 %in% names(input_gene_list),
                                msg = paste0(x1, " is not a name of a list of genes provided to input_gene_list."))

        assertthat::assert_that(x2 %in% names(input_gene_list),
                                msg = paste0(x2, " is not a name of a list of genes provided to input_gene_list."))

        assertthat::assert_that(y1 %in% names(input_gene_list),
                                msg = paste0(y1, " is not a name of a list of genes provided to input_gene_list."))

        # Retrieve metadata variables.
        variables_to_retrieve <- c(x1, x2, y1, group.by)
        # And store them as a tibble.
        scores <- sample@meta.data[, variables_to_retrieve]
        scores[["cell"]] <- rownames(scores)
        # Shuffle the cells so that we accomplish a random plotting, not sample by sample.
        scores <- tidyr::tibble(scores)

        # Compute the scores for the X axis.
        x <- unlist(sapply(seq_len(nrow(scores)), function(x) {
          score_1 <- scores[x, x1] + stats::runif(1, min=0, max=0.15)
          score_2 <- scores[x, x2] + stats::runif(1, min=0, max=0.15)
          d <- max(score_1, score_2)
          ifelse(score_1 > score_2, d, -d)
        }))

        # Compute the scores for the Y axis.
        y <- unlist(sapply(seq_len(nrow(scores)), function(x) {
          score_1 <- scores[x, x1] + stats::runif(1, min=0, max=0.15)
          score_2 <- scores[x, x2] + stats::runif(1, min=0, max=0.15)
          d <- max(score_1, score_2)
          y <- scores[x, y1] - d
          y
        }))

        names(x) <- scores[["cell"]]
        names(y) <- scores[["cell"]]

        # Define titles.
        x_lab <- ifelse(is.null(xlab), paste0(x2, "  <---->  ", x1), xlab)
        y_lab <- ifelse(is.null(ylab), y1, ylab)



        # Plot.
        df <- data.frame("set_x" = x, "set_y" = y, "group.by" = scores[[group.by]])
        p <- ggplot2::ggplot(df, mapping = ggplot2::aes(x = .data[["set_x"]],
                                                        y = .data[["set_y"]],
                                                        color = .data[["group.by"]]))

        if (isFALSE(raster)){
          p <- p +
               ggplot2::geom_point(size = pt.size)
        } else if (isTRUE(raster)){
          p <- p +
               scattermore::geom_scattermore(size = pt.size,
                                             pointsize = pt.size,
                                             pixels = c(raster.dpi, raster.dpi))
        }
        p <- p +
             ggplot2::scale_color_manual(values = colors.use) +
             ggplot2::xlab(x_lab) +
             ggplot2::ylab(y_lab) +
             ggplot2::guides(color = ggplot2::guide_legend(title = "")) +
             ggplot2::labs(title = plot.title,
                           subtitle = plot.subtitle,
                           caption = plot.caption)

        if (isTRUE(enforce_symmetry)){
          # Define limits of polots.
          lim <- max(abs(x))
          lim_x <- c(-lim, lim)
          lim <- max(abs(y))
          lim_y <- c(-lim, lim)

          p <- p +
               ggplot2::xlim(lim_x) +
               ggplot2::ylim(lim_y)

        }

    # 4-parameter plot.
    } else if (!is.null(y2) & !(is.null(x2))){
        # Check that the names provided are not repeated.
        assertthat::assert_that(sum(duplicated(c(x1, y1, x2, y2))) == 0,
                                msg = "The names of the lists to plot can not be the same.")
        # Check that the names provided match the marker genes.
        assertthat::assert_that(x1 %in% names(input_gene_list),
                                msg = paste0(x1, " is not a name of a list of genes provided to input_gene_list."))

        assertthat::assert_that(x2 %in% names(input_gene_list),
                                msg = paste0(x2, " is not a name of a list of genes provided to input_gene_list."))

        assertthat::assert_that(y1 %in% names(input_gene_list),
                                msg = paste0(y1, " is not a name of a list of genes provided to input_gene_list."))

        assertthat::assert_that(y2 %in% names(input_gene_list),
                                msg = paste0(y2, " is not a name of a list of genes provided to input_gene_list."))


        # Retrieve metadata variables to plot.
        variables_to_retrieve <- c(x1, x2, y1, y2)
        # And store them as a tibble.
        scores <- sample@meta.data[, variables_to_retrieve]
        # Shuffle the cells so that we accomplish a random plotting, not sample by sample.

        # Compute Y axis values.
        d <- apply(scores, 1, function(x){max(x[c(x1, x2)]) - max(x[c(y1, y2)])})

        # Compute X axis values.
        x <- sapply(seq_along(d), function(x) {
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
        df[["set_x"]] <- x
        df[["set_y"]] <- d
        df[["group.by"]] <- sample@meta.data[, group.by]
        p <- ggplot2::ggplot(df, mapping = ggplot2::aes(x = .data[["set_x"]],
                                                        y = .data[["set_y"]],
                                                        color = .data[["group.by"]]))

        if (isFALSE(raster)){
          p <- p +
               ggplot2::geom_point(size = pt.size)
        } else if (isTRUE(raster)){
          p <- p +
               scattermore::geom_scattermore(size = pt.size,
                                             pointsize = pt.size,
                                             pixels = c(raster.dpi, raster.dpi))
        }
        p <- p +
             ggplot2::scale_color_manual(values = colors.use) +
             ggplot2::xlab(x_lab1) +
             ggplot2::ylab(y_lab1) +
             ggplot2::labs(title = plot.title,
                           subtitle = plot.subtitle,
                           caption = plot.caption)
        suppressMessages({
          p <- p +
               ggplot2::scale_y_continuous(sec.axis = ggplot2::sec_axis(~., name = y_lab2)) +
               ggplot2::scale_x_continuous(sec.axis = ggplot2::sec_axis(~., name = x_lab2))
        })
    if (isTRUE(enforce_symmetry)){
      # Define limits of polots.
      lim_1 <- min(min(d), min(x))
      lim_2 <- max(max(d), max(x))
      value <- max(abs(c(lim_1, lim_2)))
      lim <- c(-value, value)

      suppressMessages({
        p <- p +
             ggplot2::xlim(lim) +
             ggplot2::ylim(lim) +
             ggplot2::scale_y_continuous(sec.axis = ggplot2::sec_axis(~., name = y_lab2)) +
             ggplot2::scale_x_continuous(sec.axis = ggplot2::sec_axis(~., name = x_lab2)) +
             ggplot2::coord_fixed(xlim = c(-value, value), ylim = c(-value, value))
      })

    }

    }



    # Overall formatting for the plot.
    p <- p &
         ggplot2::theme_minimal(base_size = font.size) &
         ggplot2::theme(axis.title = ggplot2::element_text(face = "bold"),
                        axis.line.y.right = ggplot2::element_line(color = "black"),
                        axis.ticks.y.right = ggplot2::element_line(color = "black"),
                        axis.line.x.top = ggplot2::element_line(color = "black"),
                        axis.ticks.x.top = ggplot2::element_line(color = "black"),
                        axis.text.x.top = ggplot2::element_text(face = "bold", color = "black"),
                        axis.text.y.right = ggplot2::element_text(face = "bold", color = "black"),
                        axis.title.x.top = ggplot2::element_text(face = "bold", color = "black"),
                        axis.title.y.right = ggplot2::element_text(face = "bold", color = "black"),
                        axis.text = ggplot2::element_text(face = "bold", color = "black"),
                        plot.title = ggplot2::element_text(face = "bold", hjust = 0, vjust = 0),
                        plot.subtitle = ggplot2::element_text(hjust = 0),
                        plot.caption = ggplot2::element_text(hjust = 1),
                        plot.title.position = "plot",
                        panel.grid = ggplot2::element_blank(),
                        text = ggplot2::element_text(family = font.type),
                        plot.caption.position = "plot",
                        legend.text = ggplot2::element_text(face = "bold"),
                        legend.position = legend.position,
                        legend.title = ggplot2::element_text(face = "bold"),
                        legend.justification = "center",
                        plot.margin = ggplot2::margin(t = 10, r = 10, b = 10, l = 10),
                        axis.ticks = ggplot2::element_line(color = "black"),
                        axis.line = ggplot2::element_line(color = "black"),
                        plot.background = ggplot2::element_rect(fill = "white", color = "white"),
                        panel.background = ggplot2::element_rect(fill = "white", color = "white"),
                        legend.background = ggplot2::element_rect(fill = "white", color = "white")) &
         ggplot2::guides(color = ggplot2::guide_legend(title = "",
                                                       ncol = legend.ncol,
                                                       nrow = legend.nrow,
                                                       byrow = legend.byrow,
                                                       override.aes = list(size = legend.icon.size)))

    # Add cell borders.
    if (isTRUE(plot_cell_borders)){
      if (isFALSE(raster)){
        base_layer <-  ggplot2::geom_point(data = df,
                                           mapping = ggplot2::aes(x = .data[["set_x"]],
                                                                  y = .data[["set_y"]]),
                                           size = pt.size * border.size,
                                           color = border.color,
                                           show.legend = FALSE)
      } else if (isTRUE(raster)){
        base_layer <-  scattermore::geom_scattermore(data = df,
                                                     mapping = ggplot2::aes(x = .data[["set_x"]],
                                                                            y = .data[["set_y"]]),
                                                     size = pt.size * border.size,
                                                     stroke = pt.size / 2,
                                                     color = border.color,
                                                     pointsize = pt.size * border.size,
                                                     pixels = c(raster.dpi, raster.dpi),
                                                     show.legend = FALSE)
      }
      p[["layers"]] <- append(base_layer, p[["layers"]])
    }

    if (isTRUE(plot_features) | isTRUE(plot_enrichment_scores)){
      if (isTRUE(plot_features)){
        assertthat::assert_that(!is.null(features),
                                msg = "Please provide features to plot.")
      }

      output_list <- list()

      # Generate a mock DimRed object for the plots.
      df.use <- df[, c("set_x", "set_y")]
      colnames(df.use) <- c("DIM_1", "DIM_2")
      sample@reductions[["test"]] <- Seurat::CreateDimReducObject(embeddings = as.matrix(df.use), assay = "SCT")

      if (isTRUE(plot_features) & isTRUE(plot_enrichment_scores)){
        features <- c(features, names(input_gene_list))
      } else if (isFALSE(plot_features) & isTRUE(plot_enrichment_scores)){
        features <- names(input_gene_list)
      }

      for (feature in features){
        p.feature <- SCpubr::do_FeaturePlot(sample = sample,
                                            features = feature,
                                            reduction = "test",
                                            plot_cell_borders = plot_cell_borders,
                                            pt.size = pt.size,
                                            legend.position = legend.position,
                                            border.size = border.size,
                                            border.color = border.color,
                                            raster = raster,
                                            raster.dpi = raster.dpi,
                                            font.type = font.type,
                                            font.size = font.size,
                                            viridis_color_map = viridis_color_map,
                                            viridis_direction = viridis_direction)

        # Add back the missing aesthetics.
        if (is.null(y2) & is.null(x2)){
          # Define titles.
          p.feature <- p.feature +
                       ggplot2::xlab(x_lab) +
                       ggplot2::ylab(y_lab)

          if (isTRUE(enforce_symmetry)){
            suppressMessages({
              p.feature <- p.feature  +
                           ggplot2::coord_fixed(xlim = lim_x, ylim = lim_y)
            })
          }
        } else if (is.null(y2) & !(is.null(x2))){
          # Define titles.
          p.feature <- p.feature +
                       ggplot2::xlab(x_lab) +
                       ggplot2::ylab(y_lab)

          if (isTRUE(enforce_symmetry)){
            suppressMessages({
              p.feature <- p.feature +
                           ggplot2::xlim(lim_x) +
                           ggplot2::ylim(lim_y)
            })

          }
        } else if (!is.null(y2) & !(is.null(x2))){
          p.feature <- p.feature +
                       ggplot2::xlab(x_lab1) +
                       ggplot2::ylab(y_lab1)

          if (isTRUE(enforce_symmetry)){
            suppressMessages({
              p.feature <- p.feature +
                           ggplot2::xlim(lim) +
                           ggplot2::ylim(lim) +
                           ggplot2::scale_y_continuous(sec.axis = ggplot2::sec_axis(~., name = y_lab2)) +
                           ggplot2::scale_x_continuous(sec.axis = ggplot2::sec_axis(~., name = x_lab2)) +
                           ggplot2::coord_fixed(xlim = c(-value, value), ylim = c(-value, value))
            })
          }
        }
        p.feature <- p.feature +
                     ggplot2::theme_minimal(base_size = font.size) &
                     ggplot2::theme(axis.title = ggplot2::element_text(face = "bold"),
                                    axis.text = ggplot2::element_text(face = "bold", color = "black"),
                                    plot.title = ggplot2::element_text(face = "bold", hjust = 0, vjust = 0),
                                    plot.subtitle = ggplot2::element_text(hjust = 0),
                                    plot.caption = ggplot2::element_text(hjust = 1),
                                    plot.title.position = "plot",
                                    panel.grid = ggplot2::element_blank(),
                                    text = ggplot2::element_text(family = font.type),
                                    plot.caption.position = "plot",
                                    legend.text = ggplot2::element_text(face = "bold"),
                                    legend.position = legend.position,
                                    legend.title = ggplot2::element_text(face = "bold"),
                                    legend.justification = "center",
                                    plot.margin = ggplot2::margin(t = 10, r = 10, b = 10, l = 10),
                                    axis.ticks = ggplot2::element_line(color = "black"),
                                    axis.line = ggplot2::element_line(color = "black"),
                                    plot.background = ggplot2::element_rect(fill = "white", color = "white"),
                                    panel.background = ggplot2::element_rect(fill = "white", color = "white"),
                                    legend.background = ggplot2::element_rect(fill = "white", color = "white"))
        output_list[[feature]] <- p.feature

      }
    }

    if (isTRUE(plot_marginal_distributions)){
      # Remove annoying warnings when violin is used as marginal distribution.
      if (marginal.type == "violin"){
        p <- suppressWarnings({ggExtra::ggMarginal(p = p,
                                                   groupColour = ifelse(isTRUE(marginal.group), TRUE, FALSE),
                                                   groupFill = ifelse(isTRUE(marginal.group), TRUE, FALSE),
                                                   type = marginal.type,
                                                   size = marginal.size)})
      } else {
        p <- ggExtra::ggMarginal(p = p,
                                 groupColour = ifelse(isTRUE(marginal.group), TRUE, FALSE),
                                 groupFill = ifelse(isTRUE(marginal.group), TRUE, FALSE),
                                 type = marginal.type,
                                 size = marginal.size)
      }

      # Transform back to ggplot2 object.
      p <- ggplotify::as.ggplot(p)

      # Fix for the plot backgrounds after applying ggMarginal.
      p[["theme"]][["plot.background"]] <- ggplot2::element_rect(fill = "white", color = "white")
      p[["theme"]][["legend.background"]] <- ggplot2::element_rect(fill = "white", color = "white")
      p[["theme"]][["panel.background"]] <- ggplot2::element_rect(fill = "white", color = "white")
    }

    # Remove axis ticks?
    if (axis.ticks == FALSE){
        p <- p +
             ggplot2::theme(axis.ticks = ggplot2::element_blank())
    }

    # Remove axis text?
    if (axis.text == FALSE){
        p <- p +
             ggplot2::theme(axis.text = ggplot2::element_blank())
    }

    if (isTRUE(plot_features) | isTRUE(plot_enrichment_scores)){
      output_list[["main"]] <- p
      return_object <- output_list
    } else if (isFALSE(plot_features) & isFALSE(plot_enrichment_scores)){
      return_object <- p
    }

    return(return_object)

}
