#' Generate Alluvial plots.
#'
#' This function is based on the \pkg{ggalluvial} package. It allows you to generate alluvial plots from a given Seurat object.
#'
#' @inheritParams doc_function
#'
#' @param first_group \strong{\code{\link[base]{character}}} | Categorical metadata variable. First group of nodes of the alluvial plot.
#' @param last_group \strong{\code{\link[base]{character}}} | Categorical metadata variable. Last group of nodes of the alluvial plot.
#' @param middle_groups \strong{\code{\link[base]{character}}} | Categorical metadata variable. Vector of groups of nodes of the alluvial plot.
#' @param colors.use \strong{\code{\link[base]{character}}} | Named list of colors corresponding to the unique values in fill.by (which defaults to last_group).
#' @param fill.by \strong{\code{\link[base]{character}}} | One of first_group, middle_groups (one of the values, if multiple mid_groups) or last_group. These values will be used to color the alluvium/flow.
#' @param use_labels \strong{\code{\link[base]{logical}}} | Whether to use labels instead of text for the stratum.
#' @param stratum.color,alluvium.color,flow.color \strong{\code{\link[base]{character}}} | Color for the border of the alluvium (and flow) and stratum.
#' @param stratum.fill \strong{\code{\link[base]{character}}} | Color to fill the stratum.
#' @param stratum.width \strong{\code{\link[base]{logical}}} | Width of the stratum.
#' @param stratum.fill.conditional \strong{\code{\link[base]{logical}}} | Whether to fill the stratum with the same colors as the alluvium/flow.
#' @param use_geom_flow \strong{\code{\link[base]{logical}}} | Whether to use \code{\link[ggalluvial]{geom_flow}} instead of \code{\link[ggalluvial]{geom_alluvium}}. Visual results might differ.
#' @param label.color \strong{\code{\link[base]{character}}} | Color for the text labels.
#' @param curve_type \strong{\code{\link[base]{character}}} | Type of curve used in \code{\link[ggalluvial]{geom_alluvium}}. One of:
#' \itemize{
#'   \item \emph{\code{linear}}.
#'   \item \emph{\code{cubic}}.
#'   \item \emph{\code{quintic}}.
#'   \item \emph{\code{sine}}.
#'   \item \emph{\code{arctangent}}.
#'   \item \emph{\code{sigmoid}}.
#'   \item \emph{\code{xspline}}.
#' }
#'
#' @return A ggplot2 object.
#' @export
#'
#' @example /man/examples/examples_do_AlluvialPlot.R
do_AlluvialPlot <- function(sample,
                            first_group,
                            last_group,
                            middle_groups = NULL,
                            colors.use = NULL,
                            plot.title = NULL,
                            plot.subtitle = NULL,
                            plot.caption = NULL,
                            font.size = 14,
                            font.type = "sans",
                            xlab = NULL,
                            ylab = "Number of cells",
                            repel = FALSE,
                            fill.by = last_group,
                            use_labels = FALSE,
                            stratum.color = "black",
                            stratum.fill = "white",
                            stratum.width = 1/3,
                            stratum.fill.conditional = FALSE,
                            use_geom_flow = FALSE,
                            alluvium.color = "white",
                            flow.color = "white",
                            flip = FALSE,
                            label.color = "black",
                            curve_type = "sigmoid",
                            use_viridis = FALSE,
                            viridis_color_map = "G",
                            viridis_direction = -1,
                            plot.grid = FALSE,
                            grid.color = "grey75",
                            grid.type = "dashed",
                            na.value = "white",
                            legend.position = "right",
                            legend.title = NULL){
  check_suggests(function_name = "do_AlluvialPlot")
  check_Seurat(sample)

  StatStratum <- ggalluvial::StatStratum
  `%>%` <- magrittr::`%>%`

  if (isTRUE(use_labels)){
    if (isTRUE(repel)){
      func_use <- ggrepel::geom_label_repel
    } else if (isFALSE(repel)){
      func_use <- ggplot2::geom_label
    }
  } else if (isFALSE(use_labels)){
    if (isTRUE(repel)){
      func_use <- ggrepel::geom_text_repel
    } else if (isFALSE(repel)){
      func_use <- ggplot2::geom_text
    }
  }

  vars.use <- c(first_group)
  for (variable in middle_groups){
    vars.use <- append(vars.use, variable)
  }
  vars.use <- append(vars.use, last_group)

  for (var in vars.use){
    assertthat::assert_that(var %in% colnames(sample@meta.data),
                            msg = "Please make sure that the variables provided to first_group, middle_groups and last_group are in sample@meta.data.")

    assertthat::assert_that(class(sample@meta.data[, var]) %in% c("character", "factor"),
                            msg = "Please make sure that the variables provided to first_group, middle_groups and last_group are either characters or factors.")
  }

  assertthat::assert_that(length(fill.by) == 1,
                          msg = "Parameter fill.by has to be a single value.")


  assertthat::assert_that(isTRUE(fill.by %in% vars.use),
                          msg = "Parameter fill.by has to be the same as one of the values in first_group, last_group or middle_groups.")
  suppressMessages({
    data <- sample@meta.data %>%
            dplyr::select(dplyr::all_of(vars.use)) %>%
            dplyr::group_by_at(vars.use) %>%
            dplyr::summarise(n = dplyr::n())
  })


  # COLORS.
  if (is.null(colors.use)){
    if (is.factor(data[[fill.by]])){
      colors.use = generate_color_scale(levels(data[[fill.by]]))
    } else {
      colors.use = generate_color_scale(sort(unique(data[[fill.by]])))
    }
  } else {
    check_colors(colors.use)
  }


  p <- prepare_ggplot_alluvial_plot(data = data,
                                    vars.use = vars.use)

  if (isTRUE(use_geom_flow)){
    p <- p +
         ggalluvial::geom_flow(mapping = ggplot2::aes(fill = .data[[fill.by]]),
                               color = flow.color)
  } else if (isFALSE(use_geom_flow)){
    p <- p +
         ggalluvial::geom_alluvium(mapping = ggplot2::aes(fill = .data[[fill.by]]),
                                                          color = alluvium.color,
                                                          curve_type = curve_type)
  }
  if (isTRUE(stratum.fill.conditional)){
    p <- p +
         ggalluvial::geom_stratum(color = stratum.color,
                                  mapping = ggplot2::aes(fill = .data[[fill.by]]),
                                  width = stratum.width)
  } else if (isFALSE(stratum.fill.conditional)){
    p <- p +
         ggalluvial::geom_stratum(color = stratum.color,
                                  fill = stratum.fill,
                                  width = stratum.width)
  }
  p <- p  +
    func_use(stat = ggalluvial::StatStratum,
             mapping = ggplot2::aes(label = ggplot2::after_stat(stratum)),
             color = label.color,
             fontface = "bold") +
    ggplot2::scale_x_discrete(limits = vars.use)

  if (isTRUE(use_viridis)){
    p <- p +
         ggplot2::scale_fill_viridis_d(option = viridis_color_map,
                                       direction = viridis_direction,
                                       na.value = na.value)
  } else if (isFALSE(use_viridis)){
    p <- p +
         ggplot2::scale_fill_manual(values = colors.use,
                                    na.value = na.value)
  }
  p <- p +
    ggplot2::xlab(xlab) +
    ggplot2::ylab(ylab) +
    ggplot2::labs(title = plot.title,
                  subtitle = plot.subtitle,
                  caption = plot.caption) +
    ggplot2::guides(fill = ggplot2::guide_legend(title = legend.title)) +
    ggplot2::theme_minimal(base_size = font.size) +
    ggplot2::theme(axis.title = ggplot2::element_text(color = "black",
                                                      face = "bold"),
                   axis.line.x = if (isFALSE(flip)){ggplot2::element_blank()} else {ggplot2::element_line(color = "black")},
                   axis.line.y = if (isFALSE(flip)){ggplot2::element_line(color = "black")} else {ggplot2::element_blank()},
                   axis.ticks.y = if (isFALSE(flip)){ggplot2::element_line(color = "black")} else {ggplot2::element_blank()},
                   axis.ticks.x = if (isFALSE(flip)){ggplot2::element_blank()} else {ggplot2::element_line(color = "black")},
                   axis.text.y = ggplot2::element_text(color = "black", face = "bold"),
                   axis.text.x = ggplot2::element_text(color = "black", face = "bold"),
                   panel.grid.major = ggplot2::element_blank(),
                   plot.title.position = "plot",
                   plot.title = ggplot2::element_text(face = "bold", hjust = 0),
                   plot.subtitle = ggplot2::element_text(hjust = 0),
                   plot.caption = ggplot2::element_text(hjust = 1),
                   panel.grid = ggplot2::element_blank(),
                   panel.grid.major.y = if (isFALSE(flip)) {if (isTRUE(plot.grid)){ggplot2::element_line(color = grid.color, linetype = grid.type)}} else if (isTRUE(flip)) {ggplot2::element_blank()},
                   panel.grid.major.x = if (isTRUE(flip)) {if (isTRUE(plot.grid)){ggplot2::element_line(color = grid.color, linetype = grid.type)}} else if (isFALSE(flip)) {ggplot2::element_blank()},
                   text = ggplot2::element_text(family = font.type),
                   plot.caption.position = "plot",
                   legend.text = ggplot2::element_text(face = "bold"),
                   legend.position = legend.position,
                   legend.title = ggplot2::element_text(face = "bold"),
                   legend.justification = "center",
                   plot.margin = ggplot2::margin(t = 10, r = 10, b = 10, l = 10),
                   plot.background = ggplot2::element_rect(fill = "white", color = "white"),
                   panel.background = ggplot2::element_rect(fill = "white", color = "white"),
                   legend.background = ggplot2::element_rect(fill = "white", color = "white"),
                   strip.text =ggplot2::element_text(color = "black", face = "bold"))

  if (isTRUE(flip)){
    p <- p + ggplot2::coord_flip()
  }
  return(p)
}
