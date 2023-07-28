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
                            viridis.palette = "G",
                            viridis.direction = -1,
                            sequential.palette = "YlGnBu",
                            sequential.direction = 1,
                            plot.grid = FALSE,
                            grid.color = "grey75",
                            grid.type = "dashed",
                            na.value = "white",
                            legend.position = "right",
                            legend.title = NULL,
                            plot.title.face = "bold",
                            plot.subtitle.face = "plain",
                            plot.caption.face = "italic",
                            axis.title.face = "bold",
                            axis.text.face = "plain",
                            legend.title.face = "bold",
                            legend.text.face = "plain"){
  # Add lengthy error messages.
  withr::local_options(.new = list("warning.length" = 8170))
  
  check_suggests(function_name = "do_AlluvialPlot")
  check_Seurat(sample)
  
  
  

  # Check logical parameters.
  logical_list <- list("use_labels" = use_labels,
                       "stratum.fill.conditional" = stratum.fill.conditional,
                       "flip" = flip,
                       "plot.grid" = plot.grid,
                       "repel" = repel,
                       "use_geom_flow" = use_geom_flow,
                       "use_viridis" = use_viridis)
  check_type(parameters = logical_list, required_type = "logical", test_function = is.logical)
  # Check numeric parameters.
  numeric_list <- list("stratum.width" = stratum.width,
                       "font.size" = font.size,
                       "viridis.direction" = viridis.direction,
                       "sequential.direction" = sequential.direction)
  check_type(parameters = numeric_list, required_type = "numeric", test_function = is.numeric)
  # Check character parameters.
  character_list <- list("first_group" = first_group,
                         "last_group" = last_group,
                         "middle_groups" = middle_groups,
                         "colors.use" = colors.use,
                         "plot.title" = plot.title,
                         "plot.subtitle" = plot.subtitle,
                         "plot.caption" = plot.caption,
                         "font.type" = font.type,
                         "xlab" = xlab,
                         "ylab" = ylab,
                         "fill.by" = fill.by,
                         "stratum.color" = stratum.color,
                         "stratum.fill" = stratum.fill,
                         "alluvium.color" = alluvium.color,
                         "flow.color" = flow.color,
                         "label.color" = label.color,
                         "curve_type" = curve_type,
                         "viridis.palette" = viridis.palette,
                         "grid.color" = grid.color,
                         "grid.type" = grid.type,
                         "na.value" = na.value,
                         "legend.position" = legend.position,
                         "legend.title" = legend.title,
                         "plot.title.face" = plot.title.face,
                         "plot.subtitle.face" = plot.subtitle.face,
                         "plot.caption.face" = plot.caption.face,
                         "axis.title.face" = axis.title.face,
                         "axis.text.face" = axis.text.face,
                         "legend.title.face" = legend.title.face,
                         "legend.text.face" = legend.text.face,
                         "sequential.palette" = sequential.palette)
  check_type(parameters = character_list, required_type = "character", test_function = is.character)

  check_parameters(parameter = font.type, parameter_name = "font.type")
  check_parameters(parameter = legend.position, parameter_name = "legend.position")
  check_parameters(plot.title.face, parameter_name = "plot.title.face")
  check_parameters(plot.subtitle.face, parameter_name = "plot.subtitle.face")
  check_parameters(plot.caption.face, parameter_name = "plot.caption.face")
  check_parameters(axis.title.face, parameter_name = "axis.title.face")
  check_parameters(axis.text.face, parameter_name = "axis.text.face")
  check_parameters(legend.title.face, parameter_name = "legend.title.face")
  check_parameters(legend.text.face, parameter_name = "legend.text.face")
  check_parameters(viridis.direction, parameter_name = "viridis.direction")
  check_parameters(sequential.direction, parameter_name = "sequential.direction")
  

  #StatStratum <- ggalluvial::StatStratum
  `%>%` <- magrittr::`%>%`
  
  colors.gradient <- compute_continuous_palette(name = ifelse(isTRUE(use_viridis), viridis.palette, sequential.palette),
                                                use_viridis = use_viridis,
                                                direction = ifelse(isTRUE(use_viridis), viridis.direction, sequential.direction),
                                                enforce_symmetry = FALSE)
  
  if (isTRUE(use_labels)){
    if (isTRUE(repel)){
      func_use <- ggrepel::geom_label_repel
    } else if (base::isFALSE(repel)){
      func_use <- ggplot2::geom_label
    }
  } else if (base::isFALSE(use_labels)){
    if (isTRUE(repel)){
      func_use <- ggrepel::geom_text_repel
    } else if (base::isFALSE(repel)){
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
                            msg = paste0(add_cross(), crayon_body("Please make sure that the variables provided to "),
                                         crayon_key("first_group"),
                                         crayon_body(", "),
                                         crayon_key("middle_groups"),
                                         crayon_body(" and "),
                                         crayon_key("last_group"),
                                         crayon_body(" are "),
                                         crayon_key("metadata variables"),
                                         crayon_body(".")))

    assertthat::assert_that(class(sample@meta.data[, var]) %in% c("character", "factor"),
                            msg = paste0(add_cross(), crayon_body("Please make sure that the variables provided to "),
                                         crayon_key("first_group"),
                                         crayon_body(", "),
                                         crayon_key("middle_groups"),
                                         crayon_body(" and "),
                                         crayon_key("last_group"),
                                         crayon_body(" are of class "),
                                         crayon_key("character"),
                                         crayon_body(" or "),
                                         crayon_key("factor"),
                                         crayon_body(".")))
  }

  assertthat::assert_that(length(fill.by) == 1,
                          msg = paste0(add_cross(), crayon_body("Please provide a single value to "),
                                       crayon_key("fill.by"),
                                       crayon_body(".")))


  assertthat::assert_that(isTRUE(fill.by %in% vars.use),
                          msg = paste0(add_cross(), crayon_body("Paramter "),
                                       crayon_key("fill.by"),
                                       crayon_body(" has to be the same as one of the values in "),
                                       crayon_key("first_group"),
                                       crayon_body(", "),
                                       crayon_key("middle_groups"),
                                       crayon_body(" and "),
                                       crayon_key("last_group"),
                                       crayon_body(".")))
  suppressMessages({
    data <- sample@meta.data %>%
            dplyr::select(dplyr::all_of(vars.use)) %>%
            dplyr::group_by_at(vars.use) %>%
            dplyr::summarise(n = dplyr::n())
  })


  # COLORS.
  if (is.null(colors.use)){
    if (is.factor(data[[fill.by]])){
      colors.use <- generate_color_scale(levels(data[[fill.by]]))
    } else {
      colors.use <- generate_color_scale(sort(unique(data[[fill.by]])))
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
  } else if (base::isFALSE(use_geom_flow)){
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
  } else if (base::isFALSE(stratum.fill.conditional)){
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

  if (is.null(colors.use)){
    p <- p + 
         ggplot2::scale_fill_gradientn(colors = colors.gradient,
                                       na.value = na.value,
                                       name = legend.title)
  } else if (base::isFALSE(use_viridis)){
    p <- p +
         ggplot2::scale_fill_manual(values = colors.use,
                                    na.value = na.value,
                                    name = legend.title)
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
                                                         face = axis.title.face),
                      axis.line.x = if (base::isFALSE(flip)){ggplot2::element_blank()} else {ggplot2::element_line(color = "black")},
                      axis.line.y = if (base::isFALSE(flip)){ggplot2::element_line(color = "black")} else {ggplot2::element_blank()},
                      axis.ticks.y = if (base::isFALSE(flip)){ggplot2::element_line(color = "black")} else {ggplot2::element_blank()},
                      axis.ticks.x = if (base::isFALSE(flip)){ggplot2::element_blank()} else {ggplot2::element_line(color = "black")},
                      axis.text.y = ggplot2::element_text(color = "black", face = axis.text.face),
                      axis.text.x = ggplot2::element_text(color = "black", face = axis.text.face),
                      panel.grid.major = ggplot2::element_blank(),
                      plot.title.position = "plot",
                      plot.title = ggplot2::element_text(face = plot.title.face, hjust = 0),
                      plot.subtitle = ggplot2::element_text(face = plot.subtitle.face, hjust = 0),
                      plot.caption = ggplot2::element_text(face = plot.caption.face, hjust = 1),
                      panel.grid = ggplot2::element_blank(),
                      panel.grid.major.y = if (base::isFALSE(flip)) {if (isTRUE(plot.grid)){ggplot2::element_line(color = grid.color, linetype = grid.type)}} else if (isTRUE(flip)) {ggplot2::element_blank()},
                      panel.grid.major.x = if (isTRUE(flip)) {if (isTRUE(plot.grid)){ggplot2::element_line(color = grid.color, linetype = grid.type)}} else if (base::isFALSE(flip)) {ggplot2::element_blank()},
                      text = ggplot2::element_text(family = font.type),
                      plot.caption.position = "plot",
                      legend.text = ggplot2::element_text(face = legend.text.face),
                      legend.position = legend.position,
                      legend.title = ggplot2::element_text(face = legend.title.face),
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
