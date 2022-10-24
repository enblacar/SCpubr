
#' Create Bar Plots.
#'
#' @inheritParams doc_function
#' @param group.by \strong{\code{\link[base]{character}}} | Metadata column to compute the counts of. Has to be either a character or factor column.
#' @param split.by \strong{\code{\link[base]{character}}} | Metadata column to split the values of group.by by. If not used, defaults to the active idents.
#' @param order \strong{\code{\link[base]{logical}}} | Whether to order the results in descending order of counts.
#' @param position \strong{\code{\link[base]{character}}} | Position function from \pkg{ggplot2}. One of:
#' \itemize{
#'   \item \emph{\code{stack}}: Set the bars side by side, displaying the total number of counts. Uses \link[ggplot2]{position_stack}.
#'   \item \emph{\code{fill}}: Set the bars on top of each other, displaying the proportion of counts from the total that each group represents. Uses \link[ggplot2]{position_fill}.
#' }
#'
#' @return A ggplot2 object containing a Bar plot.
#' @export
#'
#' @example /man/examples/examples_do_BarPlot.R
do_BarPlot <- function(sample,
                       group.by,
                       order = TRUE,
                       split.by = NULL,
                       position = "stack",
                       font.size = 14,
                       font.type = "sans",
                       legend.position = "bottom",
                       rotate_x_axis_labels = 45,
                       legend.title = NULL,
                       xlab = NULL,
                       ylab = NULL,
                       colors.use = NULL,
                       flip = FALSE,
                       plot.title = NULL,
                       plot.subtitle = NULL,
                       plot.caption = NULL,
                       plot.grid = TRUE,
                       grid.color = "grey75",
                       grid.type = "dashed") {
  check_suggests(function_name = "do_BarPlot")
  check_Seurat(sample)

  `%>%` <- magrittr::`%>%`
  `:=` <- rlang::`:=`

  # Check logical parameters.
  logical_list <- list("order" = order,
                       "flip" = flip,
                       "plot.grid" = plot.grid)
  check_type(parameters = logical_list, required_type = "logical", test_function = is.logical)
  # Check numeric parameters.
  numeric_list <- list("font.size" = font.size,
                       "rotate_x_axis_labels" = rotate_x_axis_labels)
  check_type(parameters = numeric_list, required_type = "numeric", test_function = is.numeric)
  # Check character parameters.

  character_list <- list("group.by" = group.by,
                         "split.by" = split.by,
                         "position" = position,
                         "font.type" = font.type,
                         "legend.position" = legend.position,
                         "legend.title" = legend.title,
                         "xlab" = xlab,
                         "ylab" = ylab,
                         "plot.title" = plot.title,
                         "plot.subtitle" = plot.subtitle,
                         "plot.caption" = plot.caption,
                         "grid.color" = grid.color,
                         "grid.type" = grid.type)
  # Checks
  check_type(parameters = character_list, required_type = "character", test_function = is.character)

  check_colors(grid.color, parameter_name = "grid.color")

  check_parameters(parameter = font.type, parameter_name = "font.type")
  check_parameters(parameter = legend.position, parameter_name = "legend.position")
  check_parameters(parameter = grid.type, parameter_name = "grid.type")
  check_parameters(parameter = rotate_x_axis_labels, parameter_name = "rotate_x_axis_labels")

  # Get the general table.
  assertthat::assert_that(class(sample@meta.data[, group.by]) %in% c("character", "factor"),
                          msg = "This function only works with categorical metadta variables supplied to feature.")

  if (is.null(colors.use)){
    colors.use <- generate_color_scale(names_use = if (is.factor(sample@meta.data[, group.by])) {levels(sample@meta.data[, group.by])} else {sort(unique(sample@meta.data[, group.by]))})
  } else {
    check_colors(colors.use, parameter_name = "colors.use")
    check_consistency_colors_and_names(sample = sample, colors = colors.use, grouping_variable = group.by)
  }
  data <-  sample@meta.data %>%
           tibble::as_tibble() %>%
           dplyr::select(dplyr::all_of(c(group.by, split.by))) %>%
           dplyr::mutate("{group.by}" := if(isFALSE(order)) {.data[[group.by]]} else {factor(as.character(.data[[group.by]]), levels = {sample@meta.data %>%
                                                                                                                                        tibble::as_tibble() %>%
                                                                                                                                        dplyr::select(dplyr::all_of(c(group.by, split.by))) %>%
                                                                                                                                        dplyr::group_by(.data[[group.by]]) %>%
                                                                                                                                        dplyr::summarise("n" = dplyr::n()) %>%
                                                                                                                                        dplyr::arrange(if(isFALSE(flip)){dplyr::desc(.data[["n"]])} else {.data[["n"]]}) %>%
                                                                                                                                        dplyr::pull(.data[[group.by]]) %>%
                                                                                                                                        as.character()})})
  if (is.null(split.by)){
    p <- data %>%
         ggplot2::ggplot(mapping = ggplot2::aes(x = .data[[group.by]],
                                                fill = .data[[group.by]]))
  } else {
    p <- data %>%
         ggplot2::ggplot(mapping = ggplot2::aes(x = .data[[split.by]],
                                                fill = .data[[group.by]]))
  }
  p <- p +
       ggplot2::stat_count(geom = "bar", position = position, color = "black") +
       ggplot2::xlab(if (!is.null(split.by) & is.null(xlab)) {split.by} else {ifelse(is.null(group.by), "Idents", group.by)}) +
       ggplot2::ylab(ifelse(is.null(ylab), paste0(ifelse(position == "stack", "Count", "Frequency"), " of ", group.by), ylab)) +
       ggplot2::labs(title = plot.title,
                     subtitle = plot.subtitle,
                     caption = plot.caption) +
       ggplot2::scale_fill_manual(values = colors.use) +
       ggplot2::guides(fill = ggplot2::guide_legend(title = if (!is.null(split.by) & is.null(legend.title)) {ifelse(is.null(group.by), "Idents", group.by)} else {legend.title},
                                                    title.position = "top",
                                                    title.hjust = 0.5)) +
       ggplot2::theme_minimal(base_size = font.size) +
       ggplot2::theme(axis.title = ggplot2::element_text(color = "black",
                                                         face = "bold"),
                      panel.grid.major.y = if (isFALSE(flip)) {if (isTRUE(plot.grid)){ggplot2::element_line(color = grid.color, linetype = grid.type)}} else if (isTRUE(flip)) {ggplot2::element_blank()},
                      panel.grid.major.x = if (isTRUE(flip)) {if (isTRUE(plot.grid)){ggplot2::element_line(color = grid.color, linetype = grid.type)}} else if (isFALSE(flip)) {ggplot2::element_blank()},
                      axis.line.x = if (isFALSE(flip)) {ggplot2::element_line(color = "black")} else if (isTRUE(flip)) {ggplot2::element_blank()},
                      axis.line.y = if (isTRUE(flip)) {ggplot2::element_line(color = "black")} else if (isFALSE(flip)) {ggplot2::element_blank()},
                      axis.text.x = ggplot2::element_text(color = "black",
                                                          face = "bold",
                                                          angle = get_axis_parameters(angle = rotate_x_axis_labels, flip = flip)[["angle"]],
                                                          hjust = get_axis_parameters(angle = rotate_x_axis_labels, flip = flip)[["hjust"]],
                                                          vjust = get_axis_parameters(angle = rotate_x_axis_labels, flip = flip)[["vjust"]]),
                      axis.text.y = ggplot2::element_text(color = "black", face = "bold"),
                      axis.ticks = ggplot2::element_line(color = "black"),
                      plot.title.position = "plot",
                      plot.title = ggplot2::element_text(face = "bold", hjust = 0),
                      plot.subtitle = ggplot2::element_text(hjust = 0),
                      plot.caption = ggplot2::element_text(hjust = 1),
                      panel.grid = ggplot2::element_blank(),
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
