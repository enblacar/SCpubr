
#' Create Bar Plots.
#'
#' @inheritParams doc_function
#' @param group.by \strong{\code{\link[base]{character}}} | Metadata column to compute the counts of. Has to be either a character or factor column.
#' @param split.by \strong{\code{\link[base]{character}}} | Metadata column to split the values of group.by by. If not used, defaults to the active idents.
#' @param facet.by \strong{\code{\link[base]{character}}} | Metadata column to gather the columns by. This is useful if you have other overarching metadata.
#' @param order \strong{\code{\link[base]{logical}}} | Whether to order the results in descending order of counts.
#' @param order.by \strong{\code{\link[base]{character}}} | When \strong{\code{split.by}} is used, value of \strong{\code{group.by}} to reorder the columns based on its value.
#' @param position \strong{\code{\link[base]{character}}} | Position function from \pkg{ggplot2}. Either stack or fill.
#' @param return_data \strong{\code{\link[base]{logical}}} | Returns a data.frame with the count and proportions displayed in the plot.
#' @param add.n \strong{\code{\link[base]{logical}}} | Whether to add the total counts on top of each bar.
#' @param add.n.face \strong{\code{\link[base]{character}}} | Font face of the labels added by \strong{\code{add.n}}.
#' @param add.n.size \strong{\code{\link[base]{numeric}}} | Size of the labels
#' @param add.n.expand \strong{\code{\link[base]{numeric}}} | Vector of two numerics representing the start and end of the scale. Minimum should be 0 and max should be above 1. This basically expands the Y axis so that the labels fit when \strong{\code{flip = TRUE}}.
#' \itemize{
#'   \item \emph{\code{stack}}: Set the bars side by side, displaying the total number of counts. Uses \link[ggplot2]{position_stack}.
#'   \item \emph{\code{fill}}: Set the bars on top of each other, displaying the proportion of counts from the total that each group represents. Uses \link[ggplot2]{position_fill}.
#' }
#' @param strip.text.face \strong{\code{\link[base]{character}}} | Controls the style of the font for the strip text.  One of:
#' \itemize{
#'   \item \emph{\code{plain}}: For normal text.
#'   \item \emph{\code{italic}}: For text in itallic.
#'   \item \emph{\code{bold}}: For text in bold.
#'   \item \emph{\code{bold.italic}}: For text both in itallic and bold.
#' }
#' @return A ggplot2 object containing a Bar plot.
#' @export
#'
#' @example /man/examples/examples_do_BarPlot.R
do_BarPlot <- function(sample,
                       group.by,
                       order = FALSE,
                       add.n = FALSE,
                       add.n.face = "bold",
                       add.n.expand = c(0, 1.15),
                       add.n.size = 4,
                       order.by = NULL,
                       split.by = NULL,
                       facet.by = NULL,
                       position = "stack",
                       font.size = 14,
                       font.type = "sans",
                       legend.position = "bottom",
                       legend.title = NULL,
                       legend.ncol = NULL,
                       legend.nrow = NULL,
                       legend.byrow = FALSE,
                       axis.text.x.angle = 45,
                       xlab = NULL,
                       ylab = NULL,
                       colors.use = NULL,
                       flip = FALSE,
                       plot.title = NULL,
                       plot.subtitle = NULL,
                       plot.caption = NULL,
                       plot.grid = FALSE,
                       grid.color = "grey75",
                       grid.type = "dashed",
                       plot.title.face = "bold",
                       plot.subtitle.face = "plain",
                       plot.caption.face = "italic",
                       axis.title.face = "bold",
                       axis.text.face = "plain",
                       legend.title.face = "bold",
                       legend.text.face = "plain",
                       strip.text.face = "bold",
                       return_data = FALSE) {
  # Add lengthy error messages.
  withr::local_options(.new = list("warning.length" = 8170))
  
  check_suggests(function_name = "do_BarPlot")
  check_Seurat(sample)

  `%>%` <- magrittr::`%>%`
  `:=` <- rlang::`:=`

  # Check logical parameters.
  logical_list <- list("order" = order,
                       "flip" = flip,
                       "plot.grid" = plot.grid,
                       "legend.byrow" = legend.byrow,
                       "add.n" = add.n,
                       "return_data" = return_data)
  check_type(parameters = logical_list, required_type = "logical", test_function = is.logical)
  # Check numeric parameters.
  numeric_list <- list("font.size" = font.size,
                       "axis.text.x.angle" = axis.text.x.angle,
                       "legend.ncol" = legend.ncol,
                       "legend.nrow" = legend.nrow,
                       "add.n.expand" = add.n.expand)
  check_type(parameters = numeric_list, required_type = "numeric", test_function = is.numeric)
  # Check character parameters.

  character_list <- list("group.by" = group.by,
                         "split.by" = split.by,
                         "facet.by" = facet.by,
                         "order.by" = order.by,
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
                         "grid.type" = grid.type,
                         "legend.title" = legend.title,
                         "plot.title.face" = plot.title.face,
                         "plot.subtitle.face" = plot.subtitle.face,
                         "plot.caption.face" = plot.caption.face,
                         "axis.title.face" = axis.title.face,
                         "axis.text.face" = axis.text.face,
                         "legend.title.face" = legend.title.face,
                         "legend.text.face" = legend.text.face)
  # Checks
  check_type(parameters = character_list, required_type = "character", test_function = is.character)

  check_colors(grid.color, parameter_name = "grid.color")

  check_parameters(parameter = font.type, parameter_name = "font.type")
  check_parameters(parameter = legend.position, parameter_name = "legend.position")
  check_parameters(parameter = grid.type, parameter_name = "grid.type")
  check_parameters(parameter = axis.text.x.angle, parameter_name = "axis.text.x.angle")
  check_parameters(plot.title.face, parameter_name = "plot.title.face")
  check_parameters(plot.subtitle.face, parameter_name = "plot.subtitle.face")
  check_parameters(plot.caption.face, parameter_name = "plot.caption.face")
  check_parameters(axis.title.face, parameter_name = "axis.title.face")
  check_parameters(axis.text.face, parameter_name = "axis.text.face")
  check_parameters(legend.title.face, parameter_name = "legend.title.face")
  check_parameters(legend.text.face, parameter_name = "legend.text.face")
  check_parameters(strip.text.face, parameter_name = "strip.text.face")
  # Get the general table.
  assertthat::assert_that(class(sample@meta.data[, group.by]) %in% c("character", "factor"),
                          msg = paste0(add_cross(), crayon_body("Please provide to "),
                                       crayon_key("feature"),
                                       crayon_body(" a "),
                                       crayon_key(" metadta categorical "),
                                       crayon_body(" variable.")))
  
  assertthat::assert_that(base::isFALSE(position == "fill" & is.null(split.by)),
                          msg = paste0(add_cross(),
                                       crayon_body("Please use "),
                                       crayon_key("position =  fill"),
                                       crayon_body(" alongisde "),
                                       crayon_key("split.by"),
                                       crayon_body(".")))
  
  assertthat::assert_that(base::isFALSE(position == "stack" & isTRUE(order) & !is.null(order.by)),
                          msg = paste0(add_cross(),
                                       crayon_body("Please use "),
                                       crayon_key("order.by"),
                                       crayon_body(" alongisde "),
                                       crayon_key("position = fill"),
                                       crayon_body(".")))
  
  assertthat::assert_that(base::isFALSE(position == "fill" & isTRUE(order) & is.null(order.by)),
                          msg = paste0(add_cross(),
                                       crayon_body("Please use "),
                                       crayon_key("order.by"),
                                       crayon_body(" alongisde "),
                                       crayon_key("position = fill"),
                                       crayon_body(".")))
  
  if (is.null(colors.use)){
    colors.use <- generate_color_scale(names_use = if (is.factor(sample@meta.data[, group.by])) {levels(sample@meta.data[, group.by])} else {sort(unique(sample@meta.data[, group.by]))})
  } else {
    check_colors(colors.use, parameter_name = "colors.use")
    check_consistency_colors_and_names(sample = sample, colors = colors.use, grouping_variable = group.by)
    colors.use <- colors.use[unique(sample@meta.data[, group.by])]
  }
  
  # Check group.by.
  out <- check_group_by(sample = sample,
                        group.by = group.by,
                        is.heatmap = FALSE)
  sample <- out[["sample"]]
  group.by <- out[["group.by"]]
  
  data <-  sample@meta.data %>%
           tibble::as_tibble() %>%
           dplyr::select(dplyr::all_of(c(group.by, split.by, facet.by)))
  
  if (isTRUE(order)){
    if (is.null(order.by)){
      order.use <- sample@meta.data %>%
                   tibble::as_tibble() %>%
                   dplyr::select(dplyr::all_of(c(group.by, split.by))) %>%
                   dplyr::group_by(.data[[group.by]]) %>%
                   dplyr::summarise("n" = dplyr::n()) %>%
                   dplyr::arrange(if(base::isFALSE(flip)){dplyr::desc(.data[["n"]])} else {.data[["n"]]}) %>%
                   dplyr::pull(.data[[group.by]]) %>%
                   as.character()
      data <- data %>% 
              dplyr::mutate("{group.by}" := factor(as.character(.data[[group.by]]), 
                                                   levels = order.use))
    } else {
      order.use <- sample@meta.data %>%
                   tibble::as_tibble() %>%
                   dplyr::mutate(dplyr::across(dplyr::all_of(c(group.by, split.by)), as.character),
                                 "count" = 1) %>% 
                   tidyr::complete(.data[[split.by]], .data[[group.by]], explicit = FALSE) %>% 
                   dplyr::group_by(.data[[split.by]], .data[[group.by]]) %>%
                   dplyr::summarise("n" = sum(.data$count, na.rm = TRUE),
                                    "{split.by}" := unique(.data[[split.by]])) %>% 
                   dplyr::reframe("freq" = .data$n / sum(.data$n),
                                  "{split.by}" := unique(.data[[split.by]]),
                                  "{group.by}" := unique(.data[[group.by]])) %>% 
                   dplyr::filter(.data[[group.by]] == order.by) %>% 
                   dplyr::arrange(if(base::isFALSE(flip)){dplyr::desc(.data[["freq"]])} else {.data[["freq"]]}) %>% 
                   dplyr::pull(.data[[split.by]]) %>% 
                   as.character()
      
      data <- data %>% 
              dplyr::mutate("{split.by}" := factor(as.character(.data[[split.by]]), 
                                                   levels = order.use))
    }
  }
  
  if (isTRUE(add.n)){
    assertthat::assert_that(position == "fill",
                            msg = paste0(add_cross(),
                                         crayon_body("Parameter "),
                                         crayon_key("add.n"),
                                         crayon_body(" can only be used alongside "),
                                         crayon_key("position = fill"),
                                         crayon_body(".")))
    if (is.null(split.by) & !is.null(group.by)){
      data.n <- data %>% 
                dplyr::group_by(.data[[group.by]]) %>% 
                dplyr::summarise(n = dplyr::n()) %>% 
                dplyr::mutate(n = paste0("n = ", .data$n))
    } else if (!is.null(split.by) & !is.null(group.by)){
      data.n <- data %>% 
                dplyr::group_by(.data[[split.by]]) %>% 
                dplyr::summarise(n = dplyr::n()) %>% 
                dplyr::mutate(n = paste0("n = ", .data$n))
    }
    max.char <- max(vapply(data.n$n, nchar, FUN.VALUE = integer(1)))
    data.n$n <- vapply(data.n$n, function(x){return(paste0(x, paste(rep(" ", (max.char - nchar(x))), collapse = "")))}, FUN.VALUE = character(1))
  }
  
          
  if (is.null(split.by)){
    p <- data %>%
         ggplot2::ggplot(mapping = ggplot2::aes(x = .data[[group.by]],
                                                fill = .data[[group.by]]))
  } else {
    p <- data %>%
         ggplot2::ggplot(mapping = ggplot2::aes(x = .data[[split.by]],
                                                fill = .data[[group.by]]))
  }
  
  
  if (is.null(xlab)){
    if (!is.null(group.by) & is.null(split.by)){
      xlab <- group.by
    } else if (!is.null(group.by) & !is.null(split.by)){
      xlab <- split.by
    }
  }
  
  if (is.null(ylab)){
    ylab <- ifelse(position == "stack", "Count", "Proportion")
  }
  
 if (is.null(legend.title)){
   if (position == "stack"){
     if (is.null(split.by)){
       legend.title <- NULL
     } else {
       legend.title <- split.by
     }
   } else {
     legend.title <- group.by
   }
 }

  p <- p +
       ggplot2::stat_count(geom = "bar", position = position, color = "black")
  
  if (isTRUE(add.n)){
    if (is.null(split.by) & !is.null(group.by)){
      p <- p + 
           ggplot2::geom_text(data = data.n,
                              mapping = ggplot2::aes(x = .data[[group.by]],
                                                     y = ifelse(base::isFALSE(flip), 1.03, 1.01),
                                                     label = .data$n,
                                                     fill = NULL),
                              hjust = ifelse(isTRUE(flip), 0, 0.5),
                              fontface = "plain",
                              size = add.n.size)
    } else if (!is.null(split.by) & !is.null(group.by)){
      p <- p + 
           ggplot2::geom_text(data = data.n,
                              mapping = ggplot2::aes(x = .data[[split.by]],
                                                     y = ifelse(base::isFALSE(flip), 1.03, 1.01),
                                                     label = .data$n,
                                                     fill = NULL),
                              hjust = ifelse(isTRUE(flip), 0, 0.5),
                              fontface = "plain",
                              size = add.n.size)
    }
    p <- p + 
         ggplot2::scale_y_continuous(limits = add.n.expand, labels = c("0", "0.25", "0.5", "0.75", "1"), breaks = c(0, 0.25, 0.5, 0.75, 1))

  }
  
  if (isTRUE(flip)){
    p <- p + ggplot2::coord_flip()
  }
  
  if (!is.null(facet.by)){
    if (base::isFALSE(flip)){
      p <- p + 
        ggplot2::facet_grid(cols = ggplot2::vars(.data[[facet.by]]),
                            scales = "free",
                            space = "free",
                            drop = TRUE)
    } else {
      p <- p + 
        ggplot2::facet_grid(rows = ggplot2::vars(.data[[facet.by]]),
                            scales = "free",
                            space = "free",
                            drop = TRUE)
    }
  }
  p <- p +
       ggplot2::xlab(xlab) +
       ggplot2::ylab(ylab) +
       ggplot2::labs(title = plot.title,
                     subtitle = plot.subtitle,
                     caption = plot.caption) +
       ggplot2::scale_fill_manual(values = colors.use) +
       ggplot2::guides(fill = ggplot2::guide_legend(title = legend.title,
                                                    title.position = "top",
                                                    title.hjust = 0.5,
                                                    ncol = legend.ncol,
                                                    nrow = legend.nrow,
                                                    byrow = legend.byrow)) +
       ggplot2::theme_minimal(base_size = font.size) +
       ggplot2::theme(axis.title = ggplot2::element_text(color = "black",
                                                         face = axis.title.face),
                      panel.grid.major.y = if (base::isFALSE(flip)) {if (isTRUE(plot.grid)){ggplot2::element_line(color = grid.color, linetype = grid.type)}} else if (isTRUE(flip)) {ggplot2::element_blank()},
                      panel.grid.major.x = if (isTRUE(flip)) {if (isTRUE(plot.grid)){ggplot2::element_line(color = grid.color, linetype = grid.type)}} else if (base::isFALSE(flip)) {ggplot2::element_blank()},
                      axis.line.x = if (base::isFALSE(flip)) {ggplot2::element_line(color = "black")} else if (isTRUE(flip)) {ggplot2::element_blank()},
                      axis.line.y = if (isTRUE(flip)) {ggplot2::element_line(color = "black")} else if (base::isFALSE(flip)) {ggplot2::element_blank()},
                      axis.text.x = ggplot2::element_text(color = "black",
                                                          face = axis.text.face,
                                                          angle = get_axis_parameters(angle = axis.text.x.angle, flip = flip)[["angle"]],
                                                          hjust = get_axis_parameters(angle = axis.text.x.angle, flip = flip)[["hjust"]],
                                                          vjust = get_axis_parameters(angle = axis.text.x.angle, flip = flip)[["vjust"]]),
                      axis.text.y = ggplot2::element_text(color = "black", face = axis.text.face),
                      axis.ticks = ggplot2::element_line(color = "black"),
                      plot.title.position = "plot",
                      plot.title = ggplot2::element_text(face = plot.title.face, hjust = 0),
                      plot.subtitle = ggplot2::element_text(face = plot.subtitle.face, hjust = 0),
                      plot.caption = ggplot2::element_text(face = plot.caption.face, hjust = 1),
                      panel.grid = ggplot2::element_blank(),
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
                      strip.text = ggplot2::element_text(color = "black", face = strip.text.face),
                      strip.background = ggplot2::element_blank())
 
  if (isTRUE(return_data)){
    data <- sample@meta.data %>%
            tibble::as_tibble() %>%
            dplyr::mutate(dplyr::across(dplyr::all_of(c(group.by, split.by)), as.character),
                          "count" = 1) %>% 
            tidyr::complete(.data[[split.by]], .data[[group.by]], explicit = FALSE) %>% 
            dplyr::group_by(.data[[split.by]], .data[[group.by]]) %>%
            dplyr::summarise("n" = sum(.data$count, na.rm = TRUE),
                             "{split.by}" := unique(.data[[split.by]])) %>% 
            dplyr::reframe("n" = .data$n,
                           "freq" = .data$n / sum(.data$n),
                           "{split.by}" := unique(.data[[split.by]]),
                           "{group.by}" := unique(.data[[group.by]]))
    return(list("Plot" = p,
                "Data" = data))
  } else {
    return(p)
  }
}
