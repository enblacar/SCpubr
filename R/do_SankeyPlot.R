
#' Do Sankey or Alluvial plots.
#'
#' @inheritParams doc_function
#' @param first_group \strong{\code{\link[base]{character}}} | Categorical metadata variable. First group of nodes of the sankey plot.
#' @param last_group \strong{\code{\link[base]{character}}} | Categorical metadata variable. Last group of nodes of the sankey plot.
#' @param middle_groups \strong{\code{\link[base]{character}}} | Categorical metadata variable. Vector of groups of nodes of the sankey plot.
#' @param type \strong{\code{\link[base]{character}}} | Type of plot to make. One of:
#' \itemize{
#'   \item \emph{\code{sankey}}: Generates a sankey plot.
#'   \item \emph{\code{alluvial}}: Generated an Alluvial plot, a kind of sankey plot where all groups have the same height.
#' }
#' @param width \strong{\code{\link[base]{numeric}}} | Width of the nodes.
#' @param space \strong{\code{\link[base]{numeric}}} | Vertical space between the nodes. It appears to be equal to a single cell. Use big numbers to see a difference (like, 1000 or 10000).
#' @param position \strong{\code{\link[base]{character}}} | GGplot2 position.
#' @param node.fill \strong{\code{\link[base]{character}}} | Color to fill the nodes.
#' @param node.color \strong{\code{\link[base]{character}}} | Color for the contour of the nodes.
#' @param flow.alpha \strong{\code{\link[base]{character}}} | Alpha of the connections.
#' @param flow.color \strong{\code{\link[base]{character}}} | Color for the contour of the connections.
#' @param text_size \strong{\code{\link[base]{numeric}}} | Size of the labels.
#' @param text_color \strong{\code{\link[base]{character}}} | Color of the labels.
#' @param smooth \strong{\code{\link[base]{numeric}}} | How smooth the connections are.
#' @param colors.first,colors.middle,colors.last \strong{\code{\link[base]{character}}} | Named vector of colors equal to ALL unique values in first_group, middle_groups, or last_group.
#' @param use_labels \strong{\code{\link[base]{logical}}} | Whether to use labels or text for the node names.
#' @param hjust \strong{\code{\link[base]{numeric}}} | General hjust for the labels.
#'
#' @return A ggplot2 object.
#' @export
#'
#' @example /man/examples/examples_do_SankeyPlot.R
do_SankeyPlot <- function(sample,
                          first_group,
                          last_group,
                          type = "sankey",
                          middle_groups = NULL,
                          width = 0.1,
                          space = ifelse(type == "sankey", 0.05 * ncol(sample), 0),
                          position = "identity",
                          node.fill = "white",
                          node.color = "white",
                          flow.alpha = 0.75,
                          flow.color = "black",
                          text_size = 3,
                          text_color = "black",
                          font.size = 14,
                          font.type = "sans",
                          smooth = 8,
                          use_labels = FALSE,
                          hjust = NULL,
                          colors.first = NULL,
                          colors.middle = NULL,
                          colors.last = NULL,
                          plot.title = NULL,
                          plot.subtitle = NULL,
                          plot.caption = NULL){

  # Checks for packages.
  check_suggests(function_name = "do_SankeyPlot")
  # Check if the sample provided is a Seurat object.
  check_Seurat(sample = sample)

  # Check logical parameters.
  logical_list <- list("use_labels" = use_labels)
  check_type(parameters = logical_list, required_type = "logical", test_function = is.logical)
  # Check numeric parameters.
  numeric_list <- list("width" = width,
                       "space" = space,
                       "flow.alpha" = flow.alpha,
                       "text_size" = text_size,
                       "font.size" = font.size,
                       "smooth" = smooth,
                       "hjust" = hjust)
  check_type(parameters = numeric_list, required_type = "numeric", test_function = is.numeric)
  # Check character parameters.

  character_list <- list("first_group" = first_group,
                         "last_group" = last_group,
                         "middle_groups" = middle_groups,
                         "type" = type,
                         "position" = position,
                         "node.color" = node.color,
                         "flow.color" = flow.color,
                         "text_color" = text_color,
                         "font.type" = font.type,
                         "colors.first" = colors.first,
                         "colors.middle" = colors.middle,
                         "colors.last" = colors.last,
                         "node.fill" = node.fill,
                         "plot.title" = plot.title,
                         "plot.subtitle" = plot.subtitle,
                         "plot.caption" = plot.caption)
  # Checks
  check_type(parameters = character_list, required_type = "character", test_function = is.character)

  check_colors(node.color, parameter_name = "node.color")
  check_colors(flow.color, parameter_name = "flow.color")
  check_colors(text_color, parameter_name = "text_color")

  check_parameters(parameter = font.type, parameter_name = "font.type")


  # Workaround for ggsankey bug.
  if("dplyr" %!in% (.packages())){
    stop("Due to an internal bug in ggsankey, package dplyr must be loaded using `library(dplyr)`. This behaviour will be corrected if this bug disappears. Thanks for your understanding.", call. = FALSE)
  }

  # Wrong type.
  if (type %!in% c("alluvial", "sankey")){
    stop("Please provide either sankey or alluvial to type.", call. = FALSE)
  }

  # Wrong position.
  if (position %!in% c("identity", "nudge")){
    stop("This position type has not been tested.", call. = FALSE)
  }

  # Not a metadata column.
  if (first_group %!in% colnames(sample@meta.data)){
    stop("The metadata variable for first_group is not in the metadata of the object.", call. = FALSE)
  } else {
    if (!(isTRUE(is.character(sample@meta.data[, first_group]) | is.factor(sample@meta.data[, first_group])))){
      stop("The metadata variable for first_group has to be either a character vector or a factor.", call. = FALSE)
    }
  }

  if (last_group %!in% colnames(sample@meta.data)){
    stop("The metadata variable for last_group is not in the metadata of the object.", call. = FALSE)
  } else {
    if (!(isTRUE(is.character(sample@meta.data[, last_group]) | is.factor(sample@meta.data[, last_group])))){
      stop("The metadata variable for last_group has to be either a character vector or a factor.", call. = FALSE)
    }
  }

  for (var in middle_groups){
    if (var %!in% colnames(sample@meta.data)){
      stop("The metadata variable for middle_groups is not in the metadata of the object.", call. = FALSE)
    } else {
      if (!(isTRUE(is.character(sample@meta.data[, var]) | is.factor(sample@meta.data[, var])))){
        stop("The metadata variable for middle_groups has to be either a character vector or a factor.", call. = FALSE)
      }
    }
  }


  `%>%` <- magrittr::`%>%`

  data <- suppressWarnings({sample@meta.data %>%
                            dplyr::select(dplyr::all_of(c(first_group, middle_groups, last_group))) %>%
                            tibble::rownames_to_column(var = "cell") %>%
                            dplyr::select(-.data$cell) %>%
                            ggsankey::make_long(dplyr::all_of(c(first_group, middle_groups, last_group))) %>%
                            dplyr::rowwise() %>%
                            dplyr::mutate(hjust = if(.data$x %in% middle_groups){0.5}
                                          else if (.data$x == last_group){0}
                                          else if (.data$x == first_group){1})})
  if (!is.null(hjust)){data$hjust <- hjust}

  if (!(is.null(colors.first))){
    check_colors(colors.first, parameter_name = "colors.first")
    if (sum(names(colors.first) %!in% unique(sample@meta.data[, first_group])) > 0){
      stop("Not all colors provided for the first group match the unique values for first_group.", call. = FALSE)
    }

    if (length(colors.first) != length(unique(sample@meta.data[, first_group]))){
      stop("The colors provided for the first group do not match the number of unique nodes.", call. = FALSE)
    }
  } else {
    colors.first <- viridis::viridis(n = length(unique(sample@meta.data[, first_group])), option = "G")
    names(colors.first) <- if(is.factor(sample@meta.data[, first_group])) {levels(sample@meta.data[, first_group])} else {sort(unique(sample@meta.data[, first_group]))}
  }

  if (!(is.null(colors.last))){
    check_colors(colors.last, parameter_name = "colors.last")
    if (sum(names(colors.last) %!in% unique(sample@meta.data[, last_group])) > 0){
      stop("Not all colors provided for the last group match the unique values for last_group", call. = FALSE)
    }

    if (length(colors.last) != length(unique(sample@meta.data[, last_group]))){
      stop("The colors provided for the last group do not match the number of unique nodes.", call. = FALSE)
    }
  } else{
    colors.last <- viridis::viridis(n = length(unique(sample@meta.data[, last_group])), option = "D")
    names(colors.last) <- if(is.factor(sample@meta.data[, last_group])) {levels(sample@meta.data[, last_group])} else {sort(unique(sample@meta.data[, last_group]))}
  }

  if (!(is.null(colors.middle))){
    check_colors(colors.middle, parameter_name = "colors.middle")

    unique_middle_values <- c()
    for(var in middle_groups){
      unique_middle_values <- c(unique_middle_values, if(is.factor(sample@meta.data[, var])) {levels(sample@meta.data[, var])} else {sort(unique(sample@meta.data[, var]))})
    }

    if (sum(names(colors.middle) %!in% unique_middle_values) > 0){
      stop("Not all colors provided for the middle groups match the unique values for middle_groups", call. = FALSE)
    }

    if (length(colors.middle) != length(unique_middle_values)){
      stop("The colors provided for the middle groups do not match the number of unique nodes.", call. = FALSE)
    }
  } else {
    unique_middle_values <- c()
    for(var in middle_groups){
      unique_middle_values <- c(unique_middle_values, if(is.factor(sample@meta.data[, var])) {levels(sample@meta.data[, var])} else {sort(unique(sample@meta.data[, var]))})
    }

    colors.middle <- viridis::viridis(n = length(unique_middle_values), option = "C")
    names(colors.middle) <- unique_middle_values
  }

  colors.use <- c(colors.first, colors.middle, colors.last)
  func_use <- ifelse(isTRUE(use_labels), ggsankey::geom_sankey_label, ggsankey::geom_sankey_text)

  p <- data %>%

       ggplot2::ggplot(mapping = ggplot2::aes(x = .data$x,
                                              next_x = .data$next_x,
                                              node = .data$node,
                                              next_node = .data$next_node,
                                              fill = factor(.data$node),
                                              label = .data$node,
                                              hjust = .data$hjust)) +
       ggsankey::geom_sankey(flow.alpha = flow.alpha,
                             node.color = node.color,
                             node.fill = node.fill,
                             color = flow.color,
                             width = width,
                             position = position,
                             type = type,
                             space = space) +
       func_use(size = text_size,
                color = text_color,
                fontface = "bold",
                position = position,
                type = type,
                space = space) +
       ggplot2::scale_fill_manual(values = colors.use) +
       ggplot2::xlab("") +
       ggplot2::ylab("") +
       ggplot2::labs(title = plot.title,
                     subtitle = plot.subtitle,
                     caption = plot.caption) +
       ggplot2::theme_minimal(base_size = font.size) +
       ggplot2::theme(axis.title = ggplot2::element_text(color = "black",
                                                         face = "bold"),
                      axis.line.x = ggplot2::element_blank(),
                      axis.text.x = ggplot2::element_text(color = "black",
                                                          face = "bold",
                                                          angle = 0,
                                                          hjust = 0.5,
                                                          vjust = 1),
                      axis.text.x.top = ggplot2::element_text(color = "black",
                                                              face = "bold",
                                                              angle = 0,
                                                              hjust = 0.5,
                                                              vjust = 1),
                      axis.text.y = ggplot2::element_blank(),
                      axis.ticks = ggplot2::element_blank(),
                      panel.grid.major = ggplot2::element_blank(),
                      plot.title.position = "plot",
                      plot.title = ggplot2::element_text(face = "bold", hjust = 0),
                      plot.subtitle = ggplot2::element_text(hjust = 0),
                      plot.caption = ggplot2::element_text(hjust = 1),
                      panel.grid = ggplot2::element_blank(),
                      text = ggplot2::element_text(family = font.type),
                      plot.caption.position = "plot",
                      legend.text = ggplot2::element_text(face = "bold"),
                      legend.position = "none",
                      legend.title = ggplot2::element_text(face = "bold"),
                      legend.justification = "center",
                      plot.margin = ggplot2::margin(t = 10, r = 10, b = 10, l = 10),
                      plot.background = ggplot2::element_rect(fill = "white", color = "white"),
                      panel.background = ggplot2::element_rect(fill = "white", color = "white"),
                      legend.background = ggplot2::element_rect(fill = "white", color = "white"),
                      strip.text =ggplot2::element_text(color = "black", face = "bold"))

  return(p)
}
