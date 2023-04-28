#' Generate a Chord diagram.
#'
#' @inheritParams doc_function
#' @inheritParams circlize::chordDiagram
#' @param from,to \strong{\code{\link[base]{character}}} | Categorical metadata variable to be used as origin and end points of the interactions.
#' @param big.gap \strong{\code{\link[base]{numeric}}} | Space between the groups in "from" and "to".
#' @param small.gap \strong{\code{\link[base]{numeric}}} | Space within the groups.
#' @param link.border.color \strong{\code{\link[base]{character}}} | Color for the border of the links. NA = no color.
#' @param link.border.width \strong{\code{\link[base]{numeric}}} | Width of the border line of the links.
#' @param highlight_group \strong{\code{\link[base]{character}}} | A value from from that will be used to highlight only the links coming from it.
#' @param alpha.highlight \strong{\code{\link[base]{numeric}}} | A value between 00 (double digits) and 99 to depict the alpha of the highlighted links. No transparency needs "FF"
#' @param z_index \strong{\code{\link[base]{logical}}} | Whether to bring the bigger links to the top.
#' @param self.link \strong{\code{\link[base]{numeric}}} | Behavior of the links. One of:
#' \itemize{
#'   \item \emph{\code{1}}: Prevents self linking.
#'   \item \emph{\code{2}}: Allows self linking.
#' }
#' @param directional \strong{\code{\link[base]{numeric}}} | Set the direction of the links. One of:
#' \itemize{
#'   \item \emph{\code{0}}: Non-directional data.
#'   \item \emph{\code{1}}: Links go from "from" to "to".
#'   \item \emph{\code{-1}}: Links go from "to" to "from".
#'   \item \emph{\code{2}}: Links go in both directions.
#' }
#' @param direction.type \strong{\code{\link[base]{character}}} | How to display the directions. One of:
#' \itemize{
#'   \item \emph{\code{diffHeight}}: Sets a line at the origin of the group showing to how many groups and in which proportion this group is linked to.
#'   \item \emph{\code{arrows}}: Sets the connection as arrows.
#'   \item \emph{\code{both}}: Sets up both behaviors. Use as: \code{c("diffHeight", "arrows")}.
#' }
#' @param link.arr.type \strong{\code{\link[base]{character}}} | Sets the appearance of the arrows. One of:
#' \itemize{
#'   \item \emph{\code{triangle}}: Arrow with a triangle tip at the end displayed on top of the link.
#'   \item \emph{\code{big.arrow}}: The link itself ends in a triangle shape.
#' }
#' @param scale \strong{\code{\link[base]{logical}}} | Whether to put all nodes the same width.
#' @param alignment \strong{\code{\link[base]{character}}} | How to align the diagram.  One of:
#' \itemize{
#'   \item \emph{\code{default}}: Allows \pkg{circlize} to set up the plot as it sees fit.
#'   \item \emph{\code{horizontal}}: Sets the break between "from" and "to" groups on the horizontal axis.
#'   \item \emph{\code{vertical}}: Sets the break between "from" and "to" groups on the vertical axis.
#' }
#' @param padding_labels \strong{\code{\link[base]{numeric}}} | Number of extra padding (white spaces) of the labels so that they do not overlap with the scales.
#' @param colors.from,colors.to \strong{\code{\link[SCpubr]{named_vector}}} | Named vector of colors corresponding to the unique values of "from" and "to".
#' @param ... For internal use only.
#'
#' @return A circlize plot.
#' @export
#'
#' @example /man/examples/examples_do_ChordDiagramPlot.R
do_ChordDiagramPlot <- function(sample = NULL,
                                from = NULL,
                                to = NULL,
                                colors.from = NULL,
                                colors.to = NULL,
                                big.gap = 10,
                                small.gap = 1,
                                link.border.color = NA,
                                link.border.width = 1,
                                highlight_group = NULL,
                                alpha.highlight = 25,
                                link.sort = NULL,
                                link.decreasing = TRUE,
                                z_index = FALSE,
                                self.link = 1,
                                symmetric = FALSE,
                                directional = 1,
                                direction.type = c("diffHeight", "arrows"),
                                link.arr.type = "big.arrow",
                                scale = FALSE,
                                alignment = "default",
                                annotationTrack = c("grid", "axis"),
                                padding_labels = 4,
                                ...){
  # Add lengthy error messages.
  withr::local_options(.new = list("warning.length" = 8170))

  `%>%` <- magrittr::`%>%`
  check_suggests(function_name = "do_ChordDiagramPlot")
  # Check logical parameters.
  logical_list <- list("link.decreasing" = link.decreasing,
                       "z_index" = z_index,
                       "symmetric" = symmetric,
                       "scale" = scale)
  check_type(parameters = logical_list, required_type = "logical", test_function = is.logical)
  # Check numeric parameters.
  numeric_list <- list("big.gap" = big.gap,
                       "small.gap" = small.gap,
                       "link.border.width" = link.border.width,
                       "alpha.highlight" = alpha.highlight,
                       "self.link" = self.link,
                       "directional" = directional,
                       "padding_labels" = padding_labels)
  check_type(parameters = numeric_list, required_type = "numeric", test_function = is.numeric)
  # Check character parameters.

  character_list <- list("from" = from,
                         "to" = to,
                         "highlight_group" = highlight_group,
                         "direction.type" = direction.type,
                         "link.arr.type" = link.arr.type,
                         "alignment" = alignment,
                         "annotationTrack" = annotationTrack,
                         "colors.from" = colors.from,
                         "colors.to" = colors.to)
  # Checks
  check_type(parameters = character_list, required_type = "character", test_function = is.character)
  check_parameters(parameter = direction.type, parameter_name = "direction.type")
  check_parameters(parameter = self.link, parameter_name = "self.link")
  check_parameters(parameter = directional, parameter_name = "directional")
  check_parameters(parameter = link.arr.type, parameter_name = "link.arr.type")
  check_parameters(parameter = alignment, parameter_name = "alignment")
  check_parameters(parameter = alpha.highlight, parameter_name = "alpha.highlight")
  circlize::circos.clear()

  extra_params <- list(...)

  # Internal use only
  if ("from_df" %in% names(extra_params) & "df" %in% names(extra_params)){
    data <- extra_params[["df"]]

    assertthat::assert_that("data.frame" %in% class(data),
                            msg = paste0(add_cross(), crayon_body("Please provide a "),
                                         crayon_key("data.frame"),
                                         crayon_body(" or a "),
                                         crayon_key("tibble"),
                                         crayon_body(" to the parameter "),
                                         crayon_key("df"),
                                         crayon_body(".")))

    assertthat::assert_that(ncol(data) == 3,
                            msg = paste0(add_cross(), crayon_body("If you make use of "),
                                         crayon_key("from_df"),
                                         crayon_body(" parameter, make sure that the dataframe provided to "),
                                         crayon_key("df"),
                                         crayon_body(" has at least "),
                                         crayon_key("three"),
                                         crayon_body(" columns.")))

    assertthat::assert_that(sum(colnames(data) %in% c("from", "to", "value")) == 3,
                            msg = paste0(add_cross(), crayon_body("If you make use of "),
                                         crayon_key("from_df"),
                                         crayon_body(" parameter, make sure that the dataframe provided to "),
                                         crayon_key("df"),
                                         crayon_body(" has at least "),
                                         crayon_key("three"),
                                         crayon_body(" columns named: "),
                                         crayon_key("from"),
                                         crayon_body(", "),
                                         crayon_key("to"),
                                         crayon_body(" and "),
                                         crayon_key("value"),
                                         crayon_body(".")))

    assertthat::assert_that(class(data[["from"]]) %in% c("factor", "character"),
                            msg = paste0(add_cross(), crayon_body("Make sure that the column named "),
                                         crayon_key("from"),
                                         crayon_body(" in the dataframe provided to "),
                                         crayon_key("df"),
                                         crayon_body(" is a "),
                                         crayon_key("character"),
                                         crayon_body(" or "),
                                         crayon_key("factor"),
                                         crayon_body(" column.")))

    assertthat::assert_that(class(data[["to"]]) %in% c("factor", "character"),
                            msg = paste0(add_cross(), crayon_body("Make sure that the column named "),
                                         crayon_key("to"),
                                         crayon_body(" in the dataframe provided to "),
                                         crayon_key("df"),
                                         crayon_body(" is a "),
                                         crayon_key("character"),
                                         crayon_body(" or "),
                                         crayon_key("factor"),
                                         crayon_body(" column.")))

    assertthat::assert_that("integer" %in% class(data[["value"]]),
                            msg = paste0(add_cross(), crayon_body("Make sure that the column named "),
                                         crayon_key("from"),
                                         crayon_body(" in the dataframe provided to "),
                                         crayon_key("df"),
                                         crayon_body(" is a "),
                                         crayon_key("numeric"),
                                         crayon_body(" column.")))

  }  else {
    assertthat::assert_that(!is.null(sample),
                            msg = paste0(add_cross(), crayon_body("Please provide a "),
                                         crayon_key("Seurat object"),
                                         crayon_body(" to the parameter "),
                                         crayon_key("sample"),
                                         crayon_body(".")))

    # Check if the sample provided is a Seurat object.
    check_Seurat(sample = sample)

    assertthat::assert_that(!is.null(from) | !is.null(to),
                            msg = paste0(add_cross(), crayon_body("Please provide a value to "),
                                         crayon_key("from"),
                                         crayon_body(" or "),
                                         crayon_key("to"),
                                         crayon_body(" parameters.")))

    assertthat::assert_that(from %in% colnames(sample@meta.data) | to %in% colnames(sample@meta.data),
                            msg = paste0(add_cross(), crayon_body("Parameters "),
                                         crayon_key("from"),
                                         crayon_body(" and "),
                                         crayon_key("to"),
                                         crayon_body(" need to be present in the object "),
                                         crayon_key("metadata"),
                                         crayon_body(".")))

    assertthat::assert_that(class(sample@meta.data[, from]) %in% c("factor", "character") | class(sample@meta.data[, to]) %in% c("factor", "character"),
                            msg = paste0(add_cross(), crayon_body("Parameters "),
                                         crayon_key("from"),
                                         crayon_body(" and "),
                                         crayon_key("to"),
                                         crayon_body(" need to be either a "),
                                         crayon_key("factor"),
                                         crayon_body(" or "),
                                         crayon_key("character"),
                                         crayon_body(" columns.")))

    data <- sample@meta.data %>%
            tibble::as_tibble() %>%
            dplyr::select(dplyr::all_of(c(from, to))) %>%
            dplyr::group_by(.data[[to]], .data[[from]]) %>%
            dplyr::summarize(value = dplyr::n()) %>%
            dplyr::rename("from" = dplyr::all_of(c(from)),
                          "to" = dplyr::all_of(c(to))) %>%
            dplyr::select(dplyr::all_of(c("from", "to", "value")))
  }


  max_char <- max(c(max(nchar(as.character(data[["from"]]))), max(nchar(as.character(data[["to"]]))))) + padding_labels

  if (is.factor(data[["to"]]) & is.factor(data[["from"]])){
    levels_to <- stringr::str_pad(levels(data[["to"]]), width = max_char, side = "both")
    levels_from <- stringr::str_pad(levels(data[["from"]]), width = max_char, side = "both")

    data <- data %>%
            dplyr::mutate("from" = factor(stringr::str_pad(.data[["from"]], width = max_char, side = "both"), levels = levels_from),
                          "to" = factor(stringr::str_pad(.data[["to"]], width = max_char, side = "both"), levels = levels_to))
  } else if (is.factor(data[["to"]]) & is.character(data[["from"]])){
    levels_to <- stringr::str_pad(levels(data[["to"]]), width = max_char, side = "both")
    data <- data %>%
            dplyr::mutate("from" = stringr::str_pad(.data[["from"]], width = max_char, side = "both"),
                          "to" = factor(stringr::str_pad(.data[["to"]], width = max_char, side = "both"), levels = levels_to))

  } else if (is.character(data[["to"]]) & is.factor(data[["from"]])){
    levels_from <- stringr::str_pad(levels(data[["from"]]), width = max_char, side = "both")

    data <- data %>%
            dplyr::mutate("from" = factor(stringr::str_pad(.data[["from"]], width = max_char, side = "both"), levels = levels_from),
                          "to" = stringr::str_pad(.data[["to"]], width = max_char, side = "both"))
  } else if (is.character(data[["to"]]) & is.character(data[["from"]])){
    data <- data %>%
            dplyr::mutate("from" = stringr::str_pad(.data[["from"]], width = max_char, side = "both"),
                          "to" = stringr::str_pad(.data[["to"]], width = max_char, side = "both"))
  }


  if (!(is.null(colors.from))){
    check_colors(colors.from, parameter_name = "colors.from")
    check_consistency_colors_and_names(sample = sample,
                                       colors = colors.from,
                                       grouping_variable = from)
  } else {
    if (is.factor(data[["from"]])){
      colors.from <- generate_color_scale(names_use = levels(data[["from"]]))
    } else {
      colors.from <- generate_color_scale(names_use = sort(unique(data[["from"]])))
    }
  }
  names(colors.from) <- stringr::str_pad(names(colors.from), width = max_char, side = "both")

  if (!(is.null(colors.to))){
    check_colors(colors.to, parameter_name = "colors.to")
    check_consistency_colors_and_names(sample = sample,
                                       colors = colors.to,
                                       grouping_variable = to)
  } else {
    colors.to <- viridis::viridis(n = length(unique(data[["to"]])), option = "G")
    if (is.factor(data[["to"]])){
      colors.to <- stats::setNames(colors.to, levels(data[["to"]]))
    } else {
      colors.to <- stats::setNames(colors.to, sort(unique(data[["to"]])))
    }
  }
  names(colors.to) <- stringr::str_pad(names(colors.to), width = max_char, side = "both")
  colors.use <- c(colors.from, colors.to)
  if (is.null(link.sort)){link.sort <- "default"}
  if (base::isFALSE(z_index)){link.zindex <- NULL} else {link.zindex <- rank(data[["value"]])}

  if (alignment == "vertical"){
    circlize::circos.par(start.degree = 0)
  } else if (alignment == "horizontal"){
    circlize::circos.par(start.degree = 90)
  }

  if (!is.na(link.border.color)){
    check_colors(link.border.color)
  }

  if (!(is.null(highlight_group))){
    alpha.colors <- NULL
    highlight_group <- stringr::str_pad(highlight_group, width = max_char, side = "both")
    for (color in names(colors.use)){
      name <- color
      color <- colors.use[name]
      if (nchar(color) == 7){
        if (name %!in% highlight_group){
          color <- paste0(color, alpha.highlight)
          names(color) <- name
        } else {
          names(color) <- name
        }
      } else if (nchar(color) == 9){
        if (name %!in% highlight_group){
          color <- paste0(stringr::str_sub(color, 1, 7), as.character(alpha.highlight))
          names(color) <- name
        } else {
          names(color) <- name
        }
      }
      alpha.colors <- append(alpha.colors, color)
    }
    colors.use <- alpha.colors
    rm(alpha.colors)
  }
  circlize::chordDiagram(data,
                         big.gap = big.gap,
                         small.gap = small.gap,
                         grid.col = colors.use,
                         link.border = link.border.color,
                         link.lwd = link.border.width,
                         link.sort = link.sort,
                         link.decreasing = link.decreasing,
                         link.zindex = link.zindex,
                         self.link = self.link,
                         symmetric = symmetric,
                         directional = directional,
                         direction.type = direction.type,
                         link.arr.type = link.arr.type,
                         scale = scale,
                         annotationTrack = annotationTrack,
                         preAllocateTracks = list(track.height = max(graphics::strwidth(unlist(dimnames(data))))))
  circlize::circos.track(track.index = 1,
                         panel.fun = function(x, y){circlize::circos.text(circlize::CELL_META$xcenter,
                                                                          circlize::CELL_META$ylim[[1]],
                                                                          circlize::CELL_META$sector.index,
                                                                          facing = "clockwise",
                                                                          niceFacing = TRUE,
                                                                          adj = c(-0.15, 0.5),
                                                                          font = 2)},
                         bg.border = NA)
  p <- grDevices::recordPlot()
  circlize::circos.clear()
  grDevices::dev.off()

  return(p)
}
