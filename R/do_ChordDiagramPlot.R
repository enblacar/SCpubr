#' Generate a Chord diagram.
#' @inheritParams circlize::chordDiagram
#' @param sample Seurat object.
#' @param first_group Character. Categorical metadata variable. First group to be used as origin of the interaction.
#' @param last_group Character. Categorical metadata variable. Last group to be used as the end of the interaction.
#' @param big.gap Numeric. Space between the groups in first_group and last_group.
#' @param small.gap Numeric. Space within the groups.
#' @param link.border.color Character. Color for the border of the links. NA = no color.
#' @param link.border.width Numeric. Width of the border line of the links.
#' @param highlight_group Character. A value from first_group that will be used to highlight only the links coming from it.
#' @param alpha.highlight Numeric. A value between 00 (double digits) and 99 to depict the alpha of the highlighted links. No transparency needs "FF"
#' @param z_index Logical. Whether to bring the bigger links to the top.
#' @param self.link Numeric. 1 to allow self linking and 2 to prevent them.
#' @param directional Numeric. 0 for non-directional data, 1 for first_group to last_group, -1 for last_group to first_group and 2 for bidirectional.
#' @param direction.type Character. One of diffHeight or arrows or both as a character vector.
#' @param link.arr.type Character. One of triangle or big.arrow.
#' @param scale Logical. Whether to put all nodes the same width.
#' @param alignment Character. One of default, vertical or horizontal. Default will allow circlize to put the angle as it wishes, vertical centers the symmetry on the Y axis and horizontal on the X axis.
#' @param ... Additional parameters to circlize::chordDiagram.
#'
#' @return A circlize plot.
#' @export
#'
#' @examples
#' \dontrun{
#' TBD
#' }
do_ChordDiagramPlot <- function(sample,
                                first_group,
                                last_group,
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
                                directional = 0,
                                direction.type = "diffHeight",
                                link.arr.type = "big.arrow",
                                scale = FALSE,
                                alignment = "default",
                                annotationTrack = c("name", "grid", "axis"),
                                ...){

  # Checks for packages.
  check_suggests(function_name = "do_ChordDiagramPlot")
  # Check if the sample provided is a Seurat object.
  check_Seurat(sample = sample)

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
                       "directional" = directional)
  check_type(parameters = numeric_list, required_type = "numeric", test_function = is.numeric)
  # Check character parameters.

  character_list <- list("first_group" = first_group,
                         "last_group" = last_group,
                         "highlight_group" = highlight_group,
                         "direction.type" = direction.type,
                         "link.arr.type" = link.arr.type,
                         "alignment" = alignment,
                         "annotationTrack" = annotationTrack)
  # Checks
  check_type(parameters = character_list, required_type = "character", test_function = is.character)


  circlize::circos.clear()
  data <- sample@meta.data %>%
          tibble::as_tibble() %>%
          dplyr::select(dplyr::all_of(c(first_group, last_group))) %>%
          dplyr::group_by(.data[[last_group]], .data[[first_group]]) %>%
          dplyr::summarize(value = dplyr::n()) %>%
          dplyr::rename("from" = .data[[first_group]],
                        "to" = .data[[last_group]]) %>%
          dplyr::select(dplyr::all_of(c("from", "to", "value")))

  sectors <- as.character(c(unique(data$from), unique(data$to)))
  colors.use <- SCpubr:::generate_color_scale(names_use = sectors)

  if (is.null(link.sort)){link.sort <- "default"}
  if (isFALSE(z_index)){link.zindex <- NULL} else {link.zindex <- rank(data$value)}
  if (self.link %!in% c(1, 2)){
    stop("Please set self.link as either 1 or 2.", call. = FALSE)
  }

  if (directional %!in% c(0, 1, 2, -1)){
    stop("Please set directional as either 0, 1, 2 or -1.", call. = FALSE)
  }

  for (item in direction.type){
    if (direction.type %!in% c("diffHeight", "arrows")){
      stop("Please set direction.type as either diffHeight, arrows or both.", call. = FALSE)
    }
  }

  if (link.arr.type %!in% c("big.arrow", "triangle")){
    stop("Please set link.arr.type as either big.arrow or triangle.", call. = FALSE)
  }

  if (alignment %!in% c("default", "vertical", "horizontal")){
    stop("Please set alignment as either default or vertical or horizontal.", call. = FALSE)
  } else {
    if (alignment == "vertical"){
      circlize::circos.par(start.degree = 90)
    } else if (alignment == "horizontal"){
      circlize::circos.par(start.degree = 90)
    }
  }

  if (alpha.highlight %!in% c(seq(1, 99), "FF")){
    stop("Please provide either FF or a number between 1 and 99 to alpha.highlight.", call. = FALSE)
  }

  if (!is.na(link.border.color)){
    SCpubr:::check_colors(link.border.color)
  }

  if (!(is.null(highlight_group))){
    colors.use <- stats::setNames(sapply(names(colors.use),
                                         function(x){if (nchar(colors.use[x]) == 7){
                                                       if(x %in% highlight_group) {
                                                         paste0(colors.use[x], "25")
                                                       } else {
                                                         x
                                                       }
                                                     } else if (nchar(colors.use[x]) == 9){
                                                       if(x %in% highlight_group) {
                                                         paste0(stringr::str_sub(colors.use[x], 1, 7), as.character(alpha.highlight))
                                                       } else {
                                                         colors.use[x]
                                                       }
                                                     }
                                  }), names(colors.use))
  }

  col_mat <- ifelse(data$from %in% highlight_group, "#00FF0010", "#FF000080")
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
                         ...)
  p <- grDevices::recordPlot()
  circlize::circos.clear()
  grDevices::dev.off()
  return(p)
}
