#' Generate a Chord diagram.
#' @inheritParams circlize::chordDiagram
#' @param sample Seurat object.
#' @param from Character. Categorical metadata variable. First group to be used as origin of the interaction.
#' @param to Character. Categorical metadata variable. Last group to be used as the end of the interaction.
#' @param big.gap Numeric. Space between the groups in from and to.
#' @param small.gap Numeric. Space within the groups.
#' @param link.border.color Character. Color for the border of the links. NA = no color.
#' @param link.border.width Numeric. Width of the border line of the links.
#' @param highlight_group Character. A value from from that will be used to highlight only the links coming from it.
#' @param alpha.highlight Numeric. A value between 00 (double digits) and 99 to depict the alpha of the highlighted links. No transparency needs "FF"
#' @param z_index Logical. Whether to bring the bigger links to the top.
#' @param self.link Numeric. 1 to allow self linking and 2 to prevent them.
#' @param directional Numeric. 0 for non-directional data, 1 for from to to, -1 for to to from and 2 for bidirectional.
#' @param direction.type Character. One of diffHeight or arrows or both as a character vector.
#' @param link.arr.type Character. One of triangle or big.arrow.
#' @param scale Logical. Whether to put all nodes the same width.
#' @param alignment Character. One of default, vertical or horizontal. Default will allow circlize to put the angle as it wishes, vertical centers the symmetry on the Y axis and horizontal on the X axis.
#' @param padding_labels Numeric. Number of extra padding of the labels so that they do not overlap with the scales.
#' @param colors.from,colors.to Named vector of colors corresponding to the unique values of "from" and "to".
#' @param ... For internal use only.
#'
#' @return A circlize plot.
#' @export
#'
#' @examples
#' \dontrun{
#' TBD
#' }
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

  # Checks for packages.
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

  circlize::circos.clear()

  extra_params = list(...)

  # Internal use only
  if ("from_df" %in% names(extra_params) & "df" %in% names(extra_params)){
    data <- extra_params$df
    if ("data.frame" %!in% class(data)){
      stop("Please provide a data.frame or tibble to df.", call. = FALSE)
    }
    if (ncol(data) != 3){
      stop("If you pass from_df and df to this function, make sure that the df has 3 columns.", call. = FALSE)
    } else {
      if (sum(colnames(data) %in% c("from", "to", "value")) != 3){
        stop("Please name the columns in the df as: from, to and value.", call. = FALSE)
      } else {
        if (class(data$from) %!in% c("factor", "character")){
          stop("Please make sure that the column from is either a factor or a character column.", call. = FALSE)
        }

        if (class(data$to) %!in% c("factor", "character")){
          stop("Please make sure that the column to is either a factor or a character column.", call. = FALSE)
        }

        if (class(data$value) %!in% c("integer")){
          stop("Please make sure that the column value is either an integer column.", call. = FALSE)
        }
      }
    }
  }  else {
    if (is.null(sample)){
      stop("Please provide a Seurat object.", call. = FALSE)
    } else {
      # Check if the sample provided is a Seurat object.
      check_Seurat(sample = sample)
    }

    if (is.null(from) | is.null(to)){
      stop("Please provide a value to: from or to parameters.", call. = FALSE)
    } else {
      if (from %!in% colnames(sample@meta.data) | to %!in% colnames(sample@meta.data)){
        stop("From or to parameters have to be in the object metadata.", call. = FALSE)
      }

      if (class(sample@meta.data[, from]) %!in% c("factor", "character") | class(sample@meta.data[, to]) %!in% c("factor", "character")){
        stop("Parameters from or to have to be either a factor or a character column in the object metadata.", call. = FALSE)
      }
    }
    data <- sample@meta.data %>%
            tibble::as_tibble() %>%
            dplyr::select(dplyr::all_of(c(from, to))) %>%
            dplyr::group_by(.data[[to]], .data[[from]]) %>%
            dplyr::summarize(value = dplyr::n()) %>%
            dplyr::rename("from" = .data[[from]],
                          "to" = .data[[to]]) %>%
            dplyr::select(dplyr::all_of(c("from", "to", "value")))
  }


  max_char <- max(c(max(nchar(as.character(data$from))), max(nchar(as.character(data$to))))) + padding_labels
  data <- data %>%
          dplyr::mutate("from" = stringr::str_pad(.data$from, width = max_char, side = "both"),
                        "to" = stringr::str_pad(.data$to, width = max_char, side = "both"))

  if (!(is.null(colors.from))){
    SCpubr:::check_colors(colors.from, parameter_name = "colors.from")
    SCpubr:::check_consistency_colors_and_names(sample = sample,
                                                colors = colors.from,
                                                grouping_variable = from)
  } else {
    colors.from <- SCpubr:::generate_color_scale(names_use = if(is.factor(data$from)){levels(data$from)} else {sort(unique(data$from))})
  }
  names(colors.from) <- stringr::str_pad(names(colors.from), width = max_char, side = "both")

  if (!(is.null(colors.to))){
    SCpubr:::check_colors(colors.to, parameter_name = "colors.to")
    SCpubr:::check_consistency_colors_and_names(sample = sample,
                                                colors = colors.to,
                                                grouping_variable = to)
  } else {
    colors.to <- viridis::viridis(n = length(unique(data$to)), option = "G")
    colors.to <- stats::setNames(colors.to, if(is.factor(data$to)){levels(data$to)} else {sort(unique(data$to))})
  }
  names(colors.to) <- stringr::str_pad(names(colors.to), width = max_char, side = "both")
  colors.use <- c(colors.from, colors.to)
  if (is.null(link.sort)){link.sort <- "default"}
  if (isFALSE(z_index)){link.zindex <- NULL} else {link.zindex <- rank(data$value)}
  if (self.link %!in% c(1, 2)){
    stop("Please set self.link as either 1 or 2.", call. = FALSE)
  }

  if (directional %!in% c(0, 1, 2, -1)){
    stop("Please set directional as either 0, 1, 2 or -1.", call. = FALSE)
  }

  for (item in direction.type){
    if (item %!in% c("diffHeight", "arrows")){
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
    alpha.colors <- c()
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
                         preAllocateTracks = list(track.height = max(strwidth(unlist(dimnames(data))))))
  circlize::circos.track(track.index = 1,
                         panel.fun = function(x, y){circlize::circos.text(circlize::CELL_META$xcenter,
                                                                          circlize::CELL_META$ylim[1],
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
