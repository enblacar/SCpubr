#' Generate Waffle plots to display cell group proportions.
#'
#' This function displays the proportional composition of cell groups as a
#' waffle chart, where each tile represents a fixed number of cells.
#'
#' @inheritParams doc_function
#' @param waffle.size \strong{\code{\link[base]{numeric}}} | Tile border size.
#' 
#' @return A ggplot2 object with a Waffle Plot.
#' @export
#'
#' @example man/examples/examples_do_WafflePlot.R
do_WafflePlot <- function(sample,
                          group.by,
                          waffle.size = 2,
                          flip = FALSE,
                          colors.use = NULL,
                          colorblind = FALSE,
                          na.value = "grey75",
                          font.size = 14,
                          font.type = "sans",
                          plot.title = NULL,
                          plot.subtitle = NULL,
                          plot.caption = NULL,
                          legend.title = NULL,
                          legend.ncol = NULL,
                          legend.nrow = NULL,
                          legend.byrow = FALSE,
                          legend.position = "bottom",
                          plot.title.face = "bold",
                          plot.subtitle.face = "plain",
                          plot.caption.face = "italic",
                          axis.title.face = "bold",
                          axis.text.face = "plain",
                          legend.title.face = "bold",
                          legend.text.face = "plain",
                          strip.text.face = "bold"){
  # Add lengthy error messages.
  withr::local_options(.new = list("warning.length" = 8170))

  check_suggests(function_name = "do_WafflePlot")
  check_Seurat(sample)
  
  # Define pipe operator internally.
  `%>%` <- magrittr::`%>%`
  
  # Check logical parameters
  logical_list <- list("flip" = flip,
                       "legend.byrow" = legend.byrow,
                       "colorblind" = colorblind)
  check_type(parameters = logical_list, required_type = "logical", test_function = is.logical)
  # Check numeric parameters.
  numeric_list <- list("waffle.size" = waffle.size,
                       "font.size" = font.size)
  check_type(parameters = numeric_list, required_type = "numeric", test_function = is.numeric)
  
  # Check character parameters.
  character_list <- list("group.by" = group.by,
                         "na.value" = na.value,
                         "font.type" = font.type,
                         "plot.title" = plot.title,
                         "plot.suybtitle" = plot.subtitle,
                         "plot.caption" = plot.caption,
                         "legend.title" = legend.title,
                         "plot.title.face" = plot.title.face,
                         "plot.subtitle.face" = plot.subtitle.face,
                         "plot.caption.face" = plot.caption.face,
                         "axis.title.face" = axis.title.face,
                         "axis.text.face" = axis.text.face,
                         "legend.title.face" = legend.title.face,
                         "legend.text.face" = legend.text.face,
                         "legend.position" = legend.position)
  check_type(parameters = character_list, required_type = "character", test_function = is.character)
  
  # Check group.by.
  out <- check_group_by(sample = sample,
                        group.by = group.by,
                        is.heatmap = FALSE)
  sample <- out[["sample"]]
  group.by <- out[["group.by"]]
  
  
  check_parameters(parameter = font.type, parameter_name = "font.type")
  check_parameters(parameter = legend.position, parameter_name = "legend.position")
  check_parameters(plot.title.face, parameter_name = "plot.title.face")
  check_parameters(plot.subtitle.face, parameter_name = "plot.subtitle.face")
  check_parameters(plot.caption.face, parameter_name = "plot.caption.face")
  check_parameters(axis.title.face, parameter_name = "axis.title.face")
  check_parameters(axis.text.face, parameter_name = "axis.text.face")
  check_parameters(legend.title.face, parameter_name = "legend.title.face")
  check_parameters(legend.text.face, parameter_name = "legend.text.face")
  
  

  # Check the colors provided.
  check_colors(na.value, parameter_name = "na.value")
  
  if (is.null(colors.use)){
    if (is.factor(sample@meta.data[, group.by])){
      names.use <- levels(sample@meta.data[, group.by])
    } else {
      names.use <- sort(unique(sample@meta.data[, group.by]))
    }
    
    colors.use <- generate_color_scale(names_use = names.use, colorblind = colorblind)
  } else {
    check_colors(colors.use, parameter_name = "colors.use")
    check_consistency_colors_and_names(sample = sample, colors = colors.use, grouping_variable = group.by)
    if (is.factor(sample@meta.data[, group.by])){
      colors.use <- colors.use[levels(sample@meta.data[, group.by])]
    } else {
      colors.use <- colors.use[sort(unique(sample@meta.data[, group.by]))]
    }
  }

  # Get data
  data <- sample@meta.data %>% 
          tibble::rownames_to_column(var = "cell") %>% 
          tibble::as_tibble() %>% 
          dplyr::select(dplyr::all_of(c("cell", group.by))) %>% 
          dplyr::group_by(.data[[group.by]]) %>% 
          dplyr::summarise("n" = dplyr::n()) %>% 
          dplyr::mutate("freq" = (.data$n / sum(.data$n)) * 100,
                        "Groups" = .data[[group.by]])
  
  # Add rounded percentages.
  data$Totals <- round_percent(x = data,
                               group.by = group.by)
  
  # Build waffle grid manually to avoid depending on the 'waffle' package.
  # Expand each group into tiles summing to 100 (10x10 grid).
  tile_df <- data %>%
    dplyr::mutate(tiles = .data$Totals) %>%
    dplyr::select(.data$Groups, .data$tiles) %>%
    tidyr::uncount(weights = .data$tiles, .remove = FALSE) %>%
    dplyr::group_by(.data$Groups) %>%
    dplyr::mutate(idx = dplyr::row_number()) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(global_idx = dplyr::row_number())

  # Calculate x/y positions for a 10x10 grid.
  tile_df <- tile_df %>%
    dplyr::mutate(x = (global_idx - 1) %% 10 + 1,
                  y = 10 - ((global_idx - 1) %/% 10))

  # If flip requested, swap axes.
  if (isTRUE(flip)){
    tile_df <- tile_df %>% dplyr::mutate(tmp = x, x = y, y = tmp) %>% dplyr::select(-tmp)
  }

  p <- tile_df %>%
    ggplot2::ggplot(ggplot2::aes(x = .data$x, y = .data$y, fill = .data$Groups)) +
    ggplot2::geom_tile(color = "white", linewidth = waffle.size) +
    ggplot2::geom_tile(color = "black", linewidth = 0.35, alpha = 0.25, fill = NA) +
    ggplot2::coord_fixed() +
    ggplot2::scale_fill_manual(values = colors.use, na.value = na.value) +
       # Add plot labels.
       ggplot2::labs(title = plot.title,
                     subtitle = plot.subtitle,
                     caption = ifelse(is.null(plot.caption), paste0("Grid: 100 tiles | Each: 1% | Cells: ", sum(data$n)), plot.caption)) + 
       # Customise legend.
       ggplot2::guides(fill = ggplot2::guide_legend(title = ifelse(is.null(legend.title), group.by, legend.title),
                                                    title.position = "top",
                                                    title.hjust = 0.5,
                                                    ncol = legend.ncol,
                                                    nrow = legend.nrow,
                                                    byrow = legend.byrow))  +
      # Add theme.
      ggplot2::theme_minimal(base_size = font.size) +
       # Customise theme.
      ggplot2::theme(axis.title = ggplot2::element_blank(),
                      panel.grid = ggplot2::element_blank(),
                      axis.line = ggplot2::element_blank(),
                      axis.text = ggplot2::element_blank(),
                      axis.ticks = ggplot2::element_blank(),
                      plot.title.position = "plot",
                      plot.title = ggplot2::element_text(face = plot.title.face, hjust = 0),
                      plot.subtitle = ggplot2::element_text(face = plot.subtitle.face, hjust = 0),
                      plot.caption = ggplot2::element_text(face = plot.caption.face, hjust = 1),
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
    
  # Return the plot.
  return(p)
}

