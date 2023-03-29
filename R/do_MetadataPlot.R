#' Compute a heatmap of categorical variables.
#' 
#' The main use of this function is to generate a metadata heatmap of your categorical data, 
#' normally targeted to the different patient samples one has in the Seurat object. It requires
#' that the metadata columns chosen have one and only one possible value for each of the values in 
#' group.by.
#'
#' @inheritParams doc_function
#' @param group.by \strong{\code{\link[base]{character}}} | Metadata column to use as basis for the plot.
#' @param metadata \strong{\code{\link[base]{character}}} | Metadata columns that will be used to plot the heatmap on the basis of the variable provided to group.by.
#' @param colors.use \strong{\code{\link[base]{list}}} | A named list of named vectors. The names of the list correspond to the names of the values provided to metadata and the names of the items in the named vectors correspond to the unique values of that specific metadata variable. The values are the desired colors in HEX code for the values to plot. The used are pre-defined by the pacakge but, in order to get the most out of the plot, please provide your custom set of colors for each metadata column! 
#' @param heatmap.gap \strong{\code{\link[base]{numeric}}} | Size of the gap between heatmaps in mm.
#' @param from_df \strong{\code{\link[base]{logical}}} | Whether to provide a data frame with the metadata instead.
#' @param df \strong{\code{\link[base]{data.frame}}} | Data frame containing the metadata to plot. Rows contain the unique values common to all columns (metadata variables). The columns must be named.
#' @param legend.font.size \strong{\code{\link[base]{numeric}}} | Size of the font size of the legend. NULL uses default theme font size for legend according to the strong{\code{font.size}} parameter.
#' @param legend.symbol.size \strong{\code{\link[base]{numeric}}} | Size of symbols in the legend in mm. NULL uses the default size.
#' @return A ggplot2 object.
#' @export
#'
#' @examples
#' \dontrun{
#' TBD
#' }
do_MetadataPlot <- function(sample = NULL,
                            group.by = NULL,
                            metadata = NULL,
                            from_df = FALSE,
                            df = NULL,
                            colors.use = NULL,
                            flip = TRUE,
                            heatmap.gap = 1,
                            rotate_x_axis_labels = 45,
                            legend.position = "bottom",
                            font.size = 14,
                            legend.font.size = NULL,
                            legend.symbol.size = NULL,
                            legend.ncol = NULL,
                            legend.nrow = NULL,
                            legend.byrow = FALSE,
                            na.value = "grey75",
                            font.type = "sans",
                            grid.color = "white",
                            border.color = "black"){
  check_suggests(function_name = "do_MetadataPlot")
  
  # Check logical parameters.
  logical_list <- list("flip" = flip,
                       "from_df" = from_df,
                       "legend.byrow" = legend.byrow)
  check_type(parameters = logical_list, required_type = "logical", test_function = is.logical)
  
  # Check numeric parameters.
  numeric_list <- list("heatmap.gap" = heatmap.gap,
                       "rotate_x_axis_labels" = rotate_x_axis_labels,
                       "font.size" = font.size,
                       "legend.ncol" = legend.ncol,
                       "legend.nrow" = legend.nrow)
  check_type(parameters = numeric_list, required_type = "numeric", test_function = is.numeric)
  
  # Check character parameters.
  character_list <- list("group.by" = group.by,
                         "metadata" = metadata,
                         "legend.position" = legend.position,
                         "font.type" = font.type,
                         "grid.color" = grid.color,
                         "border.color" = border.color)
  check_type(parameters = character_list, required_type = "character", test_function = is.character)
  
  check_colors(grid.color, parameter_name = "grid.color")
  check_colors(border.color, parameter_name = "border.color")
  
  `%>%` <- magrittr::`%>%`
  `:=` <- rlang::`:=`
  
  if (isFALSE(from_df)){
    check_Seurat(sample = sample)
    
    for (meta in metadata){
      assertthat::assert_that(meta %in% colnames(sample@meta.data),
                              msg = paste0(crayon_body("Metadata column "),
                                           crayon_key(meta), 
                                           crayon_body(" is not in the sample "),
                                           crayon_key("metadata"),
                                           crayon_body(". Please check.")))
    }
    
    assertthat::assert_that(!is.null(sample) & !is.null(metadata) & !is.null(group.by),
                            msg = paste0(crayon_body("If "),
                                         crayon_key("from_df = FALSE"),
                                         crayon_body(" you need to use the "),
                                         crayon_key("sample"),
                                         crayon_body(", "),
                                         crayon_key("group.by"),
                                         crayon_body(", and "),
                                         crayon_key("metadata"),
                                         crayon_body(" parameters.")))
    
    data.plot <- sample@meta.data %>% 
                 tibble::rownames_to_column(var = "cell") %>% 
                 dplyr::select(dplyr::all_of(c(group.by, metadata))) %>% 
                 dplyr::group_by(.data[[group.by]]) %>% 
                 dplyr::summarise(dplyr::across(.cols = dplyr::all_of(c(metadata)), unique))
    
    assertthat::assert_that(length(unique(data.plot %>% dplyr::pull(.data[[group.by]]))) == nrow(data.plot),
                            msg = paste0(crayon_body("Please provide only metadata column that have a "),
                                         crayon_key(" one to one assignment "),
                                         crayon_body(" to the unique values in "),
                                         crayon_key("group.by"),
                                         crayon_body(".")))
    
    data.order <-  data.plot %>% 
                   tibble::column_to_rownames(var = group.by) %>% 
                   dplyr::mutate(dplyr::across(dplyr::everything(), as.factor))
  } else {
    assertthat::assert_that(!is.null(df),
                            msg = paste0(crayon_body("If "),
                                         crayon_key("from_df = TRUE"),
                                         crayon_body(" you need to use the "),
                                         crayon_key("df"),
                                         crayon_body(" parameter.")))
    
    group.by <- "Groups"
    metadata <- colnames(df)
    data.plot <- df %>% 
                 tibble::rownames_to_column(var = group.by)
    data.order <-  data.plot %>% 
                   tibble::column_to_rownames(var = group.by) %>% 
                   dplyr::mutate(dplyr::across(dplyr::everything(), as.factor))
  }
  

  order.use <- rownames(data.order)[stats::hclust(cluster::daisy(data.order, metric = c("gower")), method = "ward.D")$order]
  
  
  list.heatmaps <- list()
  
  # Get a list of predefined colors to then compute color wheels on for each metadata variable not covered.
  colors.pool <- get_SCpubr_colors()
  counter <- 0
  for (name in metadata){
    # Colors
    colors.use.name <- colors.use[[name]]
    if (is.null(colors.use.name)){
      counter <- counter + 1
      values <- unique(data.plot %>% dplyr::pull(name))
      
      colors.use.name <- stats::setNames(do_ColorPalette(n = length(values), colors.use = colors.pool[counter]),
                                          values)
    }
    
    
    data.use <- data.plot %>% 
                dplyr::select(dplyr::all_of(c(group.by, name))) %>% 
                dplyr::mutate("{name}_fill" := factor(.data[[name]]),
                              "{name}" := .env$name,
                              "{group.by}" := factor(.data[[group.by]], levels = order.use)) %>% 
                ggplot2::ggplot(mapping = ggplot2::aes(x = if(isFALSE(flip)){.data[[group.by]]} else {.data[[name]]},
                                                       y = if(isFALSE(flip)){.data[[name]]} else {.data[[group.by]]},
                                                       fill = .data[[paste0(name, "_fill")]])) + 
                ggplot2::geom_tile(color = grid.color, linewidth = 0.5) +
                ggplot2::scale_y_discrete(expand = c(0, 0)) +
                ggplot2::scale_x_discrete(expand = c(0, 0),
                                          position = "top") +
                ggplot2::guides(y.sec = guide_axis_label_trans(~paste0(levels(.data[[name]]))),
                                x.sec = guide_axis_label_trans(~paste0(levels(.data[[group.by]])))) + 
                ggplot2::coord_equal() + 
                ggplot2::scale_fill_manual(values = colors.use.name, name = name, na.value = na.value)
    list.heatmaps[[name]] <- data.use
  }
  
  # Modify legends.
  for (name in names(list.heatmaps)){
    p <- list.heatmaps[[name]]
    p <- p + 
         ggplot2::guides(fill = ggplot2::guide_legend(legend.position = legend.position,
                                                      title.position = "top",
                                                      title.hjust = ifelse(legend.position %in% c("top", "bottom"), 0.5, 0),
                                                      ncol = legend.ncol,
                                                      nrow = legend.nrow,
                                                      byrow = legend.byrow))
    list.heatmaps[[name]] <- p
  }
  
  # Add theme
  counter <- 0
  for (name in rev(names(list.heatmaps))){
    counter <- counter + 1
    # Set axis titles.
    if (isFALSE(flip)){
      if (counter == 1){
        xlab <- NULL
        ylab <- NULL
      } else if (counter == length(metadata)){
        xlab <- group.by
        ylab <- NULL
      } else {
        xlab <- NULL
        ylab <- NULL
      }
    } else {
      if (counter == 1){
        xlab <- NULL
        ylab <- group.by
      } else {
        xlab <- NULL
        ylab <- NULL
      }
    }
    
    
    p <- list.heatmaps[[name]]
    
    axis.parameters <- handle_axis(flip = flip,
                                   group.by = rep("A", length(names(list.heatmaps))),
                                   group = name,
                                   counter = counter,
                                   rotate_x_axis_labels = rotate_x_axis_labels)
    
    p <- p +
         ggplot2::xlab(xlab) +
         ggplot2::ylab(ylab) +
         ggplot2::theme_minimal(base_size = font.size) +
         ggplot2::theme(axis.ticks.x.bottom = axis.parameters$axis.ticks.x.bottom,
                        axis.ticks.x.top = axis.parameters$axis.ticks.x.top,
                        axis.ticks.y.left = axis.parameters$axis.ticks.y.left,
                        axis.ticks.y.right = axis.parameters$axis.ticks.y.right,
                        axis.text.y.left = axis.parameters$axis.text.y.left,
                        axis.text.y.right = axis.parameters$axis.text.y.right,
                        axis.text.x.top = axis.parameters$axis.text.x.top,
                        axis.text.x.bottom = axis.parameters$axis.text.x.bottom,
                        axis.title.x.bottom = axis.parameters$axis.title.x.bottom,
                        axis.title.x.top = axis.parameters$axis.title.x.top,
                        axis.title.y.right = axis.parameters$axis.title.y.right,
                        axis.title.y.left = axis.parameters$axis.title.y.left,
                        strip.background = axis.parameters$strip.background,
                        strip.clip = axis.parameters$strip.clip,
                        strip.text = axis.parameters$strip.text,
                        legend.position = legend.position,
                        axis.line = ggplot2::element_blank(),
                        plot.title = ggplot2::element_text(face = "bold", hjust = 0),
                        plot.subtitle = ggplot2::element_text(hjust = 0),
                        plot.caption = ggplot2::element_text(hjust = 1),
                        plot.title.position = "plot",
                        panel.grid = ggplot2::element_blank(),
                        panel.grid.minor.y = ggplot2::element_line(color = "white", linewidth = 1),
                        text = ggplot2::element_text(family = font.type),
                        plot.caption.position = "plot",
                        legend.text = ggplot2::element_text(face = "bold", size = legend.font.size),
                        legend.title = ggplot2::element_text(face = "bold", size = legend.font.size),
                        legend.justification = "center",
                        plot.margin = ggplot2::margin(t = heatmap.gap, r = 0, b = 0, l = heatmap.gap, unit = "mm"),
                        panel.border = ggplot2::element_rect(fill = NA, color = border.color, linewidth = 1),
                        panel.grid.major = ggplot2::element_blank(),
                        plot.background = ggplot2::element_rect(fill = "white", color = "white"),
                        panel.background = ggplot2::element_rect(fill = "white", color = "white"),
                        legend.background = ggplot2::element_rect(fill = "white", color = "white"),
                        panel.spacing.x = ggplot2::unit(0, "cm"))
    
    if (!is.null(legend.symbol.size)){
      p <- p + ggplot2::theme(legend.key.size = ggplot2::unit(legend.symbol.size, "mm"))
    }
    
    list.heatmaps[[name]] <- p
  }
  
  if (isTRUE(flip)){
    names.use <- rev(metadata)
  } else {
    names.use <- metadata
  }
  p <- patchwork::wrap_plots(list.heatmaps[names.use],
                             ncol = if (isFALSE(flip)){1} else {NULL},
                             nrow = if(isTRUE(flip)) {1} else {NULL},
                             guides = "collect")
  p <- p +
       patchwork::plot_annotation(theme = ggplot2::theme(legend.position = legend.position,
                                                         plot.title = ggplot2::element_text(family = font.type,
                                                                                            color = "black",
                                                                                            face = "bold",
                                                                                            hjust = 0),
                                                         plot.subtitle = ggplot2::element_text(family = font.type,
                                                                                               color = "black",
                                                                                               hjust = 0),
                                                         plot.caption = ggplot2::element_text(size = font.size,
                                                                                              family = font.type,
                                                                                              color = "black",
                                                                                              hjust = 1),
                                                         plot.caption.position = "plot"))
     
  return(p)
}