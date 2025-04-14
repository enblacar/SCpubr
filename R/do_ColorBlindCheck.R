#' Generate colorblind variations of a given color palette.
#'
#' This function generate colorblind variations of a provided color palette in order to check if it is colorblind friendly. Variations are generated using colorspace package.
#'
#' @inheritParams doc_function
#' @param colors.use \strong{\code{\link[base]{character}}} | One color upon which generate the color scale. Can be a name or a HEX code.
#' @return  A character vector with the desired color scale.
#' @export
#' @example man/examples/examples_do_ColorBlindCheck.R

do_ColorBlindCheck <- function(colors.use,
                               flip = FALSE,
                               font.size = 14,
                               font.type = "sans",
                               plot.title.face = "bold",
                               plot.subtitle.face = "plain",
                               plot.caption.face = "italic",
                               axis.title.face = "bold",
                               axis.text.face = "plain",
                               legend.text.face = "plain",
                               legend.title.face = "bold",
                               grid.color = "white",
                               border.color = "black",
                               axis.text.x.angle = 45){
  
  `%>%` <- magrittr::`%>%`
  `:=` <- rlang::`:=`
  
  # Add lengthy error messages.
  withr::local_options(.new = list("warning.length" = 8170))
  
  check_suggests(function_name = "do_ColorPalette")
  # Check logical parameters.
  logical_list <- list("flip" = flip)
  check_type(parameters = logical_list, required_type = "logical", test_function = is.logical)
  
  # Check numeric parameters.
  numeric_list <- list("font.size", font.size,
                       "axis.text.x.angle" = axis.text.x.angle,
                       "font.size" = font.size)
  check_type(parameters = numeric_list, required_type = "numeric", test_function = is.numeric)
  
  # Check character parameters.
  character_list <- list("colors.use" = colors.use,
                         "font.type" = font.type,
                         "grid.color" = grid.color,
                         "border.color" = border.color,
                         "plot.title.face" = plot.title.face,
                         "plot.subtitle.face" = plot.subtitle.face,
                         "plot.caption.face" = plot.caption.face,
                         "axis.title.face" = axis.title.face,
                         "axis.text.face" = axis.text.face,
                         "legend.text.face" = legend.text.face,
                         "legend.title.face" = legend.title.face)
  check_type(parameters = character_list, required_type = "character", test_function = is.character)
  
  
  # Check that the color provided is a valid color representation.
  check_colors(colors.use, parameter_name = "colors.use")
  check_colors(grid.color, parameter_name = "grid.color")
  check_colors(border.color, parameter_name = "border.color")
  
  # Dicromatic view:
  deutan.colors <- colorspace::deutan(colors.use)  # Red-green (most common)
  protan.colors <- colorspace::protan(colors.use)  # Red-green (less common)
  tritan.colors <- colorspace::tritan(colors.use)  # Blue-yellow
  
  colors.use <- list("Normal" = colors.use,
                     "Protanopia" = protan.colors,
                     "Deuteranopia" = deutan.colors,
                     "Tritanopia" = tritan.colors)
  
  df <- as.data.frame(colors.use)
 # df <- df[rev(seq(1, length(rownames(df)))),]
  
  list.heatmaps <- list()
  metadata <- if(base::isFALSE(flip)){rev(colnames(df))} else {colnames(df)}
  group.by <- "Colors"
  
  data.plot <- df %>% 
               dplyr::mutate("{group.by}" := .data$Normal) %>% 
               tidyr::pivot_longer(cols = -"Colors",
                                   names_to = "Type",
                                   values_to = "Color") %>% 
               dplyr::mutate("{group.by}" := factor(.data[[group.by]], levels = df$Normal))
  
  # Get a list of predefined colors to then compute color wheels on for each metadata variable not covered.
  counter <- 0
  for (name in metadata){
    counter <- counter + 1
    # Colors
    colors.use.name <- df[, name]
    names(colors.use.name) <- df$Normal
    
    # Handle axis
    axis.parameters <- handle_axis(flip = flip,
                                   group.by = rep("A", length(metadata)),
                                   group = name,
                                   counter = counter,
                                   axis.text.x.angle = axis.text.x.angle,
                                   plot.title.face = plot.title.face,
                                   plot.subtitle.face = plot.subtitle.face,
                                   plot.caption.face = plot.caption.face,
                                   axis.title.face = axis.title.face,
                                   axis.text.face = axis.text.face,
                                   legend.title.face = "bold",
                                   legend.text.face = "plain")
    
    
    p <- data.plot %>% 
         dplyr::filter(.data$Type == name) %>% 
         # nocov start
         ggplot2::ggplot(mapping = ggplot2::aes(x = if(base::isFALSE(flip)){.data[[group.by]]} else {.data$Type},
                                                y = if(base::isFALSE(flip)){.data$Type} else {.data[[group.by]]},
                                                fill = .data[[group.by]])) + 
         # nocov end
         ggplot2::geom_tile(color = grid.color, linewidth = 0.5) +
         ggplot2::scale_y_discrete(expand = c(0, 0)) +
         ggplot2::scale_x_discrete(expand = c(0, 0),
                                   position = "top") +
         ggplot2::guides(y.sec = guide_axis_label_trans(~paste0(levels(.data$Type))),
                         x.sec = guide_axis_label_trans(~paste0(levels(.data[[group.by]])))) + 
         ggplot2::coord_equal() + 
         ggplot2::scale_fill_manual(values = colors.use.name, name = name, na.value = "grey75") +
         ggplot2::xlab(NULL) +
         ggplot2::ylab(NULL) +
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
                        legend.position = "none",
                        axis.line = ggplot2::element_blank(),
                        plot.title = ggplot2::element_text(face = plot.title.face, hjust = 0),
                        plot.subtitle = ggplot2::element_text(face = plot.subtitle.face, hjust = 0),
                        plot.caption = ggplot2::element_text(face = plot.caption.face, hjust = 1),
                        plot.title.position = "plot",
                        panel.grid = ggplot2::element_blank(),
                        panel.grid.minor.y = ggplot2::element_line(color = "white", linewidth = 1),
                        text = ggplot2::element_text(family = font.type),
                        plot.caption.position = "plot",
                        legend.text = ggplot2::element_text(face = "plain", size = font.size),
                        legend.title = ggplot2::element_text(face = "bold", size = font.size),
                        legend.justification = "center",
                        plot.margin = ggplot2::margin(t = 0, r = 0, b = 0, l = 0, unit = "mm"),
                        panel.border = ggplot2::element_rect(fill = NA, color = border.color, linewidth = 1),
                        panel.grid.major = ggplot2::element_blank(),
                        plot.background = ggplot2::element_rect(fill = "white", color = "white"),
                        panel.background = ggplot2::element_rect(fill = "white", color = "white"),
                        legend.background = ggplot2::element_rect(fill = "white", color = "white"),
                        legend.spacing = ggplot2::unit(0, "cm"),
                        panel.spacing.x = ggplot2::unit(0, "cm"))
    list.heatmaps[[name]] <- p
  }
  
  # Tweak Normal plot space.
  list.heatmaps[["Normal"]] <- list.heatmaps[["Normal"]] + 
                               ggplot2::theme(plot.margin = ggplot2::margin(t = 0, 
                                                                            r = if (base::isFALSE(flip)){0} else {5}, 
                                                                            b = if (base::isFALSE(flip)){5} else {0}, 
                                                                            l = 0, unit = "mm"))
  
  p <- patchwork::wrap_plots(list.heatmaps[if(base::isFALSE(flip)){rev(metadata)} else {metadata}],
                             ncol = if (base::isFALSE(flip)){1} else {NULL},
                             nrow = if(isTRUE(flip)) {1} else {NULL},
                             guides = "collect") 
      
  
  return(p)
}
