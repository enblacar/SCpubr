#' Generate color scales based on a value.
#'
#' This function is an adaptation of colortools package. As the package was removed from CRAN on 23-06-2022, this utility function came to existence in order to cover the gap. It is, on its basis,
#' an adaptation of the package into a single function. Original code, developed by Gaston Sanchez, can be found in: <https://github.com/gastonstat/colortools>
#'
#' @inheritParams doc_function
#' @param colors.use \strong{\code{\link[base]{character}}} | One color upon which generate the color scale. Can be a name or a HEX code.
#' @param n \strong{\code{\link[base]{numeric}}} | Number of colors to include in the color wheel. Use it when all other options are FALSE, otherwise, it becomes 12.
#' @param opposite \strong{\code{\link[base]{logical}}} | Return the opposing color to the one provided.
#' @param adjacent \strong{\code{\link[base]{logical}}} | Return the adjacent colors to the one provided.
#' @param triadic \strong{\code{\link[base]{logical}}} | Return the triadic combination of colors to the one provided.
#' @param split_complementary \strong{\code{\link[base]{logical}}} | Return the split complementary combination of colors to the one provided.
#' @param tetradic \strong{\code{\link[base]{logical}}} | Return the tetradic combination of colors to the one provided.
#' @param square \strong{\code{\link[base]{logical}}} | Return the square combination of colors to the one provided.
#' @param complete_output \strong{\code{\link[base]{logical}}} | Runs all the previous options and returns all the outputs as a list that contains all color vectors, all plots and a combined plot with everything.
#' @param plot \strong{\code{\link[base]{logical}}} | Whether to also return a plot displaying the values instead of a vector with the color.
#' @return  A character vector with the desired color scale.
#' @export
#' @example man/examples/examples_do_ColorPalette.R

do_ColorPalette <- function(colors.use,
                            n = 12,
                            opposite = FALSE,
                            adjacent = FALSE,
                            triadic = FALSE,
                            split_complementary = FALSE,
                            tetradic = FALSE,
                            square = FALSE,
                            complete_output = FALSE,
                            plot = FALSE,
                            font.size = 14,
                            font.type = "sans"){
  check_suggests(function_name = "do_ColorPalette")
  # Check logical parameters.
  logical_list <- list("opposite" = opposite,
                       "adjacent" = adjacent,
                       "triadic" = triadic,
                       "split_complementary" = split_complementary,
                       "tetradic" = tetradic,
                       "square" = square,
                       "complete_output" = complete_output)
  check_type(parameters = logical_list, required_type = "logical", test_function = is.logical)

  # Check numeric parameters.
  numeric_list <- list("n" = n,
                       "font.size", font.size)
  check_type(parameters = numeric_list, required_type = "numeric", test_function = is.numeric)

  # Check character parameters.
  character_list <- list("colors.use" = colors.use,
                         "font.type" = font.type)
  check_type(parameters = character_list, required_type = "character", test_function = is.character)

  # Check that the colors provided are only one.
  assertthat::assert_that(length(colors.use) == 1,
                          msg = "Please, provide a single color to colors.use.")

  # Check that the color provided is a valid color representation.
  check_colors(colors.use, parameter_name = "colors.use")

  # Check that only one option is activated.
  options_list <- c(opposite, adjacent, triadic, split_complementary, tetradic, square, complete_output)
  if (sum(options_list) > 0){
    assertthat::assert_that(sum(options_list) == 1,
                            msg = "Please select only one option to form the color scale. If you want more than one output, consider using complete_output = TRUE.")
  }

  # Check that n is actually positive.
  assertthat::assert_that(n > 0,
                          msg = "Please provide a positive integer value for 'n'.")

  # If any option is set to TRUE, pal_length is 12
  if (sum(options_list) >= 1 & n != 12){
    warning("When an color output option is selected, n parameter becomes by default 12. Please consider not using n for these purposes.", call. = FALSE)
    n <- 12
  }

  # Convert input to RGB colors: Input can be either color names, hex code.
  RGB_colors <- grDevices::col2rgb(colors.use)

  # Convert RGB values to HSV values.
  HSV_colors <- grDevices::rgb2hsv(RGB_colors)[, 1]

  # Get HSV components.
  hue <- HSV_colors[[1]] # Hue
  sat <- HSV_colors[[2]] # Saturation
  val <- HSV_colors[[3]] # Value

  # Generate a vector of hues that range a total of 1 unit, divided equally by n.
  hue_vector <- seq(hue, hue + 1, by = 1 / n)
  # Subset only the n colors.
  hue_vector <- hue_vector[1:n]
  # As this will generate hues over 1, anything over it, we deduct 1.
  hue_vector[hue_vector > 1] <- hue_vector[hue_vector > 1] - 1

  # Transform HSV values into HEX codes.
  colors <- grDevices::hsv(hue_vector, sat, val)

  # Add transparency value of the original color to the generated color scale.
  # This only works in the case the original color has a transparency value.
  if (substr(colors.use, 1, 1) == "#" && nchar(colors.use) == 9){
    alpha <- substr(colors.use, 8, 9)
    colors <- paste(colors, alpha, sep="")
  }

  # If opposite is TRUE, select the first and middle colors.
  if (isTRUE(opposite)){
    colors.mod <- colors[c(1, 7)]
    # If adjacent is TRUE, select the hues next to the original color.
  } else if (isTRUE(adjacent)){
    colors.mod <- colors[c(1, 2, 12)]
    # If triadic is TRUE, select the hues forming a triangle.
  } else if (isTRUE(triadic)){
    colors.mod <- colors[c(1, 5, 9)]
    # If split_complementary is TRUE, select the hues forming a triangle.
  } else if (isTRUE(split_complementary)){
    colors.mod <- colors[c(1, 6, 8)]
    # If tetradic is TRUE, select the hues forming a triangle.
  } else if (isTRUE(tetradic)){
    colors.mod <- colors[c(1, 3, 7, 9)]
    # If square is TRUE, select the hues forming a triangle.
  } else if (isTRUE(square)){
    colors.mod <- colors[c(1, 4, 7, 10)]
    # If complete_output is TRUE, report everything.
  } else {
    colors.mod <- colors
  }

  if (isTRUE(plot) & isFALSE(complete_output)){
    # Dummy df to plot.
    names(colors) <- colors
    df <- data.frame("values" = rep(1, n), "names" = factor(colors, levels = names(colors)))
    limits <- c(-5, 1.35)
    colors.use <- colors

    # Define name for the center of the plot.
    if (isTRUE(opposite)){
      name_center <- "Opposite"
      colors.use[!(names(colors.use) %in% colors[c(1, 7)])] <- "grey75"
      # If adjacent is TRUE, select the hues next to the original color.
    } else if (isTRUE(adjacent)){
      name_center <- "Adjacent"
      colors.use[!(names(colors.use) %in% colors[c(1, 2, 12)])] <- "grey75"
      # If triadic is TRUE, select the hues forming a triangle.
    } else if (isTRUE(triadic)){
      name_center <- "Triadic"
      colors.use[!(names(colors.use) %in% colors[c(1, 5, 9)])] <- "grey75"
      # If split_complementary is TRUE, select the hues forming a triangle.
    } else if (isTRUE(split_complementary)){
      name_center <- stringr::str_wrap("Split complementary", width = 5)
      colors.use[!(names(colors.use) %in% colors[c(1, 6, 8)])] <- "grey75"
      # If tetradic is TRUE, select the hues forming a triangle.
    } else if (isTRUE(tetradic)){
      name_center <- "Tetradic"
      colors.use[!(names(colors.use) %in% colors[c(1, 3, 7, 9)])] <- "grey75"
      # If square is TRUE, select the hues forming a triangle.
    } else if (isTRUE(square)){
      name_center <- "Square"
      colors.use[!(names(colors.use) %in% colors[c(1, 4, 7, 10)])] <- "grey75"
      # If complete_output is TRUE, report everything.
    } else {
      name_center <- "Wheel"
    }

    # Define blank labels.
    count <- 0
    if ("grey75" %in% colors.use){
      names.vector <- c()
      # Iterate over each color.
      for (name in names(colors.use)){
        if (colors.use[name] == "grey75"){
          count <- count + 1
          label.use <- paste0(rep(" ", count), collapse = "")
        } else {
          label.use <- name
        }
        names.vector <- c(names.vector, label.use)
      }
      names(colors.use) <- names.vector
      df[["names"]] <- factor(names(colors.use), levels = names(colors.use))
    }

    p <- ggplot2::ggplot(data = df, mapping = ggplot2::aes(x = .data[["names"]],
                                                           y = .data[["values"]],
                                                           fill = .data[["names"]])) +
         ggplot2::geom_col(color = "black", size = 1) +
         ggplot2::coord_polar(start = ifelse(sum(options_list) == 1,  -0.275, 0), direction = 1, clip = "off") +
         ggplot2::scale_fill_manual(values = colors.use, na.value = "grey75") +
         ggplot2::ylim(limits) +
         # Add X axis title in the center of the plot.
         ggplot2::annotate(geom = "text",
                           x = df[["names"]][[1]],
                           y = limits[[1]],
                           angle = 0,
                           hjust = 0.5,
                           vjust = 0.5,
                           label = name_center,
                           size = 8,
                           fontface = "bold") +
         ggplot2::theme_minimal(base_size = font.size) +
         ggplot2::theme(axis.title = ggplot2::element_blank(),
                        axis.ticks = ggplot2::element_blank(),
                        axis.text.y = ggplot2::element_blank(),
                        axis.text.x = ggplot2::element_text(face = "bold", color = "black"),
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
                        plot.margin = ggplot2::margin(t = 10, r = 40, b = 10, l = 40),
                        plot.background = ggplot2::element_rect(fill = "white", color = "white"),
                        panel.background = ggplot2::element_rect(fill = "white", color = "white"),
                        legend.background = ggplot2::element_rect(fill = "white", color = "white"))

  } else if (isTRUE(plot) & isTRUE(complete_output)) {
    stop("Parameter 'plot' only works when 'complete_output' is FALSE.")
  }


  # Complete output.


  # If plot = TRUE, return the plot, if not, colors. If complete_output = TRUE, return the report.
  if (isTRUE(complete_output)){
    # List of colors.
    return_colors <- list("wheel" = do_ColorPalette(colors.use = colors.use,
                                                    n = n),
                          "opposite" = do_ColorPalette(colors.use = colors.use,
                                                       n = n,
                                                       opposite = TRUE),
                          "adjacent" = do_ColorPalette(colors.use = colors.use,
                                                       n = n,
                                                       adjacent = TRUE),
                          "triadic" = do_ColorPalette(colors.use = colors.use,
                                                      n = n,
                                                      triadic = TRUE),
                          "split_complementary" = do_ColorPalette(colors.use = colors.use,
                                                                  n = n,
                                                                  split_complementary = TRUE),
                          "tetradic" = do_ColorPalette(colors.use = colors.use,
                                                       n = n,
                                                       tetradic = TRUE),
                          "square" = do_ColorPalette(colors.use = colors.use,
                                                     n = n,
                                                     square = TRUE))

    # List of plots.
    return_plots <- list("wheel" = do_ColorPalette(colors.use = colors.use,
                                                   n = n,
                                                   plot = TRUE),
                         "opposite" = do_ColorPalette(colors.use = colors.use,
                                                      n = n,
                                                      opposite = TRUE,
                                                      plot = TRUE),
                         "adjacent" = do_ColorPalette(colors.use = colors.use,
                                                      n = n,
                                                      adjacent = TRUE,
                                                      plot = TRUE),
                         "triadic" = do_ColorPalette(colors.use = colors.use,
                                                     n = n,
                                                     triadic = TRUE,
                                                     plot = TRUE),
                         "split_complementary" = do_ColorPalette(colors.use = colors.use,
                                                                 n = n,
                                                                 split_complementary = TRUE,
                                                                 plot = TRUE),
                         "tetradic" = do_ColorPalette(colors.use = colors.use,
                                                      n = n,
                                                      tetradic = TRUE,
                                                      plot = TRUE),
                         "square" = do_ColorPalette(colors.use = colors.use,
                                                    n = n,
                                                    square = TRUE,
                                                    plot = TRUE))

    layout <- "ABCD
               EFGH"

    patch <- patchwork::wrap_plots(A = return_plots[["wheel"]],
                                   B = return_plots[["opposite"]],
                                   C = return_plots[["adjacent"]],
                                   D = return_plots[["triadic"]],
                                   E = return_plots[["split_complementary"]],
                                   F = return_plots[["tetradic"]],
                                   G = return_plots[["square"]],
                                   H = patchwork::plot_spacer(),
                                   design = layout)

    # Build the output object.
    return_object <- list("colors" = return_colors,
                          "plots" = return_plots,
                          "combined_plot" = patch)

  } else {
    if (isTRUE(plot)){
      return_object <- p
    } else {
      return_object <- colors.mod
    }
  }

  return(return_object)
}
