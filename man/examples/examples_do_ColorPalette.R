# Check Suggests.
value <- SCpubr:::check_suggests(function_name = "do_ColorPalette", passive = TRUE)

if (isTRUE(value)){
  # Generate a color wheel based on a single value.
  colors <- SCpubr::do_ColorPalette(colors.use = "steelblue")
  p <- SCpubr::do_ColorPalette(colors.use = "steelblue",
                               plot = TRUE)

  # Generate a pair of opposite colors based on a given one.
  colors <- SCpubr::do_ColorPalette(colors.use = "steelblue",
                                    opposite = TRUE)
  p <- SCpubr::do_ColorPalette(colors.use = "steelblue",
                               opposite = TRUE,
                               plot = TRUE)

  # Generate a trio of adjacent colors based on a given one.
  colors <- SCpubr::do_ColorPalette(colors.use = "steelblue",
                                    adjacent = TRUE)
  p <- SCpubr::do_ColorPalette(colors.use = "steelblue",
                               adjacent = TRUE,
                               plot = TRUE)

  # Generate a trio of triadic colors based on a given one.
  colors <- SCpubr::do_ColorPalette(colors.use = "steelblue",
                                    triadic = TRUE)
  p <- SCpubr::do_ColorPalette(colors.use = "steelblue",
                               triadic = TRUE,
                               plot = TRUE)

  # Generate a trio of split complementary colors based on a given one.
  colors <- SCpubr::do_ColorPalette(colors.use = "steelblue",
                                    split_complementary = TRUE)
  p <- SCpubr::do_ColorPalette(colors.use = "steelblue",
                               split_complementary = TRUE,
                               plot = TRUE)

  # Generate a group of tetradic colors based on a given one.
  colors <- SCpubr::do_ColorPalette(colors.use = "steelblue",
                                    tetradic = TRUE)
  p <- SCpubr::do_ColorPalette(colors.use = "steelblue",
                               tetradic = TRUE,
                               plot = TRUE)

  # Generate a group of square colors based on a given one.
  colors <- SCpubr::do_ColorPalette(colors.use = "steelblue",
                                    square = TRUE)
  p <- SCpubr::do_ColorPalette(colors.use = "steelblue",
                               square = TRUE,
                               plot = TRUE)

  # Retrieve the output of all options.
  out <- SCpubr::do_ColorPalette(colors.use = "steelblue",
                                 complete_output = TRUE)
  ## Retrieve the colors.
  colors <- out$colors
  ## Retrieve the plots.
  plots <- out$plots
  ## Retrieve a combined plot with all the options.
  p <- out$combined_plot

} else if (isFALSE(value)){
  message("This function can not be used without its suggested packages.")
  message("Check out which ones are needed using `SCpubr::state_dependencies()`.")
}
