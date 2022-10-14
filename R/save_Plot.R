#' Save a plot as png, pdf and svg.
#'
#'
#' @param plot Plot to save.
#' @param figure_path \strong{\code{\link[base]{character}}} | Path where the figure will be stored.
#' @param create_path \strong{\code{\link[base]{logical}}} | Whether to create the path.
#' @param file_name \strong{\code{\link[base]{character}}} | Name of the file (without extension, it will be added automatically).
#' @param output_format \strong{\code{\link[base]{character}}} | Output format of the saved figure.  One of:
#' \itemize{
#'   \item \emph{\code{pdf}}: Saves the figure as a PDF file.
#'   \item \emph{\code{png}}: Saves the figure as a PNG file.
#'   \item \emph{\code{jpeg}}: Saves the figure as a JPEG file.
#'   \item \emph{\code{tiff}}: Saves the figure as a TIFF file.
#'   \item \emph{\code{svg}}: Saves the figure as a SVG file.
#'   \item \emph{\code{publication}}: Saves the figure as PDF, PNG and SVG files.
#'   \item \emph{\code{all}}: Saves the figure in all possible formats.
#' }
#' @param dpi \strong{\code{\link[base]{numeric}}} | Dpi to use.
#' @param width,height \strong{\code{\link[base]{numeric}}} | Width and height of the figure (inches).
#'
#' @return Nothing.
#' @export
#'
#' @example /man/examples/examples_save_Plot.R
save_Plot <- function(plot,
                      figure_path = NULL,
                      create_path = TRUE,
                      file_name = NULL,
                      dpi = 300,
                      output_format = "publication",
                      width = 8,
                      height = 8){


  # Checks for packages.
  check_suggests(function_name = "save_Plot")

  # Check logical parameters.
  logical_list <- list("create_path" = create_path)
  check_type(parameters = logical_list, required_type = "logical", test_function = is.logical)
  # Check numeric parameters.
  numeric_list <- list("dpi" = dpi,
                       "width" = width,
                       "height" = height)
  check_type(parameters = numeric_list, required_type = "numeric", test_function = is.numeric)
  # Check character parameters.
  character_list <- list("figure_path" = figure_path,
                         "file_name" = file_name)
  check_type(parameters = character_list, required_type = "character", test_function = is.character)

  # Null file name?
  if (is.null(file_name)){file_name <- "output_figure"}
  # Null figure path?
  if (is.null(figure_path)){figure_path <- paste0(".", .Platform$file.sep)}

  # Create directory.
  if (!(dir.exists(figure_path))){
    if (isTRUE(create_path)){dir.create(figure_path, recursive = TRUE)}
  }




  # Handle devices:
  output_options <- c("all", "publication", "pdf", "png", "jpeg", "svg", "tiff")

  assertthat::assert_that(sum(output_format %in% output_options) >= 1,
                          msg = "Please select a valid output format from the available options: all, publication, pdf, png, jpeg, svg, tiff")

  assertthat::assert_that(isFALSE("all" %in% output_format & "publication" %in% output_format),
                          msg = "Please select either `all` or `publication`.")

  if (output_format == "publication"){
    devices_use <- c("pdf", "png", "svg")
  } else if (output_format == "all"){
    devices_use <- c("pdf", "png", "jpeg", "svg", "tiff")
  } else {
    options <- c("pdf", "png", "jpeg", "svg", "tiff")
    devices_use <- output_format[output_format %in% options]
  }

  # is ggplot?

  if (sum(class(plot) %in% c("ggplot")) >= 1){
    # Having width = NULL and height = NULL will make the ggsave() function crash.
    for (device in devices_use){
        suppressMessages({
          ggplot2::ggsave(filename = sprintf("%s.%s", file_name, device),
                          plot = plot,
                          path = figure_path,
                          dpi = dpi,
                          width = width,
                          height = height,
                          device = device)
        })
    }
  # Is it a heatmap?
  } else if (sum(class(plot) %in% c("HeatmapList", "ComplexHeatmap")) >= 1) {
    suppressMessages({
      filename <- paste0(figure_path, "/", file_name)
      if ("png" %in% devices_use){
        grDevices::png(filename = paste0(filename, ".png"), units = "in", height = height, width = width, res = dpi)
        ComplexHeatmap::draw(plot, show_heatmap_legend = TRUE, padding = ggplot2::unit(c(20, 20, 2, 20), "mm"))
        grDevices::dev.off()
      }

      if ("pdf" %in% devices_use){
        grDevices::pdf(file = paste0(filename, ".pdf"), height = height, width = width)
        ComplexHeatmap::draw(plot, show_heatmap_legend = TRUE, padding = ggplot2::unit(c(20, 20, 2, 20), "mm"))
        grDevices::dev.off()
      }

      if ("jpeg" %in% devices_use){
        grDevices::jpeg(file = paste0(filename, ".jpeg"), units = "in", height = height, width = width, res = dpi)
        ComplexHeatmap::draw(plot, show_heatmap_legend = TRUE, padding = ggplot2::unit(c(20, 20, 2, 20), "mm"))
        grDevices::dev.off()
      }

      if ("tiff" %in% devices_use){
        grDevices::jpeg(file = paste0(filename, ".tiff"), units = "in", height = height, width = width, res = dpi)
        ComplexHeatmap::draw(plot, show_heatmap_legend = TRUE, padding = ggplot2::unit(c(20, 20, 2, 20), "mm"))
        grDevices::dev.off()
      }

      if ("svg" %in% devices_use){
        svglite::svglite(filename = paste0(filename, ".svg"), height = height, width = width)
        ComplexHeatmap::draw(plot, show_heatmap_legend = TRUE, padding = ggplot2::unit(c(20, 20, 2, 20), "mm"))
        grDevices::dev.off()
      }

    })

  } else if (sum(class(plot) %in% c("pheatmap")) >= 1){
    suppressMessages({
      filename <- paste0(figure_path, "/", file_name)
      if ("png" %in% devices_use){
        grDevices::png(filename = paste0(filename, ".png"), units = "in", height = height, width = width, res = dpi)
        print(plot)
        grDevices::dev.off()

      }

      if ("pdf" %in% devices_use){
        grDevices::pdf(file = paste0(filename, ".pdf"), height = height, width = width)
        print(plot)
        grDevices::dev.off()
      }

      if ("jpeg" %in% devices_use){
        grDevices::jpeg(file = paste0(filename, ".jpeg"), units = "in", height = height, width = width, res = dpi)
        print(plot)
        grDevices::dev.off()
      }

      if ("tiff" %in% devices_use){
        grDevices::jpeg(file = paste0(filename, ".tiff"), units = "in", height = height, width = width, res = dpi)
        print(plot)
        grDevices::dev.off()
      }

      if ("svg" %in% devices_use){
        svglite::svglite(filename = paste0(filename, ".svg"), height = height, width = width)
        print(plot)
        grDevices::dev.off()
      }
    })
  } else if (sum(class(plot) %in% c("recordedplot")) >= 1){
    suppressMessages({
      filename <- paste0(figure_path, "/", file_name)
      if ("png" %in% devices_use){
        grDevices::png(filename = paste0(filename, ".png"), units = "in", height = height, width = width, res = dpi)
        grDevices::replayPlot(plot)
        grDevices::dev.off()

      }

      if ("pdf" %in% devices_use){
        grDevices::pdf(file = paste0(filename, ".pdf"), height = height, width = width)
        grDevices::replayPlot(plot)
        grDevices::dev.off()
      }

      if ("jpeg" %in% devices_use){
        grDevices::jpeg(file = paste0(filename, ".jpeg"), units = "in", height = height, width = width, res = dpi)
        grDevices::replayPlot(plot)
        grDevices::dev.off()
      }

      if ("tiff" %in% devices_use){
        grDevices::jpeg(file = paste0(filename, ".tiff"), units = "in", height = height, width = width, res = dpi)
        grDevices::replayPlot(plot)
        grDevices::dev.off()
      }

      if ("svg" %in% devices_use){
        svglite::svglite(filename = paste0(filename, ".svg"), height = height, width = width)
        grDevices::replayPlot(plot)
        grDevices::dev.off()
      }
    })
  }
}
