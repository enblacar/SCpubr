#' Save a plot as png, pdf and svg.
#'
#'
#' @param plot Plot to save.
#' @param figure_path Path where the figure will be stored.
#' @param create_path Logical. Whether to create the path.
#' @param file_name Name of the file (without extension, it will be added automatically).
#' @param output_format Character. One of the following:
#' - pdf.
#' - png.
#' - jpeg.
#' - tiff.
#' - svg
#' - publication: Includes pdf, png and svg formats.
#' - all: Includes all formats.
#' @param dpi Dpi to use.
#' @param width,height Width and height of the figure (inches).
#'
#' @return Nothing.
#' @export
#'
#' @example /man/examples/examples_savePlot.R
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
  if (is.null(file_name)){file_name <- stringr::str_replace_all(Sys.time(), " ", "_")}
  # Null figure path?
  if (is.null(figure_path)){figure_path <- paste0(".", .Platform$file.sep)}

  # Create directory.
  if (!(dir.exists(figure_path))){
    if (isTRUE(create_path)){dir.create(figure_path, recursive = T)}
  }




  # Handle devices:
  output_options <- c("all", "publication", "pdf", "png", "jpeg", "svg", "tiff")
  if (sum(output_format %in% output_options) < 1){stop("Please select a valid output format from the available options.", call. = F)}
  if ("all" %in% output_format & "publication" %in% output_format){stop("Please select either `all` or `publication`.", call. = F)}
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
      filename <- paste0(figure_path, file_name)
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
      filename <- paste0(figure_path, file_name)
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
  }
}
