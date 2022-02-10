#' Checks for Suggests.
#'
#' @param function_name Function to check.
#' @noRd
#' @return None
#' @examples
#' \dontrun{
#' TBD
#' }
check_suggests <- function(function_name){

  pkg_list <- list("do_BarPlot" = c("Seurat", "colortools", "dplyr", "ggplot2", "ggpubr", "purrr", "utils", "rlang"),
                   "do_ButterflyPlot" = c("Seurat", "tidyr", "pbapply", "dplyr", "ggplot2", "ggpubr", "viridis", "purrr", "utils", "rlang"),
                   "do_DimPlot" = c("grDevices", "colortools", "Seurat", "ggpubr", "ggplot2", "patchwork"),
                   "do_DotPlot" = c("Seurat", "ggplot2", "ggpubr"),
                   "do_FeaturePlot" = c("Seurat", "viridis", "ggplot2", "ggpubr", "patchwork"),
                   "do_NebulosaPlot" = c("Seurat", "ggplot2", "ggpubr", "Nebulosa", "patchwork"),
                   "do_PTEA" = c("Seurat", "stringr", "pbapply", "Matrix", "dplyr", "tidyr", "stats", "purrr", "utils", "rlang"),
                   "do_RankPlot" = c("Seurat", "ggplot2", "ggpubr", "viridis", "colortools", "ggbeeswarm"),
                   "do_VlnPlot" = c("Seurat", "ggplot2", "ggpubr", "scales"))
  pkgs <- pkg_list[[function_name]]
  for (pkg in pkgs){
    if(!requireNamespace(pkg, quietly = T)){
      stop(paste0("Package ", pkg, " must be installed to use ", function_name, "."))
    }
  }
}



#' Internal check for colors.
#'
#' Adapted from: https://stackoverflow.com/a/13290832.
#
#' @param colors Vector of colors.
#' @param parameter_name The name of the parameter for which we are testing the colors.
#' @return
#' @noRd
#' @examples
#' \dontrun{
#' TBD
#' }
check_colors <- function(colors, parameter_name = "") {
  check <- sapply(colors, function(color) {
              tryCatch(is.matrix(grDevices::col2rgb(colors)),
             error = function(e) FALSE)
  })
  # Check for cols.highlight.
  if (sum(check) != length(colors)){
    stop(paste0("One/several of the value/s for ", parameter_name, " is not a valid color representation. Please check whether it is an accepted R name or a HEX code."))
  }
}

#' Internal check for named colors and unique values of the grouping variable.
#'
#' @param sample Seurat object.
#' @param colors Named vector of colors.
#' @param groping_variable Metadata variable in sample to obtain its unique values.
#' @return
#' @noRd
#' @examples
#' \dontrun{
#' TBD
#' }
check_consistency_colors_and_names <- function(sample, colors, groping_variable){
  if (sum(names(colors) %in% unique(sample[[]][, groping_variable])) != length(unique(sample[[]][, groping_variable]))){
    stop('The names of the color vector provided to "colors.highlight" do not entirely match the unique values in "split.by" parameter.')
  }
  if (length(colors) != length(unique(sample[[]][, groping_variable]))){
    stop('The number of values provided to "colors.split" is lower than the unique values in "split.by" parameter.')
  }
}

#' Generate custom color scale.
#'
#' @param names_use Vector of the names that will go alongside the color scale.
#' @return
#' @noRd
#' @examples
#' \dontrun{
#' TBD
#' }
generate_color_scale <- function(names_use){
  # Generate a vector of colors equal to the number of identities in the sample.
  colors <- colortools::setColors("#457b9d", length(names_use))
  names(colors) <- names_use
  return(colors)
}
