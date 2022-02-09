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
#'
#' @param colors Vector of colors.
#'
#' @return
#' @noRd
#' @examples
#' \dontrun{
#' TBD
#' }
check_colors <- function(colors) {
  check <- sapply(colors, function(color) {
              tryCatch(is.matrix(grDevices::col2rgb(colors)),
             error = function(e) FALSE)
  })
  return(check)
}
