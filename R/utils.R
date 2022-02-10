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


#' Compute the max and min value of a variable provided to FeaturePlot.
#'
#' @param sample Seurat object.
#' @param feature Feature to plot.
#' @param assay Assay used.
#' @param reduction Reduction used.
#'
#' @return
#' @noRd
#' @examples
#' \dontrun{
#' TBD
#' }
compute_scale_limits <- function(sample, feature, assay, reduction){
  if (feature %in% rownames(sample)){
    scale.begin <- min(sample@assays[[assay]]@data[feature,])
    scale.end <- max(sample@assays[[assay]]@data[feature,])
  } else if (feature %in% colnames(sample@meta.data)){
    scale.begin <- min(sample@meta.data[, feature])
    scale.end <- max(sample@meta.data[, feature])
  } else if (feature %in% colnames(sample@reductions[[reduction]][[]])){
    scale.begin <- min(sample@reductions[[reduction]][[]][, feature])
    scale.end <- max(sample@reductions[[reduction]][[]][, feature])
  }
  return(list("scale.begin" = scale.begin,
              "scale.end" = scale.end))
}


#' Check if the feature to plot is in the Seurat object.
#'
#' @param sample Seurat object.
#' @param features Feature to plot.
#' @param reduction Reduction used.
#'
#' @return
#' @noRd
#' @examples
#' \dontrun{
#' TBD
#' }
check_feature <- function(sample, features, reduction){
  for (feature in features){
    check <- 0
    if (!(feature %in% rownames(sample))){
      check <- check + 1
    }
    if (!(feature %in% colnames(sample@meta.data))){
      check <- check + 1
    }
    if (!(feature %in% colnames(sample@reductions[[reduction]][[]]))){
      check <- check + 1
    }
    if (check == 3) {
      stop(paste0("The requested feature (", feature, ") could not be found:\n", "    - Not matching any gene name (rownames of the provided object).\n",
                  "    - Not matching any metadata column (in sample@meta.data).\n", "    - Not part of the dimension names in the selected reduction (", reduction, ")."))
    }
  }
}

#' Check if the identity provided is in the current Seurat identities.
#'
#' @param sample Seurat object.
#' @param identities Identities to test.
#'
#' @return
#' @noRd
#' @examples
#' \dontrun{
#' TBD
#' }
check_identity <- function(sample, identities){
  for (identity in identities){
    if (!(identity %in% levels(sample))){
      stop(paste0("Could not find provided identity (", identity, ") in the current active identities of the object.\n Try running 'levels(your_seurat_object)' and see whether any typos were introduced."))
    }
  }
}
