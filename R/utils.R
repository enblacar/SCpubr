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
                   "do_FeaturePlot" = c("Seurat", "viridis", "ggplot2", "ggpubr", "patchwork", "scales"),
                   "do_NebulosaPlot" = c("Seurat", "ggplot2", "ggpubr", "Nebulosa", "patchwork"),
                   "do_PTEA" = c("Seurat", "stringr", "pbapply", "Matrix", "dplyr", "tidyr", "stats", "purrr", "utils", "rlang"),
                   "do_RankPlot" = c("Seurat", "ggplot2", "ggpubr", "viridis", "colortools", "ggbeeswarm"),
                   "do_VlnPlot" = c("Seurat", "ggplot2", "ggpubr"))
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
    stop(paste0("The value/s for ", parameter_name, " is/are not a valid color representation. Please check whether it is an accepted R name or a HEX code."))
  }
}

#' Internal check for named colors and unique values of the grouping variable.
#'
#' @param sample Seurat object.
#' @param colors Named vector of colors.
#' @param grouping_variable Metadata variable in sample to obtain its unique values.
#' @return
#' @noRd
#' @examples
#' \dontrun{
#' TBD
#' }
check_consistency_colors_and_names <- function(sample, colors, grouping_variable = NULL){
  if (is.null(grouping_variable)){
    check_values <- levels(sample)
  } else {
    check_values <- unique(sample@meta.data[, grouping_variable])
  }
  if (sum(names(colors) %in% check_values) != length(check_values)){
    stop('The names of the colors in the vector provided do not match the number of unique values in the selected grouping variable (levels(object), group.by or split.by).')
  }
  if (length(colors) != length(check_values)){
    stop('The number of colors provided is lower than the unique values in the selected grouping variable (levels(object), group.by or split.by).')
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
compute_scale_limits <- function(sample, feature, assay = NULL, reduction = NULL){
  if (is.null(assay)){
    assay <- Seurat::DefaultAssay(sample)
  }
  if (is.null(reduction)){
    dim_colnames <- c()
    for(red in Seurat::Reductions(object = sample)){
      if (feature %in% colnames(sample@reductions[[red]][[]])){
        reduction <- red
      }
    }
  }

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
#' @param dump_reduction_names Whether to return the reduction colnames.
#'
#' @return
#' @noRd
#' @examples
#' \dontrun{
#' TBD
#' }
check_feature <- function(sample, features, dump_reduction_names = FALSE, enforce_check = NULL, enforce_parameter = NULL){
  if (is.list(features)){
    features <- unlist(features)
  }
  check_enforcers <- list()
  for (feature in features){
    check <- 0
    if (!(feature %in% rownames(sample))){
      check <- check + 1
      check_enforcers[["gene"]] <- FALSE
    }
    check_enforcers[["gene"]] <- TRUE
    if (!(feature %in% colnames(sample@meta.data))){
      check <- check + 1
      check_enforcers[["metadata"]] <- FALSE
    }
    check_enforcers[["metadata"]] <- TRUE
    dim_colnames <- c()
    for(red in Seurat::Reductions(object = sample)){
      dim_colnames <- c(dim_colnames, colnames(sample@reductions[[red]][[]]))
    }
    if (!(feature %in% dim_colnames)){
      check <- check + 1
      check_enforcers[["reductions"]] <- FALSE
    }
    check_enforcers[["reductions"]] <- TRUE
    if (check == 3) {
      stop(paste0("The requested feature (", feature, ") could not be found:\n", "    - Not matching any gene name (rownames of the provided object).\n",
                  "    - Not matching any metadata column (in sample@meta.data).\n", "    - Not part of the dimension names in any of the following reductions: ", Seurat::Reductions(object = sample), "."))
    }
  }
  if (!(enforce_check %in% names(check_enforcers))){
    stop("The variable enforcer is not in the current list of checked variable types.")
  } else {
    if (isFALSE(check_enforcers[[enforce_check]])){
      stop("The provided feature (", enforce_parameter, " = ", feature, ") not found in ", enforce_check, ".")
    }
  }
  if (dump_reduction_names == TRUE){return(dim_colnames)}
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

#' Check the reduction provided and set it up.
#'
#' @param sample Seurat sample.
#' @param reduction Reduction.
#'
#' @return
#' @noRd
#' @examples
#' \dontrun{
#' TBD
#' }
check_and_set_reduction <- function(sample, reduction){
  # Check if the object has a reduction computed.
  if (length(Seurat::Reductions(sample)) == 0){stop("This object has no reductions computed!")}
  # If no reduction was provided by the user.
  if (is.null(reduction)){
    # Select umap if computed.
    if ("umap" %in% Seurat::Reductions(sample)){
      reduction <- "umap"
    } else {
      # Select the last computed one.
      reduction <- Seurat::Reductions(sample)[length(Seurat::Reductions(sample))]
    }
  # If the user provided a value for reduction.
  } else if (!(is.null(reduction))){
    # Check if the provided reduction is in the list.
    if (!(reduction %in% Seurat::Reductions(sample))){stop("The provided reduction could not be found in the object: ", reduction)}
  }
  return(reduction)
}

#' Check the provided dimensions and set them up.
#'
#' @param sample Seurat object.
#' @param reduction Provided reduction.
#' @param dims Provided dimensions.
#'
#' @return
#' @noRd
#' @examples
#' \dontrun{
#' TBD
#' }
check_and_set_dimensions <- function(sample, reduction, dims){
  # Check that the dimensions is a 2 item vector.
  if (!(is.null(dims)) & length(dims) != 2){
    stop("Provided dimensions need to be a 2-item vector.")
  }
  # Check that the dimensions are integers.
  null_check <- is.null(dims[1]) & is.null(dims[2])
  integer_check <- is.numeric(dims[1]) & is.numeric(dims[1])
  if (!(is.null(dims)) & integer_check == FALSE){
    stop("Provied dimensions need to be numerics.")
  }
  # Check that the dimensions are in the requested embedding.
  aval_dims <- length(colnames(Seurat::Embeddings(sample[[reduction]])))
  if (!(is.null(dims))){
    if (!(dims[1] %in% seq_len(aval_dims)) | !(dims[2] %in% seq_len(aval_dims))){
      stop("Dimension could not be found in the following reduction: ", reduction, ".")
    }
  }
  # Check that at least 2 dimensions are present.
  if (aval_dims < 2){
    stop("There are less than 2 dimensions in the requested reduction: ", reduction, ".")
  }
  # If no dimensions were provided, fall back to first and second.
  if (is.null(dims)){
    dims <- c(1, 2)
  }
  return(dims)
}

#' Check and set the provided assay.
#'
#' @param sample Seurat object.
#' @param assay Provided assay.
#'
#' @return
#' @noRd
#' @examples
#' \dontrun{
#' TBD
#' }
check_and_set_assay <- function(sample, assay){
  # If assay is null, set it to the active one.
  if (is.null(assay)){
    assay <- Seurat::DefaultAssay(sample)
  } else {
    # Check if the assay is a character.
    if (!(is.character(assay))){
      stop("The value for assay has to be a character.")
    }
    # Check that at least one assay is computed.
    if (length(Seurat::Assays(sample)) == 0){
      stop("There must be at least one computed assay in the object.")
    }
    # Check that the assay is in the available assays.
    aval_assays <- Seurat::Assays(sample)
    if (!(assay %in% aval_assays)){
      stop("The following assay could not be found: ", assay, ".")
    }
  }
  # Set up the assay the user has defined.
  if (assay != Seurat::DefaultAssay(sample)){
    Seurat::DefaultAssay(sample) <- assay
  }
  return(list("sample" = sample,
              "assay" = assay))
}


#' Check a parameter for a given class.
#'
#' @param parameters List of named parameters to test.
#' @param required_type Name of the required class.
#' @param test_function Testing function.
#'
#' @return
#' @noRd
#' @examples
#' \dontrun{
#' TBD
#' }
check_type <- function(parameters, required_type, test_function){
  for(parameter_name in names(parameters)){
    parameter <- parameters[[parameter_name]]
    if (!(is.null(parameter)) & !(test_function(parameter))){
      stop("Parameter ", parameter_name, " needs to be a ", required_type, ".")
    }
  }
}

#' Check the slots.
#'
#' @param slot Slot provided.
#'
#' @return
#' @noRd
#' @examples
#' \dontrun{
#' TBD
#' }
check_and_set_slot <- function(slot){
  if (is.null(slot)){
    slot <- "data"
  } else if (!(slot %in% c("counts", "data", "scale.data"))){
    stop("Only one of these 3 options can be passed to slot parameter: counts, data, scale.data.")
  }
  return(slot)
}

#' Check if a value is in the range of the values.
#'
#' @param sample Seurat object.
#' @param feature Feature to plot.
#' @param assay Assay used.
#' @param reduction Reduction used.
#' @param value Value to check.
#' @param value_name Name of the value.
#'
#' @return
#' @noRd
#' @examples
#' \dontrun{
#' TBD
#' }
check_limits <- function(sample, feature, value_name, value, assay = NULL, reduction = NULL){
  limits <- compute_scale_limits(sample = sample, feature = feature, assay = assay, reduction = reduction)

  if (!(limits[["scale.begin"]] <= value & limits[["scale.end"]] >= value)){
    stop("The value provided for ", value_name, " (", value, ") is not in the range of the feature (", feature, "), which is: Min: ", limits[["scale.begin"]], ", Max: ", limits[["scale.end"]], ".")
  }
}

#' Compute the order of the plotted bars for do_BarPlot.
#'
#' @param sample Seurat object.
#' @param feature Feature to plot.
#' @param group.by Feature to group the output by.
#' @param order.by Unique value in group.by to reorder labels in descending order.
#'
#' @return
#' @noRd
#' @examples
#' \dontrun{
#' TBD
#' }
compute_factor_levels <- function(sample, feature, group.by = NULL, order.by = NULL, position = NULL){
  `%>%` <- purrr::`%>%`
  if (is.null(order.by) & !(is.null(group.by))){
    if (is.null(position)){stop("Position parameter needs to be provided.")}
    if (position == "fill"){
      factor_levels <- as.character(rev(sort(unique(sample@meta.data[, feature]))))
    } else if (position == "stack"){
      factor_levels <- sample@meta.data %>% # Obtain metadata
        dplyr::select(!!rlang::sym(feature)) %>% # Select the feature and group.by columns.
        dplyr::group_by(!!rlang::sym(feature)) %>% # Group by feature first and then by group.by.
        dplyr::summarise(n = dplyr::n()) %>% # Compute the summarized counts by feature.
        dplyr::arrange(dplyr::desc(.data$n)) %>% # Pass on the values on group.by to a new variable that will store the X axis values.
        dplyr::pull(!!rlang::sym(feature)) %>%
        as.character()
    }
  } else if (!(is.null(order.by)) & !(is.null(group.by))){
    if (position == "fill"){
      # Obtain the order of the groups in feature (Y axis) according to one of the values in the X axis.
      factor_levels <- sample@meta.data %>% # Obtain metadata
        dplyr::select(!!rlang::sym(feature), !!rlang::sym(group.by)) %>% # Select the feature and group.by columns.
        dplyr::group_by(!!rlang::sym(group.by), !!rlang::sym(feature)) %>% # Group by feature first and then by group.by.
        dplyr::summarise(n = dplyr::n()) %>% # Compute the summarized counts by feature.
        dplyr::mutate(x_value = !!rlang::sym(group.by)) %>% # Pass on the values on group.by to a new variable that will store the X axis values.
        dplyr::filter(.data$x_value == order.by) %>% # Filter only the values equal to the value we want to reorder the bars.
        # Compute the total number of cells for each unique group in feature.
        dplyr::mutate(num_cells = {sample@meta.data %>% # Obtain metadata.
            dplyr::select(!!rlang::sym(feature), !!rlang::sym(group.by)) %>%  # Select the feature to plot.
            dplyr::group_by(!!rlang::sym(feature)) %>% # Group the values by feature.
            dplyr::summarise(n = dplyr::n()) %>% # Compute the total counts.
            # This line basically removes any row for which the value of order.by is 0. This avoids mismatches.
            dplyr::filter(!!rlang::sym(feature) %in% unique(sample@meta.data[, c(group.by, feature)][sample@meta.data[, c(group.by, feature)][, group.by] == order.by, ][, feature])) %>%
            dplyr::pull(.data$n)}) %>% # Retrieve the number of cells.
        dplyr::mutate(frac = .data$n/.data$num_cells) %>% # Compute the fraction that represent n out of the total number of cells in the group.
        dplyr::arrange(dplyr::desc(.data$frac)) %>% # Arrange it in descending order.
        dplyr::pull(!!rlang::sym(feature))  %>%
        as.character()
    } else if (position == "stack"){
      # Obtain the order of the groups in feature (Y axis) according to one of the values in the X axis.
      factor_levels <- sample@meta.data %>% # Obtain metadata
        dplyr::select(!!rlang::sym(feature), !!rlang::sym(group.by)) %>% # Select the feature and group.by columns.
        dplyr::group_by(!!rlang::sym(group.by), !!rlang::sym(feature)) %>% # Group by feature first and then by group.by.
        dplyr::summarise(n = dplyr::n()) %>% # Compute the summarized counts by feature.
        dplyr::mutate(x_value = !!rlang::sym(group.by)) %>% # Pass on the values on group.by to a new variable that will store the X axis values.
        dplyr::filter(.data$x_value == order.by) %>%
        dplyr::arrange(dplyr::desc(.data$n)) %>% # Arrange it in descending order.
        dplyr::pull(!!rlang::sym(feature))  %>%
        as.character()
    }

    # Retrieve the total number of unique values.
    total_levels <- unique(sample[[]][, feature])
    # If some are missing, add them back.
    if (length(factor_levels) != length(total_levels)){
      factor_levels <- c(factor_levels, total_levels[!(total_levels %in% factor_levels)])
    }
    factor_levels <- rev(factor_levels)
  } else if (is.null(order.by) & is.null(group.by)){
    if (position == "fill"){
      factor_levels = as.character(rev(sort(unique(sample@meta.data[, feature]))))
    } else if (position == "stack"){
      factor_levels <- sample@meta.data %>% # Obtain metadata
        dplyr::select(!!rlang::sym(feature)) %>% # Select the feature and group.by columns.
        dplyr::group_by(!!rlang::sym(feature)) %>% # Group by feature first and then by group.by.
        dplyr::summarise(n = dplyr::n()) %>%
        dplyr::arrange(dplyr::desc(.data$n)) %>%
        dplyr::pull(!!rlang::sym(feature)) %>%
        as.character()
    }
  }
  return(factor_levels)
}
