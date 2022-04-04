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

  pkg_list <- list("do_BarPlot" = c("Seurat", "colortools", "dplyr", "ggplot2", "ggpubr", "purrr", "rlang", "ggrepel"),
                   "do_ButterflyPlot" = c("Seurat", "tidyr", "pbapply", "dplyr", "ggplot2", "ggpubr", "viridis", "purrr", "rlang"),
                   "do_DimPlot" = c("colortools", "Seurat", "ggpubr", "ggplot2", "patchwork"),
                   "do_DotPlot" = c("Seurat", "ggplot2", "ggpubr"),
                   "do_FeaturePlot" = c("Seurat", "viridis", "ggplot2", "ggpubr", "patchwork", "scales"),
                   "do_NebulosaPlot" = c("Seurat", "ggplot2", "ggpubr", "Nebulosa", "patchwork"),
                   "do_PTEA" = c("Seurat", "stringr", "pbapply", "Matrix", "dplyr", "tidyr", "purrr", "rlang"),
                   "do_RankPlot" = c("Seurat", "ggplot2", "ggpubr", "viridis", "colortools", "ggbeeswarm"),
                   "do_VlnPlot" = c("Seurat", "ggplot2", "ggpubr"),
                   "savePlot" = c("ggplot2", "ComplexHeatmap", "grDevices", "svglite"))
  pkgs <- pkg_list[[function_name]]
  for (pkg in pkgs){
    if(!requireNamespace(pkg, quietly = T)){
      stop(paste0("Package ", pkg, " must be installed to use ", function_name, "."), call. = F)
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
    stop(paste0("The value/s for ", parameter_name, " is/are not a valid color representation. Please check whether it is an accepted R name or a HEX code."), call. = F)
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
  # Remove NAs.
  check_values <- check_values[!(is.na(check_values))]

  # Remove values that are not in the vector.
  if (sum(names(colors) %in% check_values) == length(check_values) & length(names(colors)) > length(check_values)){
    colors <- colors[names(colors) %in% check_values]
  }

  if (sum(names(colors) %in% check_values) != length(check_values)){
    stop('The names of provided colors does not match the number of unique values in the selected grouping variable (levels(object), group.by or split.by).', call. = F)
  }
  if (length(colors) != length(check_values)){
    stop('The number of provided colors is lower than the unique values in the selected grouping variable (levels(object), group.by or split.by).', call. = F)
  }
  return(colors)
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
  colors <- colortools::setColors("#0084a9", length(names_use))
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
    if (is.factor(sample@meta.data[, feature])){
      sample@meta.data[, feature] <- as.character(sample@meta.data[, feature])
    }
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
#' @param permissive Throw a warning or directly stops if the feature is not found.
#'
#' @return
#' @noRd
#' @examples
#' \dontrun{
#' TBD
#' }
check_feature <- function(sample, features, permissive = FALSE, dump_reduction_names = FALSE, enforce_check = NULL, enforce_parameter = NULL){
  if (is.list(features)){
    features_check <- unlist(features)
  } else {
    features_check <- features
  }
  check_enforcers <- list() # Store the results of the checks.
  not_found_features <- c() # Store the features not found.
  # Check each of the features.
  for (feature in features_check){
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
      not_found_features <- c(not_found_features, feature)
    }
  }
  # Return the error logs if there were features not found.
  if (length(not_found_features) > 0){
    if (isTRUE(permissive)){
      # Stop if neither of the features are found.
      if (length(unlist(not_found_features)) == length(unlist(features))){
        stop("Neither of the provided features are found.", call. = F)
      }
      warning(paste0("The requested features (",
                     not_found_features,
                     ") could not be found:\n",
                     "    - Not matching any gene name (rownames of the provided object).\n",
                     "    - Not matching any metadata column (in sample@meta.data).\n",
                     "    - Not part of the dimension names in any of the following reductions: ",
                     paste(Seurat::Reductions(object = sample), collapse = ", "),
                     "."), call. = F)
      features_out <- remove_not_found_features(features = features, not_found_features = not_found_features)

    } else if (isFALSE(permissive)){
      stop(paste0("The requested features (",
                  not_found_features,
                  ") could not be found:\n",
                  "    - Not matching any gene name (rownames of the provided object).\n",
                  "    - Not matching any metadata column (in sample@meta.data).\n",
                  "    - Not part of the dimension names in any of the following reductions: ",
                  paste(Seurat::Reductions(object = sample), collapse = ", "),
                  "."), call. = F)
    }
  } else {
    features_out <- features
  }
  # If we are enforcing a given check (i.e: the feature being in the metadata).
  if (!(is.null(enforce_check))){
    if (!(enforce_check %in% names(check_enforcers))){
      stop("The variable enforcer is not in the current list of checked variable types.", call. = F)
    } else {
      if (isFALSE(check_enforcers[[enforce_check]])){
        stop("The provided feature (", enforce_parameter, " = ", feature, ") not found in ", enforce_check, ".", call. = F)
      }
    }
  }

  # Return options.
  if (isTRUE(dump_reduction_names) & isFALSE(permissive)){return(dim_colnames)}
  if (isTRUE(permissive) & isFALSE(dump_reduction_names)){return(features_out)}
  if (isTRUE(dump_reduction_names) & isTRUE(permissive)){return(list("features" = features_out, "reduction_names" = dim_colnames))}
}

#' Remove not found features
#'
#' @param features Features to check.
#' @param not_found_features Features to exclude.
#'
#' @return
#' @noRd
#' @examples
#' \dontrun{
#' TBD
#' }
remove_not_found_features <- function(features, not_found_features){
  if (is.character(features)){
    features_out <- features[!(features %in% not_found_features)]
  } else if (is.list(features)){
    features_out <- list()
    for (list_name in names(features)){
      genes <- features[[list_name]]
      genes_out <- genes[!(genes %in% not_found_features)]
      features_out[[list_name]] <- genes_out
    }
  }
  return(features_out)
}

#' Remove duplicated features.
#'
#' @param features Features to check.
#'
#' @return
#' @noRd
#' @examples
#' \dontrun{
#' TBD
#' }
remove_duplicated_features <- function(features){
  if (is.character(features)){
    check <- sum(duplicated(features))
    if (check > 0){
      warning("Found duplicated features (", paste(features[duplicated(features)], collapse = ", "), "). Excluding them from the analysis.", call. = F)
      features <- features[!(duplicated(features))]
    }
  } else if (is.list(features)){
    features_out <- list()
    all_genes <- c() # Will update with the genes as they iterate to check duplicates.
    for (list_name in names(features)){
      genes <- features[[list_name]]
      # Remove genes duplicated within the list.
      if (sum(duplicated(genes)) > 0){
        warning("Found duplicated features (", paste(genes[duplicated(genes)], collapse = ", "), ") in the list '", list_name, "'. Excluding them from the analysis.", call. = F)
      }
      genes <- genes[!(duplicated(genes))]
      # Remove genes duplicated in the vector of all genes.
      duplicated_features <- genes[genes %in% all_genes]
      all_genes <- c(all_genes, genes[!(genes %in% all_genes)])
      genes <- genes[!(genes %in% duplicated_features)]
      if (length(duplicated_features) > 0){
        warning("Found duplicated features (", paste(duplicated_features, collapse = ", "), ") in list '", list_name, "' with regard to lists. Excluding them from the analysis.", call. = F)
      }
      features_out[[list_name]] <- genes
    }
    features <- features_out
  }
  return(features)
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
      stop(paste0("Could not find provided identity (", identity, ") in the current active identities of the object.\n Try running 'levels(your_seurat_object)' and see whether any typos were introduced."), call. = F)
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
  if (length(Seurat::Reductions(sample)) == 0){stop("This object has no reductions computed!", call. = F)}
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
    if (!(reduction %in% Seurat::Reductions(sample))){stop("The provided reduction could not be found in the object: ", reduction, call. = F)}
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
    stop("Provided dimensions need to be a 2-item vector.", call. = F)
  }
  # Check that the dimensions are integers.
  null_check <- is.null(dims[1]) & is.null(dims[2])
  integer_check <- is.numeric(dims[1]) & is.numeric(dims[1])
  if (!(is.null(dims)) & integer_check == FALSE){
    stop("Provied dimensions need to be numerics.", call. = F)
  }
  # Check that the dimensions are in the requested embedding.
  aval_dims <- length(colnames(Seurat::Embeddings(sample[[reduction]])))
  if (!(is.null(dims))){
    if (!(dims[1] %in% seq_len(aval_dims)) | !(dims[2] %in% seq_len(aval_dims))){
      stop("Dimension could not be found in the following reduction: ", reduction, ".", call. = F)
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
      stop("The value for assay has to be a character.", call. = F)
    }
    # Check that at least one assay is computed.
    if (length(Seurat::Assays(sample)) == 0){
      stop("There must be at least one computed assay in the object.", call. = F)
    }
    # Check that the assay is in the available assays.
    aval_assays <- Seurat::Assays(sample)
    if (!(assay %in% aval_assays)){
      stop("The following assay could not be found: ", assay, ".", call. = F)
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
    # Get each individual parameter from the list.
    parameter <- parameters[[parameter_name]]
    # Cases in which the user has to provide a vector.
    for (item in parameter){
      if (!(is.null(item))){
        if (!(is.na(item)) & !(test_function(item))){
          stop("Parameter ", parameter_name, " needs to be a ", required_type, ".", call. = F)
        }
      }
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
    stop("Only one of these 3 options can be passed to slot parameter: counts, data, scale.data.", call. = F)
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
    stop("The value provided for ", value_name, " (", value, ") is not in the range of the feature (", feature, "), which is: Min: ", limits[["scale.begin"]], ", Max: ", limits[["scale.end"]], ".", call. = F)
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
    if (is.null(position)){stop("Position parameter needs to be provided.", call. = F)}
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

#' Check viridis color map.
#'
#' @param viridis_color_map Viridis color map provided.
#' @param verbose Verbosity choice.
#'
#' @return
#' @noRd
#' @examples
#' \dontrun{
#' TBD
#' }
check_viridis_color_map <- function(viridis_color_map, verbose){
  viridis_options <- c("A", "B", "C", "D", "E", "F", "G", "H", "magma", "inferno", "plasma", "viridis", "cividis", "rocket", "mako", "turbo")
  if (!(viridis_color_map %in% viridis_options)){stop("The option provided to viridis_color_map is not an accepted option.\nPossible options: ", paste(viridis_options, collapse = ", "), call. = FALSE)}
  if (verbose){
    if (viridis_color_map %in% c("H", "turbo")){warning("The selected option is not the most adequate for a continuous color scale.", call. = F)}
  }
}


#' State SCpubr current function dependencies.
#'
#' @param func_name Name of an exported function from SCpubr. If NULL, return all functions.
#' @return None
#' @export
#'
#' @examples
#' \dontrun{
#' TBD
#' }
state_dependencies <- function(func_name = NULL){
  pkg_list <- list("do_BarPlot" = c("Seurat", "colortools", "dplyr", "ggplot2", "ggpubr", "purrr", "rlang"),
                   #"do_ButterflyPlot" = c("Seurat", "tidyr", "pbapply", "dplyr", "ggplot2", "ggpubr", "viridis", "purrr", "rlang"),
                   "do_DimPlot" = c("colortools", "Seurat", "ggpubr", "ggplot2", "patchwork"),
                   "do_DotPlot" = c("Seurat", "ggplot2", "ggpubr"),
                   "do_FeaturePlot" = c("Seurat", "viridis", "ggplot2", "ggpubr", "patchwork", "scales"),
                   "do_NebulosaPlot" = c("Seurat", "ggplot2", "ggpubr", "Nebulosa", "patchwork"),
                   #"do_PTEA" = c("Seurat", "stringr", "pbapply", "Matrix", "dplyr", "tidyr", "purrr", "rlang"),
                   "do_RankPlot" = c("Seurat", "ggplot2", "ggpubr", "viridis", "colortools", "ggbeeswarm"),
                   "do_VlnPlot" = c("Seurat", "ggplot2", "ggpubr"))

  cran_packages <- c("colortools",
                     "dplyr",
                     "ggbeeswarm",
                     "ggplot2",
                     "ggpubr",
                     "Matrix",
                     "patchwork",
                     "purrr",
                     "rlang",
                     "scales",
                     "Seurat",
                     "stringr",
                     "tidyr",
                     "viridis")

  bioconductor_packages <- c("Nebulosa")

  func_list <- sort(names(pkg_list))
  if (!(is.null(func_name))){
    for (func in func_name){
      if (!(func %in% func_list)){
        stop("Function name provided (", func, ") not part of SCpubr current functions.", call. = F)
      }
    }
    func_list <- func_name
  }

  message("\n---LIST OF PACKAGE DEPENDENCIES---\n")
  for (func in func_list){
    packages <- pkg_list[[func]]
    cran_packages_individual <- sort(packages[packages %in% cran_packages])
    bioconductor_packages_individual <- sort(packages[packages %in% bioconductor_packages])
      message("Dependencies for ", func, ":")
      if (length(cran_packages_individual >= 1)){message("  CRAN packages: ", paste(cran_packages_individual, collapse = ", "))}
      if (length(bioconductor_packages_individual >= 1)){message("  Bioconductor packages: ", paste(bioconductor_packages_individual, collapse = ", "))}
      message("")
  }
}

#' Check length of parameters compared to features.
#'
#' @param vector_of_parameters Vector of parameters to test.
#' @param vector_of_features  Vector of features to test against.
#' @param parameters_name Name of the parameters variable.
#' @param features_name Name of the features variable.

#' @return
#' @noRd
#' @examples
#' \dontrun{
#' TBD
#' }
check_length <- function(vector_of_parameters, vector_of_features, parameters_name, features_name){
  if (length(vector_of_parameters) != length(vector_of_features)){
    stop("Length of ", parameters_name, " not equal to ", features_name, ".", call. = F)
  }
}


#' Return a SC count matrix
#'
#' @return
#' @noRd
#' @examples
#' \dontrun{
#' TBD
#' }
use_dataset <- function(){
  # We want this function to be completely silent.
  suppressWarnings({
    sample <- CHETAH::headneck_ref
    sample <- Seurat::as.Seurat(sample, counts = "counts", data = NULL)
    sample <- suppressMessages(SeuratObject::RenameAssays(sample, originalexp = "RNA"))
    sample <- Seurat::PercentageFeatureSet(sample, pattern = "^MT-", col.name = "percent.mt")
    # Compute QC.
    mask1 <- sample$nCount_RNA >= 1000
    mask2 <- sample$nFeature_RNA >= 500
    mask3 <- sample$percent.mt <= 20
    mask <- mask1 & mask2 & mask3
    sample <- sample[, mask]
    # Normalize.
    sample <- Seurat::SCTransform(sample, verbose = FALSE)

    # Dimensional reduction.
    sample <- Seurat::RunPCA(sample, verbose = FALSE)
    sample <- Seurat::RunUMAP(sample, dims = 1:30, verbose = FALSE)
    # Find clusters.
    sample <- Seurat::FindNeighbors(sample, dims = 1:30, verbose = FALSE)
    sample <- Seurat::FindClusters(sample, resolution = 0.5, verbose = FALSE)
  })

  return(sample)
}

#' Add viridis color scale while suppressing the warning that comes with adding a second scale.
#'
#' @param p GGplot2 plot.
#' @param num_plots Number of plots.
#' @param function_use Coloring function to use.
#' @param scale Name of the scale. Either fill or color.
#' @param limits Whether to put limits.
#'
#' @return
#' @noRd
#' @examples
#' \dontrun{
#' TBD
#' }
add_scale <- function(p, scale, function_use, num_plots = 1, limits = NULL){
  if (scale == "color"){scale <- "colour"}
  # Compute the number of plots in this object (maybe a more efficient solution exists).
  if (num_plots == 1){
    # Find the index in which the scale is stored.
    # Adapted from: https://stackoverflow.com/a/46003178
    x <- which(sapply(p$scales$scales, function(x) scale %in% x$aesthetics))
    # Remove it.
    p$scales$scales[[x]] <- NULL
  } else {
    for (i in seq(1, num_plots)){
      # Find the index in which the scale is stored.
      # Adapted from: https://stackoverflow.com/a/46003178
      x <- which(sapply(p[[i]]$scales$scales, function(x) scale %in% x$aesthetics))
      # Remove it.
      p[[i]]$scales$scales[[x]] <- NULL
    }
  }
  # Add the scale and now it will now show up a warning since we removed the previous scale.
  p <- p & function_use
  return(p)
}


#' Compute a heatmap matrix out of expression data.
#'
#' @param sample Seurat object.
#' @param features Genes to retrieve the values from.
#' @param group.by Grouping variable to average the expression values by.
#' @param scale_features Whether to re-scale the features.
#' @param from_gene_expression Whether to retrieve values from gene expression.
#' @param from_metadata Whether to retrieve values from metadata.
#' @param feature_cutoff Whether to use a cutoff.
#' @param assay Assay to use.
#'
#' @return
#' @noRd
#' @examples
#' \dontrun{
#' TBD
#' }
produce_heatmap_matrix <- function(sample,
                                   features,
                                   group.by,
                                   scale_features,
                                   from_gene_expression = F,
                                   from_metadata = F,
                                   feature_cutoff = NULL,
                                   assay){
  `%>%`<- purrr::`%>%`
  if (isTRUE(from_gene_expression)){
    # Subset sample.
    sample <- sample[features, ]

    sample <- Seurat::ScaleData(sample, features = features, assay = assay)

    # Add scaled values as metadata.
    for (feature in features){
      sample@meta.data[, feature] <- sample@assays[[assay]]@scale.data[feature, ]
    }

    # Get the summarized Z-scores per gene per cluster.
    data <- sample@meta.data %>%
      dplyr::select(!!rlang::sym(group.by), dplyr::all_of(features)) %>%
      dplyr::group_by(!!rlang::sym(group.by)) %>%
      dplyr::summarise_at(dplyr::vars(dplyr::all_of(features)), list(mean))



    # Apply a cutoff if needed.
    if (!(is.null(feature_cutoff))){
      # Turn it into long format.
      data.long <- data %>%
                   tidyr::pivot_longer(cols = !(!!rlang::sym(group.by)), names_to = "gene") %>%
                   dplyr::arrange(dplyr::desc(.data$value))

      genes.keep <- data.long %>% dplyr::filter(.data$value >= feature_cutoff) %>% dplyr::pull(.data$gene)
      data.long <- data.long %>% dplyr::filter(.data$gene %in% unique(genes.keep)) %>% dplyr::arrange(.data$gene)

      # Turn it back to wide as this is the input for Hetmaps.
      data <- data.long %>%
              tidyr::pivot_wider(names_from = .data$gene, values_from = .data$value)
    }

    # Final formatting.
    data <- as.data.frame(data)
    rownames(data) <- data[, group.by]
    data[, group.by] <- NULL
    data <- as.matrix(data)
  }

  return(data)
}


#' Compute the data frame of the annotation for barplot annotation in heatmaps.
#'
#' @param sample Seurat object.
#' @param group.by Variable to group by.
#' @param annotation Annotation variable to use.
#'
#' @return
#' @noRd
#' @examples
#' \dontrun{
#' TBD
#' }
compute_barplot_annotation <- function(sample,
                                       group.by,
                                       annotation){
  `%>%`<- purrr::`%>%`
  # Compute column/row annotation. Obtain the percentage of a group per variable.
  annotation <- sample@meta.data %>%
                dplyr::select(!!rlang::sym(group.by), !!rlang::sym(annotation)) %>%
                dplyr::mutate(cluster = !!rlang::sym(group.by)) %>%
                dplyr::mutate(subgroup = !!rlang::sym(annotation)) %>%
                dplyr::select(.data$cluster, .data$subgroup) %>%
                dplyr::group_by(.data$cluster, .data$subgroup) %>%
                dplyr::summarise(n = dplyr::n()) %>%
                dplyr::mutate(freq = .data$n / sum(.data$n)) %>%
                dplyr::select(.data$cluster, .data$subgroup, .data$freq) %>%
                tidyr::pivot_wider(values_from = .data$freq, names_from = .data$subgroup)
  return(annotation)
}


#' Compute bar annotations for Heatmaps.
#'
#' @param data Data to plot.
#' @param annotation_use Data frame with the annotation values.
#' @param annotation_row Whether to plot row annotation.
#' @param annotation_column Wheter to plot column annotation.
#' @param cluster_columns,cluster_rows Logical. Whether to cluster the rows or the columns of the heatmap.
#' @param colors.annotation Vector of named colors to use it in the annotation.
#' @param fontsize General fontsize of the plot.
#' @param annotation_name Name for the legend annotation.
#'
#' @return
#' @noRd
#' @examples
#' \dontrun{
#' TBD
#' }
compute_bar_annotation <- function(data,
                                   annotation_use,
                                   annotation_row = FALSE,
                                   annotation_column = FALSE,
                                   cluster_columns,
                                   cluster_rows,
                                   colors.annotation,
                                   fontsize,
                                   annotation_name){
  `%>%`<- purrr::`%>%`
  if (isTRUE(annotation_row)){
    # Transform cluster column into level to rearrange it.
    annotation_use <- annotation_use %>%
                      dplyr::mutate(cluster = factor(.data$cluster, levels = rownames(data))) %>%
                      dplyr::arrange(.data$cluster)
    # Compute a test heatmap to retrieve row/cluster ordering if needed.
    return_list <- heatmap_inner(data = data,
                                 legend_name = "Z-Scores",
                                 cluster_columns = cluster_columns,
                                 cluster_rows = cluster_rows)
    h <- ComplexHeatmap::draw(return_list[["heatmap"]])
    order.rows.use <- rownames(data)[ComplexHeatmap::row_order(h)]
    order.cols.use <- colnames(data)[ComplexHeatmap::column_order(h)]

    # Row annotation.
    annotation_use <- as.data.frame(annotation_use)
    rownames(annotation_use) <- annotation_use$cluster
    annotation_use$cluster <- NULL
    annotation_use[is.na(annotation_use)] <- 0

    # Apply the reordering of the dataset according to the user's input.
    if (isTRUE(cluster_columns) & isFALSE(cluster_rows)){
      data <- data[, order.cols.use]
    } else if (isFALSE(cluster_columns) & isTRUE(cluster_rows)){
      data <- data[order.rows.use, ]
    } else if (isTRUE(cluster_columns) & isTRUE(cluster_rows)){
      data <- data[order.rows.use, order.cols.use]
    } else if (isFALSE(cluster_columns) & isFALSE(cluster_rows)){
      data <- data
    }
    # Reorder annotation.
    # Columns are reordered based on the color vector provided.
    # Rows are reordered if the rows were clusters.
    annotation_use <- annotation_use[, names(colors.annotation), drop = F]
    if (isTRUE(cluster_rows)){
      annotation_use <- annotation_use[order.rows.use, ]
    }

    # Compute row annotation legend.
    row_anno <- ComplexHeatmap::rowAnnotation(" " = ComplexHeatmap::anno_barplot(as.matrix(annotation_use),
                                                                                 gp = grid::gpar(fill = colors.annotation),
                                                                                 which = "row",
                                                                                 bar_width = 1,
                                                                                 height = grid::unit(1, "cm")))
    # Compute annotation legend.
    lgd_anno_row = ComplexHeatmap::Legend(labels = names(colors.annotation),
                                          labels_gp = grid::gpar(fontsize = fontsize,
                                                                 fontface = "bold"),
                                          title_gp = grid::gpar(fontsize = fontsize,
                                                                fontface = "bold"),
                                          legend_gp = grid::gpar(fill = colors.annotation),
                                          title = annotation_name)
    col_anno <- NULL
    lgd_anno_col <- NULL
  } else if (isTRUE(annotation_column)){
    # TBD
  }

  output_list <- list("data" = data,
                      "annotation_object_row" = row_anno,
                      "annotation_object_column" = col_anno,
                      "legend_object_row" = lgd_anno_row,
                      "legend_object_column" = lgd_anno_col)
  return(output_list)
}


#' Inner helper of \link[SCpubr]{do_Heatmap}
#'
#' @param data Matrix ready to be plotted. Use it alongside from_matrix.
#' @param legend_name Name of the general legend.
#' @param data_range One of:
#' - "both": Will compute a color scale equally balanced to both sides. Use when the values to plot are positive and negative.
#' - "only_pos": Will compute a color scale based only on the positive values. Will take the positive end of colors.use as well. Use when the values to plot are only positive.
#' - "only_neg": Will compute a color scale based only on the negative values. Will take the negative end of colors.use as well. Use when the values to plot are only negative
#' @param colors.use Vector of 2 colors defining a gradient. White color will be inserted in the middle.
#' @param grid_color Color for the grid.
#' @param fontsize General fontsize of the plot.
#' @param cell_size Size of each of the cells in the heatmap.
#' @param row_names_side,column_names_side Where to place the column or row names. "top", "bottom", "left", "right".
#' @param cluster_columns,cluster_rows Logical. Whether to cluster the rows or the columns of the heatmap.
#' @param border Logical. Whether to draw the border of the heatmap.
#' @param row_dendogram,column_dendogram Logical. Whether to plot row and column dendograms.
#' @param row_annotation,column_annotation Logical. Whether to place the annotation in the rows or the columns.
#' @param row_annotation_side,column_annotation_side Where to place the annotation. "top", "bottom", "left", "right".
#'
#' @return
#' @noRd
#' @examples
#' \dontrun{
#' TBD
#' }
heatmap_inner <- function(data,
                          legend_name = "Values",
                          data_range = "both",
                          colors.use = NULL,
                          grid_color = "grey50",
                          fontsize = 12,
                          cell_size = 5,
                          row_names_side = "left",
                          column_names_side = "bottom",
                          cluster_columns = TRUE,
                          cluster_rows = TRUE,
                          border = TRUE,
                          row_dendogram = FALSE,
                          column_dendogram = FALSE,
                          row_annotation = NULL,
                          row_annotation_side = "right",
                          column_annotation = NULL,
                          column_annotation_side = "top"){
  `%>%`<- purrr::`%>%`
  min_value <- min(data)
  max_value <- max(data)
  abs_value <- max(c(abs(min_value), abs(max_value)))

  if (is.null(colors.use)){
    colors.use <- c("#023f73", "white", "#7a0213")
  } else {
    colors.use <- c(colors.use[1], "white", colors.use[2])
  }
  if (data_range == "both"){
    breaks <-  round(c(-abs_value, (-abs_value / 2) , 0, (abs_value / 2), abs_value), 1)
    labels <- as.character(breaks)
    colors.use <- grDevices::colorRampPalette(colors.use)(length(breaks))
    names(colors.use) <- labels
    col_fun <- circlize::colorRamp2(breaks = breaks, colors = colors.use)
  } else if (data_range == "only_neg"){
    breaks <-  round(c(-abs_value, (-abs_value / 2) , 0), 1)
    labels <- as.character(breaks)
    colors.use <- grDevices::colorRampPalette(colors.use[c(1, 2)])(length(breaks))
    names(colors.use) <- labels
    col_fun <- circlize::colorRamp2(breaks = breaks, colors = colors.use[c(1, 2)])
  } else if (data_range == "only_pos"){
    breaks <-  round(c(0, (abs_value / 2), abs_value), 1)
    labels <- as.character(breaks)
    colors.use <- grDevices::colorRampPalette(colors.use[c(2, 3)])(length(breaks))
    names(colors.use) <- labels
    col_fun <- circlize::colorRamp2(breaks = breaks, colors = colors.use[c(2, 3)])
  }
  lgd = ComplexHeatmap::Legend(at = breaks,
                               labels = labels,
                               col_fun = col_fun,
                               title = legend_name,
                               break_dist = rep(1, length(breaks) - 1),
                               labels_gp = grid::gpar(fontsize = fontsize,
                                                      fontface = "bold"),
                               title_gp = grid::gpar(fontsize = fontsize,
                                                     fontface = "bold"))


  if (!(is.null(row_annotation))){
    if (row_annotation_side == "right"){
      right_annotation <- row_annotation
      left_annotation <- NULL
    } else {
      right_annotation <- NULL
      left_annotation <- row_annotation
    }
  } else {
    right_annotation <- NULL
    left_annotation <- NULL
  }

  if (!(is.null(column_annotation))){
    if (column_annotation_side == "top"){
      top_annotation <- column_annotation
      bottom_annotation <- NULL
    } else {
      top_annotation <- NULL
      bottom_annotation <- column_annotation
    }
  } else {
    top_annotation <- NULL
    bottom_annotation <- NULL
  }


  h <- ComplexHeatmap::Heatmap(matrix = data,
                               name = legend_name,
                               col = col_fun,
                               show_heatmap_legend = FALSE,
                               cluster_rows = cluster_rows,
                               cluster_columns = cluster_columns,
                               show_row_dend = row_dendogram,
                               show_column_dend = column_dendogram,
                               top_annotation = top_annotation,
                               bottom_annotation = bottom_annotation,
                               right_annotation = right_annotation,
                               left_annotation = left_annotation,
                               width = ncol(data)*grid::unit(cell_size, "mm"),
                               height = nrow(data)*grid::unit(cell_size, "mm"),
                               column_names_gp = grid::gpar(fontsize = fontsize,
                                                            fontface = "bold"),
                               row_names_gp = grid::gpar(fontsize = fontsize,
                                                         fontface = "bold"),
                               row_names_side = row_names_side,
                               column_names_side = column_names_side,
                               column_title_gp = grid::gpar(fontsize = fontsize,
                                                            fontface = "bold"),
                               border = border,
                               rect_gp = grid::gpar(col= grid_color),
                               cell_fun = function(j, i, x, y, w, h, fill) {
                                 grid::grid.rect(x, y, w, h, gp = grid::gpar(alpha = 0.25))
                               })

  return_list <- list("heatmap" = h,
                      "legend" = lgd)

  return(return_list)
}



