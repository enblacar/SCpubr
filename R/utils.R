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

  pkg_list <- list("do_BarPlot" = c("Seurat", "colorspace", "dplyr", "ggplot2", "purrr", "rlang", "ggrepel", "ggtext"),
                   "do_CellularStatesPlot" = c("Seurat", "tidyr", "pbapply", "dplyr", "ggplot2", "viridis", "purrr", "rlang", "ggExtra", "ggplotify", "scattermore"),
                   "do_DimPlot" = c("colorspace", "Seurat", "ggplot2", "patchwork", "ggtext", "ggplotify", "scattermore"),
                   "do_DotPlot" = c("Seurat", "ggplot2", "ggtext"),
                   "do_FeaturePlot" = c("Seurat", "viridis", "ggplot2", "patchwork", "scales", "ggtext", "scattermore"),
                   "do_NebulosaPlot" = c("Seurat", "ggplot2", "Nebulosa", "patchwork", "ggtext"),
                   "do_PTEA" = c("Seurat", "stringr", "pbapply", "Matrix", "dplyr", "tidyr", "purrr", "rlang"),
                   "do_BeeSwarmPlot" = c("Seurat", "ggplot2", "viridis", "colorspace", "ggbeeswarm", "ggrastr", "ggtext"),
                   "do_VlnPlot" = c("Seurat", "ggplot2", "ggtext"),
                   "save_Plot" = c("ggplot2", "ComplexHeatmap", "grDevices", "svglite"),
                   "do_TermEnrichmentPlot" = c("ggplot2", "enrichR", "stringr", "dplyr", "patchwork", "forcats", "ggtext"),
                   "do_EnrichmentHeatmap" = c("ggplot2", "stringr", "dplyr", "patchwork", "purrr", "ComplexHeatmap", "Seurat", "rlang", "grDevices", "circlize", "grid"),
                   "do_CorrelationPlot" = c("ComplexHeatmap", "purrr", "Seurat", "rlang", "ggplot2", "patchwork", "dplyr", "grDevices", "ComplexHeatmap", "circlize", "grid"),
                   "do_LigandReceptorPlot" = c("stringr", "Seurat", "liana", "dplyr", "rlang", "tibble", "tidyr", "ggplot2", "ggtext", "purrr"),
                   "do_CopyNumberVariantPlot" = c("purrr", "dplyr", "tibble", "ggplot2", "ggdist", "rlang", "ggtext"),
                   "do_PseudotimePlot" = c("monocle3", "purrr", "ggplot2", "dplyr", "ggdist", "ggtext", "patchwork"),
                   "do_GeyserPlot" = c("purrr", "Seurat", "dplyr", "tibble", "ggplot2", "ggdist", "ggtext"),
                   "do_TFActivityPlot" = c("ComplexHeatmap", "purrr", "dplyr", "tidyr", "tibble", "Seurat", "stats", "ggplot2", "grDevices", "rlang"),
                   "do_PathwayActivityPlot" = c("ComplexHeatmap", "purrr", "dplyr", "tidyr", "tibble", "Seurat", "stats", "ggplot2", "grDevices", "rlang"),
                   "do_GroupwiseDEPlot" = c("ComplexHeatmap", "purrr", "dplyr", "tidyr", "tibble", "Seurat", "grDevices", "rlang", "plyr"),
                   "testing" = c("Does_not_exist"))
  # The function is not in the current list of possibilities.
  if (!(function_name %in% names(pkg_list))){
    stop(paste0(function_name, " is not an accepted function name."), call. = FALSE)
  }
  pkgs <- pkg_list[[function_name]]
  for (pkg in pkgs){
    if(!requireNamespace(pkg, quietly = T)){
      stop(paste0("Package ", pkg, " must be installed to use ", function_name, "."), call. = F)
    }
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
  pkg_list <- list("do_BarPlot" = c("Seurat", "colorspace", "dplyr", "ggplot2", "purrr", "rlang", "ggrepel", "ggtext"),
                   "do_CellularStatesPlot" = c("Seurat", "tidyr", "pbapply", "dplyr", "ggplot2", "viridis", "purrr", "rlang", "ggExtra", "ggplotify", "scattermore"),
                   "do_DimPlot" = c("colorspace", "Seurat", "ggplot2", "patchwork", "ggtext", "ggplotify", "scattermore"),
                   "do_DotPlot" = c("Seurat", "ggplot2", "ggtext"),
                   "do_FeaturePlot" = c("Seurat", "viridis", "ggplot2", "patchwork", "scales", "ggtext", "scattermore"),
                   "do_NebulosaPlot" = c("Seurat", "ggplot2", "Nebulosa", "patchwork", "ggtext"),
                   "do_PTEA" = c("Seurat", "stringr", "pbapply", "Matrix", "dplyr", "tidyr", "purrr", "rlang"),
                   "do_BeeSwarmPlot" = c("Seurat", "ggplot2", "viridis", "colorspace", "ggbeeswarm", "ggrastr", "ggtext"),
                   "do_VlnPlot" = c("Seurat", "ggplot2", "ggtext"),
                   "save_Plot" = c("ggplot2", "ComplexHeatmap", "grDevices", "svglite"),
                   "do_TermEnrichmentPlot" = c("ggplot2", "enrichR", "stringr", "dplyr", "patchwork", "forcats", "ggtext"),
                   "do_EnrichmentHeatmap" = c("ggplot2", "stringr", "dplyr", "patchwork", "purrr", "ComplexHeatmap", "Seurat", "rlang", "grDevices", "circlize", "grid"),
                   "do_CorrelationPlot" = c("ComplexHeatmap", "purrr", "Seurat", "rlang", "ggplot2", "patchwork", "dplyr", "grDevices", "ComplexHeatmap", "circlize", "grid"),
                   "do_LigandReceptorPlot" = c("stringr", "Seurat", "liana", "dplyr", "rlang", "tibble", "tidyr", "ggplot2", "ggtext", "purrr"),
                   "do_CopyNumberVariantPlot" = c("purrr", "dplyr", "tibble", "ggplot2", "ggdist", "rlang", "ggtext"),
                   "do_PseudotimePlot" = c("monocle3", "purrr", "ggplot2", "dplyr", "ggdist", "ggtext", "patchwork"),
                   "do_GeyserPlot" = c("purrr", "Seurat", "dplyr", "tibble", "ggplot2", "ggdist", "ggtext"),
                   "do_TFActivityPlot" = c("ComplexHeatmap", "purrr", "dplyr", "tidyr", "tibble", "Seurat", "stats", "ggplot2", "grDevices", "rlang"),
                   "do_PathwayActivityPlot" = c("ComplexHeatmap", "purrr", "dplyr", "tidyr", "tibble", "Seurat", "stats", "ggplot2", "grDevices", "rlang"),
                   "do_GroupwiseDEPlot" = c("ComplexHeatmap", "purrr", "dplyr", "tidyr", "tibble", "Seurat", "grDevices", "rlang"))
  # The function is not in the current list of possibilities.
  if (!(is.null(func_name))){
    for (func in func_name){
      if (!(func %in% names(pkg_list))){
        stop(paste0(func_name, " is not an accepted function name."), call. = FALSE)
      }
    }
  }

  cran_packages <- c("colorspace",
                     "dplyr",
                     "enrichR",
                     "forcats",
                     "ggbeeswarm",
                     "ggdist",
                     "ggplot2",
                     "ggExtra",
                     "ggplotify",
                     "ggrepel",
                     "ggtext",
                     "Matrix",
                     "patchwork",
                     "purrr",
                     "rlang",
                     "scales",
                     "scattermore",
                     "Seurat",
                     "stringr",
                     "svglite",
                     "tibble",
                     "tidyr",
                     "viridis")

  bioconductor_packages <- c("infercnv",
                             "Nebulosa")

  github_packages <- c("liana")

  func_list <- sort(names(pkg_list))
  if (!(is.null(func_name))){
    func_list <- func_name
  } else {
    func_list <- names(pkg_list)
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



#' Check for Seurat class.
#'
#' @param sample Seurat object.
#'
#' @noRd
#' @return None
#'
#' @examples
#' \dontrun{
#' TBD
#' }
check_Seurat <- function(sample){
  if (isFALSE("Seurat" %in% class(sample))){
    stop("Object provided is not a Seurat object.", call. = F)
  }
}

#' Internal check for colors.
#'
#' Adapted from: https://stackoverflow.com/a/13290832.
#
#' @param colors Vector of colors.
#' @param parameter_name The name of the parameter for which we are testing the colors.
#' @return None
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
#' @return None
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

  if (length(colors) != length(check_values)){
    stop('The number of provided colors is lower than the unique values in the selected grouping variable (levels(object), group.by or split.by).', call. = F)
  }

  if (sum(names(colors) %in% check_values) != length(check_values)){
    stop('The names of provided colors does not match the number of unique values in the selected grouping variable (levels(object), group.by or split.by).', call. = F)
  }

  return(colors)
}

#' Generate custom color scale.
#'
#' @param names_use Vector of the names that will go alongside the color scale.
#' @return None
#' @noRd
#' @examples
#' \dontrun{
#' TBD
#' }
generate_color_scale <- function(names_use){
  # Generate a vector of colors equal to the number of identities in the sample.
  colors <- colorspace::qualitative_hcl(length(names_use), palette = "Dark 3")
  colors <- grDevices::col2rgb(colors)
  colors <- grDevices::rgb2hsv(colors)
  colors["v",] <- colors["v", ] - 0.1
  colors["s",] <- colors["s", ] + 0.2
  colors["s",][colors["s",] > 1] <- 1
  colors <- grDevices::hsv(h = colors["h", ],
                           s = colors["s",],
                           v = colors["v",],
                           alpha = 1)
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
#' @return None
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

#' Check if a value is in the range of the values.
#'
#' @param sample Seurat object.
#' @param feature Feature to plot.
#' @param assay Assay used.
#' @param reduction Reduction used.
#' @param value Value to check.
#' @param value_name Name of the value.
#'
#' @return None
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

#' Check if the feature to plot is in the Seurat object.
#'
#' @param sample Seurat object.
#' @param features Feature to plot.
#' @param dump_reduction_names Whether to return the reduction colnames.
#' @param permissive Throw a warning or directly stops if the feature is not found.
#'
#' @return None
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
    } else {
      check_enforcers[["gene"]] <- TRUE
    }

    if (!(feature %in% colnames(sample@meta.data))){
      check <- check + 1
      check_enforcers[["metadata"]] <- FALSE
    } else {
      check_enforcers[["metadata"]] <- TRUE
    }

    dim_colnames <- c()
    for(red in Seurat::Reductions(object = sample)){
      dim_colnames <- c(dim_colnames, colnames(sample@reductions[[red]][[]]))
    }
    if (!(feature %in% dim_colnames)){
      check <- check + 1
      check_enforcers[["reductions"]] <- FALSE
    } else {
      check_enforcers[["reductions"]] <- TRUE
    }

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
                     ".\n\n"), call. = F)
      features_out <- remove_not_found_features(features = features, not_found_features = not_found_features)

    } else if (isFALSE(permissive)){
      stop(paste0("The requested features (",
                  not_found_features,
                  ") could not be found:\n",
                  "    - Not matching any gene name (rownames of the provided object).\n",
                  "    - Not matching any metadata column (in sample@meta.data).\n",
                  "    - Not part of the dimension names in any of the following reductions: ",
                  paste(Seurat::Reductions(object = sample), collapse = ", "),
                  ".\n\n"), call. = F)
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
#' @return None
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
#' @return None
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
#' @return None
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
#' @return None
#' @noRd
#' @examples
#' \dontrun{
#' TBD
#' }
check_and_set_reduction <- function(sample, reduction = NULL){
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
#' @return None
#' @noRd
#' @examples
#' \dontrun{
#' TBD
#' }
check_and_set_dimensions <- function(sample, reduction = NULL, dims = NULL){
  # Check that the dimensions is a 2 item vector.
  if (!(is.null(dims)) & length(dims) != 2){
    stop("Provided dimensions need to be a 2-item vector.", call. = F)
  }

  # If reduction is null, select the last computed one.
  if (is.null(reduction)){
    reduction <- Seurat::Reductions(sample)[length(Seurat::Reductions(sample))]
  }

  # Check that at least 2 dimensions are present.
  aval_dims <- length(colnames(Seurat::Embeddings(sample[[reduction]])))
  if (aval_dims < 2){
    stop("There are less than 2 dimensions in the requested reduction: ", reduction, ".")
  }

  # Check that the dimensions are integers.
  null_check <- is.null(dims[1]) & is.null(dims[2])
  integer_check <- is.numeric(dims[1]) & is.numeric(dims[1])
  if (!(is.null(dims)) & integer_check == FALSE){
    stop("Provied dimensions need to be numerics.", call. = F)
  }
  # Check that the dimensions are in the requested embedding.
  if (!(is.null(dims))){
    if (!(dims[1] %in% seq_len(aval_dims)) | !(dims[2] %in% seq_len(aval_dims))){
      stop("Dimension could not be found in the following reduction: ", reduction, ".", call. = F)
    }
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
#' @return None
#' @noRd
#' @examples
#' \dontrun{
#' TBD
#' }
check_and_set_assay <- function(sample, assay = NULL){
  # Check that at least one assay is computed.
  if (length(Seurat::Assays(sample)) == 0){
    stop("There must be at least one computed assay in the object.", call. = F)
  }
  # If assay is null, set it to the active one.
  if (is.null(assay)){
    assay <- Seurat::DefaultAssay(sample)
  } else {
    # Check if the assay is a character.
    if (!(is.character(assay))){
      stop("The value for assay has to be a character.", call. = F)
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
#' @return None
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
    # Check if the parameter is not NULL already.
    if (!(is.null(parameter))){
      # For each parameter in the vector.
      for (item in parameter){
        # If not null.
        if (!(is.null(item))){
          # If not NA, if the testing function fails, report it.
          if (sum(!(is.na(item))) > 0){
            if (sum(!(test_function(item))) > 0){
              stop("Parameter ", parameter_name, " needs to be a ", required_type, ".", call. = F)
            }
          }
        }
      }
    }
  }
}

#' Check the slots.
#'
#' @param slot Slot provided.
#'
#' @return None
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


#' Compute the order of the plotted bars for do_BarPlot.
#'
#' @param sample Seurat object.
#' @param feature Feature to plot.
#' @param group.by Feature to group the output by.
#' @param order.by Unique value in group.by to reorder labels in descending order.
#'
#' @return None
#' @noRd
#' @examples
#' \dontrun{
#' TBD
#' }
compute_factor_levels <- function(sample, feature, position, group.by = NULL, order.by = NULL){
  `%>%` <- purrr::`%>%`
  if (!(position %in% c("stack", "fill"))){stop("Position needs to be either stack or fill.", call. = F)}
  if (is.null(order.by) & !(is.null(group.by))){
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
#' @return None
#' @noRd
#' @examples
#' \dontrun{
#' TBD
#' }
check_viridis_color_map <- function(viridis_color_map, verbose = F){
  viridis_options <- c("A", "B", "C", "D", "E", "F", "G", "H", "magma", "inferno", "plasma", "viridis", "cividis", "rocket", "mako", "turbo")
  if (!(viridis_color_map %in% viridis_options)){stop("The option provided to viridis_color_map is not an accepted option.\nPossible options: ", paste(viridis_options, collapse = ", "), call. = FALSE)}
  if (verbose){
    if (viridis_color_map %in% c("H", "turbo")){warning("The selected option is not the most adequate for a continuous color scale.", call. = F)}
  }
}




#' Check length of parameters compared to features.
#'
#' @param vector_of_parameters Vector of parameters to test.
#' @param vector_of_features  Vector of features to test against.
#' @param parameters_name Name of the parameters variable.
#' @param features_name Name of the features variable.

#' @return None
#' @noRd
#' @examples
#' \dontrun{
#' TBD
#' }
check_length <- function(vector_of_parameters,
                         vector_of_features,
                         parameters_name,
                         features_name){
  if (length(vector_of_parameters) != length(vector_of_features)){
    stop("Length of ", parameters_name, " not equal to ", features_name, ".", call. = F)
  }
}


#' Return a SC count matrix
#'
#' @return None
#' @noRd
#' @examples
#' \dontrun{
#' TBD
#' }
use_dataset <- function(n_cells = 180){
  # We want this function to be completely silent.
  suppressWarnings({
    test_list <- get0("test_list", envir = asNamespace("SCpubr"))
    genes <- test_list$genes
    values <- seq(0, 15, 0.1)
    counts <- matrix(ncol = n_cells, nrow = length(genes))
    cols <- c()
    for (i in seq(1, n_cells)){
      cts <- sample(values, size = length(genes), replace = T, prob = c(0.66, rep((0.34 / 150), length(values) - 1)))
      counts[, i] <- cts
      cols <- c(cols, paste0("Cell_", i))
    }
    rownames(counts) <- genes
    colnames(counts) <- cols
    sample <- Seurat::CreateSeuratObject(counts)
    sample <- Seurat::PercentageFeatureSet(sample, pattern = "^MT-", col.name = "percent.mt")
    # Compute QC.
    mask1 <- sample$nCount_RNA >= 1000
    mask2 <- sample$nFeature_RNA >= 500
    mask3 <- sample$percent.mt <= 20
    mask <- mask1 & mask2 & mask3
    sample <- sample[, mask]
    # Normalize.
    sample <- suppressWarnings({Seurat::SCTransform(sample, verbose = FALSE)})

    # Dimensional reduction.
    sample <- Seurat::RunPCA(sample, verbose = FALSE)
    sample <- Seurat::RunUMAP(sample, dims = 1:30, verbose = FALSE)
    # Find clusters.
    sample <- Seurat::FindNeighbors(sample, dims = 1:30, verbose = FALSE)
    sample <- Seurat::FindClusters(sample, resolution = 0.5, verbose = FALSE)
    sample$seurat_clusters <- as.character(sample$seurat_clusters)
    sample$seurat_clusters[1:20] <- "0"
    sample$seurat_clusters[21:40] <- "1"
    sample$seurat_clusters[41:60] <- "2"
    sample$seurat_clusters[61:80] <- "3"
    sample$seurat_clusters[81:100] <- "4"
    sample$seurat_clusters[101:120] <- "5"
    sample$seurat_clusters[121:140] <- "6"
    sample$seurat_clusters[141:160] <- "7"
    sample$seurat_clusters[161:180] <- "8"
    Seurat::Idents(sample) <- sample$seurat_clusters
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
#' @return None
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




#' Compute the data frame of the annotation for barplot annotation in heatmaps.
#'
#' @param sample Seurat object.
#' @param group.by Variable to group by.
#' @param annotation Annotation variable to use.
#'
#' @return None
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




#' Inner helper for heatmaps
#'
#' @param data Matrix ready to be plotted. Use it alongside from_matrix.
#' @param legend_name Name of the general legend.
#' @param data_range One of:
#' - "both": Will compute a color scale equally balanced to both sides. Use when the values to plot are positive and negative.
#' - "only_pos": Will compute a color scale based only on the positive values. Will take the positive end of colors.use as well. Use when the values to plot are only positive.
#' - "only_neg": Will compute a color scale based only on the negative values. Will take the negative end of colors.use as well. Use when the values to plot are only negative
#' @param colors.use Vector of 2 colors defining a gradient. White color will be inserted in the middle.
#' @param grid_color Color for the grid.
#' @param range.data Numeric. Min or max value (data_range = "only_pos" or "only_neg") or vector of min and max (data_range = "both") that will determine the span of the color scale.
#' @param outlier.data Logical. Whether there is outlier data to take into account.
#' @param fontsize General fontsize of the plot.
#' @param cell_size Size of each of the cells in the heatmap.
#' @param row_names_side,column_names_side Where to place the column or row names. "top", "bottom", "left", "right".
#' @param cluster_columns,cluster_rows Logical. Whether to cluster the rows or the columns of the heatmap.
#' @param border Logical. Whether to draw the border of the heatmap.
#' @param row_dendogram,column_dendogram Logical. Whether to plot row and column dendograms.
#' @param row_annotation,column_annotation Annotation objects.
#' @param row_annotation_side,column_annotation_side Where to place the annotation. "top", "bottom", "left", "right".
#' @param row_title,column_title Titles for the axes.
#' @param column_title_side,row_title_side Side for the titles.
#' @param column_title_rotation,row_title_rotation Angle of rotation of the titles.
#' @param row_names_rot,column_names_rot Angle of rotation of the text.
#' @param legend.framecolor Color of the lines of the box in the legend.
#' @param legend.length,legend.width Length and width of the legend. Will adjust automatically depending on legend side.
#' @param na.value Color for NAs
#' @param use_viridis Logical. Whether to use viridis color palettes.
#' @param viridis_color_map Character. Palette to use.
#' @param viridis_direction Numeric. Direction of the scale.
#' @param zeros_are_white Logical.
#' @return None
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
                          range.data = NULL,
                          outlier.data = FALSE,
                          outlier.up.color = "#4b010b",
                          outlier.down.color = "#02294b",
                          outlier.up.label = NULL,
                          outlier.down.label = NULL,
                          round_value_outlier = 2,
                          column_title = NULL,
                          row_title = NULL,
                          row_names_side = "left",
                          column_names_side = "bottom",
                          cluster_columns = TRUE,
                          cluster_rows = TRUE,
                          border = TRUE,
                          legend.position = "bottom",
                          legend.length = 20,
                          legend.width = 1,
                          legend.framecolor = "grey50",
                          row_dendogram = FALSE,
                          column_dendogram = FALSE,
                          column_title_side = "top",
                          row_title_rotation = 90,
                          column_title_rotation = 0,
                          row_names_rot = 0,
                          column_names_rot = 90,
                          row_title_side = "left",
                          row_annotation = NULL,
                          row_annotation_side = "right",
                          column_annotation = NULL,
                          column_annotation_side = "top",
                          na.value = "grey75",
                          use_viridis = FALSE,
                          viridis_color_map = "D",
                          viridis_direction = 1,
                          zeros_are_white = FALSE,
                          symmetrical_scale = FALSE){
  `%>%`<- purrr::`%>%`


  if (legend.position %in% c("top", "bottom")){
    legend_width <- grid::unit(legend.length, "mm")
    legend_height <- NULL
    grid_height <- grid::unit(legend.width, "mm")
    grid_width <- grid::unit(4, "mm")
    direction <- "horizontal"
    title_position <- "topcenter"
  } else if (legend.position %in% c("left", "right")){
    grid_width <- grid::unit(legend.width, "mm")
    legend_height <- grid::unit(legend.length, "mm")
    legend_width <- NULL
    grid_height <- grid::unit(4, "mm")
    direction <- "vertical"
    title_position <- "topleft"
  }

  if (!is.null(range.data)){

    if (data_range == "both"){
      if (isTRUE(symmetrical_scale)){
        abs_value <- max(abs(range.data))
        q100 <- abs(abs_value)
        q0 <- -abs(abs_value)
      } else {
        q0 <- range.data[1]
        q100 <- range.data[2]
      }
    } else if (data_range == "only_pos"){
      abs_value <- abs(range.data)
      q0 <- 0
      q100 <- abs_value
    } else if (data_range == "only_neg"){
      abs_value <- abs(range.data)
      q100 <- 0
      q0 <- abs_value
    }
  } else {
    q0 <- min(data)
    q100 <- max(data)
    abs_value <- max(c(abs(q0), abs(q100)))
  }

  q50 <- mean(c(q0, q100))
  q25 <- mean(c(q0, q50))
  q75 <- mean(c(q50, q100))

  # Checks.
  if (data_range == "only_neg" & q0 >= 0){
    stop("There are no negative values in the matrix.")
  }

  if (data_range == "only_pos" & q100 < 0){
    stop("There are no positive values in the matrix.")
  }


  if (is.null(colors.use)){
    colors.use <- c("#023f73", "white", "#7a0213")
  } else {
    colors.use <- c(colors.use[1], "white", colors.use[2])
  }
  if (data_range == "both"){
    if (isTRUE(symmetrical_scale)){
      breaks <-  round(c(-abs_value, (-abs_value / 2) , 0, (abs_value / 2), abs_value), 1)
      counter <- 0
      while (sum(duplicated(breaks)) > 0){
        counter <- counter + 1
        breaks <-  round(c(-abs_value, (-abs_value / 2) , 0, (abs_value / 2), abs_value), 1 + counter)
      }
    } else if (isFALSE(symmetrical_scale)){
      breaks <-  round(c(q0, q25, q50, q75, q100), 1)
      counter <- 0
      while (sum(duplicated(breaks)) > 0){
        counter <- counter + 1
        breaks <-  round(c(q0, q25, q50, q75, q100), 1 + counter)
      }
    }
    labels <- as.character(breaks)
    colors.use <- grDevices::colorRampPalette(colors.use)(length(breaks))
    if (isTRUE(outlier.data) & !is.null(range.data)){
      breaks <- c(-abs_value - 0.00001, breaks, abs_value + 0.00001)
      colors.use <- c(outlier.down.color, colors.use, outlier.up.color)
      labels <- c(if(is.null(outlier.down.label)){paste0("< ", -round(abs_value, round_value_outlier))} else {outlier.down.label},
                  labels,
                  if(is.null(outlier.up.label)){paste0("> ", -round(abs_value, round_value_outlier))} else {outlier.up.label})
    }

    names(colors.use) <- labels
  } else if (data_range == "only_neg"){
    if (isTRUE(zeros_are_white)){
      breaks <-  round(c(-abs_value, (-abs_value * 0.75), (-abs_value * 0.5), (-abs_value * 0.25), (-abs_value * 0.01), 0), 1)
      counter <- 0
      while (sum(duplicated(breaks)) > 0){
        counter <- counter + 1
        breaks <-  round(c(-abs_value, (-abs_value * 0.75), (-abs_value * 0.5), (-abs_value * 0.25), (-abs_value * 0.01), 0), 1 + counter)
      }
    } else {
      breaks <-  round(c(-abs_value, (-abs_value * 0.75), (-abs_value * 0.5), (-abs_value * 0.25), 0), 1)
      counter <- 0
      while (sum(duplicated(breaks)) > 0){
        counter <- counter + 1
        breaks <-  round(c(-abs_value, (-abs_value * 0.75), (-abs_value * 0.5), (-abs_value * 0.25), 0), 1 + counter)
      }
    }
    labels <- as.character(breaks)
    colors.use <- grDevices::colorRampPalette(colors.use[c(1, 2)])(length(breaks))
    if (isTRUE(outlier.data) & !is.null(range.data)){
      breaks <- c(-abs_value - 0.00001, breaks)
      colors.use <- c(outlier.down.color, colors.use)
      labels <- c(if(is.null(outlier.down.label)){paste0("< ", -round(abs_value, round_value_outlier))} else {outlier.down.label},
                  labels)
    }
    names(colors.use) <- labels
  } else if (data_range == "only_pos"){
    if (isTRUE(zeros_are_white)){
      breaks <-  round(c(0, (abs_value * 0.01), (abs_value * 0.25), (abs_value * 0.5), (abs_value * 0.75), abs_value), 1)
      counter <- 0
      while (sum(duplicated(breaks)) > 0){
        counter <- counter + 1
        breaks <-  round(c(0, (abs_value * 0.01), (abs_value * 0.25), (abs_value * 0.5), (abs_value * 0.75), abs_value), 1 + counter)
      }
    } else {
      breaks <-  round(c(0, (abs_value * 0.25), (abs_value * 0.5), (abs_value * 0.75), abs_value), 1)
      counter <- 0
      while (sum(duplicated(breaks)) > 0){
        counter <- counter + 1
        breaks <-  round(c(0, (abs_value * 0.25), (abs_value * 0.5), (abs_value * 0.75), abs_value), 1 + counter)
      }
    }
    labels <- as.character(breaks)
    colors.use <- grDevices::colorRampPalette(colors.use[c(2, 3)])(length(breaks))
    if (isTRUE(outlier.data) & !is.null(range.data)){
      breaks <- c(breaks, abs_value + 0.00001)
      colors.use <- c(colors.use, outlier.up.color)
      labels <- c(labels,
                  if(is.null(outlier.up.label)){paste0("> ", -round(abs_value, round_value_outlier))} else {outlier.up.label})
    }
    names(colors.use) <- labels
  }

  if (isTRUE(use_viridis)){
    if (isTRUE(zeros_are_white) & data_range %in% c("only_pos", "only_neg")){
      col_fun <- circlize::colorRamp2(breaks = breaks, colors = c("white", viridis::viridis(n = length(breaks) - 1,
                                                                                            option = viridis_color_map,
                                                                                            direction = viridis_direction)))
    } else {
      col_fun <- circlize::colorRamp2(breaks = breaks, colors = viridis::viridis(n = length(breaks),
                                                                                 option = viridis_color_map,
                                                                                 direction = viridis_direction))
    }
  } else {
    col_fun <- circlize::colorRamp2(breaks = breaks, colors = colors.use)
  }



  lgd = ComplexHeatmap::Legend(at = breaks,
                               labels = labels,
                               col_fun = col_fun,
                               title = legend_name,
                               direction = direction,
                               legend_height = legend_height,
                               legend_width = legend_width,
                               grid_width = grid_width,
                               grid_height = grid_height,
                               border = legend.framecolor,
                               title_position = title_position,
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
                               na_col = na.value,
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
                               column_title = column_title,
                               column_title_side = column_title_side,
                               row_title_side = row_title_side,
                               row_title = row_title,
                               column_title_rot = column_title_rotation,
                               row_title_rot = row_title_rotation,
                               column_names_rot = column_names_rot,
                               row_names_rot = row_names_rot,
                               column_title_gp = grid::gpar(fontsize = fontsize,
                                                            fontface = "bold"),
                               row_title_gp = grid::gpar(fontsize = fontsize,
                                                            fontface = "bold"),
                               border = border,
                               rect_gp = grid::gpar(col= grid_color),
                               cell_fun = function(j, i, x, y, w, h, fill) {
                                 grid::grid.rect(x, y, w, h, gp = grid::gpar(alpha = 0))
                               },
                               column_names_centered = F,
                               row_names_centered = F)

  return_list <- list("heatmap" = h,
                      "legend" = lgd)

  return(return_list)
}


#' Modify a string to wrap it around the middle point.
#'
#' @param string_to_modify
#'
#' @return None
#' @noRd
#' @examples
#' \dontrun{
#' TBD
#' }
modify_string <- function(string_to_modify){
  words <- stringr::str_split(string_to_modify, " ")[[1]]
  num_words <- length(words)
  middle_point <- round(num_words / 2, 0)
  string_to_modify <- paste(paste(words[1:middle_point], collapse = " "), "\n",
                            paste(words[(middle_point + 1):num_words], collapse = " "))
  return(string_to_modify)
}


#' Compute Enrichment scores using Seurat::AddModuleScore()
#'
#' @param sample  Seurat object.
#' @param list_genes  Named list of genes to compute enrichment for.
#' @param verbose  Verbose output.
#' @param nbin Number of bins.
#' @param ctrl Number of control genes.
#'
#' @return None
#' @noRd
#' @examples
#' \dontrun{
#' TBD
#' }
compute_enrichment_scores <- function(sample, list_genes, verbose = F, nbin = 24, ctrl = 100){
  if (!is.list(list_genes) & is.character(list_genes)){
    list_genes <- list("Input" = list_genes)
  }
  for (celltype in names(list_genes)){
    list_markers <- list(list_genes[[celltype]])

    # Compute Seurat AddModuleScore as well.
    if (verbose){
      sample <- Seurat::AddModuleScore(sample,
                                       list_markers,
                                       name = celltype,
                                       search = TRUE,
                                       verbose = T,
                                       nbin = nbin,
                                       ctrl = ctrl)
    } else {
      sample <- suppressMessages(suppressWarnings(Seurat::AddModuleScore(sample,
                                                                         list_markers,
                                                                         name = celltype,
                                                                         search = TRUE,
                                                                         verbose = F,
                                                                         nbin = nbin,
                                                                         ctrl = ctrl)))
    }


    # Retrieve the scores.
    col_name <- stringr::str_replace_all(paste0(celltype, "1"), " ", ".")
    col_name <- stringr::str_replace_all(col_name, "-", ".")
    col_name <- stringr::str_replace_all(col_name, "\\+", ".")

    # Modify the name that Seurat::AddModuleScore gives by default.
    sample@meta.data[, celltype] <- sample@meta.data[, col_name]
    # Remove old metadata.
    sample@meta.data[, col_name] <- NULL
  }
  return(sample)
}


#' Modify the aspect of the legend.
#'
#' @param p Plot.
#' @param legend.aes Character. Either color or fill.
#' @param legend.type Character. Type of legend to display. One of: normal, colorbar, colorsteps.
#' @param legend.position Position of the legend in the plot. Will only work if legend is set to TRUE.
#' @param legend.framewidth,legend.tickwidth Width of the lines of the box in the legend.
#' @param legend.framecolor,legend.tickcolor Color of the lines of the box in the legend.
#' @param legend.length,legend.width Length and width of the legend. Will adjust automatically depending on legend side.
#' @param legend.title Character. Title for the legend.
#'
#' @return None
#' @noRd
#' @examples
#' \dontrun{
#' TBD
#' }
modify_continuous_legend <- function(p,
                                     legend.aes,
                                     legend.type,
                                     legend.position,
                                     legend.length,
                                     legend.width,
                                     legend.framecolor,
                                     legend.tickcolor,
                                     legend.tickwidth,
                                     legend.framewidth,
                                     legend.title = NULL){
  # Define legend parameters. Width and height values will change depending on the legend orientation.
  if (legend.position %in% c("top", "bottom")){
    legend.barwidth <- legend.length
    legend.barheight <- legend.width
  } else if (legend.position %in% c("left", "right")){
    legend.barwidth <- legend.width
    legend.barheight <- legend.length
  }

  legend.title <- if (is.null(legend.title)){ggplot2::waiver()} else {legend.title}

  if (legend.aes == "color" | legend.aes == "colour"){
    if (legend.type == "normal"){
      p <- p +
        ggplot2::guides(color = ggplot2::guide_colorbar(title = legend.title,
                                                        title.position = "top",
                                                        title.hjust = 0.5))
    } else if (legend.type == "colorbar"){
      p <- p +
        ggplot2::guides(color = ggplot2::guide_colorbar(title = legend.title,
                                                        title.position = "top",
                                                        barwidth = legend.barwidth,
                                                        barheight = legend.barheight,
                                                        title.hjust = 0.5,
                                                        ticks.linewidth = legend.tickwidth,
                                                        frame.linewidth = legend.framewidth,
                                                        frame.colour = legend.framecolor,
                                                        ticks.colour = legend.tickcolor))
    } else if (legend.type == "colorsteps"){
      p <- p +
        ggplot2::guides(color = ggplot2::guide_colorsteps(title = legend.title,
                                                          title.position = "top",
                                                          barwidth = legend.barwidth,
                                                          barheight = legend.barheight,
                                                          title.hjust = 0.5,
                                                          ticks.linewidth = legend.tickwidth,
                                                          frame.linewidth = legend.framewidth,
                                                          frame.colour = legend.framecolor,
                                                          ticks.colour = legend.tickcolor))
    }
  } else if (legend.aes == "fill"){
    if (legend.type == "normal"){
      p <- p +
        ggplot2::guides(fill = ggplot2::guide_colorbar(title = legend.title,
                                                       title.position = "top",
                                                        title.hjust = 0.5))
    } else if (legend.type == "colorbar"){
      p <- p +
        ggplot2::guides(fill = ggplot2::guide_colorbar(title = legend.title,
                                                       title.position = "top",
                                                        barwidth = legend.barwidth,
                                                        barheight = legend.barheight,
                                                        title.hjust = 0.5,
                                                        ticks.linewidth = legend.tickwidth,
                                                        frame.linewidth = legend.framewidth,
                                                        frame.colour = legend.framecolor,
                                                        ticks.colour = legend.tickcolor))
    } else if (legend.type == "colorsteps"){
      p <- p +
        ggplot2::guides(fill = ggplot2::guide_colorsteps(title = legend.title,
                                                         title.position = "top",
                                                          barwidth = legend.barwidth,
                                                          barheight = legend.barheight,
                                                          title.hjust = 0.5,
                                                          ticks.linewidth = legend.tickwidth,
                                                          frame.linewidth = legend.framewidth,
                                                          frame.colour = legend.framecolor,
                                                          ticks.colour = legend.tickcolor))
    }
  }

  return(p)
}
