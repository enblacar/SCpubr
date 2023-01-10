#' Compute an overview of the GO terms associated with the input list of genes.
#'
#' @inheritParams doc_function
#' @param levels.use \strong{\code{\link[base]{numeric}}} | Vector of numerics corresponding to the GO ontology levels to plot. If NULL will compute all recursively until there are no results.
#' @param reverse.levels \strong{\code{\link[base]{logical}}} | Whether to place the higher levels first when computing the joint heatmap.
#' @param colors.use \strong{\code{\link[base]{character}}} | Vector of 2 colors to use in the heatmap. The first will correspond to the empty values and the second one to the genes present in the terms.
#'
#' @return A list containing all the matrices for the respective GO levels and all the individual and combined heatmaps.
#' @export
#'
#' @example /man/examples/examples_do_GroupedGOTermPlot.R
do_GroupedGOTermPlot <- function(genes,
                                 org.db,
                                 levels.use = NULL,
                                 GO_ontology = "BP",
                                 min.overlap = NULL,
                                 flip = TRUE,
                                 legend.position = "right",
                                 heatmap_gap = 0.5,
                                 cluster_rows = TRUE,
                                 cluster_cols = TRUE,
                                 cell_size = 8,
                                 reverse.levels = TRUE,
                                 colors.use = c("grey90", "#29353d"),
                                 rotate_x_axis_labels = 45,
                                 font.size = 10,
                                 verbose = FALSE){
  `%>%` <- magrittr::`%>%`

  check_suggests(function_name = "do_GroupedGOTermPlot")

  # Check logical parameters.
  logical_list <- list("flip" = flip,
                       "cluster_cols" = cluster_cols,
                       "cluster_rows" = cluster_rows)
  check_type(parameters = logical_list, required_type = "logical", test_function = is.logical)
  # Check numeric parameters.
  numeric_list <- list("font.size" = font.size,
                       "rotate_x_axis_labels" = rotate_x_axis_labels,
                       "cell_size" = cell_size,
                       "heatmap_gap" = heatmap_gap,
                       "min.overlap" = min.overlap,
                       "levels.use" = levels.use)
  check_type(parameters = numeric_list, required_type = "numeric", test_function = is.numeric)
  # Check character parameters.
  character_list <- list("legend.position" = legend.position,
                         "GO_ontology" = GO_ontology,
                         "colors.use" = colors.use,
                         "genes" = genes)
  if (is.null(min.overlap)){
    if (length(genes) < 4){
      min.overlap <- 1
    } else {
      min.overlap <- 3
    }
  }

  assertthat::assert_that(min.overlap >= 1,
                          msg = "Please provide a positive value higher or equal to 1 to min.overlap.")

  assertthat::assert_that("OrgDb" %in% class(org.db),
                          msg = "Please provide a valid OrgDb object to org.db parameter.")

  check_colors(colors.use)

  check_parameters(parameter = legend.position, parameter_name = "legend.position")
  check_parameters(parameter = rotate_x_axis_labels, parameter_name = "rotate_x_axis_labels")
  check_parameters(parameter = GO_ontology, parameter_name = "GO_ontology")



  # Compute the GO terms associated with each gene.
  matrices <- do_GroupedGO_matrices(genes = genes,
                                    org.db = org.db,
                                    levels.use = levels.use,
                                    ontologies = GO_ontology,
                                    min.overlap = min.overlap,
                                    verbose = verbose)
  # Compute associated plots.
  plots <- do_GroupedGO_analysis_heatmaps(result = matrices,
                                          genes = genes,
                                          flip = flip,
                                          levels.use = levels.use,
                                          legend.position = legend.position,
                                          heatmap_gap = heatmap_gap,
                                          ontologies = GO_ontology,
                                          cluster_rows = cluster_rows,
                                          cluster_cols = cluster_cols,
                                          cell_size = cell_size,
                                          reverse.levels = reverse.levels,
                                          colors.use = colors.use,
                                          rotate_x_axis_labels = rotate_x_axis_labels,
                                          font.size = font.size,
                                          verbose = verbose)

  # Return output.
  return_object <- list("Matrices" = matrices,
                        "Plots" = plots)

  return(return_object)
}




