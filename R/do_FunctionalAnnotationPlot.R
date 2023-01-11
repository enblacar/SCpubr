#' Compute functional annotation plots using GO or KEGG ontologies
#'
#'
#' @inheritParams doc_function
#' @param organism  \strong{\code{\link[base]{character}}} | Supported KEGG organism.
#' @param database \strong{\code{\link[base]{character}}} | Database to run the analysis on. One of:
#' \itemize{
#'   \item \emph{\code{GO}}.
#'   \item \emph{\code{KEGG}}.
#' }
#' @param p.adjust.cutoff \strong{\code{\link[base]{numeric}}} | Significance cutoff used to filter non-significant terms.
#' @param pAdjustMethod \strong{\code{\link[base]{character}}} | Method to adjust for multiple testing.  One of:
#' \itemize{
#'   \item \emph{\code{holm}}.
#'   \item \emph{\code{hochberg}}.
#'   \item \emph{\code{hommel}}.
#'   \item \emph{\code{bonferroni}}.
#'   \item \emph{\code{BH}}.
#'   \item \emph{\code{BY}}.
#'   \item \emph{\code{fdr}}.
#'   \item \emph{\code{none}}.
#' }
#' @param minGSSize \strong{\code{\link[base]{numeric}}} | Minimal size of genes annotated by Ontology term for testing.
#' @param maxGSSize \strong{\code{\link[base]{numeric}}} | Maximal size of genes annotated for testing.
#' @param showCategory \strong{\code{\link[base]{numeric}}} | Number of enriched terms to display in the output tree plot.
#' @param nWords \strong{\code{\link[base]{numeric}}} | The number of words in the cluster tags in the tree plot.
#' @param nCluster \strong{\code{\link[base]{numeric}}} | The number of clusters to group the resulting terms in the tree plot. Suggested value is above 2, as two can lead to some errors.
#'
#' @return A list containing a heatmap of the presence/absence of the genes in the enriched term, as well as a bar plot, dot plot and tree plot of the enriched terms.
#' @export
#'
#' @example /man/examples/examples_do_FunctionalAnnotationPlot.R
do_FunctionalAnnotationPlot <- function(genes,
                                        org.db,
                                        organism = "hsa",
                                        database = "GO",
                                        GO_ontology = "BP",
                                        min.overlap = if(length(genes) <= 4){1} else {3},
                                        p.adjust.cutoff = 0.05,
                                        pAdjustMethod = "BH",
                                        minGSSize = 10,
                                        maxGSSize = 500,
                                        cluster_cols = TRUE,
                                        cluster_rows = TRUE,
                                        cell_size = 8,
                                        heatmap_gap = 0.5,
                                        font.size = 10,
                                        font.type = "sans",
                                        rotate_x_axis_labels = 45,
                                        xlab = NULL,
                                        ylab = NULL,
                                        plot.title = NULL,
                                        plot.subtitle = NULL,
                                        plot.caption = NULL,
                                        plot.grid = TRUE,
                                        grid.color = "grey75",
                                        grid.type = "dashed",
                                        flip = TRUE,
                                        legend.type = "colorbar",
                                        legend.position = "bottom",
                                        legend.framewidth = 0.5,
                                        legend.tickwidth = 0.5,
                                        legend.length = 20,
                                        legend.width = 1,
                                        legend.framecolor = "grey50",
                                        legend.tickcolor = "white",
                                        heatmap.legend.length = 75,
                                        heatmap.legend.width = 5,
                                        heatmap.legend.framecolor = "black",
                                        viridis_color_map = "G",
                                        viridis_direction = -1,
                                        showCategory = 30,
                                        nWords = 4,
                                        nCluster = 5){
  `%>%` <- magrittr::`%>%`

  check_suggests(function_name = "do_FunctionalAnnotationPlot")

  # Check logical parameters.
  logical_list <- list("flip" = flip,
                       "cluster_cols" = cluster_cols,
                       "cluster_rows" = cluster_rows,
                       "plot.grid" = plot.grid)
  check_type(parameters = logical_list, required_type = "logical", test_function = is.logical)
  # Check numeric parameters.
  numeric_list <- list("font.size" = font.size,
                       "legend.framewidth" = legend.framewidth,
                       "legend.tickwidth" = legend.tickwidth,
                       "legend.length" = legend.length,
                       "legend.width" = legend.width,
                       "viridis_direction" = viridis_direction,
                       "rotate_x_axis_labels" = rotate_x_axis_labels,
                       "heatmap.legend.length" = heatmap.legend.length,
                       "heatmap.legend.width" = heatmap.legend.width,
                       "showCategory" = showCategory,
                       "nWords" = nWords,
                       "nCluster" = nCluster,
                       "cell_size" = cell_size,
                       "heatmap_gap" = heatmap_gap,
                       "min.overlap" = min.overlap,
                       "minGSSize" = minGSSize,
                       "maxGSSize" = maxGSSize,
                       "p.adjust.cutoff" = p.adjust.cutoff)
  check_type(parameters = numeric_list, required_type = "numeric", test_function = is.numeric)
  # Check character parameters.
  character_list <- list("legend.position" = legend.position,
                         "plot.title" = plot.title,
                         "plot.subtitle" = plot.subtitle,
                         "plot.caption" = plot.caption,
                         "xlab" = xlab,
                         "ylab" = ylab,
                         "legend.framecolor" = legend.framecolor,
                         "legend.tickcolor" = legend.tickcolor,
                         "legend.type" = legend.type,
                         "font.type" = font.type,
                         "viridis_color_map" = viridis_color_map,
                         "grid.color" = grid.color,
                         "grid.type" = grid.type,
                         "organism" = organism,
                         "database" = database,
                         "GO_ontology" = GO_ontology,
                         "pAdjustMethod" = pAdjustMethod,
                         "genes" = genes)
  check_type(parameters = character_list, required_type = "character", test_function = is.character)


  assertthat::assert_that(min.overlap >= 1,
                          msg = "Please provide a positive value higher or equal to 1 to min.overlap.")

  assertthat::assert_that("OrgDb" %in% class(org.db),
                          msg = "Please provide a valid OrgDb object to org.db parameter.")

  check_parameters(parameter = font.type, parameter_name = "font.type")
  check_parameters(parameter = legend.type, parameter_name = "legend.type")
  check_parameters(parameter = legend.position, parameter_name = "legend.position")
  check_parameters(parameter = viridis_direction, parameter_name = "viridis_direction")
  check_parameters(parameter = viridis_color_map, parameter_name = "viridis_color_map")
  check_parameters(parameter = grid.type, parameter_name = "grid.type")
  check_parameters(parameter = rotate_x_axis_labels, parameter_name = "rotate_x_axis_labels")
  check_parameters(parameter = database, parameter_name = "database")
  check_parameters(parameter = GO_ontology, parameter_name = "GO_ontology")
  check_parameters(parameter = pAdjustMethod, parameter_name = "pAdjustMethod")

  check_colors(legend.framecolor, parameter_name = "legend.framecolor")
  check_colors(legend.tickcolor, parameter_name = "legend.tickcolor")
  check_colors(grid.color, parameter_name = "grid.color")

  # Convert genes to ENTREZIDs.
  suppressMessages({
    suppressWarnings({
      conversion <-clusterProfiler::bitr(genes, fromType = "SYMBOL",
                                         toType = c("ENTREZID"),
                                         OrgDb = org.db)
    })
  })


  if (database == "GO"){
    # Compute enriched GO terms.
    result <- clusterProfiler::enrichGO(gene = conversion$ENTREZID,
                                        OrgDb = org.db,
                                        ont = GO_ontology,
                                        pvalueCutoff = 0.05,
                                        qvalueCutoff = 0.05,
                                        pAdjustMethod = pAdjustMethod,
                                        minGSSize = minGSSize,
                                        maxGSSize = maxGSSize,
                                        readable = TRUE,
                                        pool = FALSE)

    if (is.null(result)){
      # nocov start
      return_obj <- "No gene could be mapped, try another database."
      return(return_obj)
      # nocov end
    }

  } else if (database == "KEGG"){
    suppressMessages({
      result <- clusterProfiler::enrichKEGG(gene = conversion$ENTREZID,
                                            organism = organism,
                                            pvalueCutoff = 0.05,
                                            qvalueCutoff = 0.05,
                                            pAdjustMethod = pAdjustMethod,
                                            minGSSize = minGSSize,
                                            maxGSSize = maxGSSize)
    })

    if (is.null(result)){
      # nocov start
      return_obj <- "No gene could be mapped, try another database."
      return(return_obj)
      # nocov end
    }

    geneID_column <- c()
    for (input in result@result$geneID){
      genes.use <- stringr::str_split(input, pattern = "/")[[1]]
      suppressMessages({
        conversion <- clusterProfiler::bitr(genes.use, fromType = "ENTREZID",
                                            toType = c("SYMBOL"),
                                            OrgDb = org.db)
      })
      geneID_column <- c(geneID_column, paste(conversion$SYMBOL, collapse = "/"))
    }
    result@result$geneID <- geneID_column
  }

  # Filter out non-significant results
  result@result <-  result@result %>%
                    dplyr::arrange(dplyr::desc(.data$Count)) %>%
                    dplyr::filter(.data$p.adjust <= p.adjust.cutoff,
                                  .data$Count >= min.overlap)

  if (nrow(result@result) == 0){
    output.list <- "With current presets of p.adjust.cutoff and min.overlap, no enriched terms surpass the cutoffs."
  } else {
    # Retrieve enriched GO Terms heatmap.
    h.enriched <- do_EnrichedTermMatrix(genes = genes,
                                        result = result,
                                        flip = flip,
                                        cluster_cols = cluster_cols,
                                        cluster_rows = cluster_rows,
                                        cell_size = cell_size,
                                        heatmap_gap = heatmap_gap,
                                        heatmap.legend.length = 75,
                                        heatmap.legend.width = 5,
                                        heatmap.legend.framecolor = "black",
                                        legend_gap = 1,
                                        font.size = font.size,
                                        legend.position = legend.position)

    p.barplot <- do_EnrichedTermBarPlot(result = result,
                                        font.size = font.size,
                                        font.type = font.type,
                                        rotate_x_axis_labels = rotate_x_axis_labels,
                                        xlab = xlab,
                                        ylab = ylab,
                                        plot.title = plot.title,
                                        plot.subtitle = plot.subtitle,
                                        plot.caption = plot.caption,
                                        plot.grid = plot.grid,
                                        grid.color = grid.color,
                                        grid.type = grid.type,
                                        flip = flip,
                                        legend.type = legend.type,
                                        legend.position = legend.position,
                                        legend.framewidth = legend.framewidth,
                                        legend.tickwidth = legend.tickwidth,
                                        legend.length = legend.length,
                                        legend.width = legend.width,
                                        legend.framecolor = legend.framecolor,
                                        legend.tickcolor = legend.tickcolor,
                                        viridis_color_map = viridis_color_map,
                                        viridis_direction = viridis_direction)

    p.dotplot <- do_EnrichedTermDotPlot(result = result,
                                        font.size = font.size,
                                        font.type = font.type,
                                        rotate_x_axis_labels = rotate_x_axis_labels,
                                        xlab = xlab,
                                        ylab = ylab,
                                        plot.title = plot.title,
                                        plot.subtitle = plot.subtitle,
                                        plot.caption = plot.caption,
                                        plot.grid = plot.grid,
                                        grid.color = grid.color,
                                        grid.type = grid.type,
                                        flip = flip,
                                        legend.type = legend.type,
                                        legend.position = legend.position,
                                        legend.framewidth = legend.framewidth,
                                        legend.tickwidth = legend.tickwidth,
                                        legend.length = legend.length,
                                        legend.width = legend.width,
                                        legend.framecolor = legend.framecolor,
                                        legend.tickcolor = legend.tickcolor,
                                        viridis_color_map = viridis_color_map,
                                        viridis_direction = viridis_direction)
    # nocov start
    if (nrow(result@result) < nCluster){
      nCluster <- nrow(result@result) - 1
    }
    # nocov end
    p.tree <- do_EnrichedTermTreePlot(result = result,
                                      legend.type = legend.type,
                                      legend.position = legend.position,
                                      legend.framewidth = legend.framewidth,
                                      legend.tickwidth = legend.tickwidth,
                                      legend.length = legend.length,
                                      legend.width = legend.width,
                                      legend.framecolor = legend.framecolor,
                                      legend.tickcolor = legend.tickcolor,
                                      viridis_color_map = viridis_color_map,
                                      viridis_direction = viridis_direction,
                                      font.size = font.size,
                                      font.type = font.type,
                                      plot.title = plot.title,
                                      plot.subtitle = plot.subtitle,
                                      plot.caption = plot.caption,
                                      showCategory = showCategory,
                                      nWords = nWords,
                                      nCluster = nCluster)

    output.list <- list("Heatmap" = h.enriched,
                        "BarPlot" = p.barplot,
                        "DotPlot" = p.dotplot,
                        "TreePlot" = p.tree)

  }


  return(output.list)
}
