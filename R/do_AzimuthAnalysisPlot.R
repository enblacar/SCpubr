#' Generate a report from a Seurat object that has been mapped to a reference using Azimuth.
#'
#' This functions takes as input a Seurat object that has undergone reference mapping using Azimuth and
#' generates a set of visualizations and a combined output from the results of such mapping. The user can also
#' provide the reference object used for the mapping (if available) to produce a more complete output.
#'
#' @inheritParams doc_function
#'
#' @param annotation.labels \strong{\code{\link[base]{character}}} | Metadata column that stores the inferred annotation from Azimuth.
#' @param annotation.scoring \strong{\code{\link[base]{character}}} | Metadata column that stores the annotation scoring from Azimuth.
#' @param mapping.scoring \strong{\code{\link[base]{character}}} | Metadata column that stores the mapping scoring from Azimuth.
#' @param annotation.cutoff \strong{\code{\link[base]{numeric}}} | Value from 0 to 1 to use as cutoff to assign the labels to the object. This is used in conjunction with mapping.cutoff.
#' @param mapping.cutoff \strong{\code{\link[base]{numeric}}} | Value from 0 to 1 to use as cutoff to assign the labels to the object. This is used in conjunction with annotation.cutoff.
#' @param ref.obj  \strong{\code{\link[SeuratObject]{Seurat}}} | Seurat object used for reference mapping. Providing this object will add an extra plot with the UMAP of the reference and add its silhouette to the UMAP in which the original cells are showed in the context of the UMAP embedding of the reference object.
#' @param ref.reduction \strong{\code{\link[base]{character}}} | Name of the reduction embedding used to plot the UMAP in the reference object.
#' @return A list containing multiple plots.
#' @export
#'
#' @example /man/examples/examples_do_AzimuthAnalysisPlot.R
do_AzimuthAnalysisPlot <- function(sample,
                                   annotation.labels,
                                   annotation.scoring,
                                   mapping.scoring = "mapping.score",
                                   annotation.cutoff = 0.75,
                                   mapping.cutoff = 0,
                                   group.by = NULL,
                                   ref.obj = NULL,
                                   ref.reduction = "ref.umap",
                                   raster = FALSE,
                                   pt.size = if (isTRUE(raster)) {4} else {1},
                                   raster.dpi = 2048,
                                   border.size = if (isTRUE(raster)) {1.25} else {1.5},
                                   border.color = "black",
                                   na.value = "grey75",
                                   font.size = 14,
                                   font.type = "sans",
                                   colors.use = NULL,
                                   label = TRUE,
                                   legend.position = "bottom",
                                   viridis_color_map = "G",
                                   viridis_direction = 1){



  check_suggests(function_name = "do_AzimuthAnalysisPlot")
  # Check if the sample provided is a Seurat object.
  check_Seurat(sample = sample)
  if (!(is.null(ref.obj))) {check_Seurat(sample = sample)}

  # Check logical parameters.
  logical_list <- list("raster" = raster)
  check_type(parameters = logical_list, required_type = "logical", test_function = is.logical)
  # Check numeric parameters.
  numeric_list <- list("annotation.cutoff" = annotation.cutoff,
                       "pt.size" = pt.size,
                       "raster.dpi" = raster.dpi,
                       "border.size" = border.size,
                       "font.size" = font.size,
                       "viridis_direction" = viridis_direction)
  check_type(parameters = numeric_list, required_type = "numeric", test_function = is.numeric)
  # Check character parameters.
  character_list <- list("annotation.labels" = annotation.labels,
                         "group.by" = group.by,
                         "annotation.scoring" = annotation.scoring,
                         "mapping.scoring" = mapping.scoring,
                         "colors.use" = colors.use,
                         "font.type" = font.type,
                         "border.color" = border.color,
                         "na.value" = na.value,
                         "ref.reduction" = ref.reduction,
                         "viridis_color_map" = viridis_color_map)
  check_type(parameters = character_list, required_type = "character", test_function = is.character)

  check_colors(na.value, parameter_name = "na.value")
  check_colors(border.color, parameter_name = "border.color")

  check_parameters(parameter = font.type, parameter_name = "font.type")

  output_list <- list()


  `%>%` <- magrittr::`%>%`


  # Assign group.by to a metadata variable.
  if (is.null(group.by)){
    sample@meta.data[, "Groups"] <- sample@active.ident
    group.by <- "Groups"
  }

  if (!is.null(colors.use)){
    check_colors(colors.use)
    check_consistency_colors_and_names(sample = sample,
                                       colors = colors.use,
                                       grouping_variable = group.by)
  }

  # Prepare sample metadata and apply desired cutoffs.
  sample@meta.data <- sample@meta.data %>%
                      tibble::rownames_to_column(var = "cell") %>%
                      dplyr::left_join(by = "cell",
                                       y = {sample@meta.data %>%
                                            tibble::rownames_to_column(var = "cell") %>%
                                            dplyr::select(dplyr::all_of(c("cell", annotation.labels, annotation.scoring, mapping.scoring))) %>%
                                            dplyr::rowwise() %>%
                                            dplyr::mutate("inferred_annotation" = if (.data[[annotation.scoring]] >= annotation.cutoff &
                                                                                      .data[[mapping.scoring]] >= mapping.cutoff){.data[[annotation.labels]]} else {NA}) %>%
                                            dplyr::select(c("cell", "inferred_annotation"))}) %>%
                      tibble::column_to_rownames(var = "cell")

  # DimPlot with the original cells mapped to the reference dataset.
  p.umap.ref <- do_DimPlot(sample = sample,
                           group.by = group.by,
                           reduction = "ref.umap",
                           label = label,
                           legend.position = if (isTRUE(label)){"none"} else {legend.position},
                           na.value = na.value,
                           raster = raster,
                           raster.dpi = raster.dpi,
                           plot.title = "Cells in reference UMAP embedding",
                           pt.size = pt.size,
                           font.size = font.size,
                           font.type = font.type,
                           colors.use = colors.use)

  # If the user has provided a UMAP of the reference, add its silhouette.
  if (!is.null(ref.obj)){
    suppressMessages({
      p.ref <- do_DimPlot(sample = ref.obj,
                          label = label,
                          reduction = ref.reduction,
                          legend.position = if (isTRUE(label)){"none"} else {legend.position},
                          plot.title = "Reference UMAP",
                          na.value = na.value,
                          raster = raster,
                          raster.dpi = raster.dpi,
                          pt.size = pt.size,
                          font.size = font.size,
                          font.type = font.type)
    })


    data <- ggplot2::ggplot_build(p.ref)


    if (isFALSE(raster)){
      base_layer <- ggplot2::geom_point(data = data$data[[1]], mapping = ggplot2::aes(x = .data$x,
                                                                                      y = .data$y),
                                        colour = na.value,
                                        size = pt.size,
                                        show.legend = FALSE)
      p.umap.ref$layers <- append(base_layer, p.umap.ref$layers)

      base_layer <- ggplot2::geom_point(data = data$data[[1]], mapping = ggplot2::aes(x = .data$x,
                                                                                      y = .data$y),
                                        colour = border.color,
                                        size = pt.size * border.size,
                                        show.legend = FALSE)
      p.umap.ref$layers <- append(base_layer, p.umap.ref$layers)
    } else if (isTRUE(raster)){
      base_layer <- scattermore::geom_scattermore(data = data$data[[1]],
                                                  mapping = ggplot2::aes(x = .data$x,
                                                                         y = .data$y),
                                                  color = na.value,
                                                  size = pt.size,
                                                  stroke = pt.size / 2,
                                                  show.legend = FALSE,
                                                  pointsize = pt.size,
                                                  pixels = c(raster.dpi, raster.dpi))
      p.umap.ref$layers <- append(base_layer, p.umap.ref$layers)

      base_layer <- scattermore::geom_scattermore(data = data$data[[1]],
                                                  mapping = ggplot2::aes(x = .data$x,
                                                                         y = .data$y),
                                                  color = border.color,
                                                  size = pt.size * border.size,
                                                  stroke = pt.size / 2,
                                                  show.legend = FALSE,
                                                  pointsize = pt.size * border.size,
                                                  pixels = c(raster.dpi, raster.dpi))
      p.umap.ref$layers <- append(base_layer, p.umap.ref$layers)
    }
  }


  # DimPlot with the  inferred annotation.
  p.umap <- do_DimPlot(sample = sample,
                       group.by = "inferred_annotation",
                       label = label,
                       legend.position = if (isTRUE(label)){"none"} else {legend.position},
                       na.value = na.value,
                       raster = raster,
                       raster.dpi = raster.dpi,
                       plot.title = "Inferred annotation",
                       pt.size = pt.size,
                       font.size = font.size,
                       font.type = font.type)

  # DimPlot with the original annotation.
  p.umap.cluster <- do_DimPlot(sample = sample,
                               group.by = group.by,
                               label = label,
                               legend.position = if (isTRUE(label)){"none"} else {legend.position},
                               na.value = na.value,
                               raster = raster,
                               raster.dpi = raster.dpi,
                               plot.title = "Original annotation",
                               pt.size = pt.size,
                               font.size = font.size,
                               font.type = font.type,
                               colors.use = colors.use)

  # BarPlot with the proportion of inferred identities per original cluster.
  p.barplot <- do_BarPlot(sample = sample,
                          group.by = "inferred_annotation",
                          legend.position = legend.position,
                          split.by = group.by,
                          position = "fill",
                          plot.grid = FALSE,
                          plot.title = "Proportion of inferred identities per original cluster",
                          font.size = font.size,
                          font.type = font.type)

  # BarPlot with the proportion of individual datasets per original cluster.
  p.barplot2 <- do_BarPlot(sample = sample,
                           group.by = "orig.ident",
                           split.by = group.by,
                           position = "fill",
                           legend.position = legend.position,
                           plot.grid = FALSE,
                           plot.title = "Proportion of individual datasets per original cluster",
                           font.size = font.size,
                           font.type = font.type)

  # FeaturePlot with the prediction scores.
  p.prediction <- do_FeaturePlot(sample = sample,
                                 features = annotation.scoring,
                                 legend.title = "Azimuth annotation scores",
                                 na.value = na.value,
                                 raster = raster,
                                 legend.position = legend.position,
                                 raster.dpi = raster.dpi,
                                 plot.title = "Azimuth annotation scores",
                                 pt.size = pt.size,
                                 font.size = font.size,
                                 font.type = font.type,
                                 viridis_color_map = viridis_color_map,
                                 viridis_direction = viridis_direction)

  # FeaturePlot with the mapping scores.
  p.mapping <- do_FeaturePlot(sample = sample,
                              features = mapping.scoring,
                              legend.position = legend.position,
                              legend.title = "Azimuth mapping scores",
                              plot.title = "Azimuth mapping scores",
                              raster = raster,
                              raster.dpi = raster.dpi,
                              na.value = na.value,
                              pt.size = pt.size,
                              font.size = font.size,
                              font.type = font.type,
                              viridis_color_map = viridis_color_map,
                              viridis_direction = viridis_direction)


  # Generate a combined report.
  layout <- "AABBCC
             AABBCC
             DDEEFF
             DDEEFF
             GGGHHH
             GGGHHH"

  p.combined.portrait <- patchwork::wrap_plots(A = if (!is.null(ref.obj)) {p.ref} else {patchwork::plot_spacer()},
                                               B = p.umap.cluster,
                                               C = p.umap.ref,
                                               D = p.umap,
                                               E = p.prediction,
                                               F = p.mapping,
                                               G = p.barplot2,
                                               H = p.barplot,
                                               design = layout)

  layout <- "AABBCCGGG
             AABBCCGGG
             DDEEFFHHH
             DDEEFFHHH"

  p.combined.landscape <- patchwork::wrap_plots(A = if (!is.null(ref.obj)) {p.ref} else {patchwork::plot_spacer()},
                                                B = p.umap.cluster,
                                                C = p.umap.ref,
                                                D = p.umap,
                                                E = p.prediction,
                                                F = p.mapping,
                                                G = p.barplot2,
                                                H = p.barplot,
                                                design = layout)

  # Generate the output list.
  output_list[["mapping_scores"]] <- p.mapping
  output_list[["umap_in_reference"]] <- p.umap.ref
  output_list[["inferred_annotation"]] <- sample@meta.data[, "inferred_annotation", drop = F]
  output_list[["umap_prediction"]] <- p.umap
  output_list[["umap_clusters"]] <- p.umap.cluster
  output_list[["barplot_pred"]] <- p.barplot
  output_list[["barplot_orig"]] <- p.barplot2
  output_list[["annotation_scores"]] <- p.prediction
  output_list[["mapping_scores"]] <- p.mapping
  output_list[["report_portrait"]] <- p.combined.portrait
  output_list[["report_landscape"]] <- p.combined.landscape

  return(output_list)
}
