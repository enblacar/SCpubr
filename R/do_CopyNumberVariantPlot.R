#' Display CNV scores from inferCNV as Feature Plots.
#'
#'
#' @inheritParams doc_function
#' @param infercnv_object \strong{\code{\link[infercnv]{infercnv}}} | Output inferCNV object run on the same Seurat object.
#' @param using_metacells \strong{\code{\link[base]{logical}}} | Whether inferCNV was run using metacells or not.
#' @param metacell_mapping \strong{\code{\link[SCpubr]{named_vector}}} | Vector or cell - metacell mapping.
#' @param chromosome_focus \strong{\code{\link[base]{character}}} | Region stating which chromosome to plot. Eg: 1p, 19q. NULL will plot all regions.
#' @param chromosome_locations \strong{\code{\link[tibble]{tibble}}} | Tibble containing the chromosome regions to use. Can be obtained using \code{utils::data("human_chr_locations", package = "SCpubr")}.
#'
#' @return A list containing Feature Plots for different chromosome regions and corresponding dot plots by groups..
#' @export
#'
#' @example man/examples/examples_do_CopyNumberVariantPlot.R
do_CopyNumberVariantPlot <- function(sample,
                                     infercnv_object,
                                     chromosome_locations,
                                     group.by = NULL,
                                     using_metacells = FALSE,
                                     metacell_mapping = NULL,
                                     chromosome_focus = NULL,
                                     legend.type = "colorbar",
                                     legend.position = "bottom",
                                     legend.length = 20,
                                     legend.width = 1,
                                     legend.framewidth = 0.5,
                                     legend.tickwidth = 0.5,
                                     legend.framecolor = "grey50",
                                     legend.tickcolor = "white",
                                     font.size = 14,
                                     pt.size = 1,
                                     font.type = "sans",
                                     border.size = 2,
                                     border.color = "black",
                                     rotate_x_axis_labels = 45,
                                     plot_cell_borders = TRUE,
                                     enforce_symmetry = TRUE,
                                     legend.title = NULL,
                                     na.value = "grey75",
                                     viridis.palette = "G",
                                     viridis.direction = 1,
                                     verbose = FALSE,
                                     min.cutoff = NA,
                                     max.cutoff = NA,
                                     number.breaks = 5,
                                     diverging.palette = "RdBu",
                                     sequential.palette = "YlGnBu",
                                     sequential.direction = -1,
                                     use_viridis = TRUE){


  # Check logical parameters.
  logical_list <- list("using_metacells" = using_metacells,
                       "enforce_symmetry" = enforce_symmetry,
                       "plot_cell_borders" = plot_cell_borders,
                       "use_viridis" = use_viridis)
  check_type(parameters = logical_list, required_type = "logical", test_function = is.logical)
  # Check numeric parameters.
  numeric_list <- list("font.size" = font.size,
                       "legend.length" = legend.length,
                       "legend.width" = legend.width,
                       "legend.framewidth" = legend.framewidth,
                       "legend.tickwidth" = legend.tickwidth,
                       "pt.size" = pt.size,
                       "border.size" = border.size,
                       "viridis.direction" = viridis.direction,
                       "rotate_x_axis_labels" = rotate_x_axis_labels,
                       "min.cutoff" = min.cutoff,
                       "max.cutoff" = max.cutoff,
                       "number.breaks" = number.breaks,
                       "sequential.direction" = sequential.direction)
  check_type(parameters = numeric_list, required_type = "numeric", test_function = is.numeric)
  # Check character parameters.
  character_list <- list("group.by" = group.by,
                         "legend.type" = legend.type,
                         "legend.position" = legend.position,
                         "legend.framecolor" = legend.framecolor,
                         "legend.tickcolor" = legend.tickcolor,
                         "font.type" = font.type,
                         "border.color" = border.color,
                         "legend.title" = legend.title,
                         "viridis.palette" = viridis.palette,
                         "diverging.palette" = diverging.palette,
                         "sequential.palette" = sequential.palette)
  check_type(parameters = character_list, required_type = "character", test_function = is.character)


  `:=` <- rlang::`:=`
  `%>%` <- magrittr::`%>%`
  check_colors(legend.framecolor, parameter_name = "legend.framecolor")
  check_colors(legend.tickcolor, parameter_name = "legend.tickcolor")
  check_colors(na.value, parameter_name = "na.value")

  check_parameters(parameter = font.type, parameter_name = "font.type")
  check_parameters(parameter = legend.type, parameter_name = "legend.type")
  check_parameters(parameter = legend.position, parameter_name = "legend.position")
  check_parameters(parameter = viridis.palette, parameter_name = "viridis.palette")
  check_parameters(parameter = viridis.direction, parameter_name = "viridis.direction")
  check_parameters(parameter = rotate_x_axis_labels, parameter_name = "rotate_x_axis_labels")
  check_parameters(parameter = number.breaks, parameter_name = "number.breaks")
  check_parameters(parameter = diverging.palette, parameter_name = "diverging.palette")
  check_parameters(parameter = sequential.palette, parameter_name = "sequential.palette")
  check_parameters(parameter = sequential.direction, parameter_name = "sequential.direction")

  if (is.null(chromosome_focus)){
    chromosome_list <- c(as.character(seq(1, 22)))
  } else {
    chromosome_list <- chromosome_focus
  }

  # Fix for group.by.
  if (is.null(group.by)){
    aggr_entities <- levels(sample)
    sample@meta.data[, "Groups"] <- sample@active.ident
    group.by <- "Groups"
  }


  # Retrieve the genes.
  genes <- infercnv_object@gene_order

  # Retrieve chr 1p start and end coordinates.
  chr_locations <- chromosome_locations

  # This list will contain all the outputs.
  return_list <- list()
  events_list <- c()
  for (chromosome in chromosome_list){
    # Retrieve chr locations of the chromosome.
    locations <- chr_locations %>%
                 dplyr::filter(.data[["chr"]] == chromosome)

    for (chr_arm in c("p", "q", "whole")){
      if (chr_arm != "whole"){
        # Retrieve the start.
        start <- locations %>%
                 dplyr::filter(.data[["arm"]] == chr_arm) %>%
                 dplyr::pull(start)
        # Retrieve the end.
        end <- locations %>%
               dplyr::filter(.data[["arm"]] == chr_arm) %>%
               dplyr::pull(end)

        # Retrieve the genes present in the chromosome arm.
        genes_use <- rownames(genes %>%
                              dplyr::filter(.data[["chr"]] == paste0("chr", chromosome),
                                            stop <= end))
      } else {
        # Retrieve the start.
        start <- locations %>%
                 dplyr::filter(.data[["arm"]] == "p") %>%
                 dplyr::pull(start)
        # Retrieve the end.
        end <- locations %>%
               dplyr::filter(.data[["arm"]] == "q") %>%
               dplyr::pull(end)

        # Retrieve the genes present in the chromosome arm.
        genes_use <- rownames(genes %>%
                              dplyr::filter(.data[["chr"]] == paste0("chr", chromosome)))
      }

      # Retrieve the CNV scores from the inferCNV object.
      CNV_scores <- infercnv_object@expr.data

      # Make it at least 2 genes in the object (this will only be applicable in VERY LOW QUALITY DATASETS)
      if (sum(genes_use %in% rownames(CNV_scores)) > 1){
        # Filter the scores for only the genes in the chromosome arm.
        CNV_scores <- CNV_scores[genes_use[genes_use %in% rownames(CNV_scores)], ]

        scores_name <- if (chr_arm != "whole"){paste0(chromosome, chr_arm)} else {chromosome}
        CNV_scores_final <- tibble::tibble(!!scores_name := colMeans(CNV_scores),
                                           "cells" = colnames(CNV_scores))


        # If metacells were used.
        if (isTRUE(using_metacells)){
          sample[["metacell_mapping"]] <- metacell_mapping
          sample@meta.data <- sample@meta.data %>%
                              dplyr::mutate("cells" = colnames(sample)) %>%
                              dplyr::left_join(y = {sample@meta.data %>%
                                                    dplyr::select(dplyr::all_of(c("metacell_mapping"))) %>%
                                                    tibble::rownames_to_column(var = "cells") %>%
                                                    dplyr::left_join(y = {CNV_scores_final %>% dplyr::rename("metacell_mapping" = dplyr::all_of(c("cells")))},
                                                                     by = "metacell_mapping") %>%
                                                    dplyr::select(-dplyr::all_of(c("metacell_mapping")))},
                                               by = "cells") %>%
                              tibble::column_to_rownames(var = "cells")
          # If no metacells were used.
        } else if (isFALSE(using_metacells)){
          sample@meta.data[, scores_name] <- CNV_scores_final[, scores_name]
        }
        events_list <- append(events_list, scores_name)
      } else {
        #nocov start
        if(isTRUE(verbose)){message(paste0("Your sample has only one gene in ", chromosome, chr_arm, ". Skipping this chromosome arm."))}
        #nocov end
      }
    }
  }

  sample@meta.data <- sample@meta.data %>%
                      dplyr::mutate("group" = .data[[group.by]],
                                    "cells" = colnames(sample))
  for (event in events_list){
    # Fix the legend.
    if (is.null(legend.title)){
      legend.title.use <- paste0(event, " scores")
    } else {
      legend.title.use <- legend.title
    }

    p <- do_GeyserPlot(sample = sample,
                       assay = NULL,
                       slot = NULL,
                       features = event,
                       group.by = group.by,
                       pt.size = pt.size,
                       border.size = border.size,
                       enforce_symmetry = enforce_symmetry,
                       order_by_mean = FALSE,
                       plot_cell_borders = plot_cell_borders,
                       font.size = font.size,
                       font.type = font.type,
                       legend.position = legend.position,
                       legend.type = legend.type,
                       legend.framecolor = legend.framecolor,
                       legend.tickcolor = legend.tickcolor,
                       legend.framewidth = legend.framewidth,
                       legend.tickwidth = legend.tickwidth,
                       legend.length = legend.length,
                       legend.width = legend.width,
                       xlab = group.by,
                       ylab = paste0(event, " score"),
                       legend.title = legend.title.use,
                       rotate_x_axis_labels = rotate_x_axis_labels,
                       viridis.palette = viridis.palette,
                       viridis.direction = viridis.direction,
                       min.cutoff = min.cutoff,
                       max.cutoff = max.cutoff,
                       number.breaks = number.breaks,
                       diverging.palette = diverging.palette,
                       use_viridis = use_viridis,
                       sequential.palette = sequential.palette,
                       sequential.direction = sequential.direction)


    p.f <- do_FeaturePlot(sample = sample,
                          features = event,
                          plot_cell_borders = plot_cell_borders,
                          pt.size = pt.size,
                          border.size = border.size,
                          font.size = font.size,
                          font.type = font.type,
                          legend.position = legend.position,
                          legend.type = legend.type,
                          legend.framecolor = legend.framecolor,
                          legend.tickcolor = legend.tickcolor,
                          legend.framewidth = legend.framewidth,
                          legend.tickwidth = legend.tickwidth,
                          legend.length = legend.length,
                          legend.width = legend.width,
                          viridis.palette = viridis.palette,
                          viridis.direction = viridis.direction,
                          legend.title = legend.title.use,
                          min.cutoff = min.cutoff,
                          max.cutoff = max.cutoff,
                          number.breaks = number.breaks,
                          diverging.palette = diverging.palette,
                          use_viridis = use_viridis,
                          sequential.palette = sequential.palette,
                          sequential.direction = sequential.direction)

    if (isTRUE(enforce_symmetry)){
      limits <- max(abs(c(min(sample@meta.data[, event]),
                        max(sample@meta.data[, event]))))

      scale_limit <- c((1 - (limits - 1)), limits)

      scale.use <- ggplot2::scale_color_gradientn(colors = RColorBrewer::brewer.pal(n = 11, name = diverging.palette) %>% rev(),
                                                  limits = c(-limits, limits),
                                                  na.value = na.value)
      suppressMessages({
        p <- p +
             ggplot2::ylim(scale_limit) +
             ggplot2::geom_hline(mapping = ggplot2::aes(yintercept = 1),
                                 color = "grey50",
                                 linetype = "dashed",
                                 linewidth = 1) +
             ggplot2::scale_color_gradientn(colors = RColorBrewer::brewer.pal(n = 11, name = diverging.palette) %>% rev(),
                                            limits = scale_limit,
                                            na.value = na.value)
        p.f <- p.f +
               ggplot2::scale_color_gradientn(colors = RColorBrewer::brewer.pal(n = 11, name = diverging.palette) %>% rev(),
                                              limits = scale_limit,
                                              na.value = na.value)
      })
    }



    return_list[[paste0(event, "_umap")]] <- p.f
    return_list[[paste0(event, "_geyser")]] <- p
  }
  return(return_list)
}
