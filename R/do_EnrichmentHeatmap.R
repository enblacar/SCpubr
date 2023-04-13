#' Create enrichment scores heatmaps.
#'
#' This function computes the enrichment scores for the cells using \link[Seurat]{AddModuleScore} and then aggregates the scores by the metadata variables provided by the user and displays it as a heatmap, computed by \link[ComplexHeatmap]{Heatmap}.
#'
#' @inheritParams doc_function
#' @param enforce_symmetry \strong{\code{\link[base]{logical}}} | Whether the geyser and feature plot has a symmetrical color scale.
#' @param ncores \strong{\code{\link[base]{numeric}}} | Number of cores used to run UCell scoring.
#' @param storeRanks \strong{\code{\link[base]{logical}}} | Whether to store the ranks for faster UCell scoring computations. Might require large amounts of RAM.
#' @param return_object \strong{\code{\link[base]{logical}}} | Return the Seurat object with the enrichment scores stored.
#' @param geneset.order \strong{\code{\link[base]{character}}} | Should the gene sets be ordered in a specific way? Provide it as a vector of characters with the same names as the names of the gene sets.
#' @param group.order \strong{\code{\link[SCpubr]{named_list}}} | Should the groups in theheatmaps be ordered in a specific way? Provide it as a named list (as many lists as values in \strong{\code{group.by}}) with the order for each of the elements in the groups.
#' @return A ggplot2 object.
#' @export
#'
#' @example /man/examples/examples_do_EnrichmentHeatmap.R
do_EnrichmentHeatmap <- function(sample,
                                 input_gene_list,
                                 geneset.order = NULL,
                                 group.order = NULL,
                                 assay = NULL,
                                 slot = NULL,
                                 reduction = NULL,
                                 group.by = NULL,
                                 verbose = FALSE,
                                 na.value = "grey75",
                                 legend.position = "bottom",
                                 use_viridis = FALSE,
                                 viridis.palette = "G",
                                 viridis.direction = 1,
                                 legend.framewidth = 0.5,
                                 legend.tickwidth = 0.5,
                                 legend.length = 20,
                                 legend.width = 1,
                                 legend.framecolor = "grey50",
                                 legend.tickcolor = "white",
                                 legend.type = "colorbar",
                                 font.size = 14,
                                 font.type = "sans",
                                 axis.text.x.angle = 45,
                                 enforce_symmetry = FALSE,
                                 nbin = 24,
                                 ctrl = 100,
                                 flavor = "Seurat",
                                 legend.title = NULL,
                                 ncores = 1,
                                 storeRanks = TRUE,
                                 min.cutoff = NA,
                                 max.cutoff = NA,
                                 pt.size = 1,
                                 plot_cell_borders = TRUE,
                                 border.size = 2,
                                 return_object = FALSE,
                                 number.breaks = 5,
                                 sequential.palette = "YlGnBu",
                                 diverging.palette = "RdBu",
                                 sequential.direction = 1,
                                 flip = FALSE,
                                 grid.color = "white",
                                 border.color = "black",
                                 plot.title.face = "bold",
                                 plot.subtitle.face = "plain",
                                 plot.caption.face = "italic",
                                 axis.title.face = "bold",
                                 axis.text.face = "bold",
                                 legend.title.face = "bold",
                                 legend.text.face = "plain"){
  # Get defaults user warning length.
  length.use <- getOption("warning.length")
  
  # Restore the warning length on exit.
  on.exit(options(warning.length = length.use))
  
  # Set warning length to maximum.
  options(warning.length = 8170)
  
  check_suggests(function_name = "do_EnrichmentHeatmap")
  # Check if the sample provided is a Seurat object.
  check_Seurat(sample = sample)

  # Check logical parameters.
  logical_list <- list("use_viridis" = use_viridis,
                       "enforce_symmetry" = enforce_symmetry,
                       "plot_cell_borders" = plot_cell_borders,
                       "flip" = flip)
  check_type(parameters = logical_list, required_type = "logical", test_function = is.logical)
  # Check numeric parameters.
  numeric_list <- list("viridis.direction" = viridis.direction,
                       "nbin" = nbin,
                       "ctrl" = ctrl,
                       "ncores" = ncores,
                       "pt.size" = pt.size,
                       "border.size" = border.size,
                       "font.size" = font.size,
                       "legend.width" = legend.width,
                       "legend.length" = legend.length,
                       "legend.framewidth" = legend.framewidth,
                       "legend.tickwidth" = legend.tickwidth,
                       "viridis.direction" = viridis.direction,
                       "axis.text.x.angle" = axis.text.x.angle,
                       "min.cutoff" = min.cutoff,
                       "max.cutoff" = max.cutoff,
                       "number.breaks" = number.breaks,
                       "sequential.direction" = sequential.direction)
  check_type(parameters = numeric_list, required_type = "numeric", test_function = is.numeric)
  # Check character parameters.
  character_list <- list("input_gene_list" = input_gene_list,
                         "legend.title" = legend.title,
                         "legend.position" = legend.position,
                         "legend.framecolor" = legend.framecolor,
                         "font.type" = font.type,
                         "group.by" = group.by,
                         "na.value" = na.value,
                         "legend.position" = legend.position,
                         "viridis.palette" = viridis.palette,
                         "flavor" = flavor,
                         "sequential.palette" = sequential.palette,
                         "diverging.palette" = diverging.palette,
                         "grid.color" = grid.color,
                         "border.color" = border.color,
                         "plot.title.face" = plot.title.face,
                         "plot.subtitle.face" = plot.subtitle.face,
                         "plot.caption.face" = plot.caption.face,
                         "axis.title.face" = axis.title.face,
                         "axis.text.face" = axis.text.face,
                         "legend.title.face" = legend.title.face,
                         "legend.text.face" = legend.text.face)
  check_type(parameters = character_list, required_type = "character", test_function = is.character)

  check_colors(na.value, parameter_name = "na.value")
  check_colors(legend.framecolor, parameter_name = "legend.framecolor")
  check_colors(legend.tickcolor, parameter_name = "legend.tickcolor")
  check_colors(grid.color, parameter_name = "grid.color")
  check_colors(border.color, parameter_name = "border.color")


  check_parameters(parameter = viridis.direction, parameter_name = "viridis.direction")
  check_parameters(parameter = sequential.direction, parameter_name = "sequential.direction")
  check_parameters(parameter = diverging.palette, parameter_name = "diverging.palette")
  check_parameters(parameter = sequential.palette, parameter_name = "sequential.palette")
  check_parameters(parameter = font.type, parameter_name = "font.type")
  check_parameters(parameter = legend.type, parameter_name = "legend.type")
  check_parameters(parameter = legend.position, parameter_name = "legend.position")
  check_parameters(parameter = flavor, parameter_name = "flavor")
  check_parameters(parameter = number.breaks, parameter_name = "number.breaks")
  check_parameters(plot.title.face, parameter_name = "plot.title.face")
  check_parameters(plot.subtitle.face, parameter_name = "plot.subtitle.face")
  check_parameters(plot.caption.face, parameter_name = "plot.caption.face")
  check_parameters(axis.title.face, parameter_name = "axis.title.face")
  check_parameters(axis.text.face, parameter_name = "axis.text.face")
  check_parameters(legend.title.face, parameter_name = "legend.title.face")
  check_parameters(legend.text.face, parameter_name = "legend.text.face")
  
  
  `%>%` <- magrittr::`%>%`
  
  
  
  if (!(is.null(assay)) & flavor == "UCell"){
    warning(paste0(add_warning(), crayon_body("When using "),
                   crayon_key("flavor = UCell"),
                   crayon_body(" do not use the "),
                   crayon_key("assay"),
                   crayon_body(" parameter.\nInstead, make sure that the "),
                   crayon_key("assay"),
                   crayon_body(" you want to compute the scores with is set as the "),
                   crayon_key("default"),
                   crayon_body(" assay. Setting it to "),
                   crayon_key("NULL"),
                   crayon_body(".")), call. = FALSE)
  }
  
  if (!(is.null(slot)) & flavor == "Seurat"){
    warning(paste0(add_warning(), crayon_body("When using "),
                   crayon_key("flavor = Seurat"),
                   crayon_body(" do not use the "),
                   crayon_key("slot"),
                   crayon_body(" parameter.\nThis is determiend by default in "),
                   crayon_key("Seurat"),
                   crayon_body(". Setting it to "),
                   crayon_key("NULL"),
                   crayon_body(".")), call. = FALSE)
  }
  
  if (is.null(assay)){assay <- check_and_set_assay(sample)$assay}
  if (is.null(slot)){slot <- check_and_set_slot(slot)}
  
  if (is.character(input_gene_list)){
    stop(paste0(add_cross(),
                crayon_body("You have provided a string of genes to "),
                crayon_key("input_gene_list"),
                crayon_body(". Please provide a "),
                crayon_key("named list"),
                crayon_body(" instead.")), call. = FALSE)
    
  }
  
  
  if (!is.null(geneset.order)){
    assertthat::assert_that(sum(geneset.order %in% names(input_gene_list)) == length(names(input_gene_list)),
                            msg = paste0(add_cross(), crayon_body("The names provided to "),
                                         crayon_key("geneset.order"),
                                         crayon_body(" do not match the names of the gene sets in "),
                                         crayon_key("input_gene_list"),
                                         crayon_body(".")))
  }
  
  if (is.null(legend.title)){
    if (flavor == "AUCell") {
      legend.title <- "AUC"
    } else if (flavor == "UCell"){
      legend.title <- "UCell score"
    } else if (flavor == "Seurat"){
      legend.title <- "Enrichment"
    }
  }
  
  input_list <- input_gene_list
  assertthat::assert_that(!is.null(names(input_list)),
                          msg = paste0(add_cross(), crayon_body("Please provide a "),
                                       crayon_key("named list"),
                                       crayon_body(" to "),
                                       crayon_key("input_gene_list"),
                                       crayon_body(".")))
  
  if (length(unlist(stringr::str_match_all(names(input_list), "_"))) > 0){
    warning(paste0(add_warning(), crayon_body("Found "),
                   crayon_key("underscores (_)"),
                   crayon_body(" in the name of the gene sets provided. Replacing them with "),
                   crayon_key("dots (.)"),
                   crayon_body(" to avoid conflicts when generating the Seurat assay.")), call. = FALSE)
    names.use <- stringr::str_replace_all(names(input_list), "_", ".")
    names(input_list) <- names.use
    
    if (!is.null(geneset.order)){
      geneset.order <- stringr::str_replace_all(names(geneset.order), "_", ".")
    }
  }
  

  # Compute the enrichment scores.
  sample <- compute_enrichment_scores(sample = sample,
                                      input_gene_list = input_list,
                                      verbose = verbose,
                                      nbin = nbin,
                                      ctrl = ctrl,
                                      flavor = flavor,
                                      ncores = ncores,
                                      storeRanks = storeRanks,
                                      assay = if (flavor == "UCell"){NULL} else {assay},
                                      slot = if (flavor == "Seurat"){NULL} else {slot})

  out.list <- list()
  if (is.null(group.by)){
    assertthat::assert_that(!("Groups" %in% colnames(sample@meta.data)),
                            msg = paste0(add_cross(), crayon_body("Please, make sure you provide a value for "),
                                         crayon_key("group.by")))

    aggr_entities <- levels(sample)
    sample@meta.data[, "Groups"] <- sample@active.ident
    group.by <- "Groups"
  }

  for (g in group.by){
    assertthat::assert_that(g %in% colnames(sample@meta.data),
                            msg = paste0(add_cross(), crayon_body("The value "),
                                         crayon_key(g),
                                         crayon_body(" in "),
                                         crayon_key("group.by"),
                                         crayon_body(" is not present in the sample "),
                                         crayon_key(" metadata.")))

    assertthat::assert_that(class(sample@meta.data[, g]) %in% c("character", "factor"),
                            msg = paste0(add_cross(), crayon_body("The value "),
                                         crayon_key(g),
                                         crayon_body(" in "),
                                         crayon_key("group.by"),
                                         crayon_body(" is not a "),
                                         crayon_key(" character"),
                                         crayon_body( "or "),
                                         crayon_key("factor"),
                                         crayon_body(" column in the sample"),
                                         crayon_key("metadata"),
                                         crayon_body(".")))
  }



  matrix.list <- list()
  for (group in group.by){
    # Extract activities from object as a long dataframe
    suppressMessages({
      sample$group.by <- sample@meta.data[, group]

      df <- sample@meta.data %>%
            dplyr::select(dplyr::all_of(c("group.by", names(input_list)))) %>%
            tidyr::pivot_longer(cols = -c("group.by"),
                                names_to = "gene_list",
                                values_to = "enrichment") %>%
            dplyr::group_by(.data$group.by, .data$gene_list) %>%
            dplyr::summarise(mean = mean(.data$enrichment, na.rm = TRUE))

      df.order <- df
    })

    matrix.list[[group]][["df"]] <- df
    matrix.list[[group]][["df.order"]] <- df.order
  }




  counter <- 0
  for (group in group.by){
    counter <- counter + 1

    df <- matrix.list[[group]][["df"]]
    df.order <- matrix.list[[group]][["df.order"]]

    # Transform to wide to retrieve the hclust.
    df.order <- df.order %>%
                tidyr::pivot_wider(id_cols = "group.by",
                                   names_from = 'gene_list',
                                   values_from = 'mean') %>%
                tibble::column_to_rownames("group.by") %>%
                as.matrix()

    df.order[is.na(df.order)] <- 0
    if(length(rownames(df.order)) == 1){
      row_order <- rownames(df.order)[1]
    } else {
      row_order <- rownames(df.order)[stats::hclust(stats::dist(df.order, method = "euclidean"), method = "ward.D")$order]
    }
    if (counter == 1){
      if (length(colnames(df.order)) == 1){
        col_order <- colnames(df.order)[1]
      } else {
        col_order <- colnames(df.order)[stats::hclust(stats::dist(t(df.order), method = "euclidean"), method = "ward.D")$order]
      }
    }
    if (!is.null(group.order) & (group %in% names(group.order))){
      group.order.use <- group.order[[group]]
    } else {
      group.order.use <- row_order
    }
    data <- df %>%
            dplyr::mutate("gene_list" = factor(.data$gene_list, levels = if (is.null(geneset.order)){rev(col_order)} else {geneset.order}),
                          "group.by" = factor(.data$group.by, levels = group.order.use))

    if (!is.na(min.cutoff)){
      data <- data %>%
              dplyr::mutate("mean" = ifelse(.data$mean < min.cutoff, min.cutoff, .data$mean))
    }

    if (!is.na(max.cutoff)){
      data <- data %>%
              dplyr::mutate("mean" = ifelse(.data$mean > max.cutoff, max.cutoff, .data$mean))
    }
    matrix.list[[group]] <- NULL
    matrix.list[[group]][["data"]] <- data
  }



  # Compute limits.
  min.vector <- c()
  max.vector <- c()

  for (group in group.by){
    data <- matrix.list[[group]][["data"]]

    min.vector <- append(min.vector, min(data$mean, na.rm = TRUE))
    max.vector <- append(max.vector, max(data$mean, na.rm = TRUE))
  }

  # Get the absolute limits of the datasets.
  limits <- c(min(min.vector),
              max(max.vector))


  # Compute overarching scales for all heatmaps.
  scale.setup <- compute_scales(sample = sample,
                                feature = " ",
                                assay = "SCT",
                                reduction = NULL,
                                slot = "scale.data",
                                number.breaks = number.breaks,
                                min.cutoff = min.cutoff,
                                max.cutoff = max.cutoff,
                                flavor = "Seurat",
                                enforce_symmetry = enforce_symmetry,
                                from_data = TRUE,
                                limits.use = limits)


  # Plot individual heatmaps.
  counter <- 0
  list.heatmaps <- list()
  for (group in group.by){
    counter <- counter + 1
    data <- matrix.list[[group]][["data"]]

    p <- data %>% 
         # nocov start
         ggplot2::ggplot(mapping = ggplot2::aes(x = if(isFALSE(flip)){.data$gene_list} else {.data$group.by},
                                                y = if(isFALSE(flip)){.data$group.by} else {.data$gene_list},
                                                fill = .data$mean)) +
         # nocov end
         ggplot2::geom_tile(color = grid.color, linewidth = 0.5) +
         ggplot2::scale_y_discrete(expand = c(0, 0)) +
         ggplot2::scale_x_discrete(expand = c(0, 0),
                                   position = "top") +
         ggplot2::guides(y.sec = guide_axis_label_trans(~paste0(levels(.data$group.by))),
                         x.sec = guide_axis_label_trans(~paste0(levels(.data$gene)))) +
         ggplot2::coord_equal()


    if (isTRUE(enforce_symmetry)){
      p <- p +
        ggplot2::scale_fill_gradientn(colors = RColorBrewer::brewer.pal(n = 11, name = diverging.palette) %>% rev(),
                                      na.value = na.value,
                                      name = legend.title,
                                      breaks = scale.setup$breaks,
                                      labels = scale.setup$labels,
                                      limits = scale.setup$limits)
    } else {
      if (isTRUE(use_viridis)){
        p <- p +
          ggplot2::scale_fill_viridis_c(direction = viridis.direction,
                                        option = viridis.palette,
                                        na.value = na.value,
                                        breaks = scale.setup$breaks,
                                        labels = scale.setup$labels,
                                        limits = scale.setup$limits,
                                        name = legend.title)
      } else {
        p <- p +
          ggplot2::scale_fill_gradientn(colors = if(sequential.direction == 1){RColorBrewer::brewer.pal(n = 9, name = sequential.palette)[2:9]} else {rev(RColorBrewer::brewer.pal(n = 9, name = sequential.palette)[2:9])},
                                        na.value = na.value,
                                        name = legend.title,
                                        breaks = scale.setup$breaks,
                                        labels = scale.setup$labels,
                                        limits = scale.setup$limits)
      }
    }


    p <- modify_continuous_legend(p = p,
                                  legend.title = legend.title,
                                  legend.aes = "fill",
                                  legend.type = legend.type,
                                  legend.position = legend.position,
                                  legend.length = legend.length,
                                  legend.width = legend.width,
                                  legend.framecolor = legend.framecolor,
                                  legend.tickcolor = legend.tickcolor,
                                  legend.framewidth = legend.framewidth,
                                  legend.tickwidth = legend.tickwidth)
    # Set axis titles.
    if(isFALSE(flip)){
      if (counter == 1){
        if (length(group.by) > 1){
          xlab <- NULL
        } else {
          xlab <- "Gene set"
        }
        
        ylab <- group
      } else {
        if (length(group.by) > 1){
          if (counter == length(group.by)){
            xlab <- "Gene set"
          } else {
            xlab <- NULL
          }
        } else {
          xlab <- NULL
        }
        ylab <- group
      }
    } else {
      if (counter == 1){
        ylab <- "Gene set"
        xlab <- group
      } else {
        xlab <- group
        ylab <- NULL
      }
    }
    


    axis.parameters <- handle_axis(flip = flip,
                                   group.by = group.by,
                                   group = group,
                                   counter = counter,
                                   axis.text.x.angle = axis.text.x.angle,
                                   plot.title.face = plot.title.face,
                                   plot.subtitle.face = plot.subtitle.face,
                                   plot.caption.face = plot.caption.face,
                                   axis.title.face = axis.title.face,
                                   axis.text.face = axis.text.face,
                                   legend.title.face = legend.title.face,
                                   legend.text.face = legend.text.face)



    # Set theme
    p <- p +
         ggplot2::xlab(xlab) +
         ggplot2::ylab(ylab) +
         ggplot2::theme_minimal(base_size = font.size) +
         ggplot2::theme(axis.ticks.x.bottom = axis.parameters$axis.ticks.x.bottom,
                        axis.ticks.x.top = axis.parameters$axis.ticks.x.top,
                        axis.ticks.y.left = axis.parameters$axis.ticks.y.left,
                        axis.ticks.y.right = axis.parameters$axis.ticks.y.right,
                        axis.text.y.left = axis.parameters$axis.text.y.left,
                        axis.text.y.right = axis.parameters$axis.text.y.right,
                        axis.text.x.top = axis.parameters$axis.text.x.top,
                        axis.text.x.bottom = axis.parameters$axis.text.x.bottom,
                        axis.title.x.bottom = axis.parameters$axis.title.x.bottom,
                        axis.title.x.top = axis.parameters$axis.title.x.top,
                        axis.title.y.right = axis.parameters$axis.title.y.right,
                        axis.title.y.left = axis.parameters$axis.title.y.left,
                        axis.line = ggplot2::element_blank(),
                        plot.title = ggplot2::element_text(face = plot.title.face, hjust = 0),
                        plot.subtitle = ggplot2::element_text(face = plot.subtitle.face, hjust = 0),
                        plot.caption = ggplot2::element_text(face = plot.caption.face, hjust = 1),
                        plot.title.position = "plot",
                        panel.grid = ggplot2::element_blank(),
                        panel.grid.minor.y = ggplot2::element_line(color = "white", linewidth = 1),
                        text = ggplot2::element_text(family = font.type),
                        plot.caption.position = "plot",
                        legend.text = ggplot2::element_text(face = legend.text.face),
                        legend.title = ggplot2::element_text(face = legend.title.face),
                        legend.justification = "center",
                        plot.margin = ggplot2::margin(t = 0, r = 0, b = 0, l = 0),
                        panel.border = ggplot2::element_rect(fill = NA, color = border.color, linewidth = 1),
                        panel.grid.major = ggplot2::element_blank(),
                        legend.position = legend.position,
                        plot.background = ggplot2::element_rect(fill = "white", color = "white"),
                        panel.background = ggplot2::element_rect(fill = "white", color = "white"),
                        legend.background = ggplot2::element_rect(fill = "white", color = "white"))

    list.heatmaps[[group]] <- p
  }

  # Plot the combined plot
  input <- if(isFALSE(flip)){list.heatmaps[rev(group.by)]}else{list.heatmaps[group.by]}
  p <- patchwork::wrap_plots(input,
                             ncol = if (isFALSE(flip)){1} else {NULL},
                             nrow = if(isTRUE(flip)) {1} else {NULL},
                             guides = "collect")
  p <- p +
       patchwork::plot_annotation(theme = ggplot2::theme(legend.position = legend.position,
                                                         plot.title = ggplot2::element_text(size = font.size,
                                                                                            family = font.type,
                                                                                            color = "black",
                                                                                            face = plot.title.face,
                                                                                            hjust = 0),
                                                         plot.subtitle = ggplot2::element_text(size = font.size,
                                                                                               face = plot.subtitle.face,
                                                                                               family = font.type,
                                                                                               color = "black",
                                                                                               hjust = 0),
                                                         plot.caption = ggplot2::element_text(size = font.size,
                                                                                              face = plot.caption.face,
                                                                                              family = font.type,
                                                                                              color = "black",
                                                                                              hjust = 1),
                                                         plot.caption.position = "plot"))




  out.list[["Heatmap"]] <- p

  if (isTRUE(return_object)){
    # Generate a Seurat assay.
    sample[["Enrichment"]] <- sample@meta.data %>% 
                              dplyr::select(dplyr::all_of(names(input_list))) %>% 
                              t() %>% 
                              as.data.frame() %>% 
                              Seurat::CreateAssayObject(.)
    
    sample@meta.data <- sample@meta.data %>% 
                        dplyr::select(-dplyr::all_of(names(input_list)))
    
    sample@assays$Enrichment@key <- "Enrichment_"
    Seurat::DefaultAssay(sample) <- "Enrichment"
    
    out.list[["Object"]] <- sample
  }

  if (isFALSE(return_object)){
    return_me <- out.list[["Heatmap"]]
  } else {
    return_me <- out.list
  }

  return(return_me)
}
