#' BeeSwarm plot.
#'
#' @inheritParams doc_function
#' @param feature_to_rank \strong{\code{\link[base]{character}}} | Feature for which the cells are going to be ranked. Ideal case is that this feature is stored as a metadata column.
#' @param continuous_feature \strong{\code{\link[base]{logical}}} | Is the feature to rank and color for continuous? I.e: an enrichment score.
#' @param order \strong{\code{\link[base]{logical}}} | Whether to reorder the groups based on the median of the ranking.
#' @param remove_x_axis,remove_y_axis \strong{\code{\link[base]{logical}}} | Remove X axis labels and ticks from the plot.
#' @return  A ggplot2 object containing a Bee Swarm plot.
#' @export
#'
#' @example /man/examples/examples_do_BeeSwarmPlot.R
do_BeeSwarmPlot <- function(sample,
                            feature_to_rank,
                            group.by = NULL,
                            assay = NULL,
                            reduction = NULL,
                            slot = NULL,
                            continuous_feature = FALSE,
                            order = FALSE,
                            colors.use = NULL,
                            legend.title = NULL,
                            legend.type = "colorbar",
                            legend.position = if (isTRUE(continuous_feature)) {"bottom"} else {"none"},
                            legend.framewidth = 0.5,
                            legend.tickwidth = 0.5,
                            legend.length = 20,
                            legend.width = 1,
                            legend.framecolor = "grey50",
                            legend.tickcolor = "white",
                            legend.ncol = NULL,
                            legend.icon.size = 4,
                            plot.title = NULL,
                            plot.subtitle = NULL,
                            plot.caption = NULL,
                            xlab = NULL,
                            ylab = NULL,
                            font.size = 14,
                            font.type = "sans",
                            remove_x_axis = FALSE,
                            remove_y_axis = FALSE,
                            flip = FALSE,
                            use_viridis = TRUE,
                            viridis.palette = "G",
                            viridis.direction = 1,
                            sequential.palette = "YlGnBu",
                            sequential.direction = -1,
                            verbose = TRUE,
                            raster = FALSE,
                            raster.dpi = 300,
                            plot_cell_borders = TRUE,
                            border.size = 1.5,
                            border.color = "black",
                            pt.size = 2,
                            min.cutoff = NA,
                            max.cutoff = NA,
                            na.value = "grey75",
                            number.breaks = 5,
                            plot.title.face = "bold",
                            plot.subtitle.face = "plain",
                            plot.caption.face = "italic",
                            axis.title.face = "bold",
                            axis.text.face = "plain",
                            legend.title.face = "bold",
                            legend.text.face = "plain"){
  # Add lengthy error messages.
  withr::local_options(.new = list("warning.length" = 8170))
  
  check_suggests(function_name = "do_BeeSwarmPlot")
  `%>%` <- magrittr::`%>%`
  # Check ggbeeswarm version:
  # nocov start
  if(utils::packageVersion("ggbeeswarm") < "0.7.1"){
    warning(paste0(add_warning(), crayon_body("Due to recent updates in ggbeeswarm package, some internal interaction with ggplot2 have changed. Please update ggbeeswarm and ggplot2 to ensure correct plotting."), call. = FALSE))
  }
  # nocov end

  # Check if the sample provided is a Seurat object.
  check_Seurat(sample = sample)
  # Check the assay.
  out <- check_and_set_assay(sample, assay = assay)
  sample <- out[["sample"]]
  assay <- out[["assay"]]
  # Check the reduction.
  reduction <- check_and_set_reduction(sample = sample, reduction = reduction)
  # Check logical parameters.
  logical_list <- list("continuous_feature" = continuous_feature,
                       "remove_x_axis" = remove_x_axis,
                       "remove_y_axis" = remove_y_axis,
                       "flip" = flip,
                       "verbose" = verbose,
                       "raster" = raster,
                       "plot_cell_borders" = plot_cell_borders,
                       "use_viridis" = use_viridis,
                       "order" = order)
  check_type(parameters = logical_list, required_type = "logical", test_function = is.logical)
  # Check numeric parameters.
  numeric_list <- list("font.size" = font.size,
                       "raster.dpi" = raster.dpi,
                       "legend.framewidth" = legend.framewidth,
                       "legend.tickwidth" = legend.tickwidth,
                       "legend.length" = legend.length,
                       "legend.width" = legend.width,
                       "pt.size" = pt.size,
                       "border.size" = border.size,
                       "min.cutoff" = min.cutoff,
                       "max.cutoff" = max.cutoff,
                       "viridis.direction" = viridis.direction,
                       "legend.ncol" = legend.ncol,
                       "legend.icon.size" = legend.icon.size,
                       "number.breaks" = number.breaks,
                       "sequential.direction" = sequential.direction)
  check_type(parameters = numeric_list, required_type = "numeric", test_function = is.numeric)
  # Check character parameters.
  character_list <- list("legend.position" = legend.position,
                         "legend.title" = legend.title,
                         "plot.title" = plot.title,
                         "plot.subtitle" = plot.subtitle,
                         "plot.caption" = plot.caption,
                         "feature_to_rank" = feature_to_rank,
                         "group.by" = group.by,
                         "ylab" = ylab,
                         "xlab" = xlab,
                         "slot" = slot,
                         "viridis.palette" = viridis.palette,
                         "legend.framecolor" = legend.framecolor,
                         "legend.tickcolor" = legend.tickcolor,
                         "legend.type" = legend.type,
                         "font.type" = font.type,
                         "border.color" = border.color,
                         "sequential.palette" = sequential.palette,
                         "plot.title.face" = plot.title.face,
                         "plot.subtitle.face" = plot.subtitle.face,
                         "plot.caption.face" = plot.caption.face,
                         "axis.title.face" = axis.title.face,
                         "axis.text.face" = axis.text.face,
                         "legend.title.face" = legend.title.face,
                         "legend.text.face" = legend.text.face,
                         "na.value" = na.value)
  check_type(parameters = character_list, required_type = "character", test_function = is.character)
  # Check slot.
  slot <- check_and_set_slot(slot = slot)

  # Check the colors provided to legend.framecolor and legend.tickcolor and border color.
  check_colors(legend.framecolor, parameter_name = "legend.framecolor")
  check_colors(legend.tickcolor, parameter_name = "legend.tickcolor")
  check_colors(border.color, parameter_name = "border.color")
  check_colors(na.value, parameter_name = "na.value")

  check_parameters(parameter = font.type, parameter_name = "font.type")
  check_parameters(parameter = number.breaks, parameter_name = "number.breaks")
  check_parameters(parameter = legend.type, parameter_name = "legend.type")
  check_parameters(parameter = legend.position, parameter_name = "legend.position")
  check_parameters(parameter = viridis.palette, parameter_name = "viridis.palette")
  check_parameters(parameter = sequential.palette, parameter_name = "sequential.palette")
  check_parameters(plot.title.face, parameter_name = "plot.title.face")
  check_parameters(plot.subtitle.face, parameter_name = "plot.subtitle.face")
  check_parameters(plot.caption.face, parameter_name = "plot.caption.face")
  check_parameters(axis.title.face, parameter_name = "axis.title.face")
  check_parameters(axis.text.face, parameter_name = "axis.text.face")
  check_parameters(legend.title.face, parameter_name = "legend.title.face")
  check_parameters(legend.text.face, parameter_name = "legend.text.face")
  check_parameters(viridis.direction, parameter_name = "viridis.direction")
  check_parameters(sequential.direction, parameter_name = "sequential.direction")


  assertthat::assert_that(length(feature_to_rank) == 1,
                          msg = paste0(crayon_body("Please provide only "),
                                       crayon_key("one feature"),
                                       crayon_body(" to "),
                                       crayon_key("feature_to_rank"),
                                       crayon_body(".")))
  
  colors.gradient <- compute_continuous_palette(name = ifelse(isTRUE(use_viridis), viridis.palette, sequential.palette),
                                                use_viridis = use_viridis,
                                                direction = ifelse(isTRUE(use_viridis), viridis.direction, sequential.direction),
                                                enforce_symmetry = FALSE)
  
  # Check group.by.
  out <- check_group_by(sample = sample,
                        group.by = group.by,
                        is.heatmap = FALSE)
  sample <- out[["sample"]]
  group.by <- out[["group.by"]]
  
  dim_colnames <- check_feature(sample = sample, features = feature_to_rank, dump_reduction_names = TRUE)
  if (feature_to_rank %in% colnames(sample@meta.data)) {
    sample@meta.data$rank_me <- sample@meta.data[, feature_to_rank]
    sample@meta.data$rank <- rank(sample@meta.data$rank_me)
  } else if (feature_to_rank %in% rownames(sample)){
    sample@meta.data$rank_me  <- .GetAssayData(sample = sample, slot = slot, assay = assay)[feature_to_rank, ]
    sample@meta.data$rank <- rank(sample@meta.data$rank_me)
  } else if (feature_to_rank %in% dim_colnames){
    for(red in Seurat::Reductions(object = sample)){
      if (feature_to_rank %in% colnames(sample@reductions[[red]][[]])){
        reduction <- red
        sample@meta.data$rank_me  <- sample@reductions[[reduction]][[]][, feature_to_rank]
        sample@meta.data$rank <- rank(sample@meta.data$rank_me)
      }
    }
  }
  # Compute the ranking
  sample@meta.data$ranked_groups <- factor(sample@meta.data[, group.by], levels = sort(unique(sample@meta.data[, group.by])))
  
  if (isTRUE(order)){
    # Get median rank by group.
    order <- sample@meta.data %>% 
             dplyr::select(dplyr::all_of(c("ranked_groups", "rank"))) %>% 
             dplyr::group_by(.data$ranked_groups) %>% 
             dplyr::summarise("median" = stats::median(.data$rank, na.rm = TRUE)) %>% 
             dplyr::arrange(dplyr::desc(.data$median)) %>% 
             dplyr::pull(.data$ranked_groups) %>% 
             as.character()
    sample@meta.data$ranked_groups <- factor(sample@meta.data$ranked_groups, levels = rev(order))
  }
  
  color_by <- ifelse(continuous_feature == TRUE, "rank_me", "ranked_groups")
  

  # Compute the limits.
  if (isTRUE(continuous_feature)){
    data <- sample$rank_me
    range.data <- c(min(data, na.rm = TRUE),
                    max(data, na.rm = TRUE))

    scale.setup <- compute_scales(sample = NULL,
                                  feature = feature_to_rank,
                                  assay = NULL,
                                  reduction = NULL,
                                  slot = NULL,
                                  number.breaks = number.breaks,
                                  min.cutoff = min.cutoff,
                                  max.cutoff = max.cutoff,
                                  flavor = "Seurat",
                                  enforce_symmetry = FALSE,
                                  from_data = TRUE,
                                  limits.use = range.data)

    sample$rank_me[sample$rank_me < min.cutoff] <- min.cutoff
    sample$rank_me[sample$rank_me > max.cutoff] <- max.cutoff
  }
  
  p <- ggplot2::ggplot(sample@meta.data,
                       mapping = ggplot2::aes(x = .data[["rank"]],
                                              y = .data[["ranked_groups"]],
                                              color = !!rlang::sym(color_by)))

  # Add raster layer if desired.
  if (isTRUE(raster)){
    p <- p +
         ggrastr::geom_quasirandom_rast(raster.dpi = raster.dpi,
                                        size = pt.size)
  } else {
    p <- p +
         ggbeeswarm::geom_quasirandom(size = pt.size)
  }

  p <- p +
       ggplot2::labs(title = plot.title,
                     subtitle = plot.subtitle,
                     caption = plot.caption) +
       ggplot2::theme_minimal(base_size = font.size) +
       ggplot2::theme(plot.margin = ggplot2::margin(t = 10, r = 10, b = 10, l = 10),
                      plot.title = ggplot2::element_text(face = plot.title.face, hjust = 0),
                      plot.subtitle = ggplot2::element_text(face = plot.subtitle.face, hjust = 0),
                      plot.caption = ggplot2::element_text(face = plot.caption.face, hjust = 1),
                      panel.grid = ggplot2::element_blank(),
                      plot.title.position = "plot",
                      plot.caption.position = "plot",
                      text = ggplot2::element_text(family = font.type),
                      legend.text = ggplot2::element_text(face = legend.text.face),
                      legend.position = legend.position,
                      legend.title = ggplot2::element_text(face = legend.title.face),
                      legend.justification = "center",
                      axis.title.x = ggplot2::element_text(face = axis.title.face),
                      axis.title.y = ggplot2::element_text(face = axis.title.face, angle = 90),
                      axis.ticks.y = if(base::isFALSE(flip)){ggplot2::element_line(color = "black")} else {ggplot2::element_blank()},
                      axis.ticks.x = if(isTRUE(flip)){ggplot2::element_line(color = "black")} else {ggplot2::element_blank()},
                      axis.text = ggplot2::element_text(face = axis.text.face, color = "black"),
                      axis.line = ggplot2::element_line(color = "black"),
                      plot.background = ggplot2::element_rect(fill = "white", color = "white"),
                      panel.background = ggplot2::element_rect(fill = "white", color = "white"),
                      legend.background = ggplot2::element_rect(fill = "white", color = "white"))

  if (continuous_feature == TRUE){

    p <- p + 
         ggplot2::scale_color_gradientn(colors = colors.gradient,
                                        na.value = na.value,
                                        name = legend.title,
                                        breaks = scale.setup$breaks,
                                        labels = scale.setup$labels,
                                        limits = scale.setup$limits)

    p <- modify_continuous_legend(p = p,
                                  legend.title = legend.title,
                                  legend.aes = "color",
                                  legend.type = legend.type,
                                  legend.position = legend.position,
                                  legend.length = legend.length,
                                  legend.width = legend.width,
                                  legend.framecolor = legend.framecolor,
                                  legend.tickcolor = legend.tickcolor,
                                  legend.framewidth = legend.framewidth,
                                  legend.tickwidth = legend.tickwidth)
  } else if (continuous_feature == FALSE) {
    if (is.null(colors.use)){
      colors.use <- generate_color_scale(levels(sample))
    } else {
      colors.use <- check_consistency_colors_and_names(sample = sample, colors = colors.use, grouping_variable = group.by)
    }
    # Adapt the legend to categorical variables.
    if (is.null(legend.title)){
      legend.title <- "Groups"
    }
    p <- p +
         ggplot2::scale_color_manual(values = colors.use) +
         ggplot2::guides(color = ggplot2::guide_legend(title = legend.title,
                                                       title.position = "top",
                                                       title.hjust = 0.5,
                                                       ncol = legend.ncol,
                                                       override.aes = list(size = legend.icon.size))) +
         ggplot2::theme(legend.position = legend.position)
  }

  if (remove_x_axis == TRUE){
    p <- p +
         ggplot2::theme(axis.text.x = ggplot2::element_blank(),
                        axis.ticks.x = ggplot2::element_blank())
  }
  if (remove_y_axis == TRUE){
    p <- p +
         ggplot2::theme(axis.text.y = ggplot2::element_blank(),
                        axis.ticks.y = ggplot2::element_blank())
  }
  if (flip == TRUE){
    p <- p +
         ggplot2::coord_flip() +
         ggplot2::theme(axis.text.y = ggplot2::element_blank(),
                        axis.ticks.y = ggplot2::element_blank()) +
         ggplot2::xlab(ifelse(is.null(ylab), paste0("Ranking of ", feature_to_rank), ylab)) +
         ggplot2::ylab(if(is.null(xlab)) {group.by} else {xlab})

  } else {
    p <- p +
         ggplot2::theme(axis.text.x = ggplot2::element_blank(),
                        axis.ticks.x = ggplot2::element_blank()) +
         ggplot2::xlab(ifelse(is.null(xlab), paste0("Ranking of ", feature_to_rank), xlab)) +
         ggplot2::ylab(if(is.null(ylab)) {group.by} else {ylab})

  }

  if (isTRUE(plot_cell_borders)){
    # Generate base layer.
    if (isTRUE(raster)){
      base_layer <- ggrastr::geom_quasirandom_rast(data = sample@meta.data,
                                                   mapping = ggplot2::aes(x = .data[["rank"]],
                                                                          y = .data[["ranked_groups"]]),
                                                   raster.dpi = raster.dpi,
                                                   color = border.color,
                                                   size = pt.size * border.size,
                                                   show.legend = FALSE)
    } else if (base::isFALSE(raster)){
      base_layer <-ggbeeswarm::geom_quasirandom(data = sample@meta.data,
                                                 mapping = ggplot2::aes(x = .data[["rank"]],
                                                                        y = .data[["ranked_groups"]]),
                                                 color = border.color,
                                                 size = pt.size * border.size,
                                                 show.legend = FALSE)
    }
    p[["layers"]] <- append(base_layer, p[["layers"]])

  }
  

  return(p)

}
