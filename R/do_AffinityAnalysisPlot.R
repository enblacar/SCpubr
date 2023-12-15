#' Compute affinity of gene sets to cell populations using decoupleR.
#'
#' Major contributions to this function:
#' - \href{https://github.com/MarcElosua}{Marc Elosua Bayés} for the core concept code and idea.
#' - \href{https://github.com/paubadiam}{Pau Badia i Mompel} for the network generation.
#'
#' @inheritParams doc_function
#' @param statistic \strong{\code{\link[base]{character}}} | DecoupleR statistic to use for the analysis.
#' values in the Idents of the Seurat object are reported, assessing how specific a given gene set is for a given cell population compared to other gene sets of equal expression.
#'
#' @return A list containing different plots.
#' @export
#'
#' @example /man/examples/examples_do_AffinityAnalysisPlot.R

do_AffinityAnalysisPlot <- function(sample,
                                    input_gene_list,
                                    subsample = 2500,
                                    group.by = NULL,
                                    assay = NULL,
                                    slot = NULL,
                                    statistic = "ulm",
                                    number.breaks = 5,
                                    use_viridis = FALSE,
                                    viridis.palette = "G",
                                    viridis.direction = -1,
                                    sequential.palette = "YlGnBu",
                                    sequential.direction = 1,
                                    diverging.palette = "RdBu",
                                    diverging.direction = -1,
                                    enforce_symmetry = TRUE,
                                    legend.position = "bottom",
                                    legend.width = 1,
                                    legend.length = 20,
                                    legend.framewidth = 0.5,
                                    legend.tickwidth = 0.5,
                                    legend.framecolor = "grey50",
                                    legend.tickcolor = "white",
                                    legend.type = "colorbar",
                                    na.value = "grey75",
                                    font.size = 14,
                                    font.type = "sans",
                                    axis.text.x.angle = 45,
                                    flip = FALSE,
                                    colors.use = NULL,
                                    min.cutoff = NA,
                                    max.cutoff = NA,
                                    verbose = TRUE,
                                    return_object = FALSE,
                                    grid.color = "white",
                                    border.color = "black",
                                    flavor = "Seurat",
                                    nbin = 24,
                                    ctrl = 100,
                                    plot.title.face = "bold",
                                    plot.subtitle.face = "plain",
                                    plot.caption.face = "italic",
                                    axis.title.face = "bold",
                                    axis.text.face = "plain",
                                    legend.title.face = "bold",
                                    legend.text.face = "plain"){
  # Add lengthy error messages.
  withr::local_options(.new = list("warning.length" = 8170))

  check_suggests("do_AffinityAnalysisPlot")

  check_Seurat(sample)

  if (is.null(assay)){assay <- check_and_set_assay(sample)$assay}
  if (is.null(slot)){slot <- check_and_set_slot(slot)}

  # Check logical parameters.
  logical_list <- list("verbose" = verbose,
                       "flip" = flip,
                       "enforce_symmetry" = enforce_symmetry,
                       "use_viridis" = use_viridis)
  check_type(parameters = logical_list, required_type = "logical", test_function = is.logical)
  # Check numeric parameters.
  numeric_list <- list("font.size" = font.size,
                       "legend.length" = legend.length,
                       "legend.width" = legend.width,
                       "legend.framewidth" = legend.framewidth,
                       "legend.tickwidth" = legend.tickwidth,
                       "subsample" = subsample,
                       "viridis.direction" = viridis.direction,
                       "axis.text.x.angle" = axis.text.x.angle,
                       "min.cutoff" = min.cutoff,
                       "max.cutoff" = max.cutoff,
                       "number.breaks" = number.breaks,
                       "sequential.direction" = sequential.direction,
                       "nbin" = nbin,
                       "ctrl" = ctrl,
                       "diverging.direction" = diverging.direction)
  check_type(parameters = numeric_list, required_type = "numeric", test_function = is.numeric)
  # Check character parameters.
  character_list <- list("group.by" = group.by,
                         "assay" = assay,
                         "slot" = slot,
                         "statistic" = statistic,
                         "legend.type" = legend.type,
                         "legend.position" = legend.position,
                         "legend.framecolor" = legend.framecolor,
                         "legend.tickcolor" = legend.tickcolor,
                         "font.type" = font.type,
                         "viridis.palette" = viridis.palette,
                         "diverging.palette" = diverging.palette,
                         "sequential.palette" = sequential.palette,
                         "grid.color" = grid.color,
                         "border.color" = border.color,
                         "flavor" = flavor,
                         "plot.title.face" = plot.title.face,
                         "plot.subtitle.face" = plot.subtitle.face,
                         "plot.caption.face" = plot.caption.face,
                         "axis.title.face" = axis.title.face,
                         "axis.text.face" = axis.text.face,
                         "legend.title.face" = legend.title.face,
                         "legend.text.face" = legend.text.face,
                         "na.value" = na.value)
  check_type(parameters = character_list, required_type = "character", test_function = is.character)

  `%>%` <- magrittr::`%>%`

  check_colors(grid.color, parameter_name = "grid.color")
  check_colors(na.value, parameter_name = "na.value")
  check_colors(border.color, parameter_name = "border.color")
  check_colors(legend.tickcolor, parameter_name = "legend.tickcolor")
  check_colors(legend.framecolor, parameter_name = "legend.framecolor")

  check_parameters(parameter = font.type, parameter_name = "font.type")
  check_parameters(parameter = legend.position, parameter_name = "legend.position")
  check_parameters(plot.title.face, parameter_name = "plot.title.face")
  check_parameters(plot.subtitle.face, parameter_name = "plot.subtitle.face")
  check_parameters(plot.caption.face, parameter_name = "plot.caption.face")
  check_parameters(axis.title.face, parameter_name = "axis.title.face")
  check_parameters(axis.text.face, parameter_name = "axis.text.face")
  check_parameters(legend.title.face, parameter_name = "legend.title.face")
  check_parameters(legend.text.face, parameter_name = "legend.text.face")
  check_parameters(viridis.direction, parameter_name = "viridis.direction")
  check_parameters(sequential.direction, parameter_name = "sequential.direction")
  check_parameters(diverging.direction, parameter_name = "diverging.direction")

  # Assign a group.by if this is null.
  out <- check_group_by(sample = sample,
                        group.by = group.by,
                        is.heatmap = TRUE)
  sample <- out[["sample"]]
  group.by <- out[["group.by"]]

  if (!is.na(subsample)){
    sample <- sample[, sample(colnames(sample), subsample)]
  }

  # Generate the continuous color palette.
  if (isTRUE(enforce_symmetry)){
    colors.gradient <- compute_continuous_palette(name = diverging.palette,
                                                  use_viridis = FALSE,
                                                  direction = diverging.direction,
                                                  enforce_symmetry = enforce_symmetry)
  } else {
    colors.gradient <- compute_continuous_palette(name = ifelse(isTRUE(use_viridis), viridis.palette, sequential.palette),
                                                  use_viridis = use_viridis,
                                                  direction = ifelse(isTRUE(use_viridis), viridis.direction, sequential.direction),
                                                  enforce_symmetry = enforce_symmetry)
  }


  # Generate a network with the names of the list of genes as source and the gene sets as targets with 1 of mode of regulation.
  # Step 1: Check for underscores in the names of the gene sets.
  if (length(unlist(stringr::str_match_all(names(input_gene_list), "_"))) > 0){
    warning(paste0(add_warning(), crayon_body("Found "),
                   crayon_key("underscores (_)"),
                   crayon_body(" in the name of the gene sets provided. Replacing them with "),
                   crayon_key("dots (.)"),
                   crayon_body(" to avoid conflicts when generating the Seurat assay.")), call. = FALSE)
    names.use <- stringr::str_replace_all(names(input_gene_list), "_", ".")
    names(input_gene_list) <- names.use
  }

  # Step 2: make the lists of equal length.
  max_value <- max(unname(unlist(lapply(input_gene_list, length))))
  min_value <- min(unname(unlist(lapply(input_gene_list, length))))

  assertthat::assert_that(length(input_gene_list) >= 2,
                          msg = paste0(add_cross,
                                       crayon_body("Please make sure that the gene list you provide to "),
                                       crayon_key("input_gene_list"),
                                       crayon_body(" have at least "),
                                       crayon_key("two different"),
                                       crayon_body(" gene sets.")))

  assertthat::assert_that(min_value >= 5,
                          msg = paste0(add_cross,
                                       crayon_body("Please make sure that the gene list you provide to "),
                                       crayon_key("input_gene_list"),
                                       crayon_body(" have at least "),
                                       crayon_key("five genes"),
                                       crayon_body(" each.")))

  # Add fake genes until all lists have the same length so that it can be converted into a tibble.
  gene_list <- lapply(input_gene_list, function(x){
    if (length(x) != max_value){
      remaining <- max_value - length(x)
      x <- append(x, rep("deleteme", remaining))
      x
    } else{
      x
    }
  })

  # Generate the network as a tibble and filter out fake genes.
  network <- gene_list %>%
             tibble::as_tibble() %>%
             tidyr::pivot_longer(cols = dplyr::everything(),
                                 names_to = "source",
                                 values_to = "target") %>%
             dplyr::mutate("mor" = 1) %>%
             dplyr::filter(.data$target != "deleteme")

  # Get expression data.
  suppressWarnings({
  mat <- SeuratObject::GetAssayData(sample,
                                    assay = assay,
                                    slot = slot)
  })
  # Compute activities.
  if(isTRUE(verbose)){message(paste0(add_info(), crayon_body("Computing "),
                                     crayon_key("activities"),
                                     crayon_body("...")))}
  
  if (statistic == "ulm"){
    acts <- decoupleR::run_ulm(mat = mat, 
                               network = network)
  } else {
    acts <- decoupleR::run_wmean(mat = mat,
                                 network = network)
  }
  

  # Turn them into a matrix compatible to turn into a Seurat assay.
  acts.matrix <- acts %>%
                 dplyr::filter(.data$statistic == .env$statistic) %>%
                 tidyr::pivot_wider(id_cols = dplyr::all_of("source"),
                                    names_from = "condition",
                                    values_from = "score") %>%
                 tibble::column_to_rownames('source')

  # Generate a Seurat assay.
  assay.add <- Seurat::CreateAssayObject(acts.matrix)

  # Add the assay to the Seurat object.
  sample@assays$affinity <- assay.add
  sample@assays$affinity@key <- "affinity_"

  # Set it as default assay.
  Seurat::DefaultAssay(sample) <- "affinity"

  # Scale and center the activity data.
  sample <- Seurat::ScaleData(sample, verbose = FALSE, assay = "affinity")

  # Plotting.
  # Get the data frames per group.by value for plotting.
  list.data <- list()
  counter <- 0
  for (group in group.by){
    counter <- counter + 1
    suppressWarnings({
    data.use <- SeuratObject::GetAssayData(sample,
                             assay = "affinity",
                             slot = "scale.data") %>%
                t() %>%
                as.data.frame() %>%
                tibble::rownames_to_column(var = "cell") %>%
                dplyr::left_join(y = {sample@meta.data %>%
                                      tibble::rownames_to_column(var = "cell") %>%
                                      dplyr::select(dplyr::all_of(c("cell", group)))},
                                      by = "cell") %>%
                tidyr::pivot_longer(cols = -dplyr::all_of(c("cell", group)),
                                    names_to = "source",
                                    values_to = "score")
    })
    # Clustering based on the median across all cells.
    data.cluster <- data.use %>%
                    tidyr::pivot_wider(id_cols = dplyr::all_of(c("cell", group)),
                                       names_from = "source",
                                       values_from = "score") %>%
                    dplyr::group_by(.data[[group]]) %>%
                    dplyr::summarise(dplyr::across(.cols = dplyr::all_of(c(names(input_gene_list))),
                                                   function(x){stats::median(x, na.rm = TRUE)})) %>%
                    as.data.frame() %>%
                    tibble::column_to_rownames(var = group)

    list.data[[group]][["data"]] <- data.use
    list.data[[group]][["data.cluster"]] <- data.cluster
  }

  # Plot individual heatmaps.

  list.heatmaps <- list()
  counter <- 0
  row.order.list <- list()
  for (group in group.by){
    counter <- counter + 1

    data.use <- list.data[[group]][["data"]]
    data.cluster <- list.data[[group]][["data.cluster"]]

    # nocov start
    if (counter == 1){
      if (length(colnames(data.cluster)) == 1){
        col_order <- colnames(data.cluster)[1]
      } else {
        col_order <- colnames(data.cluster)[stats::hclust(stats::dist(t(data.cluster), method = "euclidean"), method = "ward.D")$order]
      }
    }
    # nocov end

    if(length(rownames(data.cluster)) == 1){
      row_order <- rownames(data.cluster)[1]
    } else {
      row_order <- rownames(data.cluster)[stats::hclust(stats::dist(data.cluster, method = "euclidean"), method = "ward.D")$order]
    }
    row.order.list[[group]] <- row_order

    data.use <- data.use %>%
                dplyr::group_by(.data[[group]], .data$source) %>%
                dplyr::summarise("mean" = mean(.data$score, na.rm = TRUE))

    list.data[[group]][["data.mean"]] <- data.use

    if (!is.na(min.cutoff)){
      data.use <- data.use %>%
                  dplyr::mutate("mean" = ifelse(.data$mean < min.cutoff, min.cutoff, .data$mean))
    }

    if (!is.na(max.cutoff)){
      data.use <- data.use %>%
                  dplyr::mutate("mean" = ifelse(.data$mean > max.cutoff, max.cutoff, .data$mean))
    }
    p <- data.use %>%
         dplyr::mutate("source" = factor(.data$source, levels = col_order),
                       "target" = factor(.data[[group]], levels = row_order)) %>%
         # nocov start
         ggplot2::ggplot(mapping = ggplot2::aes(x = if (isTRUE(flip)){.data$source} else {.data$target},
                                                y = if (isTRUE(flip)){.data$target} else {.data$source},
                                                fill = .data$mean)) +
         # nocov end
         ggplot2::geom_tile(color = grid.color, linewidth = 0.5, na.rm = TRUE) +
         ggplot2::scale_y_discrete(expand = c(0, 0)) +
         ggplot2::scale_x_discrete(expand = c(0, 0),
                                   position = "top") +
         # nocov start
         ggplot2::guides(y.sec = guide_axis_label_trans(~paste0(levels(if (isTRUE(flip)){.data$target} else {.data$source}))),
                         x.sec = guide_axis_label_trans(~paste0(levels(if (isTRUE(flip)){.data$source} else {.data$target})))) +
         # nocov end
         ggplot2::coord_equal()
    list.heatmaps[[group]] <- p
  }


  # Compute limits.
  min.vector <- NULL
  max.vector <- NULL

  for (group in group.by){
    data.limits <-  list.data[[group]][["data.mean"]]

    min.vector <- append(min.vector, min(data.limits$mean, na.rm = TRUE))
    max.vector <- append(max.vector, max(data.limits$mean, na.rm = TRUE))
  }

  # Get the absolute limits of the datasets.
  limits <- c(min(min.vector, na.rm = TRUE),
              max(max.vector, na.rm = TRUE))

  # Compute overarching scales for all heatmaps.
  scale.setup <- compute_scales(sample = sample,
                                feature = " ",
                                assay = assay,
                                reduction = NULL,
                                slot = slot,
                                number.breaks = number.breaks,
                                min.cutoff = min.cutoff,
                                max.cutoff = max.cutoff,
                                flavor = "Seurat",
                                enforce_symmetry = enforce_symmetry,
                                from_data = TRUE,
                                limits.use = limits)

  for (group in group.by){
    p <- list.heatmaps[[group]]

    p <- p +
         ggplot2::scale_fill_gradientn(colors = colors.gradient,
                                       na.value = na.value,
                                       name = paste0(statistic, " | Scaled and Centered"),
                                       breaks = scale.setup$breaks,
                                       labels = scale.setup$labels,
                                       limits = scale.setup$limits)

    list.heatmaps[[group]] <- p
  }

  # Modify legends.
  for (group in group.by){
    p <- list.heatmaps[[group]]

    p <- modify_continuous_legend(p = p,
                                  legend.aes = "fill",
                                  legend.type = legend.type,
                                  legend.position = legend.position,
                                  legend.length = legend.length,
                                  legend.width = legend.width,
                                  legend.framecolor = legend.framecolor,
                                  legend.tickcolor = legend.tickcolor,
                                  legend.framewidth = legend.framewidth,
                                  legend.tickwidth = legend.tickwidth)
    list.heatmaps[[group]] <- p
  }

  # Add theme
  counter <- 0
  for (group in group.by){
    counter <- counter + 1

    p <- list.heatmaps[[group]]

    # Set axis titles.
    if (isTRUE(flip)){
      if (counter == 1){
        ylab <- group
        xlab <- NULL
        if (length(group.by) == counter){
          xlab <- "Gene set"
        }
      } else {
        xlab <- "Gene set"
        ylab <- group
      }
    } else {
      if (counter == 1){
        ylab <- "Gene set"
        xlab <- group
      } else {
        ylab <- NULL
        xlab <- group
      }
    }


    p <- list.heatmaps[[group]]

    axis.parameters <- handle_axis(flip = !flip,
                                   group.by = rep("A", length(group.by)),
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
                        strip.background = axis.parameters$strip.background,
                        strip.clip = axis.parameters$strip.clip,
                        strip.text = axis.parameters$strip.text,
                        legend.position = legend.position,
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
                        plot.background = ggplot2::element_rect(fill = "white", color = "white"),
                        panel.background = ggplot2::element_rect(fill = "white", color = "white"),
                        legend.background = ggplot2::element_rect(fill = "white", color = "white"),
                        panel.spacing.x = ggplot2::unit(0, "cm"))

    list.heatmaps[[group]] <- p
  }


  if (isTRUE(flip)){
    list.heatmaps <- list.heatmaps[rev(group.by)]
  }
  p <- patchwork::wrap_plots(list.heatmaps,
                             ncol = if (base::isFALSE(flip)){NULL} else {1},
                             nrow = if(base::isFALSE(flip)){1} else {NULL},
                             guides = "collect")
  p <- p +
       patchwork::plot_annotation(theme = ggplot2::theme(legend.position = legend.position,
                                                         plot.title = ggplot2::element_text(family = font.type,
                                                                                            color = "black",
                                                                                            face = "bold",
                                                                                            hjust = 0),
                                                         plot.subtitle = ggplot2::element_text(family = font.type,
                                                                                               color = "black",
                                                                                               hjust = 0),
                                                         plot.caption = ggplot2::element_text(family = font.type,
                                                                                              color = "black",
                                                                                              hjust = 1),
                                                         plot.caption.position = "plot"))

  list.output <- list()

  list.output[["Heatmap"]] <- p


  if (isTRUE(return_object)){
    list.output[["Object"]] <- sample
  }

  if (isTRUE(return_object)){
    return_me <- list.output
  } else {
    return_me <- list.output$Heatmap
  }

  return(return_me)
}
