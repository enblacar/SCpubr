do_CellPhoneDBPlot <- function(sample = sample,
                               group.by = NULL,
                               keep_source = NULL,
                               keep_target = NULL,
                               top_interactions = 25,
                               assay = "SCT",
                               dot_border = TRUE,
                               border.color = "black",
                               rotate_x_labels = TRUE,
                               legend.position = "bottom",
                               legend.type = "colorbar",
                               legend.length = 20,
                               legend.width = 1,
                               legend.framecolor = "grey50",
                               legend.tickcolor = "white",
                               legend.framewidth = 1.5,
                               legend.tickwidth = 1.5,
                               viridis_color_map = "G",
                               verbose = TRUE){


  # Check border color.
  check_colors(border.color, parameter_name = "border.color")

  # Check viridis_color_map.
  check_viridis_color_map(viridis_color_map = viridis_color_map, verbose = verbose)

  # Check the colors provided to legend.framecolor and legend.tickcolor.
  check_colors(legend.framecolor, parameter_name = "legend.framecolor")
  check_colors(legend.tickcolor, parameter_name = "legend.tickcolor")

  # Check font.type.
  if (!(font.type %in% c("sans", "serif", "mono"))){
    stop("Please select one of the following for font.type: sans, serif, mono.", call. = F)
  }

  # Check the legend.type.
  if (!(legend.type %in% c("normal", "colorbar", "colorsteps"))){
    stop("Please select one of the following for legend.type: normal, colorbar, colorsteps.", call. = FALSE)
  }

  # Check the legend.position.
  if (!(legend.position %in% c("top", "bottom", "left", "right"))){
    stop("Please select one of the following for legend.position: top, bottom, left, right.", call. = FALSE)
  }

  # Cellphonedb
  if (is.null(group.by)){
    sample$dummy <- Seurat::Idents(sample)
  } else {
    sample$dummy <- sample@meta.data[, group.by]
  }
  group.by <- "dummy"

  sample.use <- sample[, sample$dummy %in% unique(c(keep_source, keep_target))]
  Seurat::Idents(sample.use) <- sample.use$dummy

  output.cb <- liana::liana_wrap(sce = sample.use,
                                 method = c("cellphonedb"),
                                 idents_col = group.by,
                                 verbose = verbose,
                                 assay = assay)
  output.cb <- readRDS("/b06x-isilon/b06x-g/G703/eblanco/projects/test_SC_datasets/liana_cellphoneDB.rds")

  output.cb <- output.cb %>%
    dplyr::filter(pvalue <= 0.05) %>%
    dplyr::mutate(pvalue = -log10(.data$pvalue + 0.0000000001)) %>%
    dplyr::arrange(dplyr::desc(.data$lr.mean), .data$pvalue) %>%
    tidyr::unite(c("ligand", "receptor"), col = "interaction", sep = "<span style = 'color:grey50;'> | </span>") %>%
    tidyr::unite(c("source", "target"), col = "interacting_clusters", remove = FALSE) %>%
    tidyr::unite(c("interaction", "interacting_clusters"), col = "filter_me", remove = FALSE)
  output.cb <- output.cb %>%
    dplyr::filter(.data$interaction %in% {output.cb %>% dplyr::top_n(top_interactions, dplyr::desc(.data$pvalue)) %>% dplyr::pull(.data$interaction) %>% unique()})
  if (!(is.null(keep_source))){
    output.cb <- output.cb %>%
      dplyr::filter(.data$source %in% keep_source)
  }

  if (!(is.null(keep_target))){
    output.cb <- output.cb %>%
      dplyr::filter(.data$target %in% keep_target)
  }


  # Define legend parameters. Width and height values will change depending on the legend orientation.
  if (legend.position %in% c("top", "bottom")){
    legend.barwidth <- legend.length
    legend.barheight <- legend.width
    size_title <- "Interaction specificity"
    fill.title <- "Expression Magnitude"
  } else if (legend.position %in% c("left", "right")){
    legend.barwidth <- legend.width
    legend.barheight <- legend.length
    size_title <- stringr::str_wrap("Interaction specificity", width = 10)
    fill.title <- stringr::str_wrap("Expression Magnitude", width = 10)
  }


  p <-  output.cb %>%
    ggplot2::ggplot(mapping = ggplot2::aes(x = .data$target,
                                           y = .data$interaction,
                                           fill = if(isTRUE(dot_border)){.data$lr.mean} else {NULL},
                                           size = .data$pvalue,
                                           group = .data$interacting_clusters)) +
    ggplot2::geom_point(mapping = ggplot2::aes(color = if(isTRUE(dot_border)){NULL} else {.data$lr.mean}),
                        shape = if(isTRUE(dot_border)){21} else {19}) +
    ggplot2::scale_size_continuous(name = size_title)
  if (isTRUE(dot_border)){
    p$layers[[1]]$aes_params$color <- border.color
    p <- p +
      ggplot2::scale_fill_viridis_c(option = viridis_color_map,
                                    name = fill.title)
  } else {
    p <- p +
      ggplot2::scale_color_viridis_c(option = viridis_color_map,
                                     name = fill.title)
  }
  p <- p +
    ggplot2::facet_grid(. ~ .data$source,
                        space = "free",
                        scales = "free",
                        switch = "y") +
    ggplot2::labs(title = "Source") +
    ggplot2::xlab("Target") +
    ggplot2::ylab(paste("Ligand", "|", "Receptor", sep = " ")) +
    ggplot2::guides(size = ggplot2::guide_legend(title.position = "top",
                                                 title.hjust = 0.5,
                                                 override.aes = ggplot2::aes(fill = "black"))) +
    ggplot2::theme_minimal(base_size = font.size) +
    ggplot2::theme(plot.title = ggplot2::element_text(face = "bold",
                                                      hjust = 0.5,
                                                      size = font.size),
                   plot.subtitle = ggtext::element_markdown(hjust = 0),
                   plot.caption = ggtext::element_markdown(hjust = 1),
                   plot.title.position = "panel",
                   plot.caption.position = "plot",
                   text = ggplot2::element_text(family = font.type),
                   legend.justification = "center",
                   legend.text = ggplot2::element_text(face = "bold"),
                   legend.title = ggplot2::element_text(face = "bold"),
                   legend.position = legend.position,
                   axis.title.x = ggplot2::element_text(face = "bold", hjust = 0.5),
                   axis.title.y = ggplot2::element_text(face = "bold", angle = 90),
                   axis.text.y = ggtext::element_markdown(face = "bold"),
                   axis.text = ggplot2::element_text(face = "bold", color = "black"),
                   axis.text.x = ggplot2::element_text(angle = ifelse(isTRUE(rotate_x_labels), 90, 0),
                                                       hjust = ifelse(isTRUE(rotate_x_labels), 1, 0.5),
                                                       vjust = ifelse(isTRUE(rotate_x_labels), 0.5, 1)),
                   strip.text = ggplot2::element_text(face = "bold"),
                   panel.grid = ggplot2::element_blank(),
                   plot.margin = ggplot2::margin(t = 10, r = 10, b = 10, l = 10),
                   plot.background = ggplot2::element_rect(fill = "white", color = "white"),
                   panel.background = ggplot2::element_rect(fill = "white", color = "black"),
                   legend.background = ggplot2::element_rect(fill = "white", color = "white"))

  if (legend.type == "normal"){
    if (isTRUE(dot_border)){
      p <- p +
        ggplot2::guides(fill = ggplot2::guide_colorbar(title.position = "top",
                                                       title.hjust = 0.5))
    } else {
      p <- p +
        ggplot2::guides(color = ggplot2::guide_colorbar(title.position = "top",
                                                        title.hjust = 0.5))
    }

  } else if (legend.type == "colorbar"){
    if (isTRUE(dot_border)){
      p <- p +
        ggplot2::guides(fill = ggplot2::guide_colorbar(title.position = "top",
                                                       barwidth = legend.barwidth,
                                                       barheight = legend.barheight,
                                                       title.hjust = 0.5,
                                                       ticks.linewidth = legend.tickwidth,
                                                       frame.linewidth = legend.framewidth,
                                                       frame.colour = legend.framecolor,
                                                       ticks.colour = legend.tickcolor))
    } else {
      p <- p +
        ggplot2::guides(color = ggplot2::guide_colorbar(title.position = "top",
                                                        barwidth = legend.barwidth,
                                                        barheight = legend.barheight,
                                                        title.hjust = 0.5,
                                                        ticks.linewidth = legend.tickwidth,
                                                        frame.linewidth = legend.framewidth,
                                                        frame.colour = legend.framecolor,
                                                        ticks.colour = legend.tickcolor))
    }
  } else if (legend.type == "colorsteps"){
    if (isTRUE(dot_border)){
      p <- p +
        ggplot2::guides(fill = ggplot2::guide_colorsteps(title.position = "top",
                                                         barwidth = legend.barwidth,
                                                         barheight = legend.barheight,
                                                         title.hjust = 0.5,
                                                         ticks.linewidth = legend.tickwidth,
                                                         frame.linewidth = legend.framewidth,
                                                         frame.colour = legend.framecolor,
                                                         ticks.colour = legend.tickcolor))
    } else {
      p <- p +
        ggplot2::guides(color = ggplot2::guide_colorsteps(title.position = "top",
                                                          barwidth = legend.barwidth,
                                                          barheight = legend.barheight,
                                                          title.hjust = 0.5,
                                                          ticks.linewidth = legend.tickwidth,
                                                          frame.linewidth = legend.framewidth,
                                                          frame.colour = legend.framecolor,
                                                          ticks.colour = legend.tickcolor))
    }
  }

}


