#' Generate a Geyser plot.
#'
#' A Geyser plot is a custom plot in which we plot continuous values on the Y axis grouped by a categorical value in the X. This is plotted as a dot plot, jittered so that the dots span
#' all the way to the other groups. On top of this, the mean and .66 and .95 of the data is plotted, depicting the overall distribution of the dots. The cells can, then, be colored by
#' a continuous variable (same as Y axis or different) or a categorical one (same as X axis or different).
#'
#' Special thanks to Christina Blume for coming up with the name of the plot.
#'
#' @param sample Seurat object.
#' @param features Character. Features to plot.
#' @param assay Character. Assay to use. Defaults to active assay if not.
#' @param slot Character. Slot to use. Defaults to data slot if not.
#' @param group.by Character. Metadata variable to group the data by.
#' @param split.by Character. Metadata variable to further split the plot by.
#' @param symmetrical_scale Logical. Whether you want a continuous scale that is symmetrical on both ends and the intensity of the colors in each end matches (blue to red).
#' @param scale_type Character. Either continuous or categorical.
#' @param color.by Character. Vector of colors to color the groups by if scale_type is categorical.
#' @param order_by_mean Logical. Whether to order the groups by the mean of the data (highest to lowest).
#' @param plot_cell_borders Logical. Whether to plot border around the cells.
#' @param jitter Numeric. Amount of jitter in the plot along the X axis. The lower the value, the more compacted the dots are.
#' @param pt.size Numeric. Size of the cells.
#' @param border.size Numeric. Thickness of the border around the cells.
#' @param border.color Character. Color for the border of the cells.
#' @param show_legend Logical. Whether to display the legend. Might be redundant depending on the setups.
#' @param legend.title Character. Title for the legend.
#' @param legend.type Character. Type of legend to display. One of: normal, colorbar, colorsteps.
#' @param legend.position Character. Position of the legend in the plot. Will only work if legend is set to TRUE.
#' @param legend.framewidth,legend.tickwidth Numeric. Width of the lines of the box in the legend.
#' @param legend.framecolor,legend.tickcolor Character. Color of the lines of the box in the legend.
#' @param legend.length,legend.width Numeric. Length and width of the legend. Will adjust automatically depending on legend side.
#' @param legend.ncol,legend.nrow Numeric. Number of rows/columns of a categorical legend.
#' @param legend.icon.size Numeric. Size of the dots in a categorical legend.
#' @param legend.byrow Logical. Whether to fill the legend by rows instead of columns in a categorical legend.
#' @param font.size Numeric. Overall font size of the plot.
#' @param font.type Character. Font family for the plot: sans, mono, serif.
#' @param rotate_x_axis_labels Logical. Whether to rotate X axis labels.
#' @param viridis_color_map Character. Viridis color map to use. One of: A, B, C, D, E, F, G, H.
#' @param colors.use Character. Named vector of colors to use. Has to match the unique values of group.by or color.by (if used) when scale_type is set to categorical.
#' @param na.value Character. Color for NAs.
#' @param plot.title,plot.subtitle,plot.caption Character. Title, Subtitle and caption to use in the plot.
#' @param xlab,ylab Character. Titles for the X and Y axis.
#'
#' @return Either a plot of a list of plots, depending on the number of features provided.
#' @export
#' @examples
#' \dontrun{
#' TBD
#' }
do_GeyserPlot <- function(sample,
                          features,
                          assay = NULL,
                          slot = "data",
                          group.by = NULL,
                          split.by = NULL,
                          symmetrical_scale = FALSE,
                          scale_type = "continuous",
                          color.by = NULL,
                          order_by_mean = TRUE,
                          plot_cell_borders = TRUE,
                          jitter = 0.45,
                          pt.size = 1,
                          border.size = 1.5,
                          border.color = "black",
                          legend.position = "bottom",
                          legend.width = 1,
                          legend.length = 20,
                          legend.framewidth = 1.5,
                          legend.tickwidth = 1.5,
                          legend.framecolor = "grey50",
                          legend.tickcolor = "white",
                          legend.type = "colorbar",
                          font.size = 14,
                          font.type = "sans",
                          rotate_x_axis_labels = FALSE,
                          viridis_color_map = "D",
                          colors.use = NULL,
                          na.value = "grey75",
                          legend.ncol = NULL,
                          legend.nrow = NULL,
                          legend.icon.size = 4,
                          legend.byrow = FALSE,
                          legend.title = NULL,
                          show_legend = TRUE,
                          plot.title = NULL,
                          plot.subtitle = NULL,
                          plot.caption = NULL,
                          xlab = NULL,
                          ylab = NULL){

  # Checks for packages.
  check_suggests(function_name = "do_GeyserPlot")
  # Check if the sample provided is a Seurat object.
  check_Seurat(sample = sample)

  # Check logical parameters.
  logical_list <- list("symmetrical_scale" = symmetrical_scale,
                       "order_by_mean" = order_by_mean,
                       "rotate_x_axis_labels" = rotate_x_axis_labels,
                       "plot_cell_borders" = plot_cell_borders)
  check_type(parameters = logical_list, required_type = "logical", test_function = is.logical)
  # Check numeric parameters.
  numeric_list <- list("pt.size" = pt.size,
                       "jitter" = jitter,
                       "font.size" = font.size,
                       "legend.framewidth" = legend.framewidth,
                       "legend.tickwidth" = legend.tickwidth,
                       "legend.length" = legend.length,
                       "legend.width" = legend.width,
                       "border.size" = border.size,
                       "legend.ncol" = legend.ncol,
                       "legend.nrow" = legend.nrow,
                       "legend.icon.size" = legend.icon.size)
  check_type(parameters = numeric_list, required_type = "numeric", test_function = is.numeric)
  # Check character parameters.

  character_list <- list("legend.position" = legend.position,
                         "features" = features,
                         "assay" = assay,
                         "group.by" = group.by,
                         "slot" = slot,
                         "split.by" = split.by,
                         "plot.title" = plot.title,
                         "plot.subtitle" = plot.subtitle,
                         "plot.caption" = plot.caption,
                         "scale_type" = scale_type,
                         "viridis_color_map" = viridis_color_map,
                         "legend.framecolor" = legend.framecolor,
                         "legend.tickcolor" = legend.tickcolor,
                         "legend.type" = legend.type,
                         "font.type" = font.type,
                         "border.color" = border.color,
                         "na.value" = na.value)
  # Checks
  check_type(parameters = character_list, required_type = "character", test_function = is.character)

  check_colors(border.color, parameter_name = "border.color")
  check_colors(legend.tickcolor, parameter_name = "legend.tickcolor")
  check_colors(legend.framecolor, parameter_name = "legend.framecolor")
  check_colors(na.value, parameter_name = "na.value")

  `%>%` <- purrr::`%>%`
  # Check the assay.
  out <- check_and_set_assay(sample = sample, assay = assay)
  sample <- out[["sample"]]
  assay <- out[["assay"]]
  rm(out)

  # Check if the scale is properly set.
  if (!(scale_type %in% c("categorical", "continuous"))){
    stop("Please provide one of the following to scale_type: continuous, categorical.", call. = FALSE)
  }

  # Check that split.by is in metadata variables.
  if (!(is.null(split.by))){
    if (!(split.by %in% colnames(sample@meta.data))){
      stop("The variable for split.by has to be on the metadata of the object.", call. = FALSE)
    }
  }

  # Check that group.by is in metadata variables.
  if (!(is.null(group.by))){
    if (!(group.by %in% colnames(sample@meta.data))){
      stop("The variable for group.by has to be on the metadata of the object.", call. = FALSE)
    }
  }

  # Check that jitter is in range.
  if (jitter >= 0.5 | jitter < 0){
    stop("Value for jitter has to be betwen 0 and 0.49.", call. = FALSE)
  }

  # Will contain the output.
  list.out <- list()

  # Assign group.by to a metadata variable.
  if (is.null(group.by)){
    sample@meta.data[, "dummy"] <- sample@active.ident
  } else {
    sample@meta.data[, "dummy"] <- sample@meta.data[, group.by]
  }
  group.by <- "dummy"

  # Iterate for each feature.
  for (feature in features){
    # Check the feature.
    check_feature(sample = sample,
                  features = feature)

    
    # If the user wants additional coloring, if not default to feature or group.by.
    if (isTRUE(scale_type == "continuous")){
      if (is.null(color.by)){
        color.by <- feature
      }
    } else if (isTRUE(scale_type == "categorical")){
      if (is.null(color.by)){
        color.by <- group.by
      } else {
        if (!(color.by %in% colnames(sample@meta.data))){
          stop("With a categorical scale, color.by needs to be present in sample@meta.data.", call. = FALSE)
        }
      }
    }
    
    # Get a vector of all dimensional reduction compontents.
    dim_colnames <- c()
    for(red in Seurat::Reductions(object = sample)){
      col.names <- colnames(sample@reductions[[red]][[]])
      dim_colnames <- c(dim_colnames, col.names)
      if (feature %in% col.names){
        # Get the reduction in which the feature is, if this is the case.
        reduction <- red
      }
      
      if (color.by %in% col.names){
        # Get the reduction in which the feature is, if this is the case.
        reduction_color.by <- red
      }
    }
    

    # Generate a column for the color.by parameter that will be added later on to the data dataframe.
    if (isTRUE(color.by %in% colnames(sample@meta.data))){
      color.by_column <- sample@meta.data %>%
                         dplyr::select(.data[[color.by]]) %>%
                         tibble::rownames_to_column(var = "cell") %>%
                         dplyr::rename("color.by" = .data[[color.by]])
    } else if (isTRUE(color.by %in% rownames(sample))){
      color.by_column <- Seurat::GetAssayData(object = sample,
                                              assay = assay,
                                              slot = slot)[color.by, , drop = F] %>%
                         as.matrix() %>%
                         t() %>%
                         as.data.frame() %>%
                         tibble::rownames_to_column(var = "cell") %>%
                         dplyr::rename("color.by" = .data[[color.by]])
    } else if (isTRUE(color.by %in% dim_colnames)){
      color.by_column <- sample@reductions[[reduction_color.by]][[]][, color.by, drop = FALSE] %>%
                         as.data.frame() %>%
                         tibble::rownames_to_column(var = "cell") %>%
                         dplyr::rename("color.by" = .data[[color.by]])
    }


    # Depending on where the feature is, generate a tibble accordingly.
    if (isTRUE(feature %in% colnames(sample@meta.data))){
      data <- sample@meta.data %>%
              dplyr::select(c(.data[[group.by]], .data[[feature]])) %>%
              tibble::rownames_to_column(var = "cell") %>%
              dplyr::left_join(y = color.by_column,
                               by = "cell") %>%
              tibble::as_tibble()
    } else if (isTRUE(feature %in% rownames(sample))){
      data <- Seurat::GetAssayData(object = sample,
                                   assay = assay,
                                   slot = slot)[feature, , drop = F] %>%
              as.matrix() %>%
              t() %>%
              as.data.frame() %>%
              tibble::rownames_to_column(var = "cell") %>%
              tibble::tibble() %>%
              dplyr::left_join(y = color.by_column,
                               by = "cell") %>%
              dplyr::left_join(y = {sample@meta.data %>%
                                    dplyr::select(.data[[group.by]]) %>%
                                    tibble::rownames_to_column(var = "cell")},
                               by = "cell")
    } else if (isTRUE(feature %in% dim_colnames)){
      data <- sample@reductions[[reduction]][[]][, feature, drop = FALSE] %>%
              as.data.frame() %>%
              tibble::rownames_to_column(var = "cell") %>%
              dplyr::left_join(y = color.by_column,
                               by = "cell") %>%
              tibble::tibble() %>%
              dplyr::left_join(y = {sample@meta.data %>%
                                    dplyr::select(.data[[group.by]]) %>%
                                    tibble::rownames_to_column(var = "cell")},
                                    by = "cell")
    }

    # If we also want additional split.by.
    if (!(is.null(split.by))){
      data <- data %>%
              dplyr::left_join(y = {sample@meta.data %>%
                                    dplyr::select(.data[[split.by]]) %>%
                                    tibble::rownames_to_column(var = "cell")},
                               by = "cell") %>%
              dplyr::mutate("split.by" = .data[[split.by]]) %>%
              dplyr::select(-.data[[split.by]])

    }


    # Proceed with the regular plot.
    if (isTRUE(order_by_mean)){
      data <- data %>%
              dplyr::mutate("group.by" = factor(.data[[group.by]], levels = {data %>%
                                                                             dplyr::group_by(.data[[group.by]]) %>%
                                                                             dplyr::summarise("mean" = mean(.data[[feature]])) %>%
                                                                             dplyr::arrange(dplyr::desc(.data$mean)) %>%
                                                                             dplyr::pull(.data[[group.by]]) %>%
                                                                             as.character()}),
                            "values" = .data[[feature]])
    } else if (isFALSE(order_by_mean)){
      data <- data %>%
              dplyr::mutate("group.by" = .data[[group.by]],
                            "values" = .data[[feature]])
    }

    # Get the final column names.
    cols.use <- c("values", "group.by", "color.by")
    if (!(is.null(split.by))){
      cols.use <- append(cols.use, "split.by")
    }

    data <- data %>%
            dplyr::select(dplyr::all_of(cols.use))

    # Plot

    p <- ggplot2::ggplot(data = data,
                         mapping = ggplot2::aes(x = .data$group.by,
                                                y = .data$values,
                                                color = .data$color.by))

    if (isTRUE(plot_cell_borders)){
      p <- p +
           ggplot2::geom_point(position = ggplot2::position_jitter(width = jitter,
                                                                   seed = 0),
                               size = pt.size * border.size,
                               color = border.color,
                               na.rm = TRUE)
    }

    if (isTRUE(scale_type == "continuous")){
      if (isTRUE(symmetrical_scale)){
        limits <- c(min(data[, "color.by"]),
                    max(data[, "color.by"]))
        end_value <- max(abs(limits))
        scale.use <- ggplot2::scale_color_gradientn(colors = c("#033270", "#4091C9", "#fdf0d5", "#c94040", "#65010C"),
                                                    limits = c(-end_value, end_value),
                                                    na.value = na.value)
      } else if (isFALSE(symmetrical_scale)){
        scale.use <- ggplot2::scale_color_viridis_c(option = viridis_color_map,
                                                    na.value = na.value)
      }
    } else if (isTRUE(scale_type == "categorical")){
      if (is.null(colors.use)){
        values <- data %>% dplyr::pull(.data$color.by)
        names.use <- if (is.factor(values)){levels(values)} else {sort(unique(values))}
        colors.use <- generate_color_scale(names = names.use)
      } else {
        check_colors(colors.use)
      }
      scale.use <- ggplot2::scale_color_manual(values = colors.use,
                                               na.value = na.value)
    }
    p <- p +
         ggplot2::geom_point(position = ggplot2::position_jitter(width = jitter,
                                                                 seed = 0),
                             size = pt.size,
                             na.rm = TRUE) +
         ggdist::stat_pointinterval(interval_size_range = c(2, 3),
                                    fatten_point = 1.5,
                                    interval_color = "white",
                                    point_color = "white",
                                    position = ggplot2::position_dodge(width = 1),
                                    na.rm = TRUE,
                                    show.legend = FALSE) +
         ggdist::stat_pointinterval(interval_size_range = c(1, 2),
                                    interval_color = "black",
                                    point_color = "black",
                                    position = ggplot2::position_dodge(width = 1),
                                    na.rm = TRUE,
                                    show.legend = FALSE) +
         scale.use

    if (!(is.null(split.by))){
      p <- p +
           ggplot2::facet_grid(. ~ split.by)
    }
    p <- p +
         ggplot2::labs(title = plot.title,
                       subtitle = plot.subtitle,
                       caption = plot.caption) +
         ggplot2::xlab(if (is.null(xlab)) {"Groups"} else (xlab)) +
         ggplot2::ylab(if (is.null(ylab)) {feature} else (ylab)) +
         ggplot2::theme_minimal(base_size = font.size) +
         ggplot2::theme(axis.title = ggplot2::element_text(color = "black",
                                                           face = "bold"),
                        axis.line.x = ggplot2::element_line(color = "black"),
                        axis.text.x = ggplot2::element_text(color = "black",
                                                            face = "bold",
                                                            angle = ifelse(isTRUE(rotate_x_axis_labels), 90, 0),
                                                            hjust = ifelse(isTRUE(rotate_x_axis_labels), 1, 0.5),
                                                            vjust = ifelse(isTRUE(rotate_x_axis_labels), 0.5, 1)),
                        axis.text.y = ggplot2::element_text(color = "black", face = "bold"),
                        axis.ticks = ggplot2::element_line(color = "black"),
                        panel.grid.major = ggplot2::element_blank(),
                        plot.title.position = "plot",
                        plot.title = ggtext::element_markdown(face = "bold", hjust = 0),
                        plot.subtitle = ggtext::element_markdown(hjust = 0),
                        plot.caption = ggtext::element_markdown(hjust = 1),
                        panel.grid = ggplot2::element_blank(),
                        text = ggplot2::element_text(family = font.type),
                        plot.caption.position = "plot",
                        legend.text = ggplot2::element_text(face = "bold"),
                        legend.position = legend.position,
                        legend.title = ggplot2::element_text(face = "bold"),
                        legend.justification = "center",
                        plot.margin = ggplot2::margin(t = 10, r = 10, b = 10, l = 10),
                        plot.background = ggplot2::element_rect(fill = "white", color = "white"),
                        panel.background = ggplot2::element_rect(fill = "white", color = "white"),
                        legend.background = ggplot2::element_rect(fill = "white", color = "white"),
                        strip.text =ggplot2::element_text(color = "black", face = "bold"))

    if (isTRUE(scale_type == "continuous")){
      if (is.null(legend.title)){
        legend.title <- color.by
      }
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
    } else if (isTRUE(scale_type == "categorical")){
      if (is.null(legend.title)){
        legend.title <- ""
      }
      p <- p +
           ggplot2::guides(color = ggplot2::guide_legend(title = legend.title,
                                                         ncol = legend.ncol,
                                                         nrow = legend.nrow,
                                                         byrow = legend.byrow,
                                                         override.aes = list(size = legend.icon.size),
                                                         title.position = "top",
                                                         title.hjust = 0.5))
    }

    # Remove legend if necessary.
    if (isFALSE(show_legend)){
      p <- p +
           ggplot2::theme(legend.position = "none")
    }
    list.out[[feature]] <- p
  }
  return(if (length(features) > 1) {list.out} else {p})
}