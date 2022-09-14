#' Wrapper for computing publication ready bar plots.
#'
#' @param sample  Seurat object.
#' @param features  Main variable in the bar plot - This will define the height of the bars. Must be a metadata variable stored in `object@meta.data`.
#' @param group.by  Secondary variable to group the bar plot for - This will define the total number of bars. Must be a metadata variable stored in `object@meta.data`.
#' @param order.by  One of the unique values in group.by that will be used to reorder the bars according to its descending proportion.
#' @param labels.order  Vector of labels to explicitly state the order or the bars. Better fed after knowing the order. This is meant to be used to tweak the end plot if needed.
#' @param position  Either "fill" or "stack." Position "fill" will generate a bar plot with one column and the proportions of values for each group inside, while "stack" plots the bars together.
#' @param xlab  Title for the X axis.
#' @param ylab  Title for the Y axis.
#' @param colors.use  Palette of colors to use. It must match the group.by variable in terms of length and names.
#' @param font.size Base font.size of the figure.
#' @param font.type Character. Base font for the plot. One of mono, serif or sans.
#' @param legend Whether to plot the legend.
#' @param legend.title  Logical stating whether the legend title is shown or not.
#' @param legend.title.position Character stating where to place the title of the legend.
#' @param legend.position  Position of the legend in the plot.
#' @param legend.ncol,legend.nrow  Number of columns/rows in the legend.
#' @param legend.position  Position of the legend in the plot. Will only work if legend is set to TRUE.
#' @param legend.byrow  Logical stating whether the legend is filled by row or not.
#' @param plot.title,plot.subtitle,plot.caption  Title, subtitle or caption to use in the plot.
#' @param horizontal Whether to plot the Bar plot horizontally.
#' @param verbose Use warnings.
#' @param add.summary_labels Logical. Whether to add the total number of values on top of each bar. Only works with position = stack.
#' @param add.subgroup_labels Logical. Whether to add the total number of values for each group in the bar. Only works with position = stack.
#' @param repel.subgroup_labels,repel.summary_labels Logical. Whether to repel labels to avoid overplotting. This will result in labels not being aligned anymore.
#' @param size.labels Numeric. Modify the size of the labels.
#' @param rotate_x_labels Logical. Whether to rotate X axis labels to horizontal or not. If multiple features, a vector of logical values of the same length.
#' @param return_data_matrix Logical. Whether to also output the data matrix used to generate the bar plot. This is useful to report it for supplementary data.
#' @param plot_line_guides Logical. Whether to plot line guides for position = "stack".
#' @return A ggplot2 object containing a Bar plot.
#' #
#'
#' @example /man/examples/examples_do_BarPlot.R
do_BarPlot <- function(sample,
                       features,
                       group.by = NULL,
                       labels.order = NULL,
                       order.by = NULL,
                       position = "stack",
                       xlab = NULL,
                       ylab = NULL,
                       plot.title = NULL,
                       plot.subtitle = NULL,
                       plot.caption = NULL,
                       legend = TRUE,
                       legend.position = "right",
                       legend.title.position = "top",
                       legend.title = FALSE,
                       legend.ncol = NULL,
                       legend.nrow = NULL,
                       font.size = 14,
                       font.type = "sans",
                       legend.byrow = FALSE,
                       colors.use = NULL,
                       horizontal = FALSE,
                       verbose = TRUE,
                       add.summary_labels = FALSE,
                       add.subgroup_labels = FALSE,
                       repel.subgroup_labels = FALSE,
                       repel.summary_labels = FALSE,
                       size.labels = 3,
                       rotate_x_labels = NULL,
                       return_data_matrix = FALSE,
                       plot_line_guides = TRUE){
    # Checks for packages.
    check_suggests(function_name = "do_BarPlot")
    # Check if the sample provided is a Seurat object.
    check_Seurat(sample = sample)
    # Define pipe operator internally.
    `%>%` <- purrr::`%>%`

    # Check logical parameters.
    logical_list <- list("legend" = legend,
                         "legend.title" = legend.title,
                         "legend.byrow" = legend.byrow,
                         "horizontal" = horizontal,
                         "add.subgroup_labels" = add.subgroup_labels,
                         "add.summary_labels" = add.summary_labels,
                         "repel.subgroup_labels" = repel.subgroup_labels,
                         "repel.summary_labels" = repel.summary_labels,
                         "rotate_x_labels" = rotate_x_labels,
                         "return_data_matrix" = return_data_matrix)
    check_type(parameters = logical_list, required_type = "logical", test_function = is.logical)
    # Check numeric parameters.
    numeric_list <- list("font.size" = font.size,
                         "legend.ncol" = legend.ncol,
                         "legend.nrow" = legend.nrow,
                         "size.labels" = size.labels)
    check_type(parameters = numeric_list, required_type = "numeric", test_function = is.numeric)
    # Check character parameters.
    character_list <- list("legend.position" = legend.position,
                           "features" = features,
                           "group.by" = group.by,
                           "colors.use" = colors.use,
                           "plot.title" = plot.title,
                           "plot.subtitle" = plot.subtitle,
                           "plot.caption" = plot.caption,
                           "xlab" = xlab,
                           "ylab" = ylab,
                           "order.by" = order.by,
                           "position" = position,
                           "legend.title.position" = legend.title.position,
                           "font.type" = font.type)
    check_type(parameters = character_list, required_type = "character", test_function = is.character)

    # Checks.
    if (!(position %in% c("fill", "stack"))){stop("Position '", position, "' not supported. Please use either fill or stack.")}

    if (!(is.null(rotate_x_labels))){
      if(length(features) != length(rotate_x_labels)){
        stop('Total number of rotate_x_labels values does not match the number of features provided.', call. = F)
      }
    }

    # Check font.type.
    if (!(font.type %in% c("sans", "serif", "mono"))){
      stop("Please select one of the following for font.type: sans, serif, mono.", call. = F)
    }



    counter <- 0
    list.plots <- list()
    list.data <- list()
    if (is.null(colors.use)){reset_colors.use <- TRUE} else (reset_colors.use <- FALSE)
    for (feature in features){
      if (isTRUE(reset_colors.use)){colors.use <- NULL}
      counter <- counter + 1
      # Enforce the features to be part of the metadata.
      check_feature(sample = sample, features = feature, enforce_check = "metadata", enforce_parameter = "features")
      if (!is.null(rotate_x_labels)){x_label_select <- rotate_x_labels[counter]}
      # If no color scale is provided, generate a custom one.
      if (is.null(colors.use)){
        if (is.null(group.by)){
          # Generate a color palette equal to the number of identities in the seurat object.
          names.use <- sort(unique(sample@meta.data[, feature]))
          if (is.factor(names.use)){names.use <- levels(names.use)}
          colors.use <- generate_color_scale(names_use = names.use)
        } else if (!is.null(group.by)) {
          # Generate a color palette equal to the number of unique values in group.by variable.
          names.use <- sort(unique(sample@meta.data[, group.by]))
          if (is.factor(names.use)){names.use <- levels(names.use)}
          colors.use <- generate_color_scale(names_use = names.use)
        }
      } else{
        if (is.null(group.by)){
          colors.use <- check_consistency_colors_and_names(sample = sample, colors = colors.use, grouping_variable = feature)
        } else {
          colors.use <- check_consistency_colors_and_names(sample = sample, colors = colors.use, grouping_variable = group.by)
        }
      }


      if (is.null(group.by)){
        factor_levels <- compute_factor_levels(sample = sample, feature = feature, position = position)
        if (is.null(labels.order) & position == "fill"){factor_levels <- rev(factor_levels)}
        # Reorder the colors according to the factor levels.
        colors.use <- colors.use[factor_levels]
        if (isTRUE(horizontal)){factor_levels <- rev(factor_levels)}
        data <- sample@meta.data %>%
                dplyr::select(!!rlang::sym(feature)) %>%
                dplyr::group_by(!!rlang::sym(feature)) %>%
                dplyr::summarise(n = dplyr::n()) %>%
                dplyr::arrange(dplyr::desc(.data$n)) %>%
                dplyr::mutate(x_values = as.factor(!!(rlang::sym(feature)))) %>%
                dplyr::mutate(x_values = factor(.data$x_values, levels = factor_levels))
        data.out <- data %>%
                    dplyr::select(!!(rlang::sym(feature)), .data$n)
        data.out.wide <- data.out %>%
                         tidyr::pivot_wider(names_from = !!(rlang::sym(feature)),
                                            values_from = .data$n)
        data.report <- list("long" = data.out,
                            "wide" = data.out.wide)
        list.data[[feature]] <- data.report

        p <- data %>%
             ggplot2::ggplot(mapping = ggplot2::aes(x = .data$x_values, y = .data$n, fill = .data$x_values)) +
             ggplot2::geom_bar(position = position, stat="identity", width = 1,
                               colour="black",
                               size = 1) +
             ggplot2::scale_fill_manual(values = colors.use, na.value = "grey75") +
             ggplot2::guides(fill = ggplot2::guide_legend(ncol = legend.ncol,
                                                          nrow = legend.nrow,
                                                          byrow = legend.byrow,
                                                          title.position = legend.title.position))
      } else {
        check_feature(sample = sample, features = group.by, enforce_check = "metadata", enforce_parameter = "group.by")
        # Check the order of labels.
        if (!(is.null(labels.order))){
          factor_levels <- labels.order
        } else {
          if (!is.null(order.by)){
            if (!(order.by %in% unique(sample@meta.data[, group.by]))){
              stop("Parameter order.by (", order.by, ") not present in the unique values of parameter group.by (", group.by, ").", call. = F)
            }
            factor_levels <- compute_factor_levels(sample = sample, feature = feature, group.by = group.by, order.by = order.by, position = position)
            if (is.null(labels.order) & position == "stack"){factor_levels <- rev(factor_levels)}
          } else {
            factor_levels <- compute_factor_levels(sample = sample, feature = feature, group.by = group.by, position = position)
          }
        }

        if (is.null(labels.order) & position == "fill"){factor_levels <- rev(factor_levels)}
        if (isTRUE(horizontal)){factor_levels <- rev(factor_levels)}
        # Reorder the colors according to the factor levels.
        #colors.use <- colors.use[factor_levels]

        data <- sample@meta.data %>%
                dplyr::select(!!rlang::sym(feature), !!rlang::sym(group.by)) %>%
                dplyr::group_by(!!rlang::sym(group.by), !!rlang::sym(feature)) %>%
                dplyr::summarise(n = dplyr::n()) %>%
                dplyr::arrange(dplyr::desc(.data$n)) %>%
                dplyr::mutate(x_values = as.factor(!!(rlang::sym(feature)))) %>%
                dplyr::mutate(x_values = factor(.data$x_values, levels = factor_levels))
        data.out <- data %>%
                    dplyr::select(!!(rlang::sym(feature)), !!(rlang::sym(group.by)), .data$n)
        data.out.wide <- data.out %>%
                         tidyr::pivot_wider(names_from = !!(rlang::sym(feature)),
                                            values_from = .data$n)
        data.report <- list("long" = data.out,
                            "wide" = data.out.wide)
        list.data[[feature]] <- data.report
        p <- data %>%
             ggplot2::ggplot(mapping = ggplot2::aes(x = .data$x_values, y = .data$n, fill = !!rlang::sym(group.by))) +
             ggplot2::geom_bar(position = position, stat="identity", width = 1,
                               colour="black",
                               size = 1) +
             ggplot2::scale_fill_manual(values = colors.use, na.value = "grey75") +
             ggplot2::guides(fill = ggplot2::guide_legend(ncol = legend.ncol,
                                                          nrow = legend.nrow,
                                                          byrow = legend.byrow,
                                                          title.position = legend.title.position))
      }
      # Add theme.
      p <- p &
        ggplot2::theme_minimal(base_size = font.size) &
        ggplot2::theme(axis.title = ggplot2::element_text(face = "bold"),
                       axis.text = ggplot2::element_text(face = "bold", color = "black"),
                       plot.title = ggtext::element_markdown(face = "bold", hjust = 0),
                       plot.subtitle = ggtext::element_markdown(hjust = 0),
                       plot.caption = ggtext::element_markdown(hjust = 1),
                       plot.title.position = "plot",
                       panel.grid = ggplot2::element_blank(),
                       text = ggplot2::element_text(family = font.type),
                       plot.caption.position = "plot",
                       legend.text = ggplot2::element_text(face = "bold"),
                       legend.position = legend.position,
                       legend.title = ggplot2::element_text(face = "bold"),
                       legend.justification = "center",
                       plot.margin = ggplot2::margin(t = 10, r = 40, b = 10, l = 10),
                       axis.ticks = ggplot2::element_line(color = "black"),
                       axis.line = ggplot2::element_line(color = "black"),
                       plot.background = ggplot2::element_rect(fill = "white", color = "white"),
                       panel.background = ggplot2::element_rect(fill = "white", color = "white"),
                       legend.background = ggplot2::element_rect(fill = "white", color = "white"))

      # Whether to flip the axis or not.
      if (isTRUE(horizontal)){
        p <- p &
             ggplot2::coord_flip()
        if (position == "stack" & isTRUE(plot_line_guides)){
          p <- p &
               ggplot2::theme(panel.grid.major.x = ggplot2::element_line(color = "grey75", linetype = "dashed"))

        }
      } else if (isFALSE(horizontal)){
        if (position == "stack" & isTRUE(plot_line_guides)){
          p <- p &
               ggplot2::theme(panel.grid.major.y = ggplot2::element_line(color = "grey75", linetype = "dashed"))

        }
      }


      # Add labels on top of bars.
      if (isTRUE(add.subgroup_labels) | isTRUE(add.summary_labels)){
        func_use <- ifelse(isTRUE(repel.subgroup_labels), ggrepel::geom_label_repel, ggplot2::geom_label)
        func_use_summary <- ifelse(isTRUE(repel.summary_labels), ggrepel::geom_label_repel, ggplot2::geom_label)
        if (is.null(group.by)){
          if (position == "stack"){
            if (isTRUE(add.summary_labels)){
              p <- p +
                func_use_summary(data = data,
                                 mapping = ggplot2::aes(label = .data$n,
                                                        group = .data$x_values),
                                 label.size = 1,
                                 size = size.labels,
                                 color = "black",
                                 fill = "white",
                                 vjust = 0.5,
                                 hjust = 0.5,
                                 position = ggplot2::position_stack(vjust = 1,
                                                                    reverse = ifelse(horizontal == TRUE, TRUE, FALSE)),
                                 fontface = "bold",
                                 show.legend = FALSE)
              suppressMessages(p <- p + ggplot2::scale_color_manual(values = colors.use))

            }
            if (isTRUE(add.subgroup_labels)){
              p <- p +
                   func_use_summary(data = data,
                                    mapping = ggplot2::aes(label = .data$n,
                                                           group = .data$x_values,
                                                           color = .data$x_values),
                                    label.size = 1,
                                    size = size.labels,
                                    fill = "white",
                                    vjust = 0.5,
                                    hjust = 0.5,
                                    position = ggplot2::position_stack(vjust = 1,
                                                                       reverse = ifelse(horizontal == TRUE, TRUE, FALSE)),
                                    fontface = "bold",
                                    show.legend = FALSE)
              suppressMessages(p <- p + ggplot2::scale_color_manual(values = colors.use))
            }

          } else if (position == "fill"){
            warning("add.summary_labels and add.subgroup_labels are not yet implemented with position = fill.", call. = F)
          }
        } else if (!(is.null(group.by))){
          if (position == "stack"){
            # Get the values plotted on top.
            totals <- sample@meta.data %>%
              dplyr::select(!!rlang::sym(feature)) %>%
              dplyr::group_by(!!rlang::sym(feature)) %>%
              dplyr::summarise(n = dplyr::n()) %>%
              dplyr::arrange(dplyr::desc(.data$n)) %>%
              dplyr::mutate(x_values = as.factor(!!(rlang::sym(feature)))) %>%
              dplyr::mutate(x_values = factor(.data$x_values, levels = factor_levels))

            if (isTRUE(add.subgroup_labels)){
              p <- p +
                   func_use(data = data,
                            mapping = ggplot2::aes(label = .data$n,
                                                   group = !!rlang::sym(group.by),
                                                   color = !!rlang::sym(group.by)),
                            fill = "white",
                            fontface = "bold",
                            label.size = 1,
                            size = size.labels,
                            position = ggplot2::position_stack(vjust = 0.5),
                            show.legend = FALSE)
              suppressMessages(p <- p + ggplot2::scale_color_manual(values = colors.use))
            }
            if (isTRUE(add.summary_labels)){
              p <- p +
                   func_use_summary(data = totals,
                                    mapping = ggplot2::aes(label = .data$n,
                                                           group = .data$x_values),
                                    label.size = 1,
                                    size = size.labels,
                                    color = "black",
                                    fill = "white",
                                    vjust = 0.5,
                                    hjust = 0.5,
                                    fontface = "bold",
                                    show.legend = FALSE)
            }
          } else if (position == "fill"){
            warning("labels.use is not yet implemented with position = fill.", call. = F)
          }
        }
      }

      if (!(is.null(rotate_x_labels))){
        if (isTRUE(x_label_select)){
          p <- p & ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1))
        }
      }
      # Add X axis label.
      if (!is.null(xlab)){
        p <- p & ggplot2::xlab(xlab)
      } else {
        p <- p & ggplot2::xlab("")
      }
      # Add Y axis label.
      if (!is.null(ylab)){
        p <- p & ggplot2::ylab(ylab)
      } else {
        p <- p & ggplot2::ylab("")
      }
      # Add plot title.
      if (!is.null(plot.title)){
        p <- p + ggplot2::labs(title = plot.title)
      }

      # Add plot subtitle.
      if (!is.null(plot.subtitle)){
        p <- p + ggplot2::labs(subtitle = plot.subtitle)
      }

      # Add plot caption
      if (!is.null(plot.caption)){
        p <- p + ggplot2::labs(caption = plot.caption)
      }

      # Whether to plot the legend.
      if (isFALSE(legend)){
        p <- p +
             ggplot2::theme(legend.position = "none")
      }
      # Whether to remove legend title.
      if (isFALSE(legend.title)){
        p <- p +
             ggplot2::theme(legend.title = ggplot2::element_blank())
      }
    list.plots[[counter]] <- p
    }
    # Return the plot.
    if (length(features) > 1){
      p <- patchwork::wrap_plots(list.plots)
    } else {
      p <- list.plots[[1]]
    }

    # Return also the data?
    if (isTRUE(return_data_matrix)){
      return(list("plot" = p,
                  "data" = list.data))
    } else {
      return(p)
    }
}
