#' Plot Pathway Activities from DecoupleR using Progeny prior knowledge.
#'
#' @param sample Seurat object.
#' @param activities Result of running decoupleR method with progeny regulon prior knowledge.
#' @param group.by Character. Variable to group the output by. Defaults to current identities.
#' @param split.by Character. Optional variable to split the heatmaps by as well.
#' @param plot_FeaturePlots Logical. Compute output feature plots for each of the top regulons.
#' @param plot_Heatmaps Logical. Compute output heatmap showcasing the average TF activity per regulon and group.by variable.
#' @param plot_DotPlots Logical. Compute output dotplot for each of the top regulons and group.by variable.
#' @param row_title,column_title Character. Title for the rows in the heatmap.
#' @param transpose Logical. Whether to transpose the heatmap matrix.
#' @param cluster_cols,cluster_rows Logical. Whether to cluster rows and columns in the heatmap.
#' @param row_names_rot,column_names_rot Numeric. Angle of rotation of the column and row names. Suggested: either 0 or 90.
#' @param cell_size Numeric. Size of each cell in the heatmap.
#' @param pt.size Numeric. Size of the dots in the non-heatmap plots.
#' @param plot_cell_borders Logical. Whether to plot borders around cells in non-heatmap plots.
#' @param border.size Numeric. Size of the border.
#' @param na.value Character. Color for NA values.
#' @param legend.position Character. Position of the legend in the plots. One of: top, bottom, left, right.
#' @param heatmap.legend.length,heatmap.legend.width Numeric. Width and length of the legend in the heatmap.
#' @param heatmap.legend.framecolor Character. Color of the edges and ticks of the legend in the heatmap.
#' @param legend.framewidth,legend.tickwidth Width of the lines of the box in the legend for both feature plots and dot plots.
#' @param legend.framecolor,legend.tickcolor Color of the lines of the box in the legend for both feature plots and dot plots.
#' @param legend.length,legend.width Length and width of the legend for both feature plots and dot plots. Will adjust automatically depending on legend side.
#' @param legend.type Character. Type of legend to display. One of: normal, colorbar, colorsteps.
#' @param font.size Numeric. Base font.size of the figure.
#' @param font.type Character. Base font for the plot. One of mono, serif or sans.
#' @param rotate_x_axis_labels Logical. Whether to rotate the X axis names in the dot plot.
#'
#' @return A list containing several output plots according to the user's specifications.
#' @export
#'
#' @examples
#' \dontrun{
#' TBD
#' }
do_PathwayActivityPlot <- function(sample,
                                   activities,
                                   group.by = NULL,
                                   split.by = NULL,
                                   plot_FeaturePlots = FALSE,
                                   plot_Heatmaps = TRUE,
                                   plot_DotPlots = FALSE,
                                   row_title = NULL,
                                   column_title = NULL,
                                   transpose = FALSE,
                                   cluster_cols = TRUE,
                                   cluster_rows = TRUE,
                                   row_names_rot = 0,
                                   column_names_rot = 90,
                                   cell_size = 5,
                                   pt.size = 1,
                                   plot_cell_borders = TRUE,
                                   border.size = 1.5,
                                   na.value = "grey75",
                                   legend.position = "bottom",
                                   heatmap.legend.length = 75,
                                   heatmap.legend.width = 5,
                                   heatmap.legend.framecolor = "black",
                                   legend.width = 1,
                                   legend.length = 20,
                                   legend.framewidth = 1.5,
                                   legend.tickwidth = 1.5,
                                   legend.framecolor = "grey50",
                                   legend.tickcolor = "white",
                                   legend.type = "colorbar",
                                   font.size = 14,
                                   font.type = "sans",
                                   rotate_x_axis_labels = FALSE){

  #Checks for packages.
  #check_suggests(function_name = "do_FeaturePlot")
  # Check if the sample provided is a Seurat object.
  check_Seurat(sample = sample)

  # Check logical parameters.
  logical_list <- list("plot_FeaturePlots" = plot_FeaturePlots,
                       "plot_Heatmaps" = plot_Heatmaps,
                       "transpose" = transpose,
                       "cluster_cols" = cluster_cols,
                       "cluster_rows" = cluster_rows,
                       "rotate_x_axis_labels" = rotate_x_axis_labels,
                       "plot_cell_borders" = plot_cell_borders)
  check_type(parameters = logical_list, required_type = "logical", test_function = is.logical)
  # Check numeric parameters.
  numeric_list <- list("row_names_rot" = row_names_rot,
                       "column_names_rot" = column_names_rot,
                       "cell_size" = cell_size,
                       "heatmap.legend.length" = heatmap.legend.length,
                       "heatmap.legend.width" = heatmap.legend.width,
                       "pt.size" = pt.size,
                       "border.size" = border.size,
                       "font.size" = font.size,
                       "legend.width" = legend.width,
                       "legend.length" = legend.length,
                       "legend.framewidth" = legend.framewidth,
                       "legend.tickwidth" = legend.tickwidth)
  check_type(parameters = numeric_list, required_type = "numeric", test_function = is.numeric)
  # Check character parameters.
  character_list <- list("group.by" = group.by,
                         "split.by" = split.by,
                         "na.value" = na.value,
                         "legend.position" = legend.position,
                         "legend.framecolor" = legend.framecolor,
                         "font.type" = font.type,
                         "legend.tickcolor" = legend.tickcolor,
                         "legend.type" = legend.type,
                         "legend.framecolor" = legend.framecolor)
  check_type(parameters = character_list, required_type = "character", test_function = is.character)

  `%v%` <- ComplexHeatmap::`%v%`
  `%>%` <- purrr::`%>%`


  #sample <- readRDS("/b06x-isilon/b06x-g/G703/eblanco/projects/test_SC_datasets/sc_dataset.rds")
  #activities <- readRDS("/b06x-isilon/b06x-g/G703/eblanco/projects/test_SC_datasets/progeny_scores_decoupleR.rds")

  # Retrieve prior knowledge network.
  # network <- decoupleR::get_dorothea(organism = "human",
  #                                    levels = c("A", "B", "C"))
  #
  # # Run weighted means algorithm.
  # activities <- decoupleR::run_wmean(mat = as.matrix(sample@assays$SCT@data),
  #                                    network = network,
  #                                    .source = "source",
  #                                    .targe = "target",
  #                                    .mor = "mor",
  #                                    times = 100,
  #                                    minsize = 5)
  sample[["progeny"]] <- activities %>%
                         dplyr::filter(statistic == "norm_wmean") %>%
                         tidyr::pivot_wider(id_cols = "source",
                                            names_from = "condition",
                                            values_from = "score") %>%
                         tibble::column_to_rownames("source") %>%
                         Seurat::CreateAssayObject(.)

  Seurat::DefaultAssay(sample) <- "progeny"
  # Scale data.
  sample <- Seurat::ScaleData(object = sample, verbose = FALSE)

  # Start the process.
  if (is.null(group.by)){
    sample@meta.data[, "dummy"] <- sample@active.ident
  } else {
    sample@meta.data[, "dummy"] <- sample@meta.data[, group.by]
  }
  group.by <- "dummy"

  # Extract activities from object as a long dataframe
  suppressMessages({
    df <- t(as.matrix(sample@assays$progeny@scale.data)) %>%
          as.data.frame() %>%
          tibble::rownames_to_column(var = "cell") %>%
          dplyr::left_join(y = {sample@meta.data[, group.by, drop = FALSE] %>%
              tibble::rownames_to_column(var = "cell")},
              by = "cell") %>%
          dplyr::select(-"cell") %>%
          tidyr::pivot_longer(cols = -.data[[group.by]],
                              names_to = "source",
                              values_to = "score") %>%
          dplyr::group_by(.data[[group.by]], .data$source) %>%
          dplyr::summarise(mean = mean(.data$score))
  })

  list.out <- list()

  if (isTRUE(plot_FeaturePlots)){
    list.features <- list()
    for (pathway in rownames(sample)){
      limits <- c(min(sample@assays$progeny@scale.data[pathway, ]),
                  max(sample@assays$progeny@scale.data[pathway, ]))
      end_value <- max(abs(limits))

      p <- do_FeaturePlot(sample = sample,
                          features = pathway,
                          assay = "progeny",
                          slot = "scale.data",
                          pt.size = pt.size,
                          order = FALSE,
                          border.size = border.size,
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
                          legend.width = legend.width)
      p <- add_scale(p = p,
                     scale = "color",
                     function_use = ggplot2::scale_color_gradientn(colors = c("#033270", "#4091C9", "#fdf0d5", "#c94040", "#65010C"),
                                                                   limits = c(-end_value, end_value)))
      p$guides$colour$title <- paste0(pathway, " activity")
      list.features[[pathway]] <- p
    }
    list.out[["feature_plots"]] <- list.features
  }

  if (isTRUE(plot_DotPlots)){
    list.dotplots <- list()
    for (pathway in rownames(sample)){
      limits <- c(min(sample@assays$progeny@scale.data[pathway, ]),
                  max(sample@assays$progeny@scale.data[pathway, ]))
      end_value <- max(abs(limits))
      data <- sample@assays$progeny@scale.data[pathway, , drop = F] %>%
              t() %>%
              as.data.frame() %>%
              tibble::rownames_to_column(var = "cell") %>%
              dplyr::left_join(y = {sample@meta.data[, group.by, drop = F] %>%
                                    tibble::rownames_to_column(var = "cell")},
                                    by = "cell") %>%
              dplyr::select(c(-.data$cell)) %>%
              tidyr::pivot_longer(cols = -.data[[group.by]],
                                  names_to = "source",
                                  values_to = "score")
      p <- data %>%
           dplyr::mutate("group.by" = factor(.data[[group.by]], levels = {data %>%
                                                                          dplyr::group_by(.data[[group.by]]) %>%
                                                                          dplyr::summarise("mean" = mean(.data$score)) %>%
                                                                          dplyr::arrange(dplyr::desc(.data$mean)) %>%
                                                                          dplyr::pull(.data[[group.by]]) %>%
                                                                          as.character()})) %>%
           ggplot2::ggplot(mapping = ggplot2::aes(x = .data$group.by,
                                                  y = .data$score,
                                                  color = .data$score))
      if (isTRUE(plot_cell_borders)){
        p <- p +
             ggplot2::geom_point(position = ggplot2::position_jitter(width = 0.455,
                                                                     seed = 0),
                                 size = pt.size * border.size,
                                 color = "black",
                                 na.rm = TRUE)
      }
      p <- p +
           ggplot2::geom_point(position = ggplot2::position_jitter(width = 0.45,
                                                                   seed = 0),size = pt.size,
                               na.rm = TRUE) +
           ggdist::stat_pointinterval(position = ggplot2::position_dodge(width = 1),
                                      na.rm = TRUE,
                                      color = "black") +
           ggplot2::scale_color_gradientn(colors = c("#033270", "#4091C9", "#fdf0d5", "#c94040", "#65010C"),
                                          limits = c(-end_value, end_value)) +
           ggplot2::theme_minimal(base_size = font.size) +
           ggplot2::theme(axis.title = ggplot2::element_blank(),
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
                          legend.background = ggplot2::element_rect(fill = "white", color = "white"))

      p <- modify_continuous_legend(p = p,
                                    legend.title = paste0(pathway, " activity"),
                                    legend.aes = "color",
                                    legend.type = legend.type,
                                    legend.position = legend.position,
                                    legend.length = legend.length,
                                    legend.width = legend.width,
                                    legend.framecolor = legend.framecolor,
                                    legend.tickcolor = legend.tickcolor,
                                    legend.framewidth = legend.framewidth,
                                    legend.tickwidth = legend.tickwidth)
      list.dotplots[[pathway]] <- p
    }
    list.out[["dotplots"]] <- list.dotplots
  }


  if (isTRUE(plot_Heatmaps)){
    if (is.null(split.by)){
      suppressMessages({
        data <- sample@assays$progeny@scale.data %>%
                t() %>%
                as.data.frame() %>%
                tibble::rownames_to_column(var = "cell") %>%
                dplyr::left_join(y = {sample@meta.data[, group.by, drop = FALSE] %>%
                    tibble::rownames_to_column(var = "cell")},
                    by = "cell") %>%
                dplyr::select(c(-cell)) %>%
                tidyr::pivot_longer(cols = -.data[[group.by]],
                                    names_to = "source",
                                    values_to = "score") %>%
                dplyr::group_by(.data[[group.by]], .data$source) %>%
                dplyr::summarise(mean = mean(.data$score)) %>%
                tidyr::pivot_wider(id_cols = group.by,
                                   names_from = 'source',
                                   values_from = 'mean') %>%
                tibble::column_to_rownames(group.by) %>%
                as.matrix()
      })

      row_title <- {
        if (!(is.null(row_title))){
          row_title
        } else {
          ""
        }}

      column_title <- {
        if (!(is.null(column_title))){
          column_title
        } else {
          ""
        }}

      if (isTRUE(transpose)){
        data <- t(data)
      }

      out <- heatmap_inner(data,
                           legend_name = "Pathway activity",
                           column_title = column_title,
                           row_title = row_title,
                           cluster_columns = cluster_cols,
                           cluster_rows = cluster_rows,
                           column_names_rot = column_names_rot,
                           row_names_rot = row_names_rot,
                           cell_size = cell_size,
                           na.value = na.value,
                           legend.position = legend.position,
                           legend.length = heatmap.legend.length,
                           legend.width = heatmap.legend.width,
                           legend.framecolor = heatmap.legend.framecolor)
      h <- out[["heatmap"]]
      h_legend <- out[["legend"]]
      ComplexHeatmap::ht_opt("HEATMAP_LEGEND_PADDING" = ggplot2::unit(8, "mm"))
      suppressWarnings({
        grDevices::pdf(NULL)
        h <- ComplexHeatmap::draw(h,
                                  heatmap_legend_list = h_legend,
                                  heatmap_legend_side = if (legend.position %in% c("top", "bottom")){"bottom"} else {"right"},
                                  padding = ggplot2::unit(c(5, 5, 5, 5), "mm"))
        grDevices::dev.off()
      })
    } else {
      split.values <- as.character(sort(unique(sample@meta.data %>% dplyr::pull(!!rlang::sym(split.by)))))
      list.heatmaps <- list()
      # Get the maximum range.
      range <- max(abs(top_acts_mat_wide))
      for (split.value in split.values){
        suppressMessages({
          data <- sample@assays$progeny@scale.data %>%
                  t() %>%
                  as.data.frame() %>%
                  tibble::rownames_to_column(var = "cell") %>%
                  dplyr::left_join(y = {sample@meta.data[, c(group.by, split.by)] %>%
                      tibble::rownames_to_column(var = "cell")},
                      by = "cell") %>%
                  dplyr::filter(.data[[split.by]] == split.value) %>%  # This is key.
                  dplyr::select(c(-cell, -split.by)) %>%
                  tidyr::pivot_longer(cols = -.data[[group.by]],
                                      names_to = "source",
                                      values_to = "score") %>%
                  dplyr::group_by(.data[[group.by]], .data$source) %>%
                  dplyr::summarise(mean = mean(.data$score)) %>%
                  tidyr::pivot_wider(id_cols = group.by,
                                     names_from = 'source',
                                     values_from = 'mean') %>%
                  tibble::column_to_rownames(group.by) %>%
                  as.matrix()
        })

        row_title <- {
          if (!(is.null(row_title))){
            row_title
          } else {
            ""
          }}

        column_title <- {
          if (!(is.null(column_title))){
            column_title
          } else {
            ""
          }}

        if (isTRUE(transpose)){
          data <- t(data)
        }

        if (isTRUE(split.horizontal)){
          row_title <- ""
          column_title <- split.value
        } else {
          row_title <- split.value
          column_title <- ""
        }
        out <- heatmap_inner(data,
                             legend_name = "Pathway activity",
                             column_title = column_title,
                             row_title = row_title,
                             cluster_columns = cluster_cols,
                             cluster_rows = cluster_rows,
                             column_names_rot = column_names_rot,
                             row_names_rot = row_names_rot,
                             cell_size = cell_size,
                             na.value = na.value,
                             legend.position = legend.position,
                             legend.length = heatmap.legend.length,
                             legend.width = heatmap.legend.width,
                             legend.framecolor = heatmap.legend.framecolor)
        h <- out[["heatmap"]]
        h_legend <- out[["legend"]]
        list.heatmaps[[split.value]] <- h
      }
      ComplexHeatmap::ht_opt("HEATMAP_LEGEND_PADDING" = ggplot2::unit(8, "mm"))
      suppressWarnings({
        grDevices::pdf(NULL)
        ht_list <- NULL
        for (name in names(list.heatmaps)){
          ht_list = ht_list %v% list.heatmaps[[name]]
        }
        h <- ComplexHeatmap::draw(ht_list,
                                  heatmap_legend_list = h_legend,
                                  heatmap_legend_side = if (legend.position %in% c("top", "bottom")){"bottom"} else {"right"},
                                  padding = ggplot2::unit(c(5, 5, 5, 5), "mm"))
        grDevices::dev.off()
      })
    }
    list.out[["heatmaps"]] <- list("average_scores" = h)
  }
  return(list.out)
}
