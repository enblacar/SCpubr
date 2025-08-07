#' This function is a wrapper for \link[Seurat]{DotPlot}. It provides most of its functionalities while adding extra.
#' You can
#'
#' @inheritParams doc_function
#' @param cluster \strong{\code{\link[base]{logical}}} | Whether to cluster the identities based on the expression of the features.
#' @param zscore.data \strong{\code{\link[base]{logical}}} | Whether to compute Z-scores instead of showing average expression values. This allows to see, for each gene, which group has the highest average expression, but prevents you from comparing values across genes. Can not be used with slot = "scale.data" or with split.by. 
#' @param dot.min \strong{\code{\link[base]{numeric}}} | Ranges from 0 to 100. Filter out dots whose Percent Expressed falls below this threshold. 
#'
#' @return A ggplot2 object containing a Dot Plot.
#' @export
#'
#' @example man/examples/examples_do_DotPlot.R
do_DotPlot <- function(sample,
                       features,
                       assay = NULL,
                       slot = "data",
                       group.by = NULL,
                       split.by = NULL,
                       zscore.data = FALSE,
                       min.cutoff = NA,
                       max.cutoff = NA,
                       dot.min = 5,
                       enforce_symmetry = ifelse(base::isTRUE(zscore.data), TRUE, FALSE), 
                       legend.title = NULL,
                       legend.type = "colorbar",
                       legend.position = "bottom",
                       legend.framewidth = 0.5,
                       legend.tickwidth = 0.5,
                       legend.length = 20,
                       legend.width = 1,
                       legend.framecolor = "grey50",
                       legend.tickcolor = "white",
                       legend.ncol = NULL,
                       legend.nrow = NULL,
                       legend.byrow = FALSE,
                       dot.scale = 8,
                       plot.title = NULL,
                       plot.subtitle = NULL,
                       plot.caption = NULL,
                       xlab = NULL,
                       ylab = NULL,
                       font.size = 14,
                       font.type = "sans",
                       cluster = FALSE,
                       flip = FALSE,
                       axis.text.x.angle = 45,
                       use_viridis = FALSE,
                       viridis.palette = "G",
                       viridis.direction = -1,
                       sequential.palette = "YlGnBu",
                       sequential.direction = 1,
                       diverging.palette = "RdBu",
                       diverging.direction = -1,
                       na.value = "grey75",
                       plot.grid = TRUE,
                       grid.color = "grey75",
                       grid.type = "dashed",
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
  
    check_suggests(function_name = "do_DotPlot")
    check_Seurat(sample = sample)
    # Check the assay.
    out <- check_and_set_assay(sample, assay = assay)
    sample <- out[["sample"]]
    assay <- out[["assay"]]
    
    # Check group.by.
    out <- check_group_by(sample = sample,
                          group.by = group.by,
                          is.heatmap = FALSE)
    sample <- out[["sample"]]
    group.by <- out[["group.by"]]
    
    # Check slot.
    slot <- if(is.null(slot)){"data"} else {slot}
                               
    # Check logical parameters.
    logical_list <- list("flip" = flip,
                         "cluster" = cluster,
                         "use_viridis" = use_viridis,
                         "plot.grid" = plot.grid,
                         "enforce_symmetry" = enforce_symmetry,
                         "legend.byrow" = legend.byrow,
                         "zscore.data" = zscore.data)
    check_type(parameters = logical_list, required_type = "logical", test_function = is.logical)
    # Check numeric parameters.
    numeric_list <- list("dot.scale" = dot.scale,
                         "font.size" = font.size,
                         "legend.framewidth" = legend.framewidth,
                         "legend.tickwidth" = legend.tickwidth,
                         "legend.length" = legend.length,
                         "legend.width" = legend.width,
                         "viridis.direction" = viridis.direction,
                         "axis.text.x.angle" = axis.text.x.angle,
                         "number.breaks" = number.breaks,
                         "sequential.direction" = sequential.direction,
                         "diverging.direction" = diverging.direction,
                         "min.cutoff" = min.cutoff,
                         "max.cutoff" = max.cutoff,
                         "legend.ncol" = legend.ncol,
                         "legend.nrow" = legend.nrow,
                         "dot.min" = dot.min)
    check_type(parameters = numeric_list, required_type = "numeric", test_function = is.numeric)
    # Check character parameters.
    character_list <- list("legend.position" = legend.position,
                           "plot.title" = plot.title,
                           "features" = unlist(features),
                           "xlab" = xlab,
                           "ylab" = ylab,
                           "group.by" = group.by,
                           "split.by" = split.by,
                           "legend.framecolor" = legend.framecolor,
                           "legend.tickcolor" = legend.tickcolor,
                           "legend.type" = legend.type,
                           "font.type" = font.type,
                           "viridis.palette" = viridis.palette,
                           "grid.color" = grid.color,
                           "grid.type" = grid.type,
                           "sequential.palette" = sequential.palette,
                           "diverging.palette" = diverging.palette,
                           "plot.title.face" = plot.title.face,
                           "plot.subtitle.face" = plot.subtitle.face,
                           "plot.caption.face" = plot.caption.face,
                           "axis.title.face" = axis.title.face,
                           "axis.text.face" = axis.text.face,
                           "legend.title.face" = legend.title.face,
                           "legend.text.face" = legend.text.face,
                           "legend.title" = legend.title,
                           "slot" = slot)
    check_type(parameters = character_list, required_type = "character", test_function = is.character)
    
    `%>%` <- magrittr::`%>%`
    
    # Check the features.
    features <- check_feature(sample = sample, features = features, permissive = TRUE)
    features <- remove_duplicated_features(features = features)

    check_parameters(parameter = font.type, parameter_name = "font.type")
    check_parameters(parameter = legend.type, parameter_name = "legend.type")
    check_parameters(parameter = legend.position, parameter_name = "legend.position")
    check_parameters(parameter = viridis.palette, parameter_name = "viridis.palette")
    check_parameters(parameter = sequential.palette, parameter_name = "sequential.palette")
    check_parameters(parameter = grid.type, parameter_name = "grid.type")
    check_parameters(parameter = axis.text.x.angle, parameter_name = "axis.text.x.angle")
    check_parameters(parameter = number.breaks, parameter_name = "number.breaks")
    check_parameters(plot.title.face, parameter_name = "plot.title.face")
    check_parameters(plot.subtitle.face, parameter_name = "plot.subtitle.face")
    check_parameters(plot.caption.face, parameter_name = "plot.caption.face")
    check_parameters(axis.title.face, parameter_name = "axis.title.face")
    check_parameters(axis.text.face, parameter_name = "axis.text.face")
    check_parameters(legend.title.face, parameter_name = "legend.title.face")
    check_parameters(legend.text.face, parameter_name = "legend.text.face")
    check_parameters(viridis.direction, parameter_name = "viridis.direction")
    check_parameters(sequential.direction, parameter_name = "sequential.direction")

    check_colors(legend.framecolor, parameter_name = "legend.framecolor")
    check_colors(legend.tickcolor, parameter_name = "legend.tickcolor")
    check_colors(na.value, parameter_name = "na.value")
    check_colors(grid.color, parameter_name = "grid.color")
    
    
    if (base::isTRUE(zscore.data)){
      assertthat::assert_that(base::isTRUE(enforce_symmetry),
                              msg = paste0(add_cross(), crayon_body("Please set "),
                                           crayon_key("enforce_symmetry"),
                                           crayon_body(" to "),
                                           crayon_key("TRUE"),
                                           crayon_body(" when scaling the data. This allows for a "),
                                           crayon_key("centered"),
                                           crayon_body(" color scale around "),
                                           crayon_key("0"),
                                           crayon_body(".")))
      
      assertthat::assert_that(slot == "data",
                              msg = paste0(add_cross(), crayon_body("Please set "),
                                           crayon_key("slot"),
                                           crayon_body(" to "),
                                           crayon_key('"data"'),
                                           crayon_body(" when scaling the data. Performing Z-scaling over "),
                                           crayon_key("already scaled"),
                                           crayon_body(" data is "),
                                           crayon_key("not advisable"),
                                           crayon_body(".")))
      
      colors.gradient <- compute_continuous_palette(name = diverging.palette,
                                                    use_viridis = FALSE,
                                                    direction = diverging.direction,
                                                    enforce_symmetry = TRUE)
      
      
    } else {
      
      colors.gradient <- compute_continuous_palette(name = ifelse(isTRUE(use_viridis), viridis.palette, sequential.palette),
                                                    use_viridis = use_viridis,
                                                    direction = ifelse(isTRUE(use_viridis), viridis.direction, sequential.direction),
                                                    enforce_symmetry = enforce_symmetry)
      
      center_on_value <- FALSE
      value_center <- NULL
    }
    
    if (is.list(features)){
      assertthat::assert_that(!is.null(names(features)),
                              msg = paste0(add_cross(), crayon_body("Please provide features as a "),
                                           crayon_key("named list"),
                                           crayon_body(" and not as a "),
                                           crayon_key("regular list"),
                                           crayon_body(".")))
      
      assertthat::assert_that(is.null(split.by),
                              msg = paste0(add_cross(), crayon_body("Please either provide features as a "),
                                           crayon_key("named list"),
                                           crayon_body(" or set up "),
                                           crayon_key("split by"),
                                           crayon_body(". A combination of both is not allowed.")))
    }
    
    if (!is.null(split.by)){
      assertthat::assert_that(base::isFALSE(cluster),
                              msg = paste0(add_cross(), crayon_body("Please when using "),
                                           crayon_key("split.by"),
                                           crayon_body(" set "),
                                           crayon_key("cluster"),
                                           crayon_body(" to "),
                                           crayon_key("FALSE"),
                                           crayon_body(".")))
      
      assertthat::assert_that(base::isFALSE(zscore.data),
                              msg = paste0(add_cross(), crayon_body("Please when using "),
                                           crayon_key("split.by"),
                                           crayon_body(" set "),
                                           crayon_key("zscore.data"),
                                           crayon_body(" to "),
                                           crayon_key("FALSE"),
                                           crayon_body(".")))
    }
    
    # Workaround parameter depreciation.
    # nocov start
    if (base::isTRUE(utils::packageVersion("Seurat") < "4.9.9")){
      data <- Seurat::GetAssayData(object = sample,
                                   assay = assay,
                                   slot = slot)
    } else {
      data <- SeuratObject::LayerData(object = sample,
                                      assay = assay,
                                      layer = slot)
    }
    
    # nocov end
    
    # Select features.
    if (is.list(features)){
      genes.unique <- unique(unlist(features))
      
      # Remove duplicates across lists.
      features.use <- list()
      start <- 1
      
      # Add back in the original order
      for (name in names(features)) {
        len <- length(features[[name]])
        features.use[[name]] <- genes.unique[start:(start + len - 1)]
        start <- start + len
      }
      
      # Get the length of the longest list.
      max_len <- max(lengths(features.use))
      # Get a padded list
      df.map <- as.data.frame(lapply(features.use, function(x) c(x, rep(NA, max_len - length(x)))), check.names = FALSE) %>% 
                tidyr::pivot_longer(cols = dplyr::everything(),
                                    names_to = "Name",
                                    values_to = "Gene") %>% 
                dplyr::mutate("Name" = factor(.data$Name, levels = names(features)))
      
      
      features.use <- unique(unlist(unname(features)))[!duplicated(unique(unlist(unname(features))))]
    } else {
      features.use <- features
    }
    
    selection <- c(split.by, group.by, "Gene", "Avg.Exp", "P.Exp")
    data <- data[features.use, , drop = FALSE] %>% 
            as.data.frame() %>% 
            tibble::rownames_to_column(var = "Gene") %>% 
            tidyr::pivot_longer(cols = -"Gene",
                                values_to = "Expression",
                                names_to = "Cell") %>% 
            dplyr::left_join(y = {sample@meta.data[ , c(group.by, split.by), drop = FALSE] %>% 
                                  tibble::rownames_to_column(var = "Cell")},
                             by = "Cell") %>% 
            dplyr::mutate("logical" = ifelse(.data$Expression == 0, 0, 1)) %>% 
            dplyr::group_by(dplyr::across(dplyr::all_of(c(split.by, group.by, "Gene")))) %>% 
            dplyr::summarise("Avg.Exp" = mean(.data$Expression, na.rm = TRUE),
                             "N.Exp" = sum(.data$logical),
                             "N" = dplyr::n(),
                             .groups = "drop") %>% 
            dplyr::mutate("P.Exp" = (.data$N.Exp / .data$N) * 100) %>% 
            dplyr::select(dplyr::all_of(selection))
    
    if (is.null(split.by)){
      data <- data %>% tidyr::complete(.data[[group.by]], .data$Gene, fill = list("Avg.Exp" = 0, "P.Exp" = 0))
    } else {
      data <- data %>% tidyr::complete(.data[[split.by]], .data[[group.by]], .data$Gene, fill = list("Avg.Exp" = 0, "P.Exp" = 0))
    }
            
    
    if (base::isTRUE(zscore.data)){
      selection <- c(group.by, "Gene", "Avg.Exp")
      data <- data %>% 
              dplyr::select(dplyr::all_of(selection)) %>% 
              tidyr::pivot_wider(names_from = dplyr::all_of(group.by),
                                 values_from = "Avg.Exp") %>% 
              as.data.frame() %>% 
              tibble::column_to_rownames(var = "Gene") %>% 
              t() %>% 
              scale(center = TRUE, scale = TRUE) %>% 
              t() %>% 
              as.data.frame() %>% 
              tibble::rownames_to_column(var = "Gene") %>% 
              tidyr::pivot_longer(-dplyr::all_of("Gene"),
                                  names_to = group.by,
                                  values_to = "Avg.Exp") %>% 
              dplyr::left_join(y = data %>% dplyr::select(-dplyr::all_of("Avg.Exp")),
                               by = c(group.by, "Gene"))
      
      
    }
    
    # Add gene map.
    if (is.list(features)){
      data <- data %>% dplyr::left_join(y = df.map, by = "Gene")
      
    }
    
    # Define cutoffs.
    range.data <- c(min(data[, "Avg.Exp"], na.rm = TRUE),
                    max(data[, "Avg.Exp"], na.rm = TRUE))
    
    out <- check_cutoffs(min.cutoff = min.cutoff,
                         max.cutoff = max.cutoff,
                         limits = range.data)
    range.data <- out$limits
    
    
    scale.setup <- compute_scales(sample = sample,
                                  feature = NULL,
                                  assay = assay,
                                  reduction = NULL,
                                  slot = slot,
                                  number.breaks = number.breaks,
                                  min.cutoff = min.cutoff,
                                  max.cutoff = max.cutoff,
                                  flavor = "Seurat",
                                  enforce_symmetry = enforce_symmetry,
                                  from_data = TRUE,
                                  limits.use = range.data)
    
    # Modify values
    if (!is.na(min.cutoff)){
      data$Avg.Exp <- ifelse(data$Avg.Exp <= min.cutoff, min.cutoff, data$Avg.Exp)
    }
    
    if (!is.na(max.cutoff)){
      data$Avg.Exp <- ifelse(data$Avg.Exp >= max.cutoff, max.cutoff, data$Avg.Exp)
    }
    
    
    selection <- c("Groups", "Gene", "Avg.Exp")
    data.cluster <- data %>% 
                    dplyr::ungroup() %>% 
                    dplyr::mutate("Groups" = .data[[group.by]]) %>% 
                    dplyr::select(dplyr::all_of(selection)) %>% 
                    tidyr::pivot_wider(names_from = "Groups",
                                       values_from = "Avg.Exp",
                                       values_fn = list) %>% 
                    as.data.frame() %>% 
                    tibble::column_to_rownames(var = "Gene") %>% 
                    as.matrix()
    

    # Set NAs to 0.
    data.cluster[is.na(data.cluster)] <- 0
    
    # Cluster rows.
    if(length(rownames(data.cluster)) == 1){
      row_order <- rownames(data.cluster)[1]
    } else {
      if (isTRUE(cluster)){
        row_order <- rownames(data.cluster)[stats::hclust(stats::dist(data.cluster, method = "euclidean"), method = "ward.D")$order]
      } else {
        row_order <- features.use
      }
    }
    
    # Cluster columns.
    if (length(colnames(data.cluster)) == 1){
      col_order <- colnames(data.cluster)[1]
    } else {
      if (isTRUE(cluster)){
        col_order <- colnames(data.cluster)[stats::hclust(stats::dist(t(data.cluster), method = "euclidean"), method = "ward.D")$order]
      } else {
        if (is.factor(sample@meta.data[, group.by])){
          col_order <- levels(sample@meta.data[, group.by])
        } else {
          col_order <- sort(unique(sample@meta.data[, group.by])) 
        }
      }
    }
    
    # Apply clustering.
    data <- data %>% 
            dplyr::ungroup() %>% 
            dplyr::mutate("Groups" = .data[[group.by]]) %>% 
            dplyr::mutate("Gene" = factor(.data$Gene, levels = row_order),
                          "Groups" = factor(.data$Groups, levels = rev(col_order)))
    
    
    # Define legend title
    if (is.null(legend.title)){
      legend.title <- ifelse(base::isTRUE(zscore.data), "Z-Scored | Avg. Exp.", "Avg. Exp.")
    }
    
    p <- data %>%
         ggplot2::ggplot(mapping = ggplot2::aes(x = if (base::isFALSE(flip)){.data$Gene} else {.data$Groups},
                                                y = if (base::isFALSE(flip)){.data$Groups} else {.data$Gene},
                                                fill = .data$Avg.Exp,
                                                size = .data$P.Exp)) + 
         ggplot2::geom_point(color = "black", shape = 21) +
         ggplot2::scale_size_continuous(range = c(0, dot.scale)) +
         ggplot2::scale_fill_gradientn(colors = colors.gradient,
                                       na.value = na.value,
                                       name = legend.title,
                                       breaks = scale.setup$breaks,
                                       labels = scale.setup$labels,
                                       limits = scale.setup$limits)
     
    # Facet grid.
    if (!is.null(split.by)){
      if (base::isFALSE(flip)){
        p <- p +
             ggplot2::facet_grid(cols = ggplot2::vars(.data[[split.by]]),
                                 scales = "free",
                                 space = "free")
      } else {
        p <- p +
             ggplot2::facet_grid(rows = ggplot2::vars(.data[[split.by]]),
                                 scales = "free",
                                 space = "free")
      }
    } else {
      if (is.list(features)){
        if (base::isFALSE(flip)){
          p <- p +
               ggplot2::facet_grid(cols = ggplot2::vars(.data$Name),
                                   scales = "free",
                                   space = "free")
        } else {
          p <- p +
               ggplot2::facet_grid(rows = ggplot2::vars(.data$Name),
                                   scales = "free",
                                   space = "free")
        }
      }
    }
     
    p <- p +
         ggplot2::xlab(xlab) +
         ggplot2::ylab(ylab) +
         ggplot2::labs(title = plot.title,
                       subtitle = plot.subtitle,
                       caption = plot.caption) +
         ggplot2::theme_minimal(base_size = font.size) +
         ggplot2::theme(axis.text.x = ggplot2::element_text(color = "black",
                                                            face = axis.text.face,
                                                            angle = get_axis_parameters(angle = axis.text.x.angle, flip = flip)[["angle"]],
                                                            hjust = get_axis_parameters(angle = axis.text.x.angle, flip = flip)[["hjust"]],
                                                            vjust = get_axis_parameters(angle = axis.text.x.angle, flip = flip)[["vjust"]]),
                        axis.text.y = ggplot2::element_text(face = axis.text.face, color = "black"),
                        axis.ticks = ggplot2::element_line(color = "black"),
                        axis.line = ggplot2::element_line(color = "black"),
                        axis.title = ggplot2::element_text(face = axis.title.face),
                        plot.title = ggplot2::element_text(face = plot.title.face, hjust = 0),
                        plot.subtitle = ggplot2::element_text(face = plot.subtitle.face, hjust = 0),
                        plot.caption = ggplot2::element_text(face = plot.caption.face, hjust = 1),
                        plot.title.position = "plot",
                        panel.grid = if (isTRUE(plot.grid)){ggplot2::element_line(color = grid.color, linetype = grid.type)} else {ggplot2::element_blank()},
                        text = ggplot2::element_text(family = font.type),
                        plot.caption.position = "plot",
                        legend.text = ggplot2::element_text(face = legend.text.face),
                        legend.position = legend.position,
                        legend.title = ggplot2::element_text(face = legend.title.face),
                        legend.justification = "center",
                        plot.margin = ggplot2::margin(t = 10, r = 10, b = 10, l = 10),
                        panel.grid.major = ggplot2::element_blank(),
                        plot.background = ggplot2::element_rect(fill = "white", color = "white"),
                        panel.background = ggplot2::element_rect(fill = "white", color = "white"),
                        legend.background = ggplot2::element_rect(fill = "white", color = "white"))
    
    # Add leyend modifiers.
    p <- modify_continuous_legend(p = p,
                                  # nocov start
                                  legend.title = if (is.null(legend.title)){"Avg. Expression"} else {legend.title},
                                  # nocov end
                                  legend.aes = "fill",
                                  legend.type = legend.type,
                                  legend.position = legend.position,
                                  legend.length = legend.length,
                                  legend.width = legend.width,
                                  legend.framecolor = legend.framecolor,
                                  legend.tickcolor = legend.tickcolor,
                                  legend.framewidth = legend.framewidth,
                                  legend.tickwidth = legend.tickwidth)

    # Modify size legend.
    p <- p +
         ggplot2::guides(size = ggplot2::guide_legend(title = "Percent Expressed",
                                                      title.position = "top",
                                                      title.hjust = 0.5,
                                                      ncol = legend.ncol,
                                                      nrow = legend.nrow,
                                                      byrow = legend.byrow,
                                                      override.aes = ggplot2::aes(fill = "black")))
    
    # Filter out dots with low percent expressed.
    p$data <- p$data %>% dplyr::filter(.data$P.Exp >= dot.min)

    return(p)
}
