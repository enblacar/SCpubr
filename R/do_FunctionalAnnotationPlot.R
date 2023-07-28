#' Compute functional annotation plots using GO or KEGG ontologies
#'
#'
#' @inheritParams doc_function
#' @param organism  \strong{\code{\link[base]{character}}} | Supported KEGG organism.
#' @param database \strong{\code{\link[base]{character}}} | Database to run the analysis on. One of:
#' \itemize{
#'   \item \emph{\code{GO}}.
#'   \item \emph{\code{KEGG}}.
#' }
#' @param p.adjust.cutoff \strong{\code{\link[base]{numeric}}} | Significance cutoff used to filter non-significant terms.
#' @param pAdjustMethod \strong{\code{\link[base]{character}}} | Method to adjust for multiple testing.  One of:
#' \itemize{
#'   \item \emph{\code{holm}}.
#'   \item \emph{\code{hochberg}}.
#'   \item \emph{\code{hommel}}.
#'   \item \emph{\code{bonferroni}}.
#'   \item \emph{\code{BH}}.
#'   \item \emph{\code{BY}}.
#'   \item \emph{\code{fdr}}.
#'   \item \emph{\code{none}}.
#' }
#' @param minGSSize \strong{\code{\link[base]{numeric}}} | Minimal size of genes annotated by Ontology term for testing.
#' @param maxGSSize \strong{\code{\link[base]{numeric}}} | Maximal size of genes annotated for testing.
#' @param return_matrix \strong{\code{\link[base]{logical}}} | Returns the matrices with the enriched Terms for further use.
#' @return A list containing a heatmap of the presence/absence of the genes in the enriched term, as well as a bar plot, dot plot and tree plot of the enriched terms.
#' @export
#'
#' @example /man/examples/examples_do_FunctionalAnnotationPlot.R
do_FunctionalAnnotationPlot <- function(genes,
                                        org.db,
                                        organism = "hsa",
                                        database = "GO",
                                        GO_ontology = "BP",
                                        min.overlap = NULL,
                                        p.adjust.cutoff = 0.05,
                                        pAdjustMethod = "BH",
                                        minGSSize = 10,
                                        maxGSSize = 500,
                                        font.size = 10,
                                        font.type = "sans",
                                        axis.text.x.angle = 45,
                                        xlab = NULL,
                                        ylab = NULL,
                                        plot.title = NULL,
                                        plot.subtitle = NULL,
                                        plot.caption = NULL,
                                        legend.type = "colorbar",
                                        legend.position = "bottom",
                                        legend.framewidth = 0.5,
                                        legend.tickwidth = 0.5,
                                        legend.length = 10,
                                        legend.width = 1,
                                        legend.framecolor = "grey50",
                                        legend.tickcolor = "white",
                                        number.breaks = 5,
                                        return_matrix = FALSE,
                                        grid.color = "white",
                                        border.color = "black",
                                        plot.title.face = "bold",
                                        plot.subtitle.face = "plain",
                                        plot.caption.face = "italic",
                                        axis.title.face = "bold",
                                        axis.text.face = "plain",
                                        legend.title.face = "bold",
                                        legend.text.face = "plain"){
  # Add lengthy error messages.
  withr::local_options(.new = list("warning.length" = 8170))
  
  

  check_suggests(function_name = "do_FunctionalAnnotationPlot")
  `%>%` <- magrittr::`%>%`
  # Check logical parameters.
  logical_list <- list("return_matrix" = return_matrix)
  check_type(parameters = logical_list, required_type = "logical", test_function = is.logical)
  # Check numeric parameters.
  numeric_list <- list("font.size" = font.size,
                       "legend.framewidth" = legend.framewidth,
                       "legend.tickwidth" = legend.tickwidth,
                       "legend.length" = legend.length,
                       "legend.width" = legend.width,
                       "axis.text.x.angle" = axis.text.x.angle,
                       "min.overlap" = min.overlap,
                       "minGSSize" = minGSSize,
                       "maxGSSize" = maxGSSize,
                       "p.adjust.cutoff" = p.adjust.cutoff,
                       "number.breaks" = number.breaks)
  check_type(parameters = numeric_list, required_type = "numeric", test_function = is.numeric)
  # Check character parameters.
  character_list <- list("legend.position" = legend.position,
                         "plot.title" = plot.title,
                         "plot.subtitle" = plot.subtitle,
                         "plot.caption" = plot.caption,
                         "xlab" = xlab,
                         "ylab" = ylab,
                         "legend.framecolor" = legend.framecolor,
                         "legend.tickcolor" = legend.tickcolor,
                         "legend.type" = legend.type,
                         "font.type" = font.type,
                         "organism" = organism,
                         "database" = database,
                         "GO_ontology" = GO_ontology,
                         "pAdjustMethod" = pAdjustMethod,
                         "genes" = genes,
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
  
  check_colors(grid.color, parameter_name = "grid.color")
  check_colors(border.color, parameter_name = "border.color")

  if (!is.null(min.overlap)){
    assertthat::assert_that(min.overlap >= 1,
                            msg = paste0(add_cross(), crayon_body("Please provide a value higher or equal to "),
                                         crayon_key("1"),
                                         crayon_body(" to "),
                                         crayon_key("min.overlap"),
                                         crayon_body(".")))
  }

  assertthat::assert_that("OrgDb" %in% class(org.db),
                          msg = paste0(add_cross(), crayon_body("Please provide a valid "),
                                       crayon_key("OrgDb object"),
                                       crayon_body(" to "),
                                       crayon_key("org.db"),
                                       crayon_body(".")))

  check_parameters(parameter = font.type, parameter_name = "font.type")
  check_parameters(parameter = legend.type, parameter_name = "legend.type")
  check_parameters(parameter = legend.position, parameter_name = "legend.position")
  check_parameters(parameter = axis.text.x.angle, parameter_name = "axis.text.x.angle")
  check_parameters(parameter = database, parameter_name = "database")
  check_parameters(parameter = GO_ontology, parameter_name = "GO_ontology")
  check_parameters(parameter = pAdjustMethod, parameter_name = "pAdjustMethod")
  check_parameters(plot.title.face, parameter_name = "plot.title.face")
  check_parameters(plot.subtitle.face, parameter_name = "plot.subtitle.face")
  check_parameters(plot.caption.face, parameter_name = "plot.caption.face")
  check_parameters(axis.title.face, parameter_name = "axis.title.face")
  check_parameters(axis.text.face, parameter_name = "axis.text.face")
  check_parameters(legend.title.face, parameter_name = "legend.title.face")
  check_parameters(legend.text.face, parameter_name = "legend.text.face")

  check_colors(legend.framecolor, parameter_name = "legend.framecolor")
  check_colors(legend.tickcolor, parameter_name = "legend.tickcolor")

  # Convert genes to ENTREZIDs.
  suppressMessages({
    suppressWarnings({
      conversion <-clusterProfiler::bitr(genes, fromType = "SYMBOL",
                                         toType = "ENTREZID",
                                         OrgDb = org.db)
    })
  })
  
  if (is.null(min.overlap)){
    min.overlap <- if(length(genes) <= 4){1} else {3}
  }
  
  colors.use <- c("Present" = "#1e3d59", 
                  "Absent" = "#bccbcd")
  
  if (database == "GO"){
    # Compute enriched GO terms.
    result <- clusterProfiler::enrichGO(gene = conversion$ENTREZID,
                                        OrgDb = org.db,
                                        ont = GO_ontology,
                                        pvalueCutoff = 0.05,
                                        qvalueCutoff = 0.05,
                                        pAdjustMethod = pAdjustMethod,
                                        minGSSize = minGSSize,
                                        maxGSSize = maxGSSize,
                                        readable = TRUE,
                                        pool = FALSE)

    if (is.null(result)){
      # nocov start
      return_obj <- "No gene could be mapped, try another database."
      return(return_obj)
      # nocov end
    }
  # nocov start
  } else if (database == "KEGG"){
    suppressMessages({
      result <- clusterProfiler::enrichKEGG(gene = conversion$ENTREZID,
                                            organism = organism,
                                            pvalueCutoff = 0.05,
                                            qvalueCutoff = 0.05,
                                            pAdjustMethod = pAdjustMethod,
                                            minGSSize = minGSSize,
                                            maxGSSize = maxGSSize)
    })

    if (is.null(result)){
      # nocov start
      return_obj <- "No gene could be mapped, try another database."
      return(return_obj)
      # nocov end
    }

    geneID_column <- NULL
    for (input in result@result$geneID){
      genes.use <- stringr::str_split(input, pattern = "/")[[1]]
      suppressMessages({
        conversion <- clusterProfiler::bitr(genes.use, fromType = "ENTREZID",
                                            toType = "SYMBOL",
                                            OrgDb = org.db)
      })
      geneID_column <- append(geneID_column, paste(conversion$SYMBOL, collapse = "/"))
    }
    result@result$geneID <- geneID_column
  }
  # nocov end

  # Filter out non-significant results
  result@result <-  result@result %>%
                    dplyr::arrange(dplyr::desc(.data$Count)) %>%
                    dplyr::filter(.data$p.adjust <= p.adjust.cutoff,
                                  .data$Count >= min.overlap)

  if (nrow(result@result) == 0){
    return_me <- "With current presets of p.adjust.cutoff and min.overlap, no enriched terms surpass the cutoffs."
  } else {
    
    df.presence <- data.frame(row.names = genes)
    df.presence.order <- data.frame(row.names = genes)
    
    data.term <- result@result %>%
                 dplyr::arrange(result@result, dplyr::desc(.data$Count)) %>%
                 dplyr::mutate("Description" = stringr::str_to_sentence(stringr::str_replace_all(.data$Description, "_", " ")))
    # Map presence/absence.
    for (term in data.term$Description){
      data.use <- data.term %>% dplyr::filter(.data$Description == term)
      genes.use <- data.use$geneID
      genes.use <- stringr::str_split(genes.use, pattern = "/")[[1]]
      df.presence.order[[term]] <- ifelse(rownames(df.presence.order) %in% genes.use, 1, 0)
      df.presence[[term]] <- ifelse(rownames(df.presence) %in% genes.use, "Present", "Absent")
    }
    # Clustering.
    if(length(rownames(df.presence.order)) == 1){
      row_order <- rownames(df.presence.order)[1]
    } else {
      row_order <- rownames(df.presence.order)[stats::hclust(stats::dist(df.presence.order, method = "euclidean"), method = "ward.D")$order]
    }
    # nocov start
    if (length(colnames(df.presence.order)) == 1){
      col_order <- colnames(df.presence.order)[1]
    # nocov end
    } else {
      col_order <- colnames(df.presence.order)[stats::hclust(stats::dist(t(df.presence.order), method = "euclidean"), method = "ward.D")$order]
    }
    
    p.terms <- df.presence %>% 
               tibble::rownames_to_column(var = "Gene") %>% 
               tidyr::pivot_longer(cols = -dplyr::all_of("Gene"),
                                   names_to = "Description",
                                   values_to = "Status") %>% 
               dplyr::mutate("Status" = factor(.data$Status, levels = c("Present", "Absent")),
                             "Gene" = factor(.data$Gene, levels = row_order),
                             "Description" = factor(.data$Description, levels = col_order)) %>% 
               ggplot2::ggplot(mapping = ggplot2::aes(x = .data$Gene,
                                                      y = .data$Description,
                                                      fill = .data$Status)) +
               ggplot2::geom_tile(color = grid.color, linewidth = 0.5, na.rm = TRUE) +
               ggplot2::scale_y_discrete(expand = c(0, 0)) +
               ggplot2::scale_x_discrete(expand = c(0, 0),
                                         position = "top") + 
               ggplot2::coord_equal() + 
               ggplot2::scale_fill_manual(values = colors.use) + 
               ggplot2::guides(y.sec = guide_axis_label_trans(~paste0(levels(.data$Gene))),
                               x.sec = guide_axis_label_trans(~paste0(levels(.data$Description))))
       
  

    data.use <- result@result %>%
                dplyr::arrange(result@result, dplyr::desc(.data$Count)) %>%
                dplyr::mutate("Description" = stringr::str_to_sentence(stringr::str_replace_all(.data$Description, "_", " "))) %>% 
                dplyr::select(dplyr::all_of(c("Description", "Count", "p.adjust"))) %>% 
                dplyr::mutate("Description" = factor(.data$Description, levels = col_order),
                              "Counts_categorical" = "Gene count",
                              "Pval_categorical" = "-log10(p.adjust)",
                              "-log10(p.adjust)" = -1 * log10(.data$p.adjust))
    
    # Counts.
    limits <- c(min(data.use$Count, na.rm = TRUE),
                max(data.use$Count, na.rm = TRUE))
    
    scale.setup <- compute_scales(sample = sample,
                                  feature = "Test",
                                  assay = "SCT",
                                  reduction = NULL,
                                  slot = "data",
                                  number.breaks = number.breaks,
                                  min.cutoff = NA,
                                  max.cutoff = NA,
                                  flavor = "Seurat",
                                  enforce_symmetry = FALSE,
                                  from_data = TRUE,
                                  limits.use = limits)
    
    p.counts <- data.use %>% 
                ggplot2::ggplot(mapping = ggplot2::aes(x = .data$Counts_categorical,
                                                       y = .data$Description,
                                                       fill = .data$Count)) +
                ggplot2::geom_tile(color = grid.color, linewidth = 0.5, na.rm = TRUE) +
                ggplot2::scale_y_discrete(expand = c(0, 0)) +
                ggplot2::scale_x_discrete(expand = c(0, 0),
                                          position = "top") + 
                ggplot2::coord_equal() + 
                ggplot2::guides(y.sec = guide_axis_label_trans(~paste0(levels(.data$Description))),
                                x.sec = guide_axis_label_trans(~paste0(levels(.data$Counts_categorical)))) + 
                ggplot2::scale_fill_gradientn(colors = RColorBrewer::brewer.pal(n = 9, name = "YlGnBu")[2:9],
                                              na.value = "grey75",
                                              name =  "Gene count",
                                              breaks = scale.setup$breaks,
                                              labels = scale.setup$labels,
                                              limits = scale.setup$limits)
    
    
    limits <- c(min(data.use$`-log10(p.adjust)`, na.rm = TRUE),
                max(data.use$`-log10(p.adjust)`, na.rm = TRUE))
    
    scale.setup <- compute_scales(sample = sample,
                                  feature = "Test",
                                  assay = "SCT",
                                  reduction = NULL,
                                  slot = "data",
                                  number.breaks = number.breaks,
                                  min.cutoff = NA,
                                  max.cutoff = NA,
                                  flavor = "Seurat",
                                  enforce_symmetry = FALSE,
                                  from_data = TRUE,
                                  limits.use = limits)
    p.pvalue <- data.use %>% 
                ggplot2::ggplot(mapping = ggplot2::aes(x = .data$Pval_categorical,
                                                       y = .data$Description,
                                                       fill = .data$`-log10(p.adjust)`)) +
                ggplot2::geom_tile(color = grid.color, linewidth = 0.5, na.rm = TRUE) +
                ggplot2::scale_y_discrete(expand = c(0, 0)) +
                ggplot2::scale_x_discrete(expand = c(0, 0),
                                          position = "top") + 
                ggplot2::coord_equal() + 
                ggplot2::guides(y.sec = guide_axis_label_trans(~paste0(levels(.data$Description))),
                                x.sec = guide_axis_label_trans(~paste0(levels(.data$Pval_categorical)))) + 
                ggplot2::scale_fill_gradientn(colors = RColorBrewer::brewer.pal(n = 9, name = "YlOrRd")[2:9],
                                              na.value = "grey75",
                                              name =  "-log10(p.adjust)",
                                              breaks = scale.setup$breaks,
                                              labels = scale.setup$labels,
                                              limits = scale.setup$limits)
    
    list.plots <- list("Terms" = p.terms,
                       "Counts" = p.counts,
                       "Signif" = p.pvalue)
    
    counter <- 0
    for (name in names(list.plots)){
      counter <- counter + 1
      
      if (name == "Terms"){
        xlab <- "Genes"
        ylab <- "Terms"
      } else if (name == "Counts"){
        xlab <- NULL
        ylab <- NULL
      } else if (name == "Signif"){
        xlab <- NULL
        ylab <- NULL
      }
      p.use <- list.plots[[name]]
      
      axis.parameters <- handle_axis(flip = "TRUE",
                                     group.by = rep("A", length(names(list.plots))),
                                     group = name,
                                     counter = counter,
                                     axis.text.x.angle = axis.text.x.angle,
                                     plot.title.face = plot.title.face,
                                     plot.subtitle.face = plot.subtitle.face,
                                     plot.caption.face = plot.caption.face,
                                     axis.title.face = axis.title.face,
                                     axis.text.face = axis.text.face,
                                     legend.title.face = legend.title.face,
                                     legend.text.face = legend.text.face)
      
      # Modify continuous legends.
      if (name %in% c("Counts", "Signif")){
        p.use <- modify_continuous_legend(p = p.use,
                                          legend.aes = "fill",
                                          legend.type = legend.type,
                                          legend.position = legend.position,
                                          legend.length = legend.length,
                                          legend.width = legend.width,
                                          legend.framecolor = legend.framecolor,
                                          legend.tickcolor = legend.tickcolor,
                                          legend.framewidth = legend.framewidth,
                                          legend.tickwidth = legend.tickwidth)
      }

      # Set theme
      p.use <- p.use +
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
                              legend.text = ggplot2::element_text(face = legend.text.face),
                              legend.title = ggplot2::element_text(face = legend.title.face),
                              plot.title.position = "plot",
                              panel.grid = ggplot2::element_blank(),
                              panel.grid.minor.y = ggplot2::element_line(color = "white", linewidth = 1),
                              text = ggplot2::element_text(family = font.type),
                              plot.caption.position = "plot",
                              legend.justification = "center",
                              plot.margin = ggplot2::margin(t = 0, 
                                                            r = 0, 
                                                            b = 0, 
                                                            l = 5),
                              panel.border = ggplot2::element_rect(fill = NA, color = border.color, linewidth = 1),
                              panel.grid.major = ggplot2::element_blank(),
                              plot.background = ggplot2::element_rect(fill = "white", color = "white"),
                              panel.background = ggplot2::element_rect(fill = "white", color = "white"),
                              legend.background = ggplot2::element_rect(fill = "white", color = "white"))
      
      if (name == "Terms"){
        p.use <- p.use + 
                 ggplot2::guides(fill = ggplot2::guide_legend(title.position = "top",
                                                              title.hjust = 0.5,
                                                              ncol = 1))
      }
      
      list.plots[[name]] <- p.use
    }
    
    # Join the heatmaps.
    p <- patchwork::wrap_plots(list.plots,
                               nrow = 1,
                               guides = "collect")
    p <- p +
         patchwork::plot_annotation(theme = ggplot2::theme(legend.position = legend.position,
                                                           plot.title = ggplot2::element_text(family = font.type,
                                                                                              color = "black",
                                                                                              face = plot.title.face,
                                                                                              hjust = 0),
                                                           plot.subtitle = ggplot2::element_text(family = font.type,
                                                                                                 face = plot.subtitle.face,
                                                                                                 color = "black",
                                                                                                 hjust = 0),
                                                           plot.caption = ggplot2::element_text(family = font.type,
                                                                                                face = plot.caption.face,
                                                                                                color = "black",
                                                                                                hjust = 1),
                                                           plot.caption.position = "plot"))
    
    
    if (isTRUE(return_matrix)){
      return_me <- list("Result" = result@result,
                        "Plot" = p)
    } else {
      return_me <- p
    }
  }
  return(return_me)
}
