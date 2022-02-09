#' Retrieve statistically enriched cells for a list of genes in gradient-like enrichment cases.
#'
#' Designed for 10X datasets. The following analysis aims to find a way to statistically select the cells that are truly enriched in a list of marker genes.
#' For this, a permutation analysis approach is followed, in which we compute a null distribution by shuffling the expression values
#' of the genes queried for enrichment. Therefore, this disrupts the equation for the enrichment scores, since they are, in short,
#' the result of the difference in means between the expression values of the genes queried (which we disrupted) and the control
#' selected randomly by the function. This is done until we have one million values, which would act as the null distribution.
#' Once we have a null distribution and an empirical one, we can assess how extreme our values in the empirical distribution are with
#' respect to the values in the null distribution, being the p-value the proportion of values in the null distribution that are higher
#' than the given value from the empirical distribution you are querying. These p-values are, then, adjusted for multiple testing, selecting
#' a FDR cutoff of 0.05/n, being n the total number of lists of marker genes we are going to query to the same cells and use altogether
#' to define a labeling (for instance, if we query two lists for the same tumor bulk, it would be 0.05/2). This is doing to avoid
#' over inflation of the alpha error.
#'
#' Estimated running time: 15 minutes.
#'
#' @param sample Seurat object.
#' @param markers Named list of marker genes. Can contain multiple lists, the important point here is that each list has to be named.
#' @param list.name Name of the list of markers to use.
#' @param compute_enrichment Whether to compute the enrichment scores for the requested list of genes using \link[Seurat]{AddModuleScore}.
#' @param FDR_cutoff FDR cutoff to apply (ranging from 0 to 1).
#' @param group.by Variable you want the cells to be colored for in the output DimPlot.
#' @param colors.use Vector of named HEX values to color the cells. It has to match the number of unique values in either `Seurat::Idents(sample)` or the group.by variable.
#' @param number_comparisons Number of different lists that are going to be queried to the same cells. FDR value will be divided by this.
#' @param verbose Defaults to TRUE. It will provide different print statements. Progress bars can not be suppressed by this.
#' @return A list containing the plots and the surpassed cells, together with the p-value matrix.
#' @export
#'
#' @examples
#' \dontrun{
#' TBD
#' }
do_PTEA <- function(sample,
                    markers,
                    list.name,
                    group.by = NULL,
                    colors.use = NULL,
                    compute_enrichment = FALSE,
                    FDR_cutoff = 0.05,
                    number_comparisons = 1,
                    verbose = T){


  # Define pipe operator internally.
  `%>%` <- purrr::`%>%`

  check_color_list <- !(is.null(colors.use))
  check_group.by <- !(is.null(group.by))

  # If the user has provided a color list.
  if (check_color_list){
    if (check_group.by == FALSE){
      # Reduce the color list to contain only values included in the sample.
      colors.use <- colors.use[names(colors.use) %in% levels(sample)]

      group.by_values <- levels(sample)
      if (names(colors.use) != unique(sample@meta.data[, group.by])){
        stop(paste0("The color list provided does not contain all possible unique values stored in the Seurat object identities."))
      }
    } else if (check_group.by == TRUE) {
      # Reduce the color list to contain only values included in the sample.
      colors.use <- colors.use[names(colors.use) %in% unique(sample@meta.data[, group.by])]

      if (names(colors.use) != unique(sample@meta.data[, group.by])){
        stop(paste0("The color list provided does not contain all possible unique values stored in ", group.by, " metadata variable."))
      }
    }
  }

  if (compute_enrichment == TRUE){
    if (verbose){message(paste("Computing enrichment scores:", list.name))}
    sample <- Seurat::AddModuleScore(sample, features = list(markers[[list.name]]), name = list.name)
  }

  if (verbose){message(paste("Running list:", list.name))}
  # Retrieve empirical distribution.
  enrichment_name <- paste0(stringr::str_replace_all(list.name, "-", "."))
  scores_name <- paste0(enrichment_name, "1") # Seurat::AddModuleScore adds a "1" to the given name because yes.
  test.dist <- sample[[]][, scores_name]
  names(test.dist) <- colnames(sample)

  # Compute null distribution.
  # What happens inside the replicate seems to work on a different environment level.
  #set.seed(777) # Reproducibility.

  if (verbose){message("Computing permutations.")}
  # We want a null distribution with at least 1.000.000 permutations.
  wanted_permutations <- 1000000
  number_cells <- ncol(sample)
  nruns <- trunc(wanted_permutations / number_cells) + 1
  null.dist <- pbapply::pbreplicate(nruns, {
    genes.query <- markers[[list.name]][markers[[list.name]] %in% rownames(sample)]
    # First, create a replacement object.
    sample.null <- sample
    # Get the normalized data assay (sparse matrix).
    data.use <- sample.null@assays$SCT@data
    row.order <- rownames(data.use) # We want to preserve the original order of the matrix.

    # Subset out the matrix we do not want to reshuffle.
    data.keep <- data.use[-which(rownames(data.use) %in% genes.query), ] # Remove the input genes from the matrix.

    # Get the subset for reshuffle.
    data.to.shuffle <- data.use[which(rownames(data.use) %in% genes.query), ]
    # For each gene in the list of markers.
    df.new <- list()
    for (gene in genes.query){
      # Get the scores.
      expression.scores <- data.to.shuffle[gene, ]
      # Permute the scores for all cells for that given gene.
      shuffled.scores <- sample(expression.scores, length(expression.scores))
      # As the cell names get shuffled as well, we have to change them back to the original order.
      names(shuffled.scores) <- names(expression.scores)
      df.new[[gene]] <- shuffled.scores
    }
    data.to.shuffle <- Matrix::as.matrix(t(sapply(df.new, unlist)))
    data.use <- Matrix::rBind(data.keep, data.to.shuffle) # Perform a rowbind of the two matrices (instead of directly modifying the first one, which takes forever to do)
    data.use <- data.use[row.order, ] # Reorder back the matrix.
    # Set the new matrix as the Assay data from which the enrichment scores will be computed on.
    sample.null <- Seurat::SetAssayData(sample.null, assay = "SCT", slot = "data", new.data = data.use)
    # Compute enrichment scores, which will be the null distribution of the iteration.
    sample.null <- Seurat::AddModuleScore(sample.null, list(markers[[list.name]]), name = enrichment_name)
    # Retrieve teh null distribution.
    null.dist <- sample.null[[]][, scores_name]
    # Return it.
    return(null.dist)
  })

  # Assign the cell names back to the output matrix.
  rownames(null.dist) <- colnames(sample)

  # Prepare the data for plotting.
  data.plot <- as.data.frame(null.dist) %>% # Make the matrix into tidiverse accepted object.
    dplyr::mutate(Empirical = test.dist) %>% # Add the empirical distribution to the matrix of AddModuleScore() iterations.
    tidyr::pivot_longer(dplyr::everything(), names_to = "Distribution", values_to = "Enrichment") %>%
    dplyr::mutate(Distribution = ifelse(rlang::.data$Distribution == "Empirical", "Empirical", "Null")) # Assign any name in Distribution that is not "Empirical" into "Null". By default, non-labelled df columns are named V1, V2...

  # Visualize the density plot of both distributions.
  p.dist <- data.plot %>% # Transform from wide to long format.
    ggplot2::ggplot() +
    ggplot2::geom_density(mapping = ggplot2::aes(x = rlang::.data$Enrichment, color = rlang::.data$Distribution)) +
    ggplot2::scale_color_manual(values = colortools::opposite("steelblue")) +
    ggplot2::ggtitle(paste0(list.name)) +
    ggpubr::theme_pubr(legend = "right")


  # Generate the p-values for each enrichment score in the empirical distribution.
  num_permutations <- sum(data.plot$Distribution == "Null")
  null_dist_values <- data.plot %>% # Gather the null distribution.
    dplyr::filter(rlang::.data$Distribution == "Null") %>% # Filter only the values for the NULL.
    dplyr::select(rlang::.data$Enrichment)
  null_dist_values <- null_dist_values$Enrichment

  if (verbose){message("Computing p-values.")}
  p.value.vector <- c() # Will store all the p-values and will become a column of dist.data.
  # This might also take 10 minutes, since the vector is of 1 million data points to be really exact on the null distribution side.
  p.value.vector <- unlist(pbapply::pblapply(test.dist, function(x){
    greater_values <- sum(null_dist_values > x)
    p.value <- (greater_values + 1) / (num_permutations + 1) #https://pubmed.ncbi.nlm.nih.gov/21044043/
    names(p.value) <- names(x) # Assign the cell name to the p-value.
    p.value.vector <- c(p.value.vector, p.value) # Add the p-value to the output vector.
  }))

  # Generate a reporting matrix for the given permutation test.
  dist.data <- tidyr::tibble(Cell = names(test.dist),
                             Empirical = test.dist,
                             p.value = p.value.vector)

  # FDR correction.
  FDR <- FDR_cutoff / number_comparisons

  if (verbose) {message(paste0("Using the FDR cutoff of: ", FDR))}
  fdr <- FDR  # FDR to use. Divided by 4 which is the total number of lists of markers that we are gonna test and compare to the same subset of cells (tumor bulk).
  dist.data <- dist.data %>%
    dplyr::arrange(rlang::.data$p.value) %>% # Order by ascending p-value.
    dplyr::mutate(q.value = stats::p.adjust(rlang::.data$p.value, method = "BH")) %>% # Adjust for multiple testing and produce q-values.
    dplyr::mutate(significant = ifelse(rlang::.data$q.value < fdr, TRUE, FALSE)) %>%  # Assign significance.
    dplyr::mutate(significant_corrected = ifelse(rlang::.data$q.value < fdr & rlang::.data$Empirical > 0, TRUE, FALSE)) # Assess the outliers with negative enrichment scores that surpass the cutoffs.

  # Check if we had weird cases of significant cells with negative enrichments.
  if (verbose) {
    message(paste("Number of significant results:",
                sum(dist.data$significant),
                "\nNumber of significant results with enrichment scores higher than 0:",
                sum(dist.data$significant_corrected),
                "\nNumber of outliers (significant results with enrichment scores lower than 0):", sum(dist.data$significant) - sum(dist.data$significant_corrected)))
  }
  # Visualizations.
  # UMAP coloring the cells that surpassed the FDR correction.
  surpassing_cells <- dist.data %>%
    dplyr::filter(rlang::.data$significant_corrected == TRUE) %>%
    dplyr::pull(rlang::.data$Cell)

  p.umap <- SCpubr::do_DimPlot(sample, label = T, legend = F, group.by = group.by, colors.use = colors.use) |
            SCpubr::do_DimPlot(sample, cells.highlight = surpassing_cells, legend = F, plot.title = paste("Cells that surpassed FDR correction",
                               sum(dist.data$q.value < fdr),
                               "\nFDR applied: ",
                               fdr,
                               "\nFalse positives according to FDR:",
                               trunc(sum(dist.data$q.value < fdr) * fdr),
                               "\nEnrichment score cutoff:",
                               round(min(dist.data$Empirical[sum(dist.data$q.value < fdr)]), 2))) |
            SCpubr::do_FeaturePlot(sample, features = scores_name, plot.title = paste0("Enrichment scores for list\n", list.name))


  # Prepare the output of the function.
  output <- list(empirical_data = dist.data,
                 null_data = data.plot,
                 surpassed_cells = surpassing_cells,
                 p.dist = p.dist,
                 p.umap = p.umap)
  return(output)
}
