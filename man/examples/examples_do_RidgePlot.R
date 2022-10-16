\donttest{
  # Check Suggests.
  value <- SCpubr:::check_suggests(function_name = "do_RidgePlot", passive = TRUE)

  if (isTRUE(value)){
    # Define your Seurat object.
    sample <- readRDS(system.file("extdata/seurat_dataset_example.rds", package = "SCpubr"))

    # Compute the most basic ridge plot.
    p <- SCpubr::do_RidgePlot(sample = sample,
                              feature = "nFeature_RNA")
    p

    # Use continuous color scale.
    p1 <- SCpubr::do_RidgePlot(sample = sample,
                               feature = "nFeature_RNA",
                               continuous_scale = TRUE,
                               viridis_direction = 1)

    p2 <- SCpubr::do_RidgePlot(sample = sample,
                               feature = "nFeature_RNA",
                               continuous_scale = TRUE,
                               viridis_direction = -1)

    p <- p1 / p2
    p

    # Draw quantiles of the distribution.
    p1 <- SCpubr::do_RidgePlot(sample = sample,
                               feature = "nFeature_RNA",
                               continuous_scale = TRUE,
                               compute_quantiles = TRUE,
                               compute_custom_quantiles = TRUE)

    p2 <- SCpubr::do_RidgePlot(sample = sample,
                               feature = "nFeature_RNA",
                               continuous_scale = TRUE,
                               compute_quantiles = TRUE,
                               compute_custom_quantiles = TRUE,
                               quantiles = c(0.1, 0.5, 0.75))

    p <- p1 / p2
    p

    # Draw probability tails.
    p1 <- SCpubr::do_RidgePlot(sample = sample,
                               feature = "nFeature_RNA",
                               continuous_scale = TRUE,
                               compute_quantiles = TRUE,
                               compute_distribution_tails = TRUE)

    p2 <- SCpubr::do_RidgePlot(sample = sample,
                               feature = "nFeature_RNA",
                               continuous_scale = TRUE,
                               compute_quantiles = TRUE,
                               compute_distribution_tails = TRUE,
                               prob_tails = 0.3)

    p <- p1 / p2
    p

    # Draw probability tails.
    p <- SCpubr::do_RidgePlot(sample = sample,
                              feature = "nFeature_RNA",
                              continuous_scale = TRUE,
                              compute_quantiles = TRUE,
                              color_by_probabilities = TRUE)
    p
  } else if (isFALSE(value)){
    message("This function can not be used without its suggested packages.")
    message("Check out which ones are needed using `SCpubr::state_dependencies()`.")
  }
}
