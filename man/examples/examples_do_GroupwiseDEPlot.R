\donttest{
  # Check Suggests.
  value <- SCpubr:::check_suggests(function_name = "do_GroupwiseDEPlot", passive = TRUE)

  if (isTRUE(value)){
    # Define your Seurat object.
    sample <- readRDS(system.file("extdata/seurat_dataset_example.rds", package = "SCpubr"))

    # Compute DE genes and transform to a tibble.
    de_genes <- readRDS(system.file("extdata/de_genes_example.rds", package = "SCpubr"))

    # Default output.
    p <- SCpubr::do_GroupwiseDEPlot(sample = sample,
                                    de_genes = de_genes)

    p

    # Increase the number of top DE genes by cluster.
    p <- SCpubr::do_GroupwiseDEPlot(sample = sample,
                                    de_genes = de_genes,
                                    top_genes = 2)

    p

    # Modify the row and column titles and the rotation.
    p <- SCpubr::do_GroupwiseDEPlot(sample = sample,
                                    de_genes = de_genes,
                                    column_title = "Title A",
                                    row_title_p_values = "Title B",
                                    row_title_logfc = "Title C",
                                    row_title_expression = "Title D",
                                    row_title_rot = 0)

    p

    sample$modified_orig.ident <- sample(x = c("Sample_A", "Sample_B", "Sample_C"),
                                         size = ncol(sample),
                                         replace = TRUE,
                                         prob = c(0.2, 0.7, 0.1))

    # Add more layers of mean expression with group.by.
    p <- SCpubr::do_GroupwiseDEPlot(sample = sample,
                                    de_genes = de_genes,
                                    group.by = c("seurat_clusters",
                                                 "modified_orig.ident",
                                                 "orig.ident"),
                                    row_title_expression = c("",
                                                             "Title A",
                                                             "Title B"))

    p

    # Change the viridis scales.
    p <- SCpubr::do_GroupwiseDEPlot(sample = sample,
                                    de_genes = de_genes,
                                    viridis_map_pvalues = "C",
                                    viridis_map_logfc = "E",
                                    viridis_map_expression = "D")

    p
  } else if (isFALSE(value)){
    message("This function can not be used without its suggested packages.")
    message("Check out which ones are needed using `SCpubr::state_dependencies()`.")
  }
}
