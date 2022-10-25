\donttest{
  # Check Suggests.
  value <- SCpubr:::check_suggests(function_name = "do_CellularStatesPlot", passive = TRUE)

  if (isTRUE(value)){
    # Consult the full documentation in https://enblacar.github.io/SCpubr-book/

    # Define your Seurat object.
    sample <- readRDS(system.file("extdata/seurat_dataset_example.rds", package = "SCpubr"))

    # Define some gene sets to query. It has to be a named list.
    gene_set <- list("A" = rownames(sample)[1:10],
                     "B" = rownames(sample)[11:20],
                     "C" = rownames(sample)[21:30],
                     "D" = rownames(sample)[31:40])

    # Using two variables: A scatter plot X vs Y.
    p <- SCpubr::do_CellularStatesPlot(sample = sample,
                                       input_gene_list = gene_set,
                                       x1 = "A",
                                       y1 = "B",
                                       nbin = 1,
                                       ctrl = 10)
    p


    # Using three variables. Figure from: https://www.nature.com/articles/nature20123.
    p <- SCpubr::do_CellularStatesPlot(sample = sample,
                                       input_gene_list = gene_set,
                                       x1 = "A",
                                       y1 = "B",
                                       x2 = "C",
                                       nbin = 1,
                                       ctrl = 10)
    p


    # Using four variables. Figure from: https://pubmed.ncbi.nlm.nih.gov/31327527/
    p <- SCpubr::do_CellularStatesPlot(sample = sample,
                                       input_gene_list = gene_set,
                                       x1 = "A",
                                       y1 = "C",
                                       x2 = "B",
                                       y2 = "D",
                                       nbin = 1,
                                       ctrl = 10)
    p
  } else if (isFALSE(value)){
    message("This function can not be used without its suggested packages.")
    message("Check out which ones are needed using `SCpubr::state_dependencies()`.")
  }
}

