if (base::isFALSE(dep_check[["do_RankedEnrichmentPlot"]])){
  
  testthat::test_that("do_RankedEnrichmentPlot: CRAN essentials", {
    genes <- list("A" = rownames(sample)[1:5],
                  "B" = rownames(sample)[6:10],
                  "C" = rownames(sample)[11:15])
    
    p <- SCpubr::do_RankedEnrichmentPlot(sample = sample,
                                     input_gene_list =  genes,
                                     subsample = NA,
                                     nbin = 1,
                                     ctrl = 5,
                                     reduction = "umap",
                                     dims = 1:2,
                                     verbose = FALSE)
    testthat::expect_type(p, "list")
    
    
  })
  
  testthat::test_that("do_RankedEnrichmentPlot: PASS - default", {
    testthat::skip_on_cran()
    genes <- list("A" = rownames(sample)[1:5],
                  "B" = rownames(sample)[6:10],
                  "C" = rownames(sample)[11:15])
    
    p <- SCpubr::do_RankedEnrichmentPlot(sample = sample,
                                     input_gene_list =  genes,
                                     subsample = NA,
                                     nbin = 1,
                                     ctrl = 5,
                                     reduction = "umap",
                                     dims = 1:2,
                                     return_object = TRUE,
                                     verbose = FALSE,
                                     flavor = "Seurat",
                                     use_viridis = TRUE,
                                     enforce_symmetry = FALSE)
    testthat::expect_type(p, "list")
    
    p <- SCpubr::do_RankedEnrichmentPlot(sample = sample,
                                     input_gene_list =  genes,
                                     subsample = NA,
                                     nbin = 1,
                                     ctrl = 5,
                                     reduction = "umap",
                                     dims = 1:2,
                                     return_object = TRUE,
                                     verbose = FALSE,
                                     flavor = "Seurat",
                                     use_viridis = FALSE,
                                     sequential.direction = 1,
                                     enforce_symmetry = FALSE)
    testthat::expect_type(p, "list")
    
    p <- SCpubr::do_RankedEnrichmentPlot(sample = sample,
                                     input_gene_list =  genes,
                                     subsample = NA,
                                     nbin = 1,
                                     ctrl = 5,
                                     reduction = "umap",
                                     dims = 1:2,
                                     return_object = TRUE,
                                     verbose = FALSE,
                                     flavor = "Seurat",
                                     use_viridis = FALSE,
                                     sequential.direction = -1,
                                     enforce_symmetry = FALSE)
    testthat::expect_type(p, "list")
    
    
    p <- SCpubr::do_RankedEnrichmentPlot(sample = sample,
                                     input_gene_list =  genes,
                                     subsample = 100,
                                     nbin = 1,
                                     ctrl = 5,
                                     reduction = "umap",
                                     dims = 1:2,
                                     return_object = TRUE,
                                     verbose = FALSE,
                                     flavor = "UCell",
                                     use_viridis = FALSE,
                                     enforce_symmetry = FALSE)
    testthat::expect_type(p, "list")
    
    testthat::expect_warning({SCpubr::do_RankedEnrichmentPlot(sample = sample,
                                                        input_gene_list =  genes,
                                                        subsample = 100,
                                                        nbin = 1,
                                                        ctrl = 5,
                                                        reduction = "umap",
                                                        dims = 1:2,
                                                        return_object = TRUE,
                                                        verbose = FALSE,
                                                        flavor = "UCell",
                                                        assay = "SCT",
                                                        use_viridis = FALSE,
                                                        enforce_symmetry = FALSE)})
    
    testthat::expect_warning({SCpubr::do_RankedEnrichmentPlot(sample = sample,
                                                          input_gene_list =  genes,
                                                          subsample = 100,
                                                          nbin = 1,
                                                          ctrl = 5,
                                                          reduction = "umap",
                                                          dims = 1:2,
                                                          return_object = TRUE,
                                                          verbose = FALSE,
                                                          flavor = "Seurat",
                                                          slot = "data",
                                                          use_viridis = FALSE,
                                                          enforce_symmetry = FALSE)})
    
    suppressMessages({testthat::expect_message({p <- SCpubr::do_RankedEnrichmentPlot(sample = sample,
                                                               input_gene_list =  genes,
                                                               subsample = 100,
                                                               nbin = 1,
                                                               ctrl = 5,
                                                               reduction = "umap",
                                                               dims = 1:2,
                                                               return_object = TRUE,
                                                               verbose = TRUE)})})
    testthat::expect_type(p, "list")
    
    p <- SCpubr::do_RankedEnrichmentPlot(sample = sample,
                                     input_gene_list =  genes,
                                     subsample = 100,
                                     group.by = c("orig.ident", "seurat_clusters"),
                                     colors.use = list("orig.ident" = c("Cell" = "red")),
                                     nbin = 1,
                                     ctrl = 5,
                                     reduction = "umap",
                                     dims = 1:2,
                                     return_object = TRUE,
                                     verbose = FALSE,
                                     flavor = "UCell",
                                     use_viridis = FALSE,
                                     enforce_symmetry = FALSE)
    testthat::expect_type(p, "list")
    
  })
}


