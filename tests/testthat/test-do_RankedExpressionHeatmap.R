if (base::isFALSE(dep_check[["do_RankedExpressionHeatmap"]])){
  
  testthat::test_that("do_RankedExpressionHeatmap: CRAN essentials", {
    genes <- Seurat::VariableFeatures(sample)[1:30]
    
    p <- SCpubr::do_RankedExpressionHeatmap(sample = sample,
                                            features =  genes,
                                            subsample = NA,
                                            reduction = "umap",
                                            dims = 1:2,
                                            verbose = FALSE)
    testthat::expect_type(p, "list")
    
    
  })
  
  testthat::test_that("do_RankedExpressionHeatmap: PASS - default", {
    testthat::skip_on_cran()
    genes <- Seurat::VariableFeatures(sample)[1:30]
    
    p <- SCpubr::do_RankedExpressionHeatmap(sample = sample,
                                            features =  genes,
                                            subsample = NA,
                                            reduction = "umap",
                                            dims = 1:2,
                                            return_object = TRUE,
                                            verbose = FALSE,
                                            use_viridis = TRUE,
                                            enforce_symmetry = FALSE)
    testthat::expect_type(p, "list")
    
    p <- SCpubr::do_RankedExpressionHeatmap(sample = sample,
                                            features =  genes,
                                            subsample = NA,
                                            reduction = "umap",
                                            dims = 1:2,
                                            return_object = TRUE,
                                            verbose = FALSE,
                                            use_viridis = FALSE,
                                            sequential.direction = 1,
                                            enforce_symmetry = FALSE)
    testthat::expect_type(p, "list")
    
    p <- SCpubr::do_RankedExpressionHeatmap(sample = sample,
                                            features =  genes,
                                            subsample = NA,
                                            reduction = "umap",
                                            dims = 1:2,
                                            return_object = TRUE,
                                            verbose = FALSE,
                                            use_viridis = FALSE,
                                            sequential.direction = -1,
                                            enforce_symmetry = FALSE)
    testthat::expect_type(p, "list")
    
    
    p <- SCpubr::do_RankedExpressionHeatmap(sample = sample,
                                            features =  genes,
                                            subsample = 120,
                                            reduction = "umap",
                                            dims = 1:2,
                                            return_object = TRUE,
                                            verbose = FALSE,
                                            use_viridis = FALSE,
                                            enforce_symmetry = FALSE)
    testthat::expect_type(p, "list")
    
    SCpubr::do_RankedExpressionHeatmap(sample = sample,
                                                                 features =  genes,
                                                                 subsample = 100,
                                                                 reduction = "umap",
                                                                 dims = 1:2,
                                                                 return_object = TRUE,
                                                                 verbose = FALSE,
                                                                 assay = "SCT",
                                                                 use_viridis = FALSE,
                                                                 enforce_symmetry = FALSE)
    
   SCpubr::do_RankedExpressionHeatmap(sample = sample,
                                                                 features =  genes,
                                                                 subsample = 100,
                                                                 reduction = "umap",
                                                                 dims = 1:2,
                                                                 return_object = TRUE,
                                                                 verbose = FALSE,
                                                                 slot = "data",
                                                                 use_viridis = FALSE,
                                                                 enforce_symmetry = FALSE)
    
    suppressMessages({testthat::expect_message({p <- SCpubr::do_RankedExpressionHeatmap(sample = sample,
                                                                                        features =  genes,
                                                                                        subsample = 100,
                                                                                        reduction = "umap",
                                                                                        dims = 1:2,
                                                                                        return_object = TRUE,
                                                                                        verbose = TRUE)})})
    testthat::expect_type(p, "list")
    
    p <- SCpubr::do_RankedExpressionHeatmap(sample = sample,
                                            features =  genes,
                                            subsample = 100,
                                            group.by = c("orig.ident", "seurat_clusters"),
                                            colors.use = list("orig.ident" = c("Cell" = "red")),
                                            reduction = "umap",
                                            dims = 1:2,
                                            return_object = TRUE,
                                            verbose = FALSE,
                                            use_viridis = FALSE,
                                            enforce_symmetry = FALSE)
    testthat::expect_type(p, "list")
    
  })
}


