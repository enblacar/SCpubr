if (base::isFALSE(dep_check[["do_LoadingsPlot"]])){
  
  testthat::test_that("do_LoadingsPlot: CRAN essentials", {
    genes <- list("A" = rownames(sample)[1:5],
                  "B" = rownames(sample)[6:10],
                  "C" = rownames(sample)[11:15])
    
    p <- SCpubr::do_LoadingsPlot(sample = sample,
                                 dims = 1:5)
    testthat::expect_type(p, "list")
    
    
  })
  
  testthat::test_that("do_LoadingsPlot: PASS - default", {
    testthat::skip_on_cran()
    p <- SCpubr::do_LoadingsPlot(sample = sample,
                                 dims = 1:10)
    testthat::expect_type(p, "list")
    
    p <- SCpubr::do_LoadingsPlot(sample = sample,
                                 dims = 1:10,
                                 subsample = 100)
    testthat::expect_type(p, "list")
    
    sample$test <- as.factor(sample$seurat_clusters)
    p <- SCpubr::do_LoadingsPlot(sample = sample,
                                 dims = 1:10,
                                 group.by = "test")
    testthat::expect_type(p, "list")
    
    p <- SCpubr::do_LoadingsPlot(sample = sample,
                                 dims = 1:10,
                                 min.cutoff.loadings = -0.01,
                                 max.cutoff.loadings = 0.01,
                                 min.cutoff.expression = 0,
                                 max.cutoff.expression = 0.75)
    testthat::expect_type(p, "list")
    
    p <- SCpubr::do_LoadingsPlot(sample = sample,
                                 dims = 1:10,
                                 use_viridis = TRUE,
                                 viridis.direction = 1)
    testthat::expect_type(p, "list")
    
    p <- SCpubr::do_LoadingsPlot(sample = sample,
                                 dims = 1:10,
                                 use_viridis = TRUE,
                                 viridis.direction = -1)
    testthat::expect_type(p, "list")
    
    
    p <- SCpubr::do_LoadingsPlot(sample = sample,
                                 dims = 1:10,
                                 use_viridis = FALSE,
                                 sequential.direction = 1)
    testthat::expect_type(p, "list")
    
    p <- SCpubr::do_LoadingsPlot(sample = sample,
                                 dims = 1:10,
                                 use_viridis = FALSE,
                                 sequential.direction = -1)
    testthat::expect_type(p, "list")
    
    
  })
}


