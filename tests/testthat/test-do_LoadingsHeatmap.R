if (base::isFALSE(dep_check[["do_LoadingsHeatmap"]])){
  
  testthat::test_that("do_LoadingsHeatmap: CRAN essentials", {
    genes <- list("A" = rownames(sample)[1:5],
                  "B" = rownames(sample)[6:10],
                  "C" = rownames(sample)[11:15])
    
    p <- SCpubr::do_LoadingsHeatmap(sample = sample,
                                    dims = 1:5)
    testthat::expect_true(ggplot2::is_ggplot(p))
    
    
  })
  
  testthat::test_that("do_LoadingsHeatmap: PASS - default", {
    testthat::skip_on_cran()
    p <- SCpubr::do_LoadingsHeatmap(sample = sample,
                                 dims = 1:10)
    testthat::expect_true(ggplot2::is_ggplot(p))
    
    p <- SCpubr::do_LoadingsHeatmap(sample = sample,
                                    dims = 1:10,
                                    min.cutoff.loadings = 0.5)
    testthat::expect_true(ggplot2::is_ggplot(p))
    
    p <- SCpubr::do_LoadingsHeatmap(sample = sample,
                                    dims = 1:10,
                                    max.cutoff.loadings = 0.5)
    testthat::expect_true(ggplot2::is_ggplot(p))
    
    p <- SCpubr::do_LoadingsHeatmap(sample = sample,
                                    dims = 1:10,
                                    min.cutoff.expression = 0.5)
    testthat::expect_true(ggplot2::is_ggplot(p))
    
    p <- SCpubr::do_LoadingsHeatmap(sample = sample,
                                    dims = 1:10,
                                    max.cutoff.expresion = 0.5)
    testthat::expect_true(ggplot2::is_ggplot(p))
    
    p <- SCpubr::do_LoadingsHeatmap(sample = sample,
                                 dims = 1:10,
                                 subsample = 100)
    testthat::expect_true(ggplot2::is_ggplot(p))
    
    sample$test <- as.factor(sample$seurat_clusters)
    p <- SCpubr::do_LoadingsHeatmap(sample = sample,
                                 dims = 1:10,
                                 group.by = "test")
    testthat::expect_true(ggplot2::is_ggplot(p))
    
    p <- SCpubr::do_LoadingsHeatmap(sample = sample,
                                 dims = 1:10,
                                 min.cutoff.loadings = -0.01,
                                 max.cutoff.loadings = 0.01,
                                 min.cutoff.expression = 0,
                                 max.cutoff.expression = 0.75)
    testthat::expect_true(ggplot2::is_ggplot(p))
    
    p <- SCpubr::do_LoadingsHeatmap(sample = sample,
                                 dims = 1:10,
                                 use_viridis = TRUE,
                                 viridis.direction = 1)
    testthat::expect_true(ggplot2::is_ggplot(p))
    
    p <- SCpubr::do_LoadingsHeatmap(sample = sample,
                                 dims = 1:10,
                                 use_viridis = TRUE,
                                 viridis.direction = -1)
    testthat::expect_true(ggplot2::is_ggplot(p))
    
    
    p <- SCpubr::do_LoadingsHeatmap(sample = sample,
                                 dims = 1:10,
                                 use_viridis = FALSE,
                                 sequential.direction = 1)
    testthat::expect_true(ggplot2::is_ggplot(p))
    
    p <- SCpubr::do_LoadingsHeatmap(sample = sample,
                                 dims = 1:10,
                                 use_viridis = FALSE,
                                 sequential.direction = -1)
    testthat::expect_true(ggplot2::is_ggplot(p))
    
    
  })
}


