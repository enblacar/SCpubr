if (base::isFALSE(dep_check[["do_SCExpressionHeatmap"]])){
  
  testthat::test_that("do_SCExpressionHeatmap: CRAN essentials", {
    
    p <- SCpubr::do_SCExpressionHeatmap(sample = sample,
                                        features = rownames(sample)[1:5])
    testthat::expect_true(ggplot2::is_ggplot(p))
    
    
  })
  
  testthat::test_that("do_SCExpressionHeatmap: PASS - default", {
    testthat::skip_on_cran()
    p <- SCpubr::do_SCExpressionHeatmap(sample = sample,
                                        features = rownames(sample)[1:5],
                                        features.order = rownames(sample)[c(4, 2, 1, 3, 5)],
                                        cluster = FALSE)
    testthat::expect_true(ggplot2::is_ggplot(p))
    
    p <- SCpubr::do_SCExpressionHeatmap(sample = sample,
                                        features = rownames(sample)[1:5])
    testthat::expect_true(ggplot2::is_ggplot(p))
    
    testthat::expect_warning({p <- SCpubr::do_SCExpressionHeatmap(sample = sample,
                                                                  features = c(rownames(sample)[1:5], "pepe"))})
    testthat::expect_true(ggplot2::is_ggplot(p))
    
    
    sample$test <- as.factor(sample$seurat_clusters)
    p <- SCpubr::do_SCExpressionHeatmap(sample = sample,
                                        features = rownames(sample)[1:5],
                                        group.by = "test")
    testthat::expect_true(ggplot2::is_ggplot(p))
    
    p <- SCpubr::do_SCExpressionHeatmap(sample = sample,
                                        features = rownames(sample)[1:5],
                                        subsample = 100)
    testthat::expect_true(ggplot2::is_ggplot(p))
    
    p <- SCpubr::do_SCExpressionHeatmap(sample = sample,
                                        features = rownames(sample)[1])
    testthat::expect_true(ggplot2::is_ggplot(p))
    
    p <- SCpubr::do_SCExpressionHeatmap(sample = sample,
                                        features = rownames(sample)[1:5],
                                        cluster = FALSE)
    testthat::expect_true(ggplot2::is_ggplot(p))
    
    p <- SCpubr::do_SCExpressionHeatmap(sample = sample,
                                        features = rownames(sample)[1:5],
                                        metadata = c("orig.ident", "seurat_clusters"))
    testthat::expect_true(ggplot2::is_ggplot(p))
    
    p <- SCpubr::do_SCExpressionHeatmap(sample = sample,
                                        features = rownames(sample)[1:5],
                                        metadata = c("orig.ident", "seurat_clusters"),
                                        metadata.colors = list("orig.ident" = c("Cell" = "blue")))
    testthat::expect_true(ggplot2::is_ggplot(p))
    
    p <- SCpubr::do_SCExpressionHeatmap(sample = sample,
                                        features = rownames(sample)[1:5],
                                        metadata = c("orig.ident", "seurat_clusters"),
                                        min.cutoff = 1,
                                        max.cutoff = 2)
    testthat::expect_true(ggplot2::is_ggplot(p))
    
    p <- SCpubr::do_SCExpressionHeatmap(sample = sample,
                                        features = rownames(sample)[1:5],
                                        metadata = c("orig.ident", "seurat_clusters"),
                                        min.cutoff = 1,
                                        max.cutoff = 2,
                                        proportional.size = FALSE)
    testthat::expect_true(ggplot2::is_ggplot(p))
    
    p <- SCpubr::do_SCExpressionHeatmap(sample = sample,
                                        features = rownames(sample)[1:5],
                                        metadata = c("orig.ident", "seurat_clusters"),
                                        min.cutoff = 1,
                                        max.cutoff = 2,
                                        proportional.size = FALSE,
                                        enforce_symmetry = FALSE,
                                        use_viridis = FALSE)
    testthat::expect_true(ggplot2::is_ggplot(p))
    
    p <- SCpubr::do_SCExpressionHeatmap(sample = sample,
                                        features = rownames(sample)[1:5],
                                        metadata = c("orig.ident", "seurat_clusters"),
                                        min.cutoff = 1,
                                        max.cutoff = 2,
                                        proportional.size = FALSE,
                                        enforce_symmetry = FALSE,
                                        use_viridis = TRUE)
    testthat::expect_true(ggplot2::is_ggplot(p))
    
    p <- SCpubr::do_SCExpressionHeatmap(sample = sample,
                                        features = rownames(sample)[1:5],
                                        metadata = c("orig.ident", "seurat_clusters"),
                                        min.cutoff = 1,
                                        max.cutoff = 2,
                                        proportional.size = FALSE,
                                        enforce_symmetry = TRUE,
                                        use_viridis = FALSE)
    testthat::expect_true(ggplot2::is_ggplot(p))
  })
}


