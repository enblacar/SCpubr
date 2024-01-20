if (base::isFALSE(dep_check[["do_WafflePlot"]])){
  
  testthat::test_that("do_WafflePlot: CRAN essential tests", {

    p <- SCpubr::do_WafflePlot(sample = sample,
                               group.by = "seurat_clusters")
    testthat::expect_type(p, "list")
  })
  
  testthat::test_that("do_WafflePlot: PASS - flip", {
    testthat::skip_on_cran()
    
    p <- SCpubr::do_WafflePlot(sample = sample,
                               group.by = "seurat_clusters",
                               flip = FALSE)
    testthat::expect_type(p, "list")
    
    p <- SCpubr::do_WafflePlot(sample = sample,
                               group.by = "seurat_clusters",
                               flip = TRUE)
    testthat::expect_type(p, "list")
  })
}