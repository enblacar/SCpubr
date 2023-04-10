if (isFALSE(dep_check[["do_SCExpressionHeatmap"]])){
  
  testthat::test_that("do_SCExpressionHeatmap: CRAN essentials", {
    
    p <- SCpubr::do_SCExpressionHeatmap(sample = sample,
                                        features = rownames(sample)[1:5])
    testthat::expect_type(p, "list")
    
    
  })
  
  testthat::test_that("do_SCExpressionHeatmap: PASS - default", {
    testthat::skip_on_cran()

    p <- SCpubr::do_SCExpressionHeatmap(sample = sample,
                                        features = rownames(sample)[1:5])
    testthat::expect_type(p, "list")
  })
}


