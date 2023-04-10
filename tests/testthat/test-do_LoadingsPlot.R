if (isFALSE(dep_check[["do_LoadingsPlot"]])){
  
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
  })
}


