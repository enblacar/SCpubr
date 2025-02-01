if (base::isFALSE(dep_check[["do_ColorBlindCheck"]])){
  testthat::test_that("do_ColorBlindCheck: PASS - color vectors", {
    
    p <- SCpubr::do_ColorBlindCheck(colors.use = c("red", "blue", "green"))
    testthat::expect_type(p, "list")
    
  })
}