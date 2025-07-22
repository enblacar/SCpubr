if (base::isFALSE(dep_check[["do_ColorBlindCheck"]])){
  testthat::test_that("do_ColorBlindCheck: PASS - color vectors", {
    
    p <- SCpubr::do_ColorBlindCheck(colors.use = c("red", "blue", "green"), flip = TRUE)
    testthat::expect_true(ggplot2::is_ggplot(p))
    
    p <- SCpubr::do_ColorBlindCheck(colors.use = c("red", "blue", "green"), flip = FALSE)
    testthat::expect_true(ggplot2::is_ggplot(p))
    
  })
}
