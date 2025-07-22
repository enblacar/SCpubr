if (base::isFALSE(dep_check[["do_WafflePlot"]])){
  
  testthat::test_that("do_WafflePlot: CRAN essential tests", {

    p <- SCpubr::do_WafflePlot(sample = sample,
                               group.by = "seurat_clusters")
    testthat::expect_true(ggplot2::is_ggplot(p))
  })
  
  testthat::test_that("do_WafflePlot: PASS - flip", {
    testthat::skip_on_cran()
    
    p <- SCpubr::do_WafflePlot(sample = sample,
                               group.by = "annotation",
                               flip = FALSE,
                               colors.use = c("A" = "red", "B"= "blue"))
    testthat::expect_true(ggplot2::is_ggplot(p))
    
    sample$annotation <- factor(sample$annotation)
    p <- SCpubr::do_WafflePlot(sample = sample,
                               group.by = "annotation",
                               flip = FALSE,
                               colors.use = c("A" = "red", "B"= "blue"))
    testthat::expect_true(ggplot2::is_ggplot(p))
    
    p <- SCpubr::do_WafflePlot(sample = sample,
                               group.by = "annotation",
                               flip = FALSE)
    testthat::expect_true(ggplot2::is_ggplot(p))
    
    p <- SCpubr::do_WafflePlot(sample = sample,
                               group.by = "seurat_clusters",
                               flip = TRUE)
    testthat::expect_true(ggplot2::is_ggplot(p))
  })
}
