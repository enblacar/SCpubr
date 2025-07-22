if(base::isFALSE(dep_check[["do_RidgePlot"]])){

  testthat::test_that("do_RidgePlot: CRAN essentials", {

    p <- SCpubr::do_RidgePlot(sample = sample,
                              feature = "nCount_RNA")
    testthat::expect_true(ggplot2::is_ggplot(p))

    p <- SCpubr::do_RidgePlot(sample = sample,
                              feature = "nCount_RNA",
                              continuous_scale = TRUE,
                              use_viridis = TRUE,
                              viridis.direction = 1)
    testthat::expect_true(ggplot2::is_ggplot(p))




  })

  testthat::test_that("do_RidgePlot: PASS - default", {
    testthat::skip_on_cran()
    
    p <- SCpubr::do_RidgePlot(sample = sample,
                              feature = "nCount_RNA",
                              continuous_scale = FALSE)
    testthat::expect_true(ggplot2::is_ggplot(p))
    
    p <- SCpubr::do_RidgePlot(sample = sample,
                              feature = "nCount_RNA",
                              continuous_scale = FALSE,
                              group.by = "annotation",
                              colors.use = c("A" = "red", "B" = "blue"))
    testthat::expect_true(ggplot2::is_ggplot(p))
    
    p <- SCpubr::do_RidgePlot(sample = sample,
                              feature = "nCount_RNA",
                              split.by = "annotation",
                              continuous_scale = TRUE,
                              use_viridis = TRUE,
                              viridis.direction = 1)
    testthat::expect_true(ggplot2::is_ggplot(p))
    
    p <- SCpubr::do_RidgePlot(sample = sample,
                              feature = "nCount_RNA",
                              continuous_scale = TRUE,
                              use_viridis = TRUE,
                              viridis.direction = 1)
    testthat::expect_true(ggplot2::is_ggplot(p))
    
    p <- SCpubr::do_RidgePlot(sample = sample,
                              feature = "nCount_RNA",
                              continuous_scale = TRUE,
                              use_viridis = TRUE,
                              viridis.direction = -1)
    testthat::expect_true(ggplot2::is_ggplot(p))
    
    p <- SCpubr::do_RidgePlot(sample = sample,
                              feature = "nCount_RNA",
                              continuous_scale = TRUE,
                              use_viridis = FALSE,
                              sequential.direction = 1)
    testthat::expect_true(ggplot2::is_ggplot(p))
    
    p <- SCpubr::do_RidgePlot(sample = sample,
                              feature = "nCount_RNA",
                              continuous_scale = TRUE,
                              use_viridis = FALSE,
                              sequential.direction = -1)
    testthat::expect_true(ggplot2::is_ggplot(p))
    


    p <- SCpubr::do_RidgePlot(sample = sample,
                              feature = "nCount_RNA")
    testthat::expect_true(ggplot2::is_ggplot(p))

    p <- SCpubr::do_RidgePlot(sample = sample,
                              feature = "nCount_RNA",
                              group.by = "orig.ident")
    testthat::expect_true(ggplot2::is_ggplot(p))

    sample$orig.ident <- factor(sample$orig.ident)

    p <- SCpubr::do_RidgePlot(sample = sample,
                              feature = "nCount_RNA",
                              group.by = "orig.ident")
    testthat::expect_true(ggplot2::is_ggplot(p))

    p <- SCpubr::do_RidgePlot(sample = sample,
                              feature = "nCount_RNA",
                              legend.position = "bottom")
    testthat::expect_true(ggplot2::is_ggplot(p))

    p <- SCpubr::do_RidgePlot(sample = sample,
                              feature = "nCount_RNA",
                              group.by = "orig.ident",
                              legend.position = "bottom",
                              colors.use = c("Cell" = "red"))
    testthat::expect_true(ggplot2::is_ggplot(p))
  })

  testthat::test_that("do_RidgePlot: PASS - plot.grid", {
    testthat::skip_on_cran()


    p <- SCpubr::do_RidgePlot(sample = sample,
                              feature = "nCount_RNA",
                              plot.grid = TRUE)
    testthat::expect_true(ggplot2::is_ggplot(p))
    
    p <- SCpubr::do_RidgePlot(sample = sample,
                              feature = "nCount_RNA",
                              plot.grid = TRUE,
                              flip = TRUE)
    testthat::expect_true(ggplot2::is_ggplot(p))
    
    p <- SCpubr::do_RidgePlot(sample = sample,
                              feature = "nCount_RNA",
                              plot.grid = TRUE,
                              flip = FALSE)
    testthat::expect_true(ggplot2::is_ggplot(p))

    p <- SCpubr::do_RidgePlot(sample = sample,
                              feature = "nCount_RNA",
                              plot.grid = FALSE)
    testthat::expect_true(ggplot2::is_ggplot(p))
  })


  testthat::test_that("do_RidgePlot: PASS - split.by", {
    testthat::skip_on_cran()


    p <- SCpubr::do_RidgePlot(sample = sample,
                              feature = "nCount_RNA",
                              split.by = "orig.ident")
    testthat::expect_true(ggplot2::is_ggplot(p))
  })

  testthat::test_that("do_RidgePlot: PASS - continuous scale", {
    testthat::skip_on_cran()


    p <- SCpubr::do_RidgePlot(sample = sample,
                              feature = "nCount_RNA",
                              continuous_scale = TRUE,
                              viridis.direction = 1)
    testthat::expect_true(ggplot2::is_ggplot(p))

    p <- SCpubr::do_RidgePlot(sample = sample,
                              feature = "nCount_RNA",
                              continuous_scale = TRUE,
                              viridis.direction = -1)
    testthat::expect_true(ggplot2::is_ggplot(p))
  })

  testthat::test_that("do_RidgePlot: PASS - group.by", {
    testthat::skip_on_cran()


    p <- SCpubr::do_RidgePlot(sample = sample,
                              feature = "nCount_RNA",
                              group.by = "orig.ident")
    testthat::expect_true(ggplot2::is_ggplot(p))
  })




  testthat::test_that("do_RidgePlot: PASS - flip", {
    testthat::skip_on_cran()


    p <- SCpubr::do_RidgePlot(sample = sample,
                              feature = "nFeature_RNA",
                              flip = TRUE)
    testthat::expect_true(ggplot2::is_ggplot(p))

    p <- SCpubr::do_RidgePlot(sample = sample,
                              feature = "nFeature_RNA",
                              flip = FALSE)
    testthat::expect_true(ggplot2::is_ggplot(p))
  })
}



