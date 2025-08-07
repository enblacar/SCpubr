if (isFALSE(dep_check[["do_DotPlot"]])){

  testthat::test_that("do_DotPlot: CRAN essentials", {

    p <- SCpubr::do_DotPlot(sample = sample,
                            features = "EPC1")
    testthat::expect_true(ggplot2::is_ggplot(p))

  })

  testthat::test_that("do_DotPlot: PASS - one variable", {
    testthat::skip_on_cran()


    p <- SCpubr::do_DotPlot(sample = sample,
                            features = "EPC1")
    testthat::expect_true(ggplot2::is_ggplot(p))
    
    p <- SCpubr::do_DotPlot(sample = sample,
                            features = "EPC1", flip = TRUE)
    testthat::expect_true(ggplot2::is_ggplot(p))
    
    p <- SCpubr::do_DotPlot(sample = sample,
                            features = "EPC1", flip = FALSE)
    testthat::expect_true(ggplot2::is_ggplot(p))
    
    p <- SCpubr::do_DotPlot(sample = sample,
                            features = "EPC1",
                            max.cutoff = 0.5,
                            min.cutoff = 0.4)
    testthat::expect_true(ggplot2::is_ggplot(p))
    
    p <- SCpubr::do_DotPlot(sample = sample,
                            features = "EPC1",
                            zscore.data = TRUE)
    testthat::expect_true(ggplot2::is_ggplot(p))
    
    p <- SCpubr::do_DotPlot(sample = sample,
                            features = "EPC1",
                            split.by = "annotation",
                            flip = TRUE)
    testthat::expect_true(ggplot2::is_ggplot(p))
    
    p <- SCpubr::do_DotPlot(sample = sample,
                            features = "EPC1",
                            split.by = "annotation",
                            flip = FALSE)
    testthat::expect_true(ggplot2::is_ggplot(p))
    
    
    p <- SCpubr::do_DotPlot(sample = sample,
                            features = "EPC1",
                            zscore.data = TRUE)
    testthat::expect_true(ggplot2::is_ggplot(p))
    
    p <- SCpubr::do_DotPlot(sample = sample,
                            features = list("A" = "EPC1"),
                            zscore.data = TRUE)
    testthat::expect_true(ggplot2::is_ggplot(p))
    
    p <- SCpubr::do_DotPlot(sample = sample,
                            features = list("A" = "EPC1"),
                            zscore.data = TRUE)
    testthat::expect_true(ggplot2::is_ggplot(p))

    p <- SCpubr::do_DotPlot(sample = sample,
                            features = "EPC1",
                            legend.position = "right")
    testthat::expect_true(ggplot2::is_ggplot(p))
    
    p <- SCpubr::do_DotPlot(sample = sample,
                            features = "EPC1",
                            legend.position = "top")
    testthat::expect_true(ggplot2::is_ggplot(p))
    
    p <- SCpubr::do_DotPlot(sample = sample,
                            features = "EPC1",
                            legend.position = "none")
    testthat::expect_true(ggplot2::is_ggplot(p))
    
    
    p <- SCpubr::do_DotPlot(sample = sample,
                            features = "EPC1",
                            legend.position = "none",
                            use_viridis = FALSE)
    testthat::expect_true(ggplot2::is_ggplot(p))
    
    
    p <- SCpubr::do_DotPlot(sample = sample,
                            features = "EPC1",
                            legend.position = "none",
                            use_viridis = FALSE)
    testthat::expect_true(ggplot2::is_ggplot(p))
  })

  testthat::test_that("do_DotPlot: PASS - plot grid", {
    testthat::skip_on_cran()


    p <- SCpubr::do_DotPlot(sample = sample,
                            features = "EPC1",
                            plot.grid = TRUE)
    testthat::expect_true(ggplot2::is_ggplot(p))

    p <- SCpubr::do_DotPlot(sample = sample,
                            features = "EPC1",
                            plot.grid = FALSE)
    testthat::expect_true(ggplot2::is_ggplot(p))
  })

  testthat::test_that("do_DotPlot: PASS - use_viridis", {
    testthat::skip_on_cran()


    p <- SCpubr::do_DotPlot(sample = sample,
                            features = "EPC1",
                            use_viridis = TRUE)
    testthat::expect_true(ggplot2::is_ggplot(p))
  })

  testthat::test_that("do_DotPlot: PASS - one variable legend normal", {
    testthat::skip_on_cran()


    p <- SCpubr::do_DotPlot(sample = sample,
                            features = "EPC1",
                            legend.type = "normal")
    testthat::expect_true(ggplot2::is_ggplot(p))
  })

  testthat::test_that("do_DotPlot: PASS - one variable legend colorbar", {
    testthat::skip_on_cran()


    p <- SCpubr::do_DotPlot(sample = sample,
                            features = "EPC1",
                            legend.type = "colorbar")
    testthat::expect_true(ggplot2::is_ggplot(p))
  })





  testthat::test_that("do_DotPlot: FAIL - wrong legend type", {
    testthat::skip_on_cran()


    testthat::expect_error(suppressWarnings({SCpubr::do_DotPlot(sample = sample,
                                                                features = "EPC1",
                                                                flip = TRUE,
                                                                legend.type = "wrong")}))

  })

  testthat::test_that("do_DotPlot: FAIL - wrong legend position", {
    testthat::skip_on_cran()


    testthat::expect_error(suppressWarnings({SCpubr::do_DotPlot(sample = sample,
                                                                features = "EPC1",
                                                                flip = TRUE,
                                                                legend.position = "wrong")}))

  })

  testthat::test_that("do_DotPlot: FAIL - wrong font.type", {
    testthat::skip_on_cran()


    testthat::expect_error(suppressWarnings({SCpubr::do_DotPlot(sample = sample,
                                                                features = "EPC1",
                                                                flip = TRUE,
                                                                font.type = "wrong")}))

  })

  testthat::test_that("do_DotPlot: PASS - one variable flip", {
    testthat::skip_on_cran()


    p <- SCpubr::do_DotPlot(sample = sample,
                            features = "EPC1",
                            flip = TRUE)
    testthat::expect_true(ggplot2::is_ggplot(p))
  })

  testthat::test_that("do_DotPlot: PASS - multiple features", {
    testthat::skip_on_cran()
    
    genes <- Seurat::VariableFeatures(sample)[1:10]

    p <- suppressWarnings({SCpubr::do_DotPlot(sample = sample,
                                              features = genes)})
    testthat::expect_true(ggplot2::is_ggplot(p))
  })

  testthat::test_that("do_DotPlot: PASS - multiple features flip", {
    testthat::skip_on_cran()


    genes <- Seurat::VariableFeatures(sample)[1:10]
    
    p <- suppressWarnings({SCpubr::do_DotPlot(sample = sample,
                                              features = genes,
                                              flip = TRUE)})
    testthat::expect_true(ggplot2::is_ggplot(p))
  })

  testthat::test_that("do_DotPlot: PASS - multiple features flip rotate x labels", {
    testthat::skip_on_cran()


    genes <- Seurat::VariableFeatures(sample)[1:10]
    p <- suppressWarnings({SCpubr::do_DotPlot(sample = sample,
                                              features = genes,
                                              flip = TRUE,
                                              axis.text.x.angle = 45)})
    testthat::expect_true(ggplot2::is_ggplot(p))
  })





  testthat::test_that("do_DotPlot: PASS - one variable xlab, ylab, title, subtitle, caption", {
    testthat::skip_on_cran()


    p <- SCpubr::do_DotPlot(sample = sample,
                            features = "EPC1",
                            xlab = "A",
                            ylab = "B",
                            plot.title = "C",
                            plot.subtitle = "D",
                            plot.caption = "E")
    testthat::expect_true(ggplot2::is_ggplot(p))
  })
}

