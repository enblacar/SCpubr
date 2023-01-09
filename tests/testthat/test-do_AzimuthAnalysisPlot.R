if (isFALSE(dep_check[["do_AzimuthAnalysisPlot"]])){

  testthat::test_that("do_AzimuthAnalysisPlot: CRAN essential tests", {
    p <- SCpubr::do_AzimuthAnalysisPlot(sample,
                                        annotation.labels = "annotation",
                                        annotation.scoring = "annotation.score")

    testthat::expect_type(p, "list")
  })

  testthat::test_that("do_AzimuthAnalysisPlot: PASS", {
    testthat::skip_on_cran()
    p <- SCpubr::do_AzimuthAnalysisPlot(sample,
                                        annotation.labels = "annotation",
                                        annotation.scoring = "annotation.score",
                                        ref.obj = sample,
                                        raster = TRUE)

    testthat::expect_type(p, "list")

    p <- SCpubr::do_AzimuthAnalysisPlot(sample,
                                        annotation.labels = "annotation",
                                        annotation.scoring = "annotation.score",
                                        ref.obj = sample,
                                        label = FALSE,
                                        raster = TRUE)

    testthat::expect_type(p, "list")

    p <- SCpubr::do_AzimuthAnalysisPlot(sample,
                                        annotation.labels = "annotation",
                                        annotation.scoring = "annotation.score",
                                        ref.obj = sample,
                                        label = FALSE,
                                        legend.position = "right",
                                        raster = TRUE)

    testthat::expect_type(p, "list")

    p <- SCpubr::do_AzimuthAnalysisPlot(sample,
                                        annotation.labels = "annotation",
                                        annotation.scoring = "annotation.score",
                                        ref.obj = sample,
                                        raster = FALSE)

    testthat::expect_type(p, "list")

    p <- SCpubr::do_AzimuthAnalysisPlot(sample,
                                        annotation.labels = "annotation",
                                        annotation.scoring = "annotation.score",
                                        ref.obj = sample,
                                        label = FALSE,
                                        raster = FALSE)

    testthat::expect_type(p, "list")

    p <- SCpubr::do_AzimuthAnalysisPlot(sample,
                                        annotation.labels = "annotation",
                                        annotation.scoring = "annotation.score",
                                        ref.obj = sample,
                                        label = FALSE,
                                        legend.position = "right",
                                        raster = FALSE)

    testthat::expect_type(p, "list")

    p <- SCpubr::do_AzimuthAnalysisPlot(sample,
                                        annotation.labels = "annotation",
                                        annotation.scoring = "annotation.score",
                                        ref.obj = sample,
                                        label = FALSE,
                                        legend.position = "right",
                                        raster = FALSE,
                                        colors.use = SCpubr:::generate_color_scale(names_use = unique(sample$seurat_clusters)))

    testthat::expect_type(p, "list")

  })


}
