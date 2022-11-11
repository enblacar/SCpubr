if (isFALSE(dep_check[["do_AzimuthAnalysisPlot"]])){

  testthat::test_that("do_AzimuthAnalysisPlot: CRAN essential tests", {
    sample@reductions$ref.umap <- sample@reductions$umap
    p <- SCpubr::do_AzimuthAnalysisPlot(sample,
                                        annotation.labels = "annotation",
                                        annotation.scoring = "annotation.score")

    testthat::expect_type(p, "list")
  })

  testthat::test_that("do_AzimuthAnalysisPlot: PASS", {
    testthat::skip_on_cran()
    sample@reductions$ref.umap <- sample@reductions$umap
    p <- SCpubr::do_AzimuthAnalysisPlot(sample,
                                        annotation.labels = "annotation",
                                        annotation.scoring = "annotation.score",
                                        ref.obj = sample)

    testthat::expect_type(p, "list")

    p <- SCpubr::do_AzimuthAnalysisPlot(sample,
                                        annotation.labels = "annotation",
                                        annotation.scoring = "annotation.score",
                                        ref.obj = sample,
                                        label = FALSE)

    testthat::expect_type(p, "list")

    p <- SCpubr::do_AzimuthAnalysisPlot(sample,
                                        annotation.labels = "annotation",
                                        annotation.scoring = "annotation.score",
                                        ref.obj = sample,
                                        label = FALSE,
                                        legend.position = "right")

    testthat::expect_type(p, "list")

  })


}
