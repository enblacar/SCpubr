if (isFALSE(dep_check[["do_ExpressionHeatmap"]])){

  testthat::test_that("do_ExpressionHeatmap: CRAN essential tests", {

    p <- SCpubr::do_ExpressionHeatmap(sample,
                                      features = rownames(sample)[1:5])

    testthat::expect_true("HeatmapList" %in% class(p))
  })
}


testthat::test_that("do_EnrichmentHeatmap: PASS - normal", {
  testthat::skip_on_cran()

  p <- SCpubr::do_ExpressionHeatmap(sample,
                                    features = rownames(sample)[1:5],
                                    group.by = c("orig.ident"))

  testthat::expect_true("HeatmapList" %in% class(p))

  p <- SCpubr::do_ExpressionHeatmap(sample,
                                    features = rownames(sample)[1:5],
                                    group.by = c("orig.ident", "seurat_clusters"))

  testthat::expect_true("HeatmapList" %in% class(p))
})
