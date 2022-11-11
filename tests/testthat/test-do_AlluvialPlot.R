if (isFALSE(dep_check[["do_AlluvialPlot"]])){

  testthat::test_that("do_AlluvialPlot: CRAN essential tests", {

    p <- SCpubr::do_AlluvialPlot(sample,
                                 first_group = "orig.ident",
                                 last_group = "seurat_clusters")

    testthat::expect_type(p, "list")
  })
}


testthat::test_that("do_EnrichmentHeatmap: PASS - normal", {
  testthat::skip_on_cran()

  p <- SCpubr::do_AlluvialPlot(sample,
                               first_group = "orig.ident",
                               middle_groups = "annotation",
                               last_group = "seurat_clusters")

  testthat::expect_type(p, "list")

  p <- SCpubr::do_AlluvialPlot(sample,
                               first_group = "orig.ident",
                               middle_groups = "annotation",
                               last_group = "seurat_clusters",
                               flip = TRUE)

  testthat::expect_type(p, "list")
})
