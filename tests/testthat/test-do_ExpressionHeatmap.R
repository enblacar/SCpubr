if (isFALSE(dep_check[["do_ExpressionHeatmap"]])){

  testthat::test_that("do_ExpressionHeatmap: CRAN essential tests", {

    p <- SCpubr::do_ExpressionHeatmap(sample,
                                      features = rownames(sample)[1:5])

    testthat::expect_true("HeatmapList" %in% class(p))
  })

  testthat::test_that("do_ExpressionHeatmap: PASS - normal", {
    testthat::skip_on_cran()

    p <- SCpubr::do_ExpressionHeatmap(sample,
                                      features = rownames(sample)[1:5],
                                      group.by = c("orig.ident"),
                                      flip = TRUE)

    testthat::expect_true("HeatmapList" %in% class(p))

    p <- SCpubr::do_ExpressionHeatmap(sample,
                                      features = rownames(sample)[1:5],
                                      group.by = c("orig.ident"),
                                      flip = FALSE)

    testthat::expect_true("HeatmapList" %in% class(p))


    p <- SCpubr::do_ExpressionHeatmap(sample,
                                      features = rownames(sample)[1:5],
                                      group.by = c("orig.ident"),
                                      row_title = "A",
                                      column_title = "B",
                                      flip = TRUE)

    testthat::expect_true("HeatmapList" %in% class(p))

    p <- SCpubr::do_ExpressionHeatmap(sample,
                                      features = rownames(sample)[1:5],
                                      group.by = c("orig.ident"),
                                      row_title = "A",
                                      column_title = "B",
                                      flip = FALSE)

    testthat::expect_true("HeatmapList" %in% class(p))

    p <- SCpubr::do_ExpressionHeatmap(sample,
                                      features = rownames(sample)[1:5],
                                      group.by = c("orig.ident", "seurat_clusters"),
                                      flip = TRUE)

    testthat::expect_true("HeatmapList" %in% class(p))

    p <- SCpubr::do_ExpressionHeatmap(sample,
                                      features = rownames(sample)[1:5],
                                      group.by = c("orig.ident", "seurat_clusters"),
                                      flip = FALSE)

    testthat::expect_true("HeatmapList" %in% class(p))

    p <- SCpubr::do_ExpressionHeatmap(sample,
                                      features = rownames(sample)[1:5],
                                      group.by = c("orig.ident", "seurat_clusters"),
                                      column_title = c("A", "B"),
                                      row_title = c("C", "D"),
                                      flip = FALSE)

    p <- SCpubr::do_ExpressionHeatmap(sample,
                                      features = rownames(sample)[1:5],
                                      group.by = c("orig.ident", "seurat_clusters"),
                                      column_title = c("A", "B"),
                                      row_title = c("C", "D"),
                                      flip = TRUE)

    testthat::expect_true("HeatmapList" %in% class(p))
  })

  testthat::test_that("do_ExpressionHeatmap: PASS - assay", {
    testthat::skip_on_cran()

    p <- SCpubr::do_ExpressionHeatmap(sample,
                                      features = rownames(sample)[1:5],
                                      group.by = c("orig.ident"),
                                      assay = NULL)

    testthat::expect_true("HeatmapList" %in% class(p))

    p <- SCpubr::do_ExpressionHeatmap(sample,
                                      features = rownames(sample)[1:5],
                                      group.by = c("orig.ident"),
                                      assay = "SCT")

    testthat::expect_true("HeatmapList" %in% class(p))
  })

  testthat::test_that("do_ExpressionHeatmap: PASS - legend.position", {
    testthat::skip_on_cran()

    p <- SCpubr::do_ExpressionHeatmap(sample,
                                      features = rownames(sample)[1:5],
                                      group.by = c("orig.ident"),
                                      legend.position = "bottom")

    testthat::expect_true("HeatmapList" %in% class(p))

    p <- SCpubr::do_ExpressionHeatmap(sample,
                                      features = rownames(sample)[1:5],
                                      group.by = c("orig.ident"),
                                      legend.position = "right")

    testthat::expect_true("HeatmapList" %in% class(p))
  })

  testthat::test_that("do_ExpressionHeatmap: PASS - cutoffs", {
    testthat::skip_on_cran()

    p <- SCpubr::do_ExpressionHeatmap(sample,
                                      features = rownames(sample)[1:5],
                                      group.by = c("orig.ident"),
                                      assay = NULL,
                                      min.cutoff = 0.7)

    testthat::expect_true("HeatmapList" %in% class(p))

    p <- SCpubr::do_ExpressionHeatmap(sample,
                                      features = rownames(sample)[1:5],
                                      group.by = c("orig.ident"),
                                      assay = "SCT",
                                      max.cutoff = 0.72)

    testthat::expect_true("HeatmapList" %in% class(p))

    p <- SCpubr::do_ExpressionHeatmap(sample,
                                      features = rownames(sample)[1:5],
                                      group.by = c("orig.ident"),
                                      assay = "SCT",
                                      min.cutoff = 0.7,
                                      max.cutoff = 0.72)

    testthat::expect_true("HeatmapList" %in% class(p))
  })

  testthat::test_that("do_ExpressionHeatmap: FAIL", {
    testthat::skip_on_cran()
    testthat::expect_error({SCpubr::do_ExpressionHeatmap(sample = sample,
                                                         features = c("EPC1"),
                                                         min.cutoff = -10)})

    testthat::expect_error({SCpubr::do_ExpressionHeatmap(sample = sample,
                                                         features = c("EPC1"),
                                                         max.cutoff = 200)})

    testthat::expect_error({SCpubr::do_ExpressionHeatmap(sample = sample,
                                                         features = c("EPC1"),
                                                         max.cutoff = 1,
                                                         min.cutoff = 2)})

    testthat::expect_message({SCpubr::do_ExpressionHeatmap(sample = sample,
                                                         features = list("A" = c("EPC1")))})
    testthat::expect_warning({SCpubr::do_ExpressionHeatmap(sample = sample,
                                                           features =c("EPC1", "NOTFOUND"))})

  })

}


