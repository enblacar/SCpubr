sample <- SCpubr:::test_list$sample

testthat::test_that("do_ChordDiagramPlot: PASS - default", {
  p <- SCpubr::do_ChordDiagramPlot(sample = sample,
                                   from = "seurat_clusters",
                                   to = "orig.ident")
  testthat::expect_type(p, "list")
})

testthat::test_that("do_ChordDiagramPlot: PASS - link border color", {
  p <- SCpubr::do_ChordDiagramPlot(sample = sample,
                                   from = "seurat_clusters",
                                   to = "orig.ident",
                                   link.border.color = "black")
  testthat::expect_type(p, "list")
})

testthat::test_that("do_ChordDiagramPlot: PASS - alignment", {
  p <- SCpubr::do_ChordDiagramPlot(sample = sample,
                                   from = "seurat_clusters",
                                   to = "orig.ident",
                                   alignment = "vertical")
  testthat::expect_type(p, "list")

  p <- SCpubr::do_ChordDiagramPlot(sample = sample,
                                   from = "seurat_clusters",
                                   to = "orig.ident",
                                   alignment = "horizontal")
  testthat::expect_type(p, "list")
})

testthat::test_that("do_ChordDiagramPlot: PASS - highlight group", {
  p <- SCpubr::do_ChordDiagramPlot(sample = sample,
                                   from = "seurat_clusters",
                                   to = "orig.ident",
                                   highlight_group = "0")
  testthat::expect_type(p, "list")
})

testthat::test_that("do_ChordDiagramPlot: FAILS", {
  testthat::expect_error({SCpubr::do_ChordDiagramPlot(sample = sample,
                                                      from = "seurat_clusters",
                                                      to = "orig.ident",
                                                      alignment = "wrong")})

  testthat::expect_error({SCpubr::do_ChordDiagramPlot(sample = sample,
                                                      from = "seurat_clusters",
                                                      to = "orig.ident",
                                                      link.arr.type = "wrong")})

  testthat::expect_error({SCpubr::do_ChordDiagramPlot(sample = sample,
                                                      from = "seurat_clusters",
                                                      to = "orig.ident",
                                                      highlight_group = "0",
                                                      alpha.highlight = 120)})
})
