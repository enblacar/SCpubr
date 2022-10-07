
testthat::test_that("do_BoxPlot: PASS - default", {
  sample <- SCpubr:::test_list$sample
  p <- SCpubr::do_BoxPlot(sample = sample,
                          feature = "nCount_RNA")
  testthat::expect_type(p, "list")

  p <- SCpubr::do_BoxPlot(sample = sample,
                          feature = "nCount_RNA",
                          split.by = "orig.ident")
  testthat::expect_type(p, "list")
})

testthat::test_that("do_BoxPlot: PASS - custom_grouping", {
  sample <- SCpubr:::test_list$sample
  p <- SCpubr::do_BoxPlot(sample = sample,
                          feature = "nCount_RNA",
                          group.by = "orig.ident")
  testthat::expect_type(p, "list")

  sample$orig.ident <- factor(sample$orig.ident)
  p <- SCpubr::do_BoxPlot(sample = sample,
                          feature = "nCount_RNA",
                          group.by = "orig.ident")
  testthat::expect_type(p, "list")

  p <- SCpubr::do_BoxPlot(sample = sample,
                          feature = "nCount_RNA",
                          group.by = "orig.ident",
                          colors.use = c("Cell" = "blue"))
  testthat::expect_type(p, "list")
})


testthat::test_that("do_BoxPlot: PASS - split.by", {
  sample <- SCpubr:::test_list$sample
  sample$orig.ident <- ifelse(sample$seurat_clusters == "0", "C", "B")
  p <- SCpubr::do_BoxPlot(sample = sample,
                          feature = "nCount_RNA",
                          split.by = "orig.ident")
  testthat::expect_type(p, "list")
})


testthat::test_that("do_BoxPlot: PASS - silhouette", {
  sample <- SCpubr:::test_list$sample
  p <- SCpubr::do_BoxPlot(sample = sample,
                          feature = "nCount_RNA",
                          use_silhouette = TRUE)
  testthat::expect_type(p, "list")
})

testthat::test_that("do_BoxPlot: PASS - silhouette", {
  sample <- SCpubr:::test_list$sample
  p <- SCpubr::do_BoxPlot(sample = sample,
                          feature = "nCount_RNA",
                          use_test = TRUE,
                          comparisons = list(c("0", "1")))
  testthat::expect_type(p, "list")
})

testthat::test_that("do_BoxPlot: PASS - order", {
  sample <- SCpubr:::test_list$sample
  p <- SCpubr::do_BoxPlot(sample = sample,
                          feature = "nCount_RNA",
                          order = TRUE)
  testthat::expect_type(p, "list")

  p <- SCpubr::do_BoxPlot(sample = sample,
                          feature = "nCount_RNA",
                          order = TRUE,
                          flip = TRUE)
  testthat::expect_type(p, "list")
})

testthat::test_that("do_BoxPlot: PASS - flip", {
  sample <- SCpubr:::test_list$sample
  p <- SCpubr::do_BoxPlot(sample = sample,
                          feature = "nCount_RNA",
                          flip = TRUE)
  testthat::expect_type(p, "list")
})


testthat::test_that("do_BoxPlot: FAILS ", {
  sample <- SCpubr:::test_list$sample
  testthat::expect_error({SCpubr::do_BoxPlot(sample = sample,
                                             feature = "nCount_RNA",
                                             use_test = TRUE,
                                             split.by = "orig.ident")})

  testthat::expect_error({SCpubr::do_BoxPlot(sample = sample,
                                             feature = "nCount_RNA",
                                             use_silhouette = TRUE,
                                             split.by = "orig.ident")})

  testthat::expect_error({SCpubr::do_BoxPlot(sample = sample,
                                             feature = "nCount_RNA",
                                             use_test = TRUE)})

  testthat::expect_error({SCpubr::do_BoxPlot(sample = sample,
                                             feature = "nCount_RNA",
                                             order = TRUE,
                                             split.by = "orig.ident")})
})

