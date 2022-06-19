sample <- use_dataset()
testthat::test_that("do_DimPlot: PASS - sample", {
  p <- SCpubr::do_DimPlot(sample = sample)
  testthat::expect_type(p, "list")
})

testthat::test_that("do_DimPlot: PASS - title", {
  p <- SCpubr::do_DimPlot(sample = sample,
                          plot.title = "My awesome SC data set")
  testthat::expect_type(p, "list")
})

testthat::test_that("do_DimPlot: PASS - subtitle", {
  p <- SCpubr::do_DimPlot(sample = sample,
                          plot.subtitle = "My awesome SC data set")
  testthat::expect_type(p, "list")
})

testthat::test_that("do_DimPlot: PASS - caption", {
  p <- SCpubr::do_DimPlot(sample = sample,
                          plot.caption = "My awesome SC data set")
  testthat::expect_type(p, "list")
})

testthat::test_that("do_DimPlot: PASS - sample + group.by", {
  p <- SCpubr::do_DimPlot(sample = sample, group.by = "seurat_clusters")
  testthat::expect_type(p, "list")
})

testthat::test_that("do_DimPlot: PASS - sample + split.by", {
  p <- SCpubr::do_DimPlot(sample = sample, split.by = "seurat_clusters")
  testthat::expect_type(p, "list")
})

testthat::test_that("do_DimPlot: PASS - sample + split.by + idents.keep", {
  p <- SCpubr::do_DimPlot(sample = sample, split.by = "seurat_clusters", idents.keep = c("1", "3", "5"))
  testthat::expect_type(p, "list")
})

testthat::test_that("do_DimPlot: PASS - dims", {
  p <- SCpubr::do_DimPlot(sample = sample, dims = c(1, 2))
  testthat::expect_type(p, "list")
})

testthat::test_that("do_DimPlot: PASS - legend.position", {
  p <- SCpubr::do_DimPlot(sample = sample, legend.position = "top")
  testthat::expect_type(p, "list")
})

testthat::test_that("do_DimPlot: PASS - legend.ncol", {
  p <- SCpubr::do_DimPlot(sample = sample, legend.ncol = 2)
  testthat::expect_type(p, "list")
})

testthat::test_that("do_DimPlot: PASS - legend.nrow", {
  p <- SCpubr::do_DimPlot(sample = sample, legend.nrow = 2)
  testthat::expect_type(p, "list")
})

testthat::test_that("do_DimPlot: PASS - label", {
  p <- SCpubr::do_DimPlot(sample = sample, label = T)
  testthat::expect_type(p, "list")
})

testthat::test_that("do_DimPlot: PASS - order", {
  p <- SCpubr::do_DimPlot(sample = sample, order = "5", shuffle = F)
  testthat::expect_type(p, "list")
})

testthat::test_that("do_DimPlot: PASS - colors.use", {
  p <- SCpubr::do_DimPlot(sample = sample, colors.use = c("0" = "#001219",
                                                          "1" = "#005f73",
                                                          "2" = "#0a9396",
                                                          "3" = "#94d2bd",
                                                          "4" = "#e9d8a6",
                                                          "5" = "#ee9b00",
                                                          "6" = "#ca6702",
                                                          "7" = "#bb3e03",
                                                          "8" = "#ae2012"))
  testthat::expect_type(p, "list")
})

testthat::test_that("do_DimPlot: PASS - cells.highlight", {
  p <- SCpubr::do_DimPlot(sample = sample, cells.highlight = sample(colnames(sample), 50))
  testthat::expect_type(p, "list")
})

testthat::test_that("do_DimPlot: PASS - idents.highlight", {
  p <- SCpubr::do_DimPlot(sample = sample, idents.highlight = "5")
  testthat::expect_type(p, "list")
})

testthat::test_that("do_DimPlot: PASS - cells.highlight and idents.highlight", {
  p <- SCpubr::do_DimPlot(sample = sample, cells.highlight = sample(colnames(sample), 50), idents.highlight = "2")
  testthat::expect_type(p, "list")
})

testthat::test_that("do_DimPlot: PASS - idents.keep", {
  p <- SCpubr::do_DimPlot(sample = sample, idents.keep = "5")
  testthat::expect_type(p, "list")
})
