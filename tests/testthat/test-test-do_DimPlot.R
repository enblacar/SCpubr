suppressWarnings(sample <- SCpubr:::use_dataset())

testthat::test_that("do_DimPlot: PASS - sample", {
  p <- SCpubr::do_DimPlot(sample = sample)
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
