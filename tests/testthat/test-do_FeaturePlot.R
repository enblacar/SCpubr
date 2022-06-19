sample <- use_dataset()

testthat::test_that("do_FeaturePlot: PASS - single feature", {
  p <- SCpubr::do_FeaturePlot(sample = sample,
                              features = "nCount_RNA")
  testthat::expect_type(p, "list")
})

testthat::test_that("do_FeaturePlot: PASS - multiple features", {
  p <- SCpubr::do_FeaturePlot(sample = sample,
                              features = c("nCount_RNA", "nFeature_RNA"))
  testthat::expect_type(p, "list")
})

testthat::test_that("do_FeaturePlot: PASS - title", {
  p <- SCpubr::do_FeaturePlot(sample = sample,
                              features = "nCount_RNA",
                              plot.title = "Feature Plot")
  testthat::expect_type(p, "list")
})

testthat::test_that("do_FeaturePlot: PASS - subtitle", {
  p <- SCpubr::do_FeaturePlot(sample = sample,
                              features = "nCount_RNA",
                              plot.subtitle = "Feature Plot")
  testthat::expect_type(p, "list")
})

testthat::test_that("do_FeaturePlot: PASS - caption", {
  p <- SCpubr::do_FeaturePlot(sample = sample,
                              features = "nCount_RNA",
                              plot.caption = "Feature Plot")
  testthat::expect_type(p, "list")
})

testthat::test_that("do_FeaturePlot: PASS - individual titles", {
  p <- SCpubr::do_FeaturePlot(sample = sample,
                              features = c("nCount_RNA", "nFeature_RNA"),
                              individual.titles = c("A", NA))
  testthat::expect_type(p, "list")
})

testthat::test_that("do_FeaturePlot: PASS - individual subtitles ", {
  p <- SCpubr::do_FeaturePlot(sample = sample,
                              features = c("nCount_RNA", "nFeature_RNA"),
                              individual.subtitles = c("A", NA))
  testthat::expect_type(p, "list")
})

testthat::test_that("do_FeaturePlot: PASS - individual captions", {
  p <- SCpubr::do_FeaturePlot(sample = sample,
                              features = c("nCount_RNA", "nFeature_RNA"),
                              individual.captions = c("A", NA))
  testthat::expect_type(p, "list")
})

testthat::test_that("do_FeaturePlot: PASS - dims", {
  p <- SCpubr::do_FeaturePlot(sample = sample,
                              features = "nCount_RNA",
                              dims = c(2, 1))
  testthat::expect_type(p, "list")
})

testthat::test_that("do_FeaturePlot: PASS - subset of cells", {
  cells.plot <- colnames(sample[, !(sample$seurat_clusters %in% c("2", "5", "8"))])
  p <- SCpubr::do_FeaturePlot(sample,
                              features = c("CD14"),
                              cells.highlight = cells.plot)
  testthat::expect_type(p, "list")
})

testthat::test_that("do_FeaturePlot: PASS - subset of identities", {
  p <- SCpubr::do_FeaturePlot(sample,
                              features = c("CD14"),
                              idents.highlight = c("1", "2"))
  testthat::expect_type(p, "list")
})

testthat::test_that("do_FeaturePlot: PASS - subset of cells and identities", {
  cells.plot <- colnames(sample[, !(sample$seurat_clusters %in% c("2", "5", "8"))])
  p <- SCpubr::do_FeaturePlot(sample,
                              features = c("CD14"),
                              cells.highlight = cells.plot,
                              idents.highlight = c("1", "2"))
  testthat::expect_type(p, "list")
})

testthat::test_that("do_FeaturePlot: PASS - split.by", {
  p <- SCpubr::do_FeaturePlot(sample,
                              features = c("CD14"),
                              split.by = "seurat_clusters")
  testthat::expect_type(p, "list")
})

testthat::test_that("do_FeaturePlot: PASS - split.by and split.by.idents", {
  p <- SCpubr::do_FeaturePlot(sample,
                              features = c("CD14"),
                              split.by = "seurat_clusters",
                              split.by.idents = c("1", "2"))
  testthat::expect_type(p, "list")
})

testthat::test_that("do_FeaturePlot: PASS - split.by and split.by.idents multiple features", {
  p <- SCpubr::do_FeaturePlot(sample,
                              features = c("CD14", "nCount_RNA"),
                              split.by = "seurat_clusters",
                              split.by.idents = c("1", "2"))
  testthat::expect_type(p, "list")
})

testthat::test_that("do_FeaturePlot: PASS - modify color maps", {
  p <- SCpubr::do_FeaturePlot(sample,
                              features = c("nCount_RNA"),
                              viridis_color_map = "F")
  testthat::expect_type(p, "list")
})
