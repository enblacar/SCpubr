sample <- use_dataset()
testthat::test_that("do_NebulosaPlot: PASS - single feature", {
  p <- SCpubr::do_NebulosaPlot(sample = sample,
                               features = c("CD14"))
  testthat::expect_type(p, "list")
})

testthat::test_that("do_NebulosaPlot: PASS - single feature legend normal", {
  p <- SCpubr::do_NebulosaPlot(sample = sample,
                               features = c("CD14"),
                               legend.type = "normal")
  testthat::expect_type(p, "list")
})

testthat::test_that("do_NebulosaPlot: PASS - single feature legend colorbar", {
  p <- SCpubr::do_NebulosaPlot(sample = sample,
                               features = c("CD14"),
                               legend.type = "colorbar")
  testthat::expect_type(p, "list")
})

testthat::test_that("do_NebulosaPlot: PASS - single feature legend colorsteps", {
  p <- SCpubr::do_NebulosaPlot(sample = sample,
                               features = c("CD14"),
                               legend.type = "colorsteps")
  testthat::expect_type(p, "list")
})

testthat::test_that("do_NebulosaPlot: FAIL - wrong legend type ", {
  testthat::expect_error(SCpubr::do_NebulosaPlot(sample = sample,
                                                 features = c("CD14"),
                                                 legend.type = "wrong"))
})

testthat::test_that("do_NebulosaPlot: FAIL - wrong legend position ", {
  testthat::expect_error(SCpubr::do_NebulosaPlot(sample = sample,
                                                 features = c("CD14"),
                                                 legend.position = "wrong"))
})

testthat::test_that("do_NebulosaPlot: FAIL - wrong font.type", {
  testthat::expect_error(SCpubr::do_NebulosaPlot(sample = sample,
                                                 features = c("CD14"),
                                                 font.type = "wrong"))
})

testthat::test_that("do_NebulosaPlot: PASS - single feature distinct dims", {
  p <- SCpubr::do_NebulosaPlot(sample = sample,
                               features = c("CD14"),
                               dims = c(2, 1))
  testthat::expect_type(p, "list")
})

testthat::test_that("do_FeaturePlot: PASS - diffusion", {
  test <- sample@reductions$umap[[]]
  colnames(test) <- c("DC_1", "DC_2")
  obj <- Seurat::CreateDimReducObject(test, assay = "SCT", key = "DC_")
  sample@reductions$diffusion <- obj
  p <- suppressWarnings(SCpubr::do_NebulosaPlot(sample,
                               features = c("PC_1"),
                               reduction = "diffusion"))
  testthat::expect_type(p, "list")
})

testthat::test_that("do_NebulosaPlot: PASS - several", {
  p <- SCpubr::do_NebulosaPlot(sample = sample,
                               features = c("CD14", "CD8A"))
  testthat::expect_type(p, "list")
})

testthat::test_that("do_NebulosaPlot: PASS - several, joint", {
  p <- SCpubr::do_NebulosaPlot(sample = sample,
                               features = c("CD14", "CD8A"),
                               joint = T)
  testthat::expect_type(p, "list")
})

testthat::test_that("do_NebulosaPlot: PASS - several, joint only joint", {
  p <- SCpubr::do_NebulosaPlot(sample = sample,
                               features = c("CD14", "CD8A"),
                               joint = T,
                               return_only_joint = T)
  testthat::expect_type(p, "list")
})

testthat::test_that("do_NebulosaPlot: PASS - title", {
  p <- SCpubr::do_NebulosaPlot(sample = sample,
                               features = c("CD14", "CD8A"),
                               joint = T,
                               return_only_joint = T,
                               plot.title = "Title")
  testthat::expect_type(p, "list")
})

testthat::test_that("do_NebulosaPlot: PASS - subtitle", {
  p <- SCpubr::do_NebulosaPlot(sample = sample,
                               features = c("CD14", "CD8A"),
                               joint = T,
                               return_only_joint = T,
                               plot.subtitle = "Subtitle")
  testthat::expect_type(p, "list")
})

testthat::test_that("do_NebulosaPlot: PASS - caption", {
  p <- SCpubr::do_NebulosaPlot(sample = sample,
                               features = c("CD14", "CD8A"),
                               joint = T,
                               return_only_joint = T,
                               plot.caption = "Caption")
  testthat::expect_type(p, "list")
})

testthat::test_that("do_NebulosaPlot: PASS - individual titles", {
  p <- SCpubr::do_NebulosaPlot(sample = sample,
                               features = c("CD14", "CD8A"),
                               joint = T,
                               return_only_joint = F,
                               individual.titles = c("A", NA, "B"))
  testthat::expect_type(p, "list")
})

testthat::test_that("do_NebulosaPlot: PASS - individual subtitles", {
  p <- SCpubr::do_NebulosaPlot(sample = sample,
                               features = c("CD14", "CD8A"),
                               joint = T,
                               return_only_joint = F,
                               individual.subtitles = c("A", NA, "B"))
  testthat::expect_type(p, "list")
})

testthat::test_that("do_NebulosaPlot: PASS - individual captions", {
  p <- SCpubr::do_NebulosaPlot(sample = sample,
                               features = c("CD14", "CD8A"),
                               joint = T,
                               return_only_joint = F,
                               individual.captions = c("A", NA, "B"))
  testthat::expect_type(p, "list")
})

testthat::test_that("do_NebulosaPlot: PASS - color map", {
  p <- SCpubr::do_NebulosaPlot(sample = sample,
                               features = c("CD14"),
                               viridis_color_map = "F")
  testthat::expect_type(p, "list")
})

testthat::test_that("do_NebulosaPlot: PASS - legend top", {
  p <- SCpubr::do_NebulosaPlot(sample = sample,
                               features = c("CD14"),
                               legend.position = "left")
  testthat::expect_type(p, "list")
})

testthat::test_that("do_NebulosaPlot: PASS - legend top", {
  p <- SCpubr::do_NebulosaPlot(sample = sample,
                               features = c("CD14"),
                               legend.position = "top")
  testthat::expect_type(p, "list")
})

testthat::test_that("do_NebulosaPlot: WARNING - features as list", {
  testthat::expect_warning(SCpubr::do_NebulosaPlot(sample = sample,
                                                   features = list("CD14"),
                                                   viridis_color_map = "F"))
})

testthat::test_that("do_NebulosaPlot: FAIL - individual titles + joint + return only joint", {
  testthat::expect_error(SCpubr::do_NebulosaPlot(sample = sample,
                                                 features = c("CD14", "CD8A"),
                                                 joint = T,
                                                 return_only_joint = T,
                                                 individual.titles = "A"))
})

testthat::test_that("do_NebulosaPlot: FAIL - not enough individual titles ", {
  testthat::expect_error(SCpubr::do_NebulosaPlot(sample = sample,
                                                 features = c("CD14", "CD8A"),
                                                 joint = F,
                                                 individual.titles = "A"))
})

testthat::test_that("do_NebulosaPlot: FAIL - not enough individual titles for joint ", {
  testthat::expect_error(SCpubr::do_NebulosaPlot(sample = sample,
                                                 features = c("CD14", "CD8A"),
                                                 joint = T,
                                                 individual.titles = "A"))
})

testthat::test_that("do_NebulosaPlot: PASS - no legend", {
  p <- SCpubr::do_NebulosaPlot(sample = sample,
                               features = c("CD14"),
                               legend = F)
  testthat::expect_type(p, "list")
})

testthat::test_that("do_NebulosaPlot: PASS - patchwork title, subtitle and caption", {
  p <- SCpubr::do_NebulosaPlot(sample = sample,
                               features = c("CD14", "CD8A"),
                               plot.title = "A",
                               plot.subtitle = "B",
                               plot.caption = "C")
  testthat::expect_type(p, "list")
})
