
testthat::test_that("do_NebulosaPlot: PASS - single feature", {
  sample <- SCpubr:::test_list$sample

  p <- SCpubr::do_NebulosaPlot(sample = sample,
                               features = c("CD14"))
  testthat::expect_type(p, "list")
})

testthat::test_that("do_NebulosaPlot: PASS - cell_borders", {
  sample <- SCpubr:::test_list$sample

  p <- SCpubr::do_NebulosaPlot(sample = sample, features = "CD14", plot_cell_borders = T)
  testthat::expect_type(p, "list")
  p <- suppressWarnings({SCpubr::do_NebulosaPlot(sample = sample, features = c("CD14", "PC_1"), plot_cell_borders = T)})
  testthat::expect_type(p, "list")
})

testthat::test_that("do_NebulosaPlot: PASS - single feature legend normal", {
  sample <- SCpubr:::test_list$sample

  p <- SCpubr::do_NebulosaPlot(sample = sample,
                               features = c("CD14"),
                               legend.type = "normal")
  testthat::expect_type(p, "list")
})

testthat::test_that("do_NebulosaPlot: PASS - single feature legend colorbar", {
  sample <- SCpubr:::test_list$sample

  p <- SCpubr::do_NebulosaPlot(sample = sample,
                               features = c("CD14"),
                               legend.type = "colorbar")
  testthat::expect_type(p, "list")
})

testthat::test_that("do_NebulosaPlot: PASS - single feature legend colorsteps", {
  sample <- SCpubr:::test_list$sample

  p <- SCpubr::do_NebulosaPlot(sample = sample,
                               features = c("CD14"),
                               legend.type = "colorsteps")
  testthat::expect_type(p, "list")
})

testthat::test_that("do_NebulosaPlot: FAIL - wrong legend type ", {
  sample <- SCpubr:::test_list$sample

  testthat::expect_error(SCpubr::do_NebulosaPlot(sample = sample,
                                                 features = c("CD14"),
                                                 legend.type = "wrong"))
})

testthat::test_that("do_NebulosaPlot: FAIL - wrong legend position ", {
  sample <- SCpubr:::test_list$sample

  testthat::expect_error(SCpubr::do_NebulosaPlot(sample = sample,
                                                 features = c("CD14"),
                                                 legend.position = "wrong"))
})

testthat::test_that("do_NebulosaPlot: FAIL - wrong font.type", {
  sample <- SCpubr:::test_list$sample

  testthat::expect_error(SCpubr::do_NebulosaPlot(sample = sample,
                                                 features = c("CD14"),
                                                 font.type = "wrong"))
})

testthat::test_that("do_NebulosaPlot: PASS - single feature distinct dims", {
  sample <- SCpubr:::test_list$sample

  p <- SCpubr::do_NebulosaPlot(sample = sample,
                               features = c("CD14"),
                               dims = c(2, 1))
  testthat::expect_type(p, "list")
})

testthat::test_that("do_FeaturePlot: PASS - diffusion", {
  sample <- SCpubr:::test_list$sample

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
  sample <- SCpubr:::test_list$sample

  p <- SCpubr::do_NebulosaPlot(sample = sample,
                               features = c("CD14", "CD8A"))
  testthat::expect_type(p, "list")
})

testthat::test_that("do_NebulosaPlot: PASS - several, joint", {
  sample <- SCpubr:::test_list$sample

  p <- SCpubr::do_NebulosaPlot(sample = sample,
                               features = c("CD14", "CD8A"),
                               joint = T)
  testthat::expect_type(p, "list")
})

testthat::test_that("do_NebulosaPlot: PASS - several, joint only joint", {
  sample <- SCpubr:::test_list$sample

  p <- SCpubr::do_NebulosaPlot(sample = sample,
                               features = c("CD14", "CD8A"),
                               joint = T,
                               return_only_joint = T)
  testthat::expect_type(p, "list")
})

testthat::test_that("do_NebulosaPlot: PASS - title", {
  sample <- SCpubr:::test_list$sample

  p <- SCpubr::do_NebulosaPlot(sample = sample,
                               features = c("CD14", "CD8A"),
                               joint = T,
                               return_only_joint = T,
                               plot.title = "Title")
  testthat::expect_type(p, "list")
})

testthat::test_that("do_NebulosaPlot: PASS - subtitle", {
  sample <- SCpubr:::test_list$sample

  p <- SCpubr::do_NebulosaPlot(sample = sample,
                               features = c("CD14", "CD8A"),
                               joint = T,
                               return_only_joint = T,
                               plot.subtitle = "Subtitle")
  testthat::expect_type(p, "list")
})

testthat::test_that("do_NebulosaPlot: PASS - caption", {
  sample <- SCpubr:::test_list$sample

  p <- SCpubr::do_NebulosaPlot(sample = sample,
                               features = c("CD14", "CD8A"),
                               joint = T,
                               return_only_joint = T,
                               plot.caption = "Caption")
  testthat::expect_type(p, "list")
})

testthat::test_that("do_NebulosaPlot: PASS - individual titles", {
  sample <- SCpubr:::test_list$sample

  p <- SCpubr::do_NebulosaPlot(sample = sample,
                               features = c("CD14", "CD8A"),
                               joint = T,
                               return_only_joint = F,
                               individual.titles = c("A", NA, "B"))
  testthat::expect_type(p, "list")
})

testthat::test_that("do_NebulosaPlot: PASS - individual subtitles", {
  sample <- SCpubr:::test_list$sample

  p <- SCpubr::do_NebulosaPlot(sample = sample,
                               features = c("CD14", "CD8A"),
                               joint = T,
                               return_only_joint = F,
                               individual.subtitles = c("A", NA, "B"))
  testthat::expect_type(p, "list")
})

testthat::test_that("do_NebulosaPlot: PASS - individual captions", {
  sample <- SCpubr:::test_list$sample

  p <- SCpubr::do_NebulosaPlot(sample = sample,
                               features = c("CD14", "CD8A"),
                               joint = T,
                               return_only_joint = F,
                               individual.captions = c("A", NA, "B"))
  testthat::expect_type(p, "list")
})

testthat::test_that("do_NebulosaPlot: PASS - color map", {
  sample <- SCpubr:::test_list$sample

  p <- SCpubr::do_NebulosaPlot(sample = sample,
                               features = c("CD14"),
                               viridis_color_map = "F")
  testthat::expect_type(p, "list")
})

testthat::test_that("do_NebulosaPlot: PASS - legend top", {
  sample <- SCpubr:::test_list$sample

  p <- SCpubr::do_NebulosaPlot(sample = sample,
                               features = c("CD14"),
                               legend.position = "left")
  testthat::expect_type(p, "list")
})

testthat::test_that("do_NebulosaPlot: PASS - legend top", {
  sample <- SCpubr:::test_list$sample

  p <- SCpubr::do_NebulosaPlot(sample = sample,
                               features = c("CD14"),
                               legend.position = "top")
  testthat::expect_type(p, "list")
})

testthat::test_that("do_NebulosaPlot: WARNING - features as list", {
  sample <- SCpubr:::test_list$sample

  testthat::expect_warning(SCpubr::do_NebulosaPlot(sample = sample,
                                                   features = list("CD14"),
                                                   viridis_color_map = "F"))
})

testthat::test_that("do_NebulosaPlot: FAIL - individual titles + joint + return only joint", {
  sample <- SCpubr:::test_list$sample

  testthat::expect_error(SCpubr::do_NebulosaPlot(sample = sample,
                                                 features = c("CD14", "CD8A"),
                                                 joint = T,
                                                 return_only_joint = T,
                                                 individual.titles = "A"))
})

testthat::test_that("do_NebulosaPlot: FAIL - not enough individual titles ", {
  sample <- SCpubr:::test_list$sample

  testthat::expect_error(SCpubr::do_NebulosaPlot(sample = sample,
                                                 features = c("CD14", "CD8A"),
                                                 joint = F,
                                                 individual.titles = "A"))
})

testthat::test_that("do_NebulosaPlot: FAIL - not enough individual titles for joint ", {
  sample <- SCpubr:::test_list$sample

  testthat::expect_error(SCpubr::do_NebulosaPlot(sample = sample,
                                                 features = c("CD14", "CD8A"),
                                                 joint = T,
                                                 individual.titles = "A"))
})

testthat::test_that("do_NebulosaPlot: PASS - no legend", {
  sample <- SCpubr:::test_list$sample

  p <- SCpubr::do_NebulosaPlot(sample = sample,
                               features = c("CD14"),
                               legend.position = "none")
  testthat::expect_type(p, "list")
})

testthat::test_that("do_NebulosaPlot: PASS - patchwork title, subtitle and caption", {
  sample <- SCpubr:::test_list$sample

  p <- SCpubr::do_NebulosaPlot(sample = sample,
                               features = c("CD14", "CD8A"),
                               plot.title = "A",
                               plot.subtitle = "B",
                               plot.caption = "C")
  testthat::expect_type(p, "list")
})

testthat::test_that("do_NebulosaPlot: PASS - plot axis", {
  sample <- SCpubr:::test_list$sample

  p <- SCpubr::do_NebulosaPlot(sample = sample, plot.axes = TRUE, features = "nCount_RNA")
  testthat::expect_type(p, "list")

  p <- SCpubr::do_NebulosaPlot(sample = sample, reduction = "pca", plot.axes = TRUE, features = "nCount_RNA")
  testthat::expect_type(p, "list")

  p <- SCpubr::do_NebulosaPlot(sample = sample, dims = c(2, 1), plot.axes = TRUE, features = "nCount_RNA")
  testthat::expect_type(p, "list")

  sample@reductions$diffusion <- sample@reductions$umap
  p <- SCpubr::do_NebulosaPlot(sample = sample,
                              reduction = "diffusion",
                              plot.axes = TRUE,
                              features = "nCount_RNA")
  testthat::expect_type(p, "list")
})
