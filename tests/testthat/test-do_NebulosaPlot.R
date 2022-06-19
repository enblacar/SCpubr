sample <- use_dataset()
testthat::test_that("do_NebulosaPlot: PASS - single feature", {
  p <- SCpubr::do_NebulosaPlot(sample = sample,
                               features = c("CD14"))
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
