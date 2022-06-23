sample <- use_dataset()
testthat::test_that("do_BeeSwarmPlot: PASS - categorical variable dimred component", {
  p <- SCpubr::do_BeeSwarmPlot(sample = sample,
                               feature_to_rank = "PC_1",
                               group.by = "seurat_clusters",
                               continuous_feature = F)
  testthat::expect_type(p, "list")
})

testthat::test_that("do_BeeSwarmPlot: PASS - categorical variable gene", {
  p <- SCpubr::do_BeeSwarmPlot(sample = sample,
                               feature_to_rank = "CD14",
                               group.by = "seurat_clusters",
                               continuous_feature = F)
  testthat::expect_type(p, "list")
})

testthat::test_that("do_BeeSwarmPlot: PASS - categorical variable metadata", {
  p <- SCpubr::do_BeeSwarmPlot(sample = sample,
                               feature_to_rank = "nCount_RNA",
                               group.by = "seurat_clusters",
                               continuous_feature = F)
  testthat::expect_type(p, "list")
})

testthat::test_that("do_BeeSwarmPlot: PASS - continuous variable", {
  p <- SCpubr::do_BeeSwarmPlot(sample = sample,
                               feature_to_rank = "PC_1",
                               group.by = "seurat_clusters",
                               continuous_feature = T)
  testthat::expect_type(p, "list")
})

testthat::test_that("do_BeeSwarmPlot: PASS - continuous variable viridis scale", {
  p <- SCpubr::do_BeeSwarmPlot(sample = sample,
                               feature_to_rank = "PC_1",
                               group.by = "seurat_clusters",
                               continuous_feature = T,
                               viridis_color_map = "F")
  testthat::expect_type(p, "list")
})

testthat::test_that("do_BeeSwarmPlot: PASS - continuous variable legend position = top", {
  p <- SCpubr::do_BeeSwarmPlot(sample = sample,
                               feature_to_rank = "PC_1",
                               group.by = "seurat_clusters",
                               continuous_feature = T,
                               legend.position = "top")
  testthat::expect_type(p, "list")
})

testthat::test_that("do_BeeSwarmPlot: FAIL - feature not found", {
  testthat::expect_error(SCpubr::do_BeeSwarmPlot(sample = sample,
                                                 feature_to_rank = "not_found",
                                                 group.by = "seurat_clusters",
                                                 continuous_feature = T,
                                                 viridis_color_map = "F"))
})

testthat::test_that("do_BeeSwarmPlot: PASS - raster", {
  p <- SCpubr::do_BeeSwarmPlot(sample = sample,
                               feature_to_rank = "PC_1",
                               group.by = "seurat_clusters",
                               raster = T)
  testthat::expect_type(p, "list")
})

testthat::test_that("do_BeeSwarmPlot: PASS - colors.use", {
  p <- SCpubr::do_BeeSwarmPlot(sample = sample,
                               feature_to_rank = "PC_1",
                               group.by = "seurat_clusters",
                               colors.use = c("0" = "#001219",
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

testthat::test_that("do_BeeSwarmPlot: PASS - remove x axis", {
  p <- SCpubr::do_BeeSwarmPlot(sample = sample,
                               feature_to_rank = "PC_1",
                               group.by = "seurat_clusters",
                               remove_x_axis = T)
  testthat::expect_type(p, "list")
})

testthat::test_that("do_BeeSwarmPlot: PASS - remove y axis", {
  p <- SCpubr::do_BeeSwarmPlot(sample = sample,
                               feature_to_rank = "PC_1",
                               group.by = "seurat_clusters",
                               remove_y_axis = T)
  testthat::expect_type(p, "list")
})

testthat::test_that("do_BeeSwarmPlot: PASS - flip", {
  p <- SCpubr::do_BeeSwarmPlot(sample = sample,
                               feature_to_rank = "PC_1",
                               group.by = "seurat_clusters",
                               flip = T)
  testthat::expect_type(p, "list")
})
