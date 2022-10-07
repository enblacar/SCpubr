

testthat::test_that("do_BeeSwarmPlot: PASS - categorical variable dimred component", {
  sample <- SCpubr:::test_list$sample
  p <- SCpubr::do_BeeSwarmPlot(sample = sample,
                               feature_to_rank = "PC_1",
                               group.by = "seurat_clusters",
                               continuous_feature = F)
  testthat::expect_type(p, "list")
})

testthat::test_that("do_BeeSwarmPlot: PASS - cell_borders", {
  sample <- SCpubr:::test_list$sample
  p <- SCpubr::do_BeeSwarmPlot(sample = sample, feature_to_rank = "CD14", group.by = "seurat_clusters", plot_cell_borders = T)
  testthat::expect_type(p, "list")
  p <- SCpubr::do_BeeSwarmPlot(sample = sample, feature_to_rank = "CD14", group.by = "seurat_clusters", plot_cell_borders = T, raster = T, pt.size = 1)
  testthat::expect_type(p, "list")
})

testthat::test_that("do_BeeSwarmPlot: PASS - categorical variable gene", {
  sample <- SCpubr:::test_list$sample
  p <- SCpubr::do_BeeSwarmPlot(sample = sample,
                               feature_to_rank = "CD14",
                               group.by = "seurat_clusters",
                               continuous_feature = F)
  testthat::expect_type(p, "list")
})

testthat::test_that("do_BeeSwarmPlot: PASS - legend position = right", {
  sample <- SCpubr:::test_list$sample
  p <- SCpubr::do_BeeSwarmPlot(sample = sample,
                               feature_to_rank = "CD14",
                               group.by = "seurat_clusters",
                               continuous_feature = F,
                               legend.position = "right")
  testthat::expect_type(p, "list")
})

testthat::test_that("do_BeeSwarmPlot: PASS - categorical variable metadata", {
  sample <- SCpubr:::test_list$sample
  p <- SCpubr::do_BeeSwarmPlot(sample = sample,
                               feature_to_rank = "nCount_RNA",
                               group.by = "seurat_clusters",
                               continuous_feature = F)
  testthat::expect_type(p, "list")
})

testthat::test_that("do_BeeSwarmPlot: PASS - continuous variable", {
  sample <- SCpubr:::test_list$sample
  p <- SCpubr::do_BeeSwarmPlot(sample = sample,
                               feature_to_rank = "PC_1",
                               group.by = "seurat_clusters",
                               continuous_feature = T)
  testthat::expect_type(p, "list")
})

testthat::test_that("do_BeeSwarmPlot: PASS - continuous variable legend normal", {
  sample <- SCpubr:::test_list$sample
  p <- SCpubr::do_BeeSwarmPlot(sample = sample,
                               feature_to_rank = "PC_1",
                               group.by = "seurat_clusters",
                               continuous_feature = T,
                               legend.type = "normal")
  testthat::expect_type(p, "list")
})

testthat::test_that("do_BeeSwarmPlot: PASS - continuous variable legend colorbar", {
  sample <- SCpubr:::test_list$sample
  p <- SCpubr::do_BeeSwarmPlot(sample = sample,
                               feature_to_rank = "PC_1",
                               group.by = "seurat_clusters",
                               continuous_feature = T,
                               legend.type = "colorbar")
  testthat::expect_type(p, "list")
})

testthat::test_that("do_BeeSwarmPlot: PASS - continuous variable legend colorsteps", {
  sample <- SCpubr:::test_list$sample
  p <- SCpubr::do_BeeSwarmPlot(sample = sample,
                               feature_to_rank = "PC_1",
                               group.by = "seurat_clusters",
                               continuous_feature = T,
                               legend.type = "colorsteps")
  testthat::expect_type(p, "list")
})

testthat::test_that("do_BeeSwarmPlot: FAIL - wrong legend type", {
  sample <- SCpubr:::test_list$sample
  testthat::expect_error(SCpubr::do_BeeSwarmPlot(sample = sample,
                                                 feature_to_rank = "PC_1",
                                                 group.by = "seurat_clusters",
                                                 continuous_feature = T,
                                                 legend.type = "wrong"))
})

testthat::test_that("do_BeeSwarmPlot: FAIL - more than one feature", {
  sample <- SCpubr:::test_list$sample
  testthat::expect_error(SCpubr::do_BeeSwarmPlot(sample = sample,
                                                 feature_to_rank = c("PC_1", "PC_2"),
                                                 group.by = "seurat_clusters",
                                                 continuous_feature = T))
})

testthat::test_that("do_BeeSwarmPlot: FAIL - wrong legend position", {
  sample <- SCpubr:::test_list$sample
  testthat::expect_error(SCpubr::do_BeeSwarmPlot(sample = sample,
                                                 feature_to_rank = "PC_1",
                                                 group.by = "seurat_clusters",
                                                 continuous_feature = T,
                                                 legend.position = "wrong"))
})

testthat::test_that("do_BeeSwarmPlot: FAIL - wrong font.type", {
  sample <- SCpubr:::test_list$sample
  testthat::expect_error(SCpubr::do_BeeSwarmPlot(sample = sample,
                                                 feature_to_rank = "PC_1",
                                                 group.by = "seurat_clusters",
                                                 continuous_feature = T,
                                                 font.type = "wrong"))
})

testthat::test_that("do_BeeSwarmPlot: PASS - continuous variable viridis scale", {
  sample <- SCpubr:::test_list$sample
  p <- SCpubr::do_BeeSwarmPlot(sample = sample,
                               feature_to_rank = "PC_1",
                               group.by = "seurat_clusters",
                               continuous_feature = T,
                               viridis_color_map = "F")
  testthat::expect_type(p, "list")
})

testthat::test_that("do_BeeSwarmPlot: PASS - continuous variable legend position = top", {
  sample <- SCpubr:::test_list$sample
  p <- SCpubr::do_BeeSwarmPlot(sample = sample,
                               feature_to_rank = "PC_1",
                               group.by = "seurat_clusters",
                               continuous_feature = T,
                               legend.position = "top")
  testthat::expect_type(p, "list")
})

testthat::test_that("do_BeeSwarmPlot: FAIL - feature not found", {
  sample <- SCpubr:::test_list$sample
  testthat::expect_error(SCpubr::do_BeeSwarmPlot(sample = sample,
                                                 feature_to_rank = "not_found",
                                                 group.by = "seurat_clusters",
                                                 continuous_feature = T,
                                                 viridis_color_map = "F"))
})

testthat::test_that("do_BeeSwarmPlot: PASS - raster", {
  sample <- SCpubr:::test_list$sample
  p <- SCpubr::do_BeeSwarmPlot(sample = sample,
                               feature_to_rank = "PC_1",
                               group.by = "seurat_clusters",
                               raster = T)
  testthat::expect_type(p, "list")
})

testthat::test_that("do_BeeSwarmPlot: PASS - colors.use", {
  sample <- SCpubr:::test_list$sample
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
  sample <- SCpubr:::test_list$sample
  p <- SCpubr::do_BeeSwarmPlot(sample = sample,
                               feature_to_rank = "PC_1",
                               group.by = "seurat_clusters",
                               remove_x_axis = T)
  testthat::expect_type(p, "list")
})

testthat::test_that("do_BeeSwarmPlot: PASS - remove y axis", {
  sample <- SCpubr:::test_list$sample
  p <- SCpubr::do_BeeSwarmPlot(sample = sample,
                               feature_to_rank = "PC_1",
                               group.by = "seurat_clusters",
                               remove_y_axis = T)
  testthat::expect_type(p, "list")
})

testthat::test_that("do_BeeSwarmPlot: PASS - flip", {
  sample <- SCpubr:::test_list$sample
  p <- SCpubr::do_BeeSwarmPlot(sample = sample,
                               feature_to_rank = "PC_1",
                               group.by = "seurat_clusters",
                               flip = T)
  testthat::expect_type(p, "list")
})
