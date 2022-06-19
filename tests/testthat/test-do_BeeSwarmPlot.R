sample <- use_dataset()
testthat::test_that("do_BeeSwarmPlot: PASS - categorical variable", {
  p <- SCpubr::do_BeeSwarmPlot(sample = sample,
                               feature_to_rank = "PC_1",
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
