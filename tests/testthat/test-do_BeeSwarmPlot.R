if (isFALSE(dep_check[["do_BeeSwarmPlot"]])){
  testthat::test_that("do_BeeSwarmPlot: PASS - categorical variable dimred component", {

    p <- SCpubr::do_BeeSwarmPlot(sample = sample,
                                 feature_to_rank = "PC_1",
                                 group.by = "seurat_clusters",
                                 continuous_feature = FALSE)
    testthat::expect_type(p, "list")
  })

  testthat::test_that("do_BeeSwarmPlot: PASS - cell_borders", {

    p <- SCpubr::do_BeeSwarmPlot(sample = sample, feature_to_rank = "EPC1", group.by = "seurat_clusters", plot_cell_borders = TRUE)
    testthat::expect_type(p, "list")
    p <- SCpubr::do_BeeSwarmPlot(sample = sample, feature_to_rank = "EPC1", group.by = "seurat_clusters", plot_cell_borders = TRUE, raster = TRUE, pt.size = 1)
    testthat::expect_type(p, "list")
  })

  testthat::test_that("do_BeeSwarmPlot: PASS - categorical variable gene", {

    p <- SCpubr::do_BeeSwarmPlot(sample = sample,
                                 feature_to_rank = "EPC1",
                                 group.by = "seurat_clusters",
                                 continuous_feature = FALSE)
    testthat::expect_type(p, "list")
  })

  testthat::test_that("do_BeeSwarmPlot: PASS - legend position = right", {

    p <- SCpubr::do_BeeSwarmPlot(sample = sample,
                                 feature_to_rank = "EPC1",
                                 group.by = "seurat_clusters",
                                 continuous_feature = FALSE,
                                 legend.position = "right")
    testthat::expect_type(p, "list")
  })

  testthat::test_that("do_BeeSwarmPlot: PASS - categorical variable metadata", {

    p <- SCpubr::do_BeeSwarmPlot(sample = sample,
                                 feature_to_rank = "nCount_RNA",
                                 group.by = "seurat_clusters",
                                 continuous_feature = FALSE)
    testthat::expect_type(p, "list")
  })

  testthat::test_that("do_BeeSwarmPlot: PASS - continuous variable", {

    p <- SCpubr::do_BeeSwarmPlot(sample = sample,
                                 feature_to_rank = "PC_1",
                                 group.by = "seurat_clusters",
                                 continuous_feature = TRUE)
    testthat::expect_type(p, "list")
  })

  testthat::test_that("do_BeeSwarmPlot: PASS - continuous variable legend normal", {

    p <- SCpubr::do_BeeSwarmPlot(sample = sample,
                                 feature_to_rank = "PC_1",
                                 group.by = "seurat_clusters",
                                 continuous_feature = TRUE,
                                 legend.type = "normal")
    testthat::expect_type(p, "list")
  })

  testthat::test_that("do_BeeSwarmPlot: PASS - continuous variable legend colorbar", {

    p <- SCpubr::do_BeeSwarmPlot(sample = sample,
                                 feature_to_rank = "PC_1",
                                 group.by = "seurat_clusters",
                                 continuous_feature = TRUE,
                                 legend.type = "colorbar")
    testthat::expect_type(p, "list")
  })

  testthat::test_that("do_BeeSwarmPlot: PASS - continuous variable legend colorsteps", {

    p <- SCpubr::do_BeeSwarmPlot(sample = sample,
                                 feature_to_rank = "PC_1",
                                 group.by = "seurat_clusters",
                                 continuous_feature = TRUE,
                                 legend.type = "colorsteps")
    testthat::expect_type(p, "list")
  })

  testthat::test_that("do_BeeSwarmPlot: FAIL - wrong legend type", {

    testthat::expect_error(SCpubr::do_BeeSwarmPlot(sample = sample,
                                                   feature_to_rank = "PC_1",
                                                   group.by = "seurat_clusters",
                                                   continuous_feature = TRUE,
                                                   legend.type = "wrong"))
  })

  testthat::test_that("do_BeeSwarmPlot: FAIL - more than one feature", {

    testthat::expect_error(SCpubr::do_BeeSwarmPlot(sample = sample,
                                                   feature_to_rank = c("PC_1", "PC_2"),
                                                   group.by = "seurat_clusters",
                                                   continuous_feature = TRUE))
  })

  testthat::test_that("do_BeeSwarmPlot: FAIL - wrong legend position", {

    testthat::expect_error(SCpubr::do_BeeSwarmPlot(sample = sample,
                                                   feature_to_rank = "PC_1",
                                                   group.by = "seurat_clusters",
                                                   continuous_feature = TRUE,
                                                   legend.position = "wrong"))
  })

  testthat::test_that("do_BeeSwarmPlot: FAIL - wrong font.type", {

    testthat::expect_error(SCpubr::do_BeeSwarmPlot(sample = sample,
                                                   feature_to_rank = "PC_1",
                                                   group.by = "seurat_clusters",
                                                   continuous_feature = TRUE,
                                                   font.type = "wrong"))
  })

  testthat::test_that("do_BeeSwarmPlot: PASS - continuous variable viridis scale", {

    p <- SCpubr::do_BeeSwarmPlot(sample = sample,
                                 feature_to_rank = "PC_1",
                                 group.by = "seurat_clusters",
                                 continuous_feature = TRUE,
                                 viridis_color_map = "F")
    testthat::expect_type(p, "list")
  })

  testthat::test_that("do_BeeSwarmPlot: PASS - continuous variable legend position = top", {

    p <- SCpubr::do_BeeSwarmPlot(sample = sample,
                                 feature_to_rank = "PC_1",
                                 group.by = "seurat_clusters",
                                 continuous_feature = TRUE,
                                 legend.position = "top")
    testthat::expect_type(p, "list")
  })

  testthat::test_that("do_BeeSwarmPlot: FAIL - feature not found", {

    testthat::expect_error(SCpubr::do_BeeSwarmPlot(sample = sample,
                                                   feature_to_rank = "not_found",
                                                   group.by = "seurat_clusters",
                                                   continuous_feature = TRUE,
                                                   viridis_color_map = "F"))
  })

  testthat::test_that("do_BeeSwarmPlot: PASS - raster", {

    p <- SCpubr::do_BeeSwarmPlot(sample = sample,
                                 feature_to_rank = "PC_1",
                                 group.by = "seurat_clusters",
                                 raster = TRUE)
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
                                 remove_x_axis = TRUE)
    testthat::expect_type(p, "list")
  })

  testthat::test_that("do_BeeSwarmPlot: PASS - remove y axis", {

    p <- SCpubr::do_BeeSwarmPlot(sample = sample,
                                 feature_to_rank = "PC_1",
                                 group.by = "seurat_clusters",
                                 remove_y_axis = TRUE)
    testthat::expect_type(p, "list")
  })

  testthat::test_that("do_BeeSwarmPlot: PASS - flip", {

    p <- SCpubr::do_BeeSwarmPlot(sample = sample,
                                 feature_to_rank = "PC_1",
                                 group.by = "seurat_clusters",
                                 flip = TRUE)
    testthat::expect_type(p, "list")
  })
}

