
testthat::test_that("do_GeyserPlot: PASS - default parameters", {
  sample <- SCpubr:::test_list$sample

  p <- SCpubr::do_GeyserPlot(sample = sample,
                             features = "CD14")
  testthat::expect_type(p, "list")
})

testthat::test_that("do_GeyserPlot: PASS - default parameters = symmetrical scale", {
  sample <- SCpubr:::test_list$sample

  p <- SCpubr::do_GeyserPlot(sample = sample,
                             features = "CD14",
                             enforce_symmetry = TRUE)
  testthat::expect_type(p, "list")

  p <- SCpubr::do_GeyserPlot(sample = sample,
                             features = "CD14",
                             enforce_symmetry = FALSE)
  testthat::expect_type(p, "list")
})


testthat::test_that("do_GeyserPlot: PASS - default parameters = categorical scale", {
  sample <- SCpubr:::test_list$sample

  p <- SCpubr::do_GeyserPlot(sample = sample,
                             features = "CD14",
                             scale_type = "categorical")
  testthat::expect_type(p, "list")
})

testthat::test_that("do_GeyserPlot: PASS - default parameters = color.by", {
  sample <- SCpubr:::test_list$sample

  p <- SCpubr::do_GeyserPlot(sample = sample,
                             features = "CD14",
                             scale_type = "categorical",
                             color.by = "orig.ident")
  testthat::expect_type(p, "list")

  p <- SCpubr::do_GeyserPlot(sample = sample,
                             features = "CD14",
                             scale_type = "continuous",
                             color.by = "PC_1")
  testthat::expect_type(p, "list")

  p <- SCpubr::do_GeyserPlot(sample = sample,
                             features = "PC_2",
                             scale_type = "categorical",
                             color.by = "orig.ident")
  testthat::expect_type(p, "list")

  p <- SCpubr::do_GeyserPlot(sample = sample,
                             features = "PC_2",
                             scale_type = "continuous",
                             color.by = "nCount_RNA")
  testthat::expect_type(p, "list")

  p <- SCpubr::do_GeyserPlot(sample = sample,
                             features = "nCount_RNA",
                             scale_type = "categorical",
                             color.by = "orig.ident")
  testthat::expect_type(p, "list")

  p <- SCpubr::do_GeyserPlot(sample = sample,
                             features = "nCount_RNA",
                             scale_type = "continuous",
                             color.by = "CD14")
  testthat::expect_type(p, "list")
})


testthat::test_that("do_GeyserPlot: PASS - categorical colors.use", {
  sample <- SCpubr:::test_list$sample

  p <- SCpubr::do_GeyserPlot(sample = sample,
                             features = "CD14",
                             scale_type = "categorical",
                             group.by = "orig.ident",
                             colors.use = c("Cell" = "green"))
  testthat::expect_type(p, "list")

  p <- SCpubr::do_GeyserPlot(sample = sample,
                             features = "CD14",
                             scale_type = "categorical",
                             group.by = "orig.ident",
                             colors.use = NULL)
  testthat::expect_type(p, "list")
})

testthat::test_that("do_GeyserPlot: PASS - order by mean", {
  sample <- SCpubr:::test_list$sample

  p <- SCpubr::do_GeyserPlot(sample = sample,
                             features = "CD14",
                             scale_type = "continuous",
                             order_by_mean = TRUE)
  testthat::expect_type(p, "list")

  p <- SCpubr::do_GeyserPlot(sample = sample,
                             features = "CD14",
                             scale_type = "continuous",
                             order_by_mean = FALSE)
  testthat::expect_type(p, "list")
})

testthat::test_that("do_GeyserPlot: PASS - plot cell borders", {
  sample <- SCpubr:::test_list$sample

  p <- SCpubr::do_GeyserPlot(sample = sample,
                             features = "CD14",
                             scale_type = "continuous",
                             plot_cell_borders = TRUE)
  testthat::expect_type(p, "list")

  p <- SCpubr::do_GeyserPlot(sample = sample,
                             features = "CD14",
                             scale_type = "continuous",
                             plot_cell_borders = FALSE)
  testthat::expect_type(p, "list")
})

testthat::test_that("do_GeyserPlot: PASS - split.by", {
  sample <- SCpubr:::test_list$sample

  p <- SCpubr::do_GeyserPlot(sample = sample,
                             features = "CD14",
                             scale_type = "continuous",
                             split.by = "seurat_clusters",
                             plot_cell_borders = TRUE)
  testthat::expect_type(p, "list")

  p <- SCpubr::do_GeyserPlot(sample = sample,
                             features = "CD14",
                             scale_type = "continuous",
                             split.by = "seurat_clusters",
                             plot_cell_borders = FALSE)
  testthat::expect_type(p, "list")
})


testthat::test_that("do_BarPlot: FAIL - wrong paramters", {
  sample <- SCpubr:::test_list$sample

  testthat::expect_error(SCpubr::do_GeyserPlot(sample = sample,
                                               features = "CD14",
                                               scale_type = "wrong"))
  testthat::expect_error(SCpubr::do_GeyserPlot(sample = sample,
                                               features = "CD14",
                                               split.by = "wrong"))
  testthat::expect_error(SCpubr::do_GeyserPlot(sample = sample,
                                               features = "CD14",
                                               group.by = "wrong"))
  testthat::expect_error(SCpubr::do_GeyserPlot(sample = sample,
                                               features = "CD14",
                                               jitter = 1))
  testthat::expect_error(SCpubr::do_GeyserPlot(sample = sample,
                                               features = "CD14",
                                               scale_type = "categorical",
                                               color.by = "CD14"))
})

testthat::test_that("do_GeyserPlot: PASS - show legend", {
  sample <- SCpubr:::test_list$sample

  p <- SCpubr::do_GeyserPlot(sample = sample,
                             features = "CD14")
  testthat::expect_type(p, "list")

  p <- SCpubr::do_GeyserPlot(sample = sample,
                             features = "CD14")
  testthat::expect_type(p, "list")
})

testthat::test_that("do_GeyserPlot: PASS - several features", {
  sample <- SCpubr:::test_list$sample

  p <- SCpubr::do_GeyserPlot(sample = sample,
                             features = c("CD14", "PC_1"))
  testthat::expect_type(p, "list")
  testthat::expect_length(p, 2)
})

testthat::test_that("do_GeyserPlot: PASS - color.by factor", {
  sample <- SCpubr:::test_list$sample

  sample$seurat_clusters_factor <- as.factor(sample$seurat_clusters)
  sample$seurat_clusters_character <- as.character(sample$seurat_clusters)
  p <- SCpubr::do_GeyserPlot(sample = sample,
                             features = "CD14",
                             color.by = "seurat_clusters_factor")
  testthat::expect_type(p, "list")

  p <- SCpubr::do_GeyserPlot(sample = sample,
                             features = "CD14",
                             color.by = "seurat_clusters_character")
  testthat::expect_type(p, "list")

  p <- SCpubr::do_GeyserPlot(sample = sample,
                             features = "CD14",
                             scale_type = "categorical",
                             color.by = "seurat_clusters_character")
  testthat::expect_type(p, "list")

  p <- SCpubr::do_GeyserPlot(sample = sample,
                             features = "CD14",
                             scale_type = "categorical",
                             color.by = "seurat_clusters_factor")
  testthat::expect_type(p, "list")
})



