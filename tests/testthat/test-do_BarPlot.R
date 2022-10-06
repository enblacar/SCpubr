sample <- SCpubr:::use_dataset()
testthat::test_that("do_BarPlot: PASS - one variable - stack", {
  p <- SCpubr::do_BarPlot(sample = sample,
                          group.by = "seurat_clusters",
                          position = "stack")
  testthat::expect_type(p, "list")
})

testthat::test_that("do_BarPlot: PASS - one variable - fill", {
  p <- SCpubr::do_BarPlot(sample = sample,
                          group.by = "seurat_clusters",
                          position = "fill")
  testthat::expect_type(p, "list")
})

testthat::test_that("do_BarPlot: PASS - one variable - remove guides", {
  p <- SCpubr::do_BarPlot(sample = sample,
                          group.by = "seurat_clusters",
                          position = "stack",
                          plot.grid = F)
  testthat::expect_type(p, "list")
})

testthat::test_that("do_BarPlot: PASS - one variable - fill - flip", {
  p <- SCpubr::do_BarPlot(sample = sample,
                          group.by = "seurat_clusters",
                          position = "fill",
                          flip = T)
  testthat::expect_type(p, "list")
})
sample$orig.ident <- ifelse(sample$seurat_clusters %in% c("1", "2"), "A", "B")

testthat::test_that("do_BarPlot: PASS - two variables - fill - flip", {
  p <- SCpubr::do_BarPlot(sample = sample,
                          group.by = "seurat_clusters",
                          split.by = "orig.ident",
                          position = "fill",
                          flip = T)
  testthat::expect_type(p, "list")
})

testthat::test_that("do_BarPlot: PASS - two variables - stack - flip", {
  p <- SCpubr::do_BarPlot(sample = sample,
                          group.by = "seurat_clusters",
                          split.by = "orig.ident",
                          position = "stack",
                          flip = T)
  testthat::expect_type(p, "list")
})

testthat::test_that("do_BarPlot: PASS - two variables - stack - flip - ordered", {
  p <- SCpubr::do_BarPlot(sample = sample,
                          group.by = "orig.ident",
                          split.by = "seurat_clusters",
                          position = "stack",
                          flip = T,
                          order = T)
  testthat::expect_type(p, "list")
})



testthat::test_that("do_BarPlot: FAIL - wrong position", {
  testthat::expect_error(SCpubr::do_BarPlot(sample = sample,
                                            group.by = "orig.ident",
                                            group.by = "seurat_clusters",
                                            position = "wrong"))
})

testthat::test_that("do_BarPlot: FAIL - wrong font.type", {
  testthat::expect_error(SCpubr::do_BarPlot(sample = sample,
                                            group.by = "orig.ident",
                                            group.by = "seurat_clusters",
                                            font.type = "wrong"))
})

testthat::test_that("do_BarPlot: PASS - rotate x labels", {
  p <- SCpubr::do_BarPlot(sample = sample,
                          group.by = "orig.ident",
                          rotate_x_axis_labels  = T)
  testthat::expect_type(p, "list")
})

testthat::test_that("do_BarPlot: PASS - rotate x labels", {
  sample$seurat.clusters.factor <- factor(sample$seurat_clusters)
  p <- SCpubr::do_BarPlot(sample = sample,
                          group.by = "seurat.clusters.factor")
  testthat::expect_type(p, "list")
})

testthat::test_that("do_BarPlot: PASS - rotate x labels with group.by", {
  sample$seurat.clusters.factor <- factor(sample$seurat_clusters)
  p <- SCpubr::do_BarPlot(sample = sample,
                          group.by = "seurat.clusters.factor",
                          split.by = "seurat_clusters")
  testthat::expect_type(p, "list")
})


testthat::test_that("do_BarPlot: PASS - colors.use and group.by", {
  sample$seurat.clusters.factor <- factor(sample$seurat_clusters)

  colors <- c("0" = "#001219",
              "1" = "#005f73",
              "2" = "#0a9396",
              "3" = "#94d2bd",
              "4" = "#e9d8a6",
              "5" = "#ee9b00",
              "6" = "#ca6702",
              "7" = "#bb3e03",
              "8" = "#ae2012")

  p <- SCpubr::do_BarPlot(sample = sample,
                          group.by = "seurat.clusters.factor",
                          split.by = "seurat_clusters",
                          colors.use = colors)
  testthat::expect_type(p, "list")
})

testthat::test_that("do_BarPlot: PASS - colors.use ", {
  sample$seurat.clusters.factor <- factor(sample$seurat_clusters)

  colors <- c("0" = "#001219",
              "1" = "#005f73",
              "2" = "#0a9396",
              "3" = "#94d2bd",
              "4" = "#e9d8a6",
              "5" = "#ee9b00",
              "6" = "#ca6702",
              "7" = "#bb3e03",
              "8" = "#ae2012")

  p <- SCpubr::do_BarPlot(sample = sample,
                          group.by = "seurat.clusters.factor",
                          colors.use = colors)
  testthat::expect_type(p, "list")
})


testthat::test_that("do_BarPlot: FAIL - order by not in group.by", {
  testthat::expect_error(SCpubr::do_BarPlot(sample = sample,
                                            group.by = c("orig.ident", "seurat_clusters"),
                                            group.by = "seurat_clusters",
                                            order.by = "wrong"))
})


testthat::test_that("do_BarPlot: PASS - one variable - rotate x labels", {
  p <- SCpubr::do_BarPlot(sample = sample,
                          group.by = "orig.ident",
                          position = "stack",
                          rotate_x_axis_labels = T)
  testthat::expect_type(p, "list")
})

testthat::test_that("do_BarPlot: PASS - one variable - xlab, ylab and title", {
  p <- SCpubr::do_BarPlot(sample = sample,
                          group.by = "orig.ident",
                          position = "stack",
                          xlab = "A",
                          ylab = "B",
                          plot.title = "C",
                          plot.subtitle = "D",
                          plot.caption = "E")
  testthat::expect_type(p, "list")
})

testthat::test_that("do_BarPlot: PASS - one variable - no legend", {
  p <- SCpubr::do_BarPlot(sample = sample,
                          group.by = "orig.ident",
                          position = "stack",
                          legend.position = "none")
  testthat::expect_type(p, "list")
})


testthat::test_that("do_BarPlot: PASS - group.by factor", {
  sample$factor_seurat_clusters <- factor(sample$seurat_clusters)
  p <- SCpubr::do_BarPlot(sample = sample,
                          group.by = c("seurat_clusters"),
                          split.by = "factor_seurat_clusters",
                          position = "stack",
                          legend.position = "none")
  testthat::expect_type(p, "list")
})



