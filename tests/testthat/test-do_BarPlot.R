sample <- use_dataset()
testthat::test_that("do_BarPlot: PASS - one variable - stack", {
  p <- SCpubr::do_BarPlot(sample = sample,
                          features = "seurat_clusters",
                          position = "stack")
  testthat::expect_type(p, "list")
})

testthat::test_that("do_BarPlot: PASS - one variable - fill", {
  p <- SCpubr::do_BarPlot(sample = sample,
                          features = "seurat_clusters",
                          position = "fill")
  testthat::expect_type(p, "list")
})

testthat::test_that("do_BarPlot: PASS - one variable - fill - horizontal", {
  p <- SCpubr::do_BarPlot(sample = sample,
                          features = "seurat_clusters",
                          position = "fill",
                          horizontal = T)
  testthat::expect_type(p, "list")
})
sample$orig.ident <- ifelse(sample$seurat_clusters %in% c("1", "2"), "A", "B")

testthat::test_that("do_BarPlot: PASS - two variables - fill - horizontal", {
  p <- SCpubr::do_BarPlot(sample = sample,
                          features = "seurat_clusters",
                          group.by = "orig.ident",
                          position = "fill",
                          horizontal = T)
  testthat::expect_type(p, "list")
})

testthat::test_that("do_BarPlot: PASS - two variables - stack - horizontal", {
  p <- SCpubr::do_BarPlot(sample = sample,
                          features = "seurat_clusters",
                          group.by = "orig.ident",
                          position = "stack",
                          horizontal = T)
  testthat::expect_type(p, "list")
})

testthat::test_that("do_BarPlot: PASS - two variables - stack - horizontal - ordered", {
  p <- SCpubr::do_BarPlot(sample = sample,
                          features = "orig.ident",
                          group.by = "seurat_clusters",
                          position = "stack",
                          horizontal = T,
                          order.by = "2")
  testthat::expect_type(p, "list")
})

testthat::test_that("do_BarPlot: PASS - two variables - stack - horizontal - adding summary values", {
  p <- SCpubr::do_BarPlot(sample = sample,
                          features = "orig.ident",
                          group.by = "seurat_clusters",
                          position = "stack",
                          horizontal = T,
                          add.summary_labels = T)
  testthat::expect_type(p, "list")
})

testthat::test_that("do_BarPlot: PASS - two variables - stack - horizontal - adding summary values and size", {
  p <- SCpubr::do_BarPlot(sample = sample,
                          features = "orig.ident",
                          group.by = "seurat_clusters",
                          position = "stack",
                          horizontal = T,
                          add.summary_labels = T,
                          size.labels = 2.5)
  testthat::expect_type(p, "list")
})

testthat::test_that("do_BarPlot: PASS - two variables - stack - horizontal - adding subgroup labels", {
  p <- SCpubr::do_BarPlot(sample = sample,
                          features = "orig.ident",
                          group.by = "seurat_clusters",
                          position = "stack",
                          horizontal = T,
                          add.summary_labels = T,
                          size.labels = 5,
                          add.subgroup_labels = T)
  testthat::expect_type(p, "list")
})

testthat::test_that("do_BarPlot: PASS - two variables - stack - horizontal - adding subgroup labels - repel subgroup labels", {
  p <- SCpubr::do_BarPlot(sample = sample,
                          features = "orig.ident",
                          group.by = "seurat_clusters",
                          position = "stack",
                          horizontal = T,
                          add.summary_labels = T,
                          size.labels = 5,
                          add.subgroup_labels = T,
                          repel.subgroup_labels = T)
  testthat::expect_type(p, "list")
})

testthat::test_that("do_BarPlot: PASS - two variables - stack - horizontal - adding subgroup labels - repel summary labels", {
  p <- SCpubr::do_BarPlot(sample = sample,
                          features = "orig.ident",
                          group.by = "seurat_clusters",
                          position = "stack",
                          horizontal = T,
                          add.summary_labels = T,
                          size.labels = 5,
                          add.subgroup_labels = T,
                          repel.subgroup_labels = T,
                          repel.summary_labels = T)
  testthat::expect_type(p, "list")
})
