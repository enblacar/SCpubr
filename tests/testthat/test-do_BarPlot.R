sample <- SCpubr:::use_dataset()
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


testthat::test_that("do_BarPlot: FAIL - wrong position", {
  testthat::expect_error(SCpubr::do_BarPlot(sample = sample,
                                            features = "orig.ident",
                                            group.by = "seurat_clusters",
                                            position = "wrong"))
})

testthat::test_that("do_BarPlot: FAIL - wrong number of rotate x labels", {
  testthat::expect_error(SCpubr::do_BarPlot(sample = sample,
                                            features = c("orig.ident", "seurat_clusters"),
                                            group.by = "seurat_clusters",
                                            rotate_x_labels = T,
                                            position = "fill"))
})

testthat::test_that("do_BarPlot: PASS - rotate x labels", {
  p <- SCpubr::do_BarPlot(sample = sample,
                          features = "orig.ident",
                          rotate_x_labels = T)
  testthat::expect_type(p, "list")
})

testthat::test_that("do_BarPlot: PASS - rotate x labels", {
  sample$seurat.clusters.factor <- factor(sample$seurat_clusters)
  p <- SCpubr::do_BarPlot(sample = sample,
                          features = "seurat.clusters.factor")
  testthat::expect_type(p, "list")
})

testthat::test_that("do_BarPlot: PASS - rotate x labels with group.by", {
  sample$seurat.clusters.factor <- factor(sample$seurat_clusters)
  p <- SCpubr::do_BarPlot(sample = sample,
                          features = "seurat.clusters.factor",
                          group.by = "seurat_clusters")
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
                          features = "seurat.clusters.factor",
                          group.by = "seurat_clusters",
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
                          features = "seurat.clusters.factor",
                          colors.use = colors)
  testthat::expect_type(p, "list")
})

testthat::test_that("do_BarPlot: PASS - labels.order ", {
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
                          features = "seurat.clusters.factor",
                          labels.order = rev(names(colors)))
  testthat::expect_type(p, "list")
})


testthat::test_that("do_BarPlot: FAIL - order by not in group.by", {
  testthat::expect_error(SCpubr::do_BarPlot(sample = sample,
                                            features = c("orig.ident", "seurat_clusters"),
                                            group.by = "seurat_clusters",
                                            order.by = "wrong"))
})

testthat::test_that("do_BarPlot: PASS - one variable - adding labels - stack", {
  p <- SCpubr::do_BarPlot(sample = sample,
                          features = "orig.ident",
                          position = "stack",
                          horizontal = T,
                          add.summary_labels = T,
                          size.labels = 5,
                          add.subgroup_labels = T)
  testthat::expect_type(p, "list")
})

testthat::test_that("do_BarPlot: WARNING - one variable - adding labels - fill", {
  testthat::expect_warning(SCpubr::do_BarPlot(sample = sample,
                                              features = "orig.ident",
                                              position = "fill",
                                              horizontal = T,
                                              add.summary_labels = T,
                                              size.labels = 5,
                                              add.subgroup_labels = F))
  testthat::expect_warning(SCpubr::do_BarPlot(sample = sample,
                                              features = "orig.ident",
                                              position = "fill",
                                              horizontal = T,
                                              add.summary_labels = F,
                                              size.labels = 5,
                                              add.subgroup_labels = T))
})

testthat::test_that("do_BarPlot: PASS - one variable - rotate x labels", {
  p <- SCpubr::do_BarPlot(sample = sample,
                          features = "orig.ident",
                          position = "stack",
                          rotate_x_labels = T)
  testthat::expect_type(p, "list")
})

testthat::test_that("do_BarPlot: PASS - one variable - xlab, ylab and title", {
  p <- SCpubr::do_BarPlot(sample = sample,
                          features = "orig.ident",
                          position = "stack",
                          xlab = "A",
                          ylab = "B",
                          plot.title = "C")
  testthat::expect_type(p, "list")
})

testthat::test_that("do_BarPlot: PASS - one variable - no legend", {
  p <- SCpubr::do_BarPlot(sample = sample,
                          features = "orig.ident",
                          position = "stack",
                          legend = F)
  testthat::expect_type(p, "list")
})

testthat::test_that("do_BarPlot: PASS - two features", {
  p <- SCpubr::do_BarPlot(sample = sample,
                          features = c("orig.ident", "seurat_clusters"),
                          position = "stack",
                          legend = F)
  testthat::expect_type(p, "list")
})

testthat::test_that("do_BarPlot: PASS - return data", {
  out <- SCpubr::do_BarPlot(sample = sample,
                            features = c("orig.ident"),
                            return_data_matrix = T)
  testthat::expect_type(out$plot, "list")
  testthat::expect_type(out$data$orig.ident$long, "list")
  testthat::expect_type(out$data$orig.ident$wide, "list")
})
