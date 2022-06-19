sample <- use_dataset()
testthat::test_that("do_BeeSwarmPlot: PASS - one variable", {
  p <- SCpubr::do_VlnPlot(sample = sample,
                          features = "CD14")
  testthat::expect_type(p, "list")
})

testthat::test_that("do_BeeSwarmPlot: PASS - without boxplot", {
  p <- SCpubr::do_VlnPlot(sample = sample,
                          features = "CD14",
                          plot_boxplot = F)
  testthat::expect_type(p, "list")
})

testthat::test_that("do_BeeSwarmPlot: PASS - rotate axis", {
  p <- SCpubr::do_VlnPlot(sample = sample,
                          features = "CD14",
                          rotate_x_labels = T)
  testthat::expect_type(p, "list")
})

testthat::test_that("do_BeeSwarmPlot: PASS - several features", {
  p <- SCpubr::do_VlnPlot(sample = sample,
                          features = c("CD14", "nCount_RNA"))
  testthat::expect_type(p, "list")
})

testthat::test_that("do_BeeSwarmPlot: PASS - several features plot boxplots", {
  p <- SCpubr::do_VlnPlot(sample = sample,
                          features = c("CD14", "nCount_RNA"),
                          plot_boxplot = F)
  testthat::expect_type(p, "list")
})

testthat::test_that("do_BeeSwarmPlot: PASS - several features rotate labels", {
  p <- SCpubr::do_VlnPlot(sample = sample,
                          features = c("CD14", "nCount_RNA"),
                          rotate_x_labels = c(T, F))
  testthat::expect_type(p, "list")
})

testthat::test_that("do_BeeSwarmPlot: PASS - several features ycut", {
  p <- SCpubr::do_VlnPlot(sample = sample,
                          features = c("CD14", "nCount_RNA"),
                          rotate_x_labels = c(T, F),
                          y_cut = c(2, 20000))
  testthat::expect_type(p, "list")
})

testthat::test_that("do_BeeSwarmPlot: PASS - one feature ycut", {
  p <- SCpubr::do_VlnPlot(sample = sample,
                          features = c("CD14"),
                          rotate_x_labels = T,
                          y_cut = 2)
  testthat::expect_type(p, "list")
})


testthat::test_that("do_BeeSwarmPlot: PASS - one feature line width", {
  p <- SCpubr::do_VlnPlot(sample = sample,
                          features = c("CD14"),
                          rotate_x_labels = T,
                          y_cut = 2,
                          line_width = 3)
  testthat::expect_type(p, "list")
})

testthat::test_that("do_BeeSwarmPlot: PASS - one feature boxplot width", {
  p <- SCpubr::do_VlnPlot(sample = sample,
                          features = c("CD14"),
                          rotate_x_labels = T,
                          y_cut = 2,
                          boxplot_width = 0.1)
  testthat::expect_type(p, "list")
})

testthat::test_that("do_BeeSwarmPlot: PASS - change colors", {
  p <- SCpubr::do_VlnPlot(sample = sample,
                          features = c("CD14"),
                          rotate_x_labels = T,
                          y_cut = 2,
                          boxplot_width = 0.1,
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
