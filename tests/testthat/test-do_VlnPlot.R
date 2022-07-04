sample <- use_dataset()
testthat::test_that("do_VlnPlot: PASS - one variable", {
  p <- SCpubr::do_VlnPlot(sample = sample,
                          features = "CD14")
  testthat::expect_type(p, "list")
})

testthat::test_that("do_VlnPlot: PASS - without boxplot", {
  p <- SCpubr::do_VlnPlot(sample = sample,
                          features = "CD14",
                          plot_boxplot = F)
  testthat::expect_type(p, "list")
})

testthat::test_that("do_VlnPlot: PASS - rotate axis", {
  p <- SCpubr::do_VlnPlot(sample = sample,
                          features = "CD14",
                          rotate_x_labels = T)
  testthat::expect_type(p, "list")
})

testthat::test_that("do_VlnPlot: PASS - several features", {
  p <- SCpubr::do_VlnPlot(sample = sample,
                          features = c("CD14", "nCount_RNA"))
  testthat::expect_type(p, "list")
})

testthat::test_that("do_VlnPlot: PASS - several features plot boxplots", {
  p <- SCpubr::do_VlnPlot(sample = sample,
                          features = c("CD14", "nCount_RNA"),
                          plot_boxplot = F)
  testthat::expect_type(p, "list")
})

testthat::test_that("do_VlnPlot: PASS - show_y_axis_lines", {
  p <- SCpubr::do_VlnPlot(sample = sample,
                          features = c("CD14"),
                          show_y_axis_lines = TRUE)
  testthat::expect_type(p, "list")
})

testthat::test_that("do_VlnPlot: PASS - several features rotate labels", {
  p <- SCpubr::do_VlnPlot(sample = sample,
                          features = c("CD14", "nCount_RNA"),
                          rotate_x_labels = c(T, F))
  testthat::expect_type(p, "list")
})

testthat::test_that("do_VlnPlot: PASS - several features ycut", {
  p <- SCpubr::do_VlnPlot(sample = sample,
                          features = c("CD14", "nCount_RNA"),
                          rotate_x_labels = c(T, F),
                          y_cut = c(2, 20000))
  testthat::expect_type(p, "list")
})

testthat::test_that("do_VlnPlot: PASS - one feature ycut", {
  p <- SCpubr::do_VlnPlot(sample = sample,
                          features = c("CD14"),
                          rotate_x_labels = T,
                          y_cut = 2)
  testthat::expect_type(p, "list")
})


testthat::test_that("do_VlnPlot: PASS - one feature line width", {
  p <- SCpubr::do_VlnPlot(sample = sample,
                          features = c("CD14"),
                          rotate_x_labels = T,
                          y_cut = 2,
                          line_width = 3)
  testthat::expect_type(p, "list")
})

testthat::test_that("do_VlnPlot: PASS - one feature boxplot width", {
  p <- SCpubr::do_VlnPlot(sample = sample,
                          features = c("CD14"),
                          rotate_x_labels = T,
                          y_cut = 2,
                          boxplot_width = 0.1)
  testthat::expect_type(p, "list")
})

testthat::test_that("do_VlnPlot: PASS - change colors", {
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

testthat::test_that("do_VlnPlot: WARNING - features as list", {
  testthat::expect_warning(SCpubr::do_VlnPlot(sample = sample,
                                              features = list("CD14")))
})

testthat::test_that("do_VlnPlot: FAIL - split.by", {
  testthat::expect_error(SCpubr::do_VlnPlot(sample = sample,
                                            features = c("CD14"),
                                            split.by = "orig.ident"))
})

testthat::test_that("do_VlnPlot: FAIL - different number of y_cut values", {
  testthat::expect_error(SCpubr::do_VlnPlot(sample = sample,
                                            features = c("CD14", "nCount_RNA"),
                                            y_cut = 2000))
})

testthat::test_that("do_VlnPlot: FAIL - different number of rotate_x_labels", {
  testthat::expect_error(SCpubr::do_VlnPlot(sample = sample,
                                            features = c("CD14", "nCount_RNA"),
                                           rotate_x_labels = T))
})

testthat::test_that("do_VlnPlot: FAIL - different number of individual titles", {
  testthat::expect_error(SCpubr::do_VlnPlot(sample = sample,
                                            features = c("CD14", "nCount_RNA"),
                                            individual.titles = c("A")))
})

testthat::test_that("do_VlnPlot: PASS - one variable, group by", {
  p <- SCpubr::do_VlnPlot(sample = sample,
                          features = "CD14",
                          group.by = "orig.ident")
  testthat::expect_type(p, "list")
})

testthat::test_that("do_VlnPlot: PASS - one variable, xlab y lab", {
  p <- SCpubr::do_VlnPlot(sample = sample,
                          features = "CD14",
                          xlab = "y",
                          ylab = "x")
  testthat::expect_type(p, "list")
})

testthat::test_that("do_VlnPlot: PASS - individual titles", {
  p <- SCpubr::do_VlnPlot(sample = sample,
                          c("CD14", "nCount_RNA"),
                          individual.titles = c("A", "B"))
  testthat::expect_type(p, "list")
})

testthat::test_that("do_VlnPlot: PASS - one variable, plot.title, subtitle and caption", {
  p <- SCpubr::do_VlnPlot(sample = sample,
                          features = "CD14",
                          plot.title = "A",
                          plot.subtitle = "B",
                          plot.caption = "C")
  testthat::expect_type(p, "list")
})

testthat::test_that("do_VlnPlot: PASS - multiple variables plot.title, subtitle and caption", {
  p <- SCpubr::do_VlnPlot(sample = sample,
                          c("CD14", "nCount_RNA"),
                          plot.title = "A",
                          plot.subtitle = "B",
                          plot.caption = "C")
  testthat::expect_type(p, "list")
})

testthat::test_that("do_VlnPlot: PASS - multiple variables plot.title, subtitle and caption", {
  p <- SCpubr::do_VlnPlot(sample = sample,
                          c("CD14"),
                          group.by = "orig.ident",
                          colors.use = c("Cell" = "red"))
  testthat::expect_type(p, "list")
})
