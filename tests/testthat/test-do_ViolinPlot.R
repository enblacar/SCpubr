sample <- SCpubr:::test_list$sample

testthat::test_that("do_ViolinPlot: PASS - one variable", {
  p <- SCpubr::do_ViolinPlot(sample = sample,
                          feature = "CD14")
  testthat::expect_type(p, "list")
})

testthat::test_that("do_ViolinPlot: PASS - without boxplot", {
  p <- SCpubr::do_ViolinPlot(sample = sample,
                          feature = "CD14",
                          plot_boxplot = F)
  testthat::expect_type(p, "list")
})

testthat::test_that("do_ViolinPlot: PASS - rotate axis", {
  p <- SCpubr::do_ViolinPlot(sample = sample,
                          feature = "CD14",
                          rotate_x_axis_labels = T)
  testthat::expect_type(p, "list")
})

testthat::test_that("do_ViolinPlot: PASS - plot.grid", {
  p <- SCpubr::do_ViolinPlot(sample = sample,
                          feature = c("CD14"),
                          plot.grid = TRUE)
  testthat::expect_type(p, "list")
})




testthat::test_that("do_ViolinPlot: PASS - one feature ycut", {
  p <- SCpubr::do_ViolinPlot(sample = sample,
                          feature = c("CD14"),
                          rotate_x_axis_labels = T,
                          y_cut = 2)
  testthat::expect_type(p, "list")
})


testthat::test_that("do_ViolinPlot: PASS - one feature line width", {
  p <- SCpubr::do_ViolinPlot(sample = sample,
                          feature = c("CD14"),
                          rotate_x_axis_labels = T,
                          y_cut = 2,
                          line_width = 3)
  testthat::expect_type(p, "list")
})

testthat::test_that("do_ViolinPlot: PASS - one feature boxplot width", {
  p <- SCpubr::do_ViolinPlot(sample = sample,
                          feature = c("CD14"),
                          rotate_x_axis_labels = T,
                          y_cut = 2,
                          boxplot_width = 0.1)
  testthat::expect_type(p, "list")
})

testthat::test_that("do_ViolinPlot: PASS - change colors", {
  p <- SCpubr::do_ViolinPlot(sample = sample,
                          feature = c("CD14"),
                          rotate_x_axis_labels = T,
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



testthat::test_that("do_ViolinPlot: FAIL - split.by", {
  testthat::expect_error(SCpubr::do_ViolinPlot(sample = sample,
                                            feature = c("CD14"),
                                            split.by = "orig.ident"))
})




testthat::test_that("do_ViolinPlot: PASS - one variable, group by", {
  p <- SCpubr::do_ViolinPlot(sample = sample,
                          feature = "CD14",
                          group.by = "orig.ident")
  testthat::expect_type(p, "list")
})

testthat::test_that("do_ViolinPlot: PASS - one variable, xlab y lab", {
  p <- SCpubr::do_ViolinPlot(sample = sample,
                          feature = "CD14",
                          xlab = "y",
                          ylab = "x")
  testthat::expect_type(p, "list")
})


testthat::test_that("do_ViolinPlot: PASS - one variable, plot.title, subtitle and caption", {
  p <- SCpubr::do_ViolinPlot(sample = sample,
                          feature = "CD14",
                          plot.title = "A",
                          plot.subtitle = "B",
                          plot.caption = "C")
  testthat::expect_type(p, "list")
})


testthat::test_that("do_ViolinPlot: PASS - multiple variables plot.title, subtitle and caption", {
  p <- SCpubr::do_ViolinPlot(sample = sample,
                          c("CD14"),
                          group.by = "orig.ident",
                          colors.use = c("Cell" = "red"))
  testthat::expect_type(p, "list")
})

testthat::test_that("do_ViolinPlot: FAIL - wrong font.type", {
  testthat::expect_error(SCpubr::do_ViolinPlot(sample = sample,
                                            feature = c("CD14"),
                                            font.type = "wrong"))
})
