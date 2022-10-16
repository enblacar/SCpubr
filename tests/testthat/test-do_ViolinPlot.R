if(isFALSE(dep_check[["do_ViolinPlot"]])){
  testthat::test_that("do_ViolinPlot: PASS - one variable", {



    p <- SCpubr::do_ViolinPlot(sample = sample,
                               feature = "EPC1",
                               plot.grid = TRUE)
    testthat::expect_type(p, "list")

    p <- SCpubr::do_ViolinPlot(sample = sample,
                               feature = "EPC1",
                               plot.grid = FALSE)
    testthat::expect_type(p, "list")
  })

  testthat::test_that("do_ViolinPlot: PASS - group.by", {



    p <- SCpubr::do_ViolinPlot(sample = sample,
                               feature = "EPC1",
                               plot.grid = TRUE,
                               group.by = "seurat_clusters")
    testthat::expect_type(p, "list")

    sample$seurat_clusters <- as.character(sample$seurat_clusters)
    p <- SCpubr::do_ViolinPlot(sample = sample,
                               feature = "EPC1",
                               plot.grid = FALSE,
                               group.by = "seurat_clusters")
    testthat::expect_type(p, "list")
  })


  testthat::test_that("do_ViolinPlot: PASS - without boxplot", {



    p <- SCpubr::do_ViolinPlot(sample = sample,
                               feature = "EPC1",
                               plot_boxplot = FALSE)
    testthat::expect_type(p, "list")
  })

  testthat::test_that("do_ViolinPlot: PASS - rotate axis", {



    p <- SCpubr::do_ViolinPlot(sample = sample,
                               feature = "EPC1",
                               rotate_x_axis_labels = TRUE)
    testthat::expect_type(p, "list")
  })

  testthat::test_that("do_ViolinPlot: PASS - plot.grid", {



    p <- SCpubr::do_ViolinPlot(sample = sample,
                               feature = c("EPC1"),
                               plot.grid = TRUE)
    testthat::expect_type(p, "list")
  })




  testthat::test_that("do_ViolinPlot: PASS - one feature ycut", {



    p <- SCpubr::do_ViolinPlot(sample = sample,
                               feature = c("EPC1"),
                               rotate_x_axis_labels = TRUE,
                               y_cut = 2)
    testthat::expect_type(p, "list")
  })


  testthat::test_that("do_ViolinPlot: PASS - one feature line width", {



    p <- SCpubr::do_ViolinPlot(sample = sample,
                               feature = c("EPC1"),
                               rotate_x_axis_labels = TRUE,
                               y_cut = 2,
                               line_width = 3)
    testthat::expect_type(p, "list")
  })

  testthat::test_that("do_ViolinPlot: PASS - one feature boxplot width", {



    p <- SCpubr::do_ViolinPlot(sample = sample,
                               feature = c("EPC1"),
                               rotate_x_axis_labels = TRUE,
                               y_cut = 2,
                               boxplot_width = 0.1)
    testthat::expect_type(p, "list")
  })

  testthat::test_that("do_ViolinPlot: PASS - change colors", {



    p <- SCpubr::do_ViolinPlot(sample = sample,
                               feature = c("EPC1"),
                               rotate_x_axis_labels = TRUE,
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
                                                 feature = c("EPC1"),
                                                 split.by = "orig.ident"))
  })




  testthat::test_that("do_ViolinPlot: PASS - one variable, group by", {



    p <- SCpubr::do_ViolinPlot(sample = sample,
                               feature = "EPC1",
                               group.by = "orig.ident")
    testthat::expect_type(p, "list")
  })

  testthat::test_that("do_ViolinPlot: PASS - one variable, xlab y lab", {



    p <- SCpubr::do_ViolinPlot(sample = sample,
                               feature = "EPC1",
                               xlab = "y",
                               ylab = "x")
    testthat::expect_type(p, "list")
  })


  testthat::test_that("do_ViolinPlot: PASS - one variable, plot.title, subtitle and caption", {



    p <- SCpubr::do_ViolinPlot(sample = sample,
                               feature = "EPC1",
                               plot.title = "A",
                               plot.subtitle = "B",
                               plot.caption = "C")
    testthat::expect_type(p, "list")
  })


  testthat::test_that("do_ViolinPlot: PASS - multiple variables plot.title, subtitle and caption", {



    p <- SCpubr::do_ViolinPlot(sample = sample,
                               c("EPC1"),
                               group.by = "orig.ident",
                               colors.use = c("Cell" = "red"))
    testthat::expect_type(p, "list")

    sample$orig.ident <- factor(sample$orig.ident)
    p <- SCpubr::do_ViolinPlot(sample = sample,
                               c("EPC1"),
                               group.by = "orig.ident",
                               colors.use = c("Cell" = "red"))
    testthat::expect_type(p, "list")
  })

  testthat::test_that("do_ViolinPlot: FAIL - wrong font.type", {



    testthat::expect_error(SCpubr::do_ViolinPlot(sample = sample,
                                                 feature = c("EPC1"),
                                                 font.type = "wrong"))
  })
}

