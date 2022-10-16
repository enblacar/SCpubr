if (isFALSE(dep_check[["do_PathwayActivityPlot"]])){
  testthat::test_that("do_PathwayActivityPlot: PASS - minimal input", {




    out <- SCpubr::do_PathwayActivityPlot(sample = sample,
                                          activities = progeny_activities)
    testthat::expect_type(out, "list")

    out <- SCpubr::do_PathwayActivityPlot(sample = sample,
                                          activities = progeny_activities,
                                          flip = TRUE)
    testthat::expect_type(out, "list")

    out <- SCpubr::do_PathwayActivityPlot(sample = sample,
                                          activities = progeny_activities,
                                          flip = TRUE,
                                          split.by = "orig.ident")
    testthat::expect_type(out, "list")

    out <- SCpubr::do_PathwayActivityPlot(sample = sample,
                                          activities = progeny_activities,
                                          legend.position = "right")
    testthat::expect_type(out, "list")

    out <- SCpubr::do_PathwayActivityPlot(sample = sample,
                                          activities = progeny_activities,
                                          legend.position = "right",
                                          split.by = "orig.ident")
    testthat::expect_type(out, "list")
  })

  testthat::test_that("do_PathwayActivityPlot: PASS - plot featureplots", {




    out <- SCpubr::do_PathwayActivityPlot(sample = sample,
                                          activities = progeny_activities,
                                          plot_FeaturePlots = TRUE)
    testthat::expect_type(out, "list")
    testthat::expect_length(out, 2)
  })

  testthat::test_that("do_PathwayActivityPlot: PASS - plot geysers", {




    out <- SCpubr::do_PathwayActivityPlot(sample = sample,
                                          activities = progeny_activities,
                                          plot_GeyserPlots = TRUE)
    testthat::expect_type(out, "list")
    testthat::expect_length(out, 2)

    out <- SCpubr::do_PathwayActivityPlot(sample = sample,
                                          activities = progeny_activities,
                                          plot_GeyserPlots = TRUE,
                                          geyser_color.by = "nCount_RNA")
    testthat::expect_type(out, "list")
    testthat::expect_length(out, 2)
  })

  testthat::test_that("do_PathwayActivityPlot: PASS - all", {




    out <- SCpubr::do_PathwayActivityPlot(sample = sample,
                                          activities = progeny_activities,
                                          plot_GeyserPlots = TRUE,
                                          plot_FeaturePlots = TRUE)
    testthat::expect_type(out, "list")
    testthat::expect_length(out, 3)
  })

  testthat::test_that("do_PathwayActivityPlot: PASS - all group.by", {




    out <- SCpubr::do_PathwayActivityPlot(sample = sample,
                                          activities = progeny_activities,
                                          plot_GeyserPlots = TRUE,
                                          plot_FeaturePlots = TRUE,
                                          group.by = "orig.ident")
    testthat::expect_type(out, "list")
    testthat::expect_length(out, 3)
  })

  testthat::test_that("do_PathwayActivityPlot: PASS - all split.by", {




    out <- SCpubr::do_PathwayActivityPlot(sample = sample,
                                          activities = progeny_activities,
                                          plot_GeyserPlots = TRUE,
                                          plot_FeaturePlots = TRUE,
                                          split.by = "orig.ident")
    testthat::expect_type(out, "list")
    testthat::expect_length(out, 3)
  })

  testthat::test_that("do_PathwayActivityPlot: PASS - column.title and row.title", {




    out <- SCpubr::do_PathwayActivityPlot(sample = sample,
                                          activities = progeny_activities,
                                          plot_GeyserPlots = FALSE,
                                          plot_FeaturePlots = FALSE,
                                          split.by = "orig.ident",
                                          column_title = "A",
                                          row_title = "B")
    testthat::expect_type(out, "list")

    out <- SCpubr::do_PathwayActivityPlot(sample = sample,
                                          activities = progeny_activities,
                                          plot_GeyserPlots = FALSE,
                                          plot_FeaturePlots = FALSE,
                                          group.by = "orig.ident",
                                          column_title = "A",
                                          row_title = "B")
    testthat::expect_type(out, "list")

    out <- SCpubr::do_PathwayActivityPlot(sample = sample,
                                          activities = progeny_activities,
                                          plot_GeyserPlots = FALSE,
                                          plot_FeaturePlots = FALSE,
                                          column_title = "A",
                                          row_title = "B")
    testthat::expect_type(out, "list")
  })

}
