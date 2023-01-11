if(isFALSE(dep_check[["do_TFActivityPlot"]])){

  testthat::test_that("do_TFActivityPlot: PASS - minimal input", {
    out <- SCpubr::do_TFActivityPlot(sample = sample,
                                     activities = dorothea_activities)
    testthat::expect_type(out, "list")

  })

  testthat::test_that("do_TFActivityPlot: PASS - minimal input", {
    testthat::skip_on_cran()



    out <- SCpubr::do_TFActivityPlot(sample = sample,
                                     activities = dorothea_activities)
    testthat::expect_type(out, "list")

    out <- SCpubr::do_TFActivityPlot(sample = sample,
                                     activities = dorothea_activities,
                                     flip = TRUE)
    testthat::expect_type(out, "list")

    out <- SCpubr::do_TFActivityPlot(sample = sample,
                                     activities = dorothea_activities,
                                     flip = TRUE,
                                     split.by = "orig.ident")
    testthat::expect_type(out, "list")

    out <- SCpubr::do_TFActivityPlot(sample = sample,
                                     activities = dorothea_activities,
                                     legend.position = "right")
    testthat::expect_type(out, "list")

    out <- SCpubr::do_TFActivityPlot(sample = sample,
                                     activities = dorothea_activities,
                                     legend.position = "right",
                                     split.by = "orig.ident")
    testthat::expect_type(out, "list")
  })

  testthat::test_that("do_TFActivityPlot: PASS - minimal input", {
    testthat::skip_on_cran()



    out <- SCpubr::do_TFActivityPlot(sample = sample,
                                     activities = dorothea_activities,
                                     flip = TRUE)
    testthat::expect_type(out, "list")
  })

  testthat::test_that("do_TFActivityPlot: PASS - plot featureplots", {
    testthat::skip_on_cran()



    out <- SCpubr::do_TFActivityPlot(sample = sample,
                                     activities = dorothea_activities,
                                     plot_FeaturePlots = TRUE)
    testthat::expect_type(out, "list")
    testthat::expect_length(out, 2)
  })

  testthat::test_that("do_TFActivityPlot: PASS - plot geysers", {
    testthat::skip_on_cran()



    out <- SCpubr::do_TFActivityPlot(sample = sample,
                                     activities = dorothea_activities,
                                     plot_GeyserPlots = TRUE)
    testthat::expect_type(out, "list")
    testthat::expect_length(out, 2)
  })

  testthat::test_that("do_TFActivityPlot: PASS - all", {
    testthat::skip_on_cran()



    out <- SCpubr::do_TFActivityPlot(sample = sample,
                                     activities = dorothea_activities,
                                     plot_GeyserPlots = TRUE,
                                     plot_FeaturePlots = TRUE)
    testthat::expect_type(out, "list")
    testthat::expect_length(out, 3)
  })

  testthat::test_that("do_TFActivityPlot: PASS - all group.by", {
    testthat::skip_on_cran()



    out <- SCpubr::do_TFActivityPlot(sample = sample,
                                     activities = dorothea_activities,
                                     plot_GeyserPlots = TRUE,
                                     plot_FeaturePlots = TRUE,
                                     group.by = "orig.ident")
    testthat::expect_type(out, "list")
    testthat::expect_length(out, 3)
  })

  testthat::test_that("do_TFActivityPlot: PASS - all split.by", {
    testthat::skip_on_cran()



    out <- SCpubr::do_TFActivityPlot(sample = sample,
                                     activities = dorothea_activities,
                                     plot_GeyserPlots = TRUE,
                                     plot_FeaturePlots = TRUE,
                                     split.by = "orig.ident")
    testthat::expect_type(out, "list")
    testthat::expect_length(out, 3)
  })

  testthat::test_that("do_TFActivityPlot: PASS - all split.by", {
    testthat::skip_on_cran()



    out <- SCpubr::do_TFActivityPlot(sample = sample,
                                          activities = dorothea_activities,
                                          plot_GeyserPlots = TRUE,
                                          plot_FeaturePlots = TRUE,
                                          split.by = "orig.ident")
    testthat::expect_type(out, "list")
    testthat::expect_length(out, 3)
  })

  testthat::test_that("do_TFActivityPlot: PASS - all split.by 2", {
    testthat::skip_on_cran()

    out <- SCpubr::do_TFActivityPlot(sample = sample,
                                    activities = dorothea_activities,
                                    min.cutoff = -0.1,
                                    max.cutoff = NULL,
                                    plot_FeaturePlots = TRUE)
    testthat::expect_type(out, "list")

    out <- SCpubr::do_TFActivityPlot(sample = sample,
                                     activities = dorothea_activities,
                                     min.cutoff = NULL,
                                     max.cutoff = 0.1,
                                     plot_FeaturePlots = TRUE)
    testthat::expect_type(out, "list")

    out <- SCpubr::do_TFActivityPlot(sample = sample,
                                     activities = dorothea_activities,
                                     min.cutoff = -0.1)
    testthat::expect_type(out, "list")

    out <- SCpubr::do_TFActivityPlot(sample = sample,
                                     activities = dorothea_activities,
                                     min.cutoff = -0.1,
                                     split.by = "orig.ident")
    testthat::expect_type(out, "list")


    out <- SCpubr::do_TFActivityPlot(sample = sample,
                                     activities = dorothea_activities,
                                     max.cutoff = 0.1)
    testthat::expect_type(out, "list")

    out <- SCpubr::do_TFActivityPlot(sample = sample,
                                     activities = dorothea_activities,
                                     max.cutoff = 0.1,
                                     split.by = "orig.ident")
    testthat::expect_type(out, "list")


    out <- SCpubr::do_TFActivityPlot(sample = sample,
                                     activities = dorothea_activities,
                                     max.cutoff = 0.1,
                                     min.cutoff = -0.1)
    testthat::expect_type(out, "list")

    out <- SCpubr::do_TFActivityPlot(sample = sample,
                                     activities = dorothea_activities,
                                     max.cutoff = 0.1,
                                     min.cutoff = -0.1,
                                     split.by = "orig.ident")
    testthat::expect_type(out, "list")



    out <- SCpubr::do_TFActivityPlot(sample = sample,
                                     activities = dorothea_activities,
                                     min.cutoff = NULL)
    testthat::expect_type(out, "list")

    out <- SCpubr::do_TFActivityPlot(sample = sample,
                                     activities = dorothea_activities,
                                     min.cutoff = NULL,
                                     split.by = "orig.ident")
    testthat::expect_type(out, "list")


    out <- SCpubr::do_TFActivityPlot(sample = sample,
                                     activities = dorothea_activities,
                                     max.cutoff = NULL)
    testthat::expect_type(out, "list")

    out <- SCpubr::do_TFActivityPlot(sample = sample,
                                     activities = dorothea_activities,
                                     max.cutoff = NULL,
                                     split.by = "orig.ident")
    testthat::expect_type(out, "list")


    out <- SCpubr::do_TFActivityPlot(sample = sample,
                                     activities = dorothea_activities,
                                     max.cutoff = NULL,
                                     min.cutoff = NULL)
    testthat::expect_type(out, "list")

    out <- SCpubr::do_TFActivityPlot(sample = sample,
                                     activities = dorothea_activities,
                                     max.cutoff = NULL,
                                     min.cutoff = NULL,
                                     split.by = "orig.ident")
    testthat::expect_type(out, "list")
  })

  testthat::test_that("do_TFActivityPlot: PASS - column.title and row.title", {
    testthat::skip_on_cran()



    out <- SCpubr::do_TFActivityPlot(sample = sample,
                                     activities = dorothea_activities,
                                     plot_GeyserPlots = FALSE,
                                     plot_FeaturePlots = FALSE,
                                     split.by = "orig.ident",
                                     column_title = "A",
                                     row_title = "B")
    testthat::expect_type(out, "list")

    out <- SCpubr::do_TFActivityPlot(sample = sample,
                                     activities = dorothea_activities,
                                     plot_GeyserPlots = FALSE,
                                     plot_FeaturePlots = FALSE,
                                     group.by = "orig.ident",
                                     column_title = "A",
                                     row_title = "B")
    testthat::expect_type(out, "list")

    out <- SCpubr::do_TFActivityPlot(sample = sample,
                                     activities = dorothea_activities,
                                     plot_GeyserPlots = FALSE,
                                     plot_FeaturePlots = FALSE,
                                     column_title = "A",
                                     row_title = "B")
    testthat::expect_type(out, "list")
  })

  testthat::test_that("do_PathwayActivityPlot: FAIL", {
    testthat::skip_on_cran()

    testthat::expect_error({SCpubr::do_TFActivityPlot(sample = sample,
                                                           activities = progeny_activities,
                                                           plot_GeyserPlots = FALSE,
                                                           plot_FeaturePlots = FALSE,
                                                           min.cutoff = -10)})

    testthat::expect_error({SCpubr::do_TFActivityPlot(sample = sample,
                                                           activities = progeny_activities,
                                                           plot_GeyserPlots = FALSE,
                                                           plot_FeaturePlots = FALSE,
                                                           max.cutoff = 200)})

    testthat::expect_error({SCpubr::do_TFActivityPlot(sample = sample,
                                                           activities = progeny_activities,
                                                           plot_GeyserPlots = FALSE,
                                                           plot_FeaturePlots = FALSE,
                                                           max.cutoff = 1,
                                                           min.cutoff = 2)})

  })
}

