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
                                     flip = TRUE)
    testthat::expect_type(out, "list")

    out <- SCpubr::do_TFActivityPlot(sample = sample,
                                     activities = dorothea_activities,
                                     legend.position = "right")
    testthat::expect_type(out, "list")

    out <- SCpubr::do_TFActivityPlot(sample = sample,
                                     activities = dorothea_activities,
                                     legend.position = "right")
    testthat::expect_type(out, "list")
  })

  testthat::test_that("do_TFActivityPlot: PASS - minimal input", {
    testthat::skip_on_cran()



    out <- SCpubr::do_TFActivityPlot(sample = sample,
                                     activities = dorothea_activities,
                                     flip = TRUE)
    testthat::expect_type(out, "list")
  })

  

  testthat::test_that("do_TFActivityPlot: PASS - all group.by", {
    testthat::skip_on_cran()



    out <- SCpubr::do_TFActivityPlot(sample = sample,
                                     activities = dorothea_activities,
                                     group.by = "orig.ident")
    testthat::expect_type(out, "list")
  })




  testthat::test_that("do_PathwayActivityPlot: FAIL", {
    testthat::skip_on_cran()

    testthat::expect_error({SCpubr::do_TFActivityPlot(sample = sample,
                                                           activities = progeny_activities,
                                                           min.cutoff = -10)})

    testthat::expect_error({SCpubr::do_TFActivityPlot(sample = sample,
                                                           activities = progeny_activities,
                                                           max.cutoff = 200)})

    testthat::expect_error({SCpubr::do_TFActivityPlot(sample = sample,
                                                           activities = progeny_activities,
                                                           max.cutoff = 1,
                                                           min.cutoff = 2)})

  })
}

