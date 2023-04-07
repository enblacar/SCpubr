if (isFALSE(dep_check[["do_PathwayActivityPlot"]])){

  testthat::test_that("do_PathwayActivityPlot: CRAN essentials", {

    out <- SCpubr::do_PathwayActivityPlot(sample = sample,
                                          activities = progeny_activities)
    testthat::expect_type(out, "list")
  })


  testthat::test_that("do_PathwayActivityPlot: PASS - minimal input", {
    testthat::skip_on_cran()



    out <- SCpubr::do_PathwayActivityPlot(sample = sample,
                                          activities = progeny_activities)
    testthat::expect_type(out, "list")

    out <- SCpubr::do_PathwayActivityPlot(sample = sample,
                                          activities = progeny_activities,
                                          flip = TRUE)
    testthat::expect_type(out, "list")

    out <- SCpubr::do_PathwayActivityPlot(sample = sample,
                                          activities = progeny_activities,
                                          legend.position = "right")
    testthat::expect_type(out, "list")

  })

  

  

 

  testthat::test_that("do_PathwayActivityPlot: PASS - all group.by", {
    testthat::skip_on_cran()



    out <- SCpubr::do_PathwayActivityPlot(sample = sample,
                                          activities = progeny_activities,
                                          group.by = "orig.ident")
    testthat::expect_type(out, "list")
  })


  testthat::test_that("do_PathwayActivityPlot: PASS - all split.by 2", {
    testthat::skip_on_cran()

    out <- SCpubr::do_PathwayActivityPlot(sample = sample,
                                          activities = progeny_activities,
                                          min.cutoff = -0.1,
                                          max.cutoff = NA)
    testthat::expect_type(out, "list")

    out <- SCpubr::do_PathwayActivityPlot(sample = sample,
                                          activities = progeny_activities,
                                          min.cutoff = NA,
                                          max.cutoff = 0.1)
    testthat::expect_type(out, "list")

    out <- SCpubr::do_PathwayActivityPlot(sample = sample,
                                          activities = progeny_activities,
                                          min.cutoff = -0.1)
    testthat::expect_type(out, "list")



    out <- SCpubr::do_PathwayActivityPlot(sample = sample,
                                          activities = progeny_activities,
                                          max.cutoff = 0.1)
    testthat::expect_type(out, "list")




    out <- SCpubr::do_PathwayActivityPlot(sample = sample,
                                          activities = progeny_activities,
                                          max.cutoff = 0.1,
                                          min.cutoff = -0.1)
    testthat::expect_type(out, "list")
  })



  testthat::test_that("do_PathwayActivityPlot: FAIL", {
    testthat::skip_on_cran()

    testthat::expect_error({SCpubr::do_PathwayActivityPlot(sample = sample,
                                                           activities = progeny_activities,
                                                           min.cutoff = -10)})

    testthat::expect_error({SCpubr::do_PathwayActivityPlot(sample = sample,
                                                           activities = progeny_activities,
                                                           max.cutoff = 200)})

    testthat::expect_error({SCpubr::do_PathwayActivityPlot(sample = sample,
                                                           activities = progeny_activities,
                                                           max.cutoff = 1,
                                                           min.cutoff = 2)})

  })
}
