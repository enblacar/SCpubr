if (base::isFALSE(dep_check[["do_PathwayActivityHeatmap"]])){

  testthat::test_that("do_PathwayActivityHeatmap: CRAN essentials", {

    out <- SCpubr::do_PathwayActivityHeatmap(sample = sample,
                                          activities = progeny_activities)
    testthat::expect_true(ggplot2::is_ggplot(out))
  })


  testthat::test_that("do_PathwayActivityHeatmap: PASS - minimal input", {
    testthat::skip_on_cran()


    sample$annotation <- sample(c("A", "B"), ncol(sample), replace = TRUE)
    out <- SCpubr::do_PathwayActivityHeatmap(sample = sample,
                                          activities = progeny_activities)
    testthat::expect_true(ggplot2::is_ggplot(out))
    
    out <- SCpubr::do_PathwayActivityHeatmap(sample = sample,
                                             activities = progeny_activities,
                                             values.show = TRUE,
                                             values.threshold = 0.2)
    testthat::expect_true(ggplot2::is_ggplot(out))
    
    sample$annotation <- sample(c("A", "B"), ncol(sample), replace = TRUE)
    out <- SCpubr::do_PathwayActivityHeatmap(sample = sample,
                                             activities = progeny_activities,
                                             enforce_symmetry = TRUE)
    testthat::expect_true(ggplot2::is_ggplot(out))
    
    sample$annotation <- sample(c("A", "B"), ncol(sample), replace = TRUE)
    out <- SCpubr::do_PathwayActivityHeatmap(sample = sample,
                                             activities = progeny_activities,
                                             enforce_symmetry = FALSE)
    testthat::expect_true(ggplot2::is_ggplot(out))
    
    out <- SCpubr::do_PathwayActivityHeatmap(sample = sample,
                                     activities = progeny_activities,
                                     group.by = "orig.ident")
    testthat::expect_true(ggplot2::is_ggplot(out))
    
    out <- SCpubr::do_PathwayActivityHeatmap(sample = sample,
                                     activities = progeny_activities,
                                     group.by = c("orig.ident", "seurat_clusters", "annotation"))
    testthat::expect_true(ggplot2::is_ggplot(out))
    
    out <- SCpubr::do_PathwayActivityHeatmap(sample = sample,
                                          activities = progeny_activities,
                                          split.by = "annotation")
    testthat::expect_true(ggplot2::is_ggplot(out))
    

    out <- SCpubr::do_PathwayActivityHeatmap(sample = sample,
                                          activities = progeny_activities,
                                          flip = TRUE)
    testthat::expect_true(ggplot2::is_ggplot(out))
    
    out <- SCpubr::do_PathwayActivityHeatmap(sample = sample,
                                          activities = progeny_activities,
                                          flip = FALSE)
    testthat::expect_true(ggplot2::is_ggplot(out))

    out <- SCpubr::do_PathwayActivityHeatmap(sample = sample,
                                          activities = progeny_activities,
                                          legend.position = "right")
    testthat::expect_true(ggplot2::is_ggplot(out))
    
    out <- SCpubr::do_PathwayActivityHeatmap(sample = sample,
                                          activities = progeny_activities,
                                          legend.position = "right",
                                          return_object = TRUE)
    testthat::expect_type(out, "list")

  })

  

  

 

  testthat::test_that("do_PathwayActivityHeatmap: PASS - all group.by", {
    testthat::skip_on_cran()



    out <- SCpubr::do_PathwayActivityHeatmap(sample = sample,
                                          activities = progeny_activities,
                                          group.by = "orig.ident")
    testthat::expect_true(ggplot2::is_ggplot(out))
  })


  testthat::test_that("do_PathwayActivityHeatmap: PASS - all split.by 2", {
    testthat::skip_on_cran()

    out <- SCpubr::do_PathwayActivityHeatmap(sample = sample,
                                          activities = progeny_activities,
                                          min.cutoff = -0.1,
                                          max.cutoff = NA)
    testthat::expect_true(ggplot2::is_ggplot(out))

    out <- SCpubr::do_PathwayActivityHeatmap(sample = sample,
                                          activities = progeny_activities,
                                          min.cutoff = NA,
                                          max.cutoff = 0.1)
    testthat::expect_true(ggplot2::is_ggplot(out))

    out <- SCpubr::do_PathwayActivityHeatmap(sample = sample,
                                          activities = progeny_activities,
                                          min.cutoff = -0.1)
    testthat::expect_true(ggplot2::is_ggplot(out))



    out <- SCpubr::do_PathwayActivityHeatmap(sample = sample,
                                          activities = progeny_activities,
                                          max.cutoff = 0.1)
    testthat::expect_true(ggplot2::is_ggplot(out))




    out <- SCpubr::do_PathwayActivityHeatmap(sample = sample,
                                          activities = progeny_activities,
                                          max.cutoff = 0.1,
                                          min.cutoff = -0.1)
    testthat::expect_true(ggplot2::is_ggplot(out))
  })



  testthat::test_that("do_PathwayActivityHeatmap: FAIL", {
    testthat::skip_on_cran()

    testthat::expect_error({SCpubr::do_PathwayActivityHeatmap(sample = sample,
                                                           activities = progeny_activities,
                                                           min.cutoff = -10)})

    testthat::expect_error({SCpubr::do_PathwayActivityHeatmap(sample = sample,
                                                           activities = progeny_activities,
                                                           max.cutoff = 200)})

    testthat::expect_error({SCpubr::do_PathwayActivityHeatmap(sample = sample,
                                                           activities = progeny_activities,
                                                           max.cutoff = 1,
                                                           min.cutoff = 2)})
    
    sample$annotation <- sample(c("A", "B"), ncol(sample), replace = TRUE)
    testthat::expect_error({SCpubr::do_PathwayActivityHeatmap(sample = sample,
                                                      activities = progeny_activities,
                                                      group.by = c("seurat_clusters", "orig.ident"),
                                                      split.by = "annotation")})

  })
}
