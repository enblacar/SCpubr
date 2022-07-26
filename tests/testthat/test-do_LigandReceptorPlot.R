sample <- SCpubr:::use_dataset()
liana_output <- SCpubr:::liana_output


testthat::test_that("do_LigandReceptorPlot: PASS - from output", {
  p <- SCpubr::do_LigandReceptorPlot(from_output = TRUE,
                                     liana_output = liana_output)
  testthat::expect_type(p, "list")
})

testthat::test_that("do_LigandReceptorPlot: PASS - from output different n", {
  p <- SCpubr::do_LigandReceptorPlot(from_output = TRUE,
                                     liana_output = liana_output,
                                     top_interactions = 50)
  testthat::expect_type(p, "list")
})

testthat::test_that("do_LigandReceptorPlot: PASS - from output, angle ", {
  p <- SCpubr::do_LigandReceptorPlot(from_output = TRUE,
                                     liana_output = liana_output,
                                     x_labels_angle = 0)
  testthat::expect_type(p, "list")

  p <- SCpubr::do_LigandReceptorPlot(from_output = TRUE,
                                     liana_output = liana_output,
                                     x_labels_angle = 45)
  testthat::expect_type(p, "list")

  p <- SCpubr::do_LigandReceptorPlot(from_output = TRUE,
                                     liana_output = liana_output,
                                     x_labels_angle = 90)
  testthat::expect_type(p, "list")
})

testthat::test_that("do_LigandReceptorPlot: PASS - from output flip", {
  p <- SCpubr::do_LigandReceptorPlot(from_output = TRUE,
                                     liana_output = liana_output,
                                     top_interactions = 50,
                                     flip = TRUE)
  testthat::expect_type(p, "list")
})

testthat::test_that("do_LigandReceptorPlot: PASS - from output different keep", {
  p <- SCpubr::do_LigandReceptorPlot(from_output = TRUE,
                                     liana_output = liana_output,
                                     top_interactions = 50,
                                     keep_source = "0",
                                     keep_target = "0")
  testthat::expect_type(p, "list")

  p <- SCpubr::do_LigandReceptorPlot(from_output = TRUE,
                                     liana_output = liana_output,
                                     top_interactions = 50,
                                     keep_source = "0",
                                     keep_target = "0",
                                     flip = TRUE)
  testthat::expect_type(p, "list")
})

testthat::test_that("do_LigandReceptorPlot: PASS - from output legend.position", {
  p <- SCpubr::do_LigandReceptorPlot(from_output = TRUE,
                                     liana_output = liana_output,
                                     top_interactions = 50,
                                     flip = TRUE,
                                     legend.position = "bottom")
  testthat::expect_type(p, "list")
  p <- SCpubr::do_LigandReceptorPlot(from_output = TRUE,
                                     liana_output = liana_output,
                                     top_interactions = 50,
                                     flip = TRUE,
                                     legend.position = "right")
  testthat::expect_type(p, "list")
})

testthat::test_that("do_LigandReceptorPlot: FAIL - wrong parameters", {
  testthat::expect_error({SCpubr::do_LigandReceptorPlot(from_output = TRUE,
                                                        liana_output = liana_output,
                                                        font.type = "wrong")})

  testthat::expect_error({SCpubr::do_LigandReceptorPlot(from_output = TRUE,
                                                        liana_output = liana_output,
                                                        legend.type = "wrong")})

  testthat::expect_error({SCpubr::do_LigandReceptorPlot(from_output = TRUE,
                                                        liana_output = liana_output,
                                                        x_labels_angle  = 10)})

  testthat::expect_error({SCpubr::do_LigandReceptorPlot(from_output = TRUE,
                                                        liana_output = liana_output,
                                                        font.type = "wrong")})


})
