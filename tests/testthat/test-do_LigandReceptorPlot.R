sample <- SCpubr:::use_dataset()
test_list <- SCpubr:::test_list
liana_output <- test_list$liana_output


testthat::test_that("do_LigandReceptorPlot: PASS - from output", {
  p <- SCpubr::do_LigandReceptorPlot(liana_output = liana_output)
  testthat::expect_type(p, "list")

  p <- SCpubr::do_LigandReceptorPlot(liana_output = liana_output,
                                     dot_border = FALSE)
  testthat::expect_type(p, "list")
  p <- SCpubr::do_LigandReceptorPlot(liana_output = liana_output,
                                     flip = TRUE)
  testthat::expect_type(p, "list")

  p <- SCpubr::do_LigandReceptorPlot(liana_output = liana_output,
                                     rotate_strip_text = TRUE)
  testthat::expect_type(p, "list")
  p <- SCpubr::do_LigandReceptorPlot(liana_output = liana_output,
                                     rotate_strip_text = TRUE,
                                     flip = TRUE)
  testthat::expect_type(p, "list")
  p <- SCpubr::do_LigandReceptorPlot(liana_output = liana_output,
                                     rotate_strip_text = TRUE,
                                     flip = TRUE,
                                     dot_border = FALSE)
  testthat::expect_type(p, "list")
})

testthat::test_that("do_LigandReceptorPlot: PASS - from output different n", {
  p <- SCpubr::do_LigandReceptorPlot(liana_output = liana_output,
                                     top_interactions = 50)
  testthat::expect_type(p, "list")
})

testthat::test_that("do_LigandReceptorPlot: PASS - legend.type", {
  p <- SCpubr::do_LigandReceptorPlot(liana_output = liana_output,
                                     legend.type = "normal")
  testthat::expect_type(p, "list")

  p <- SCpubr::do_LigandReceptorPlot(liana_output = liana_output,
                                     legend.type = "colorbar")
  testthat::expect_type(p, "list")

  p <- SCpubr::do_LigandReceptorPlot(liana_output = liana_output,
                                     legend.type = "colorsteps")
  testthat::expect_type(p, "list")

  p <- SCpubr::do_LigandReceptorPlot(liana_output = liana_output,
                                     legend.type = "colorsteps",
                                     dot_border = FALSE)
  testthat::expect_type(p, "list")

  p <- SCpubr::do_LigandReceptorPlot(liana_output = liana_output,
                                     legend.type = "normal",
                                     dot_border = FALSE)
  testthat::expect_type(p, "list")

  p <- SCpubr::do_LigandReceptorPlot(liana_output = liana_output,
                                     legend.type = "colorbar",
                                     dot_border = FALSE)
  testthat::expect_type(p, "list")

  p <- SCpubr::do_LigandReceptorPlot(liana_output = liana_output,
                                     legend.type = "colorsteps")
  testthat::expect_type(p, "list")
})

testthat::test_that("do_LigandReceptorPlot: PASS - split.by", {
  p <- SCpubr::do_LigandReceptorPlot(liana_output = liana_output,
                                     split.by = "ligand.complex",
                                     dot_border = FALSE)
  testthat::expect_type(p, "list")

  p <- SCpubr::do_LigandReceptorPlot(liana_output = liana_output,
                                     split.by = "receptor.complex")
  testthat::expect_type(p, "list")

  p <- SCpubr::do_LigandReceptorPlot(liana_output = liana_output,
                                     split.by = "ligand.complex",
                                     flip = TRUE)
  testthat::expect_type(p, "list")

  p <- SCpubr::do_LigandReceptorPlot(liana_output = liana_output,
                                     split.by = "receptor.complex",
                                     flip = TRUE)
  testthat::expect_type(p, "list")
})

testthat::test_that("do_LigandReceptorPlot: PASS - from output, angle ", {
  p <- SCpubr::do_LigandReceptorPlot(liana_output = liana_output,
                                     x_labels_angle = 0)
  testthat::expect_type(p, "list")

  p <- SCpubr::do_LigandReceptorPlot(liana_output = liana_output,
                                     x_labels_angle = 45)
  testthat::expect_type(p, "list")

  p <- SCpubr::do_LigandReceptorPlot(liana_output = liana_output,
                                     x_labels_angle = 90)
  testthat::expect_type(p, "list")
})

testthat::test_that("do_LigandReceptorPlot: PASS - from output flip", {
  p <- SCpubr::do_LigandReceptorPlot(liana_output = liana_output,
                                     top_interactions = 50,
                                     flip = TRUE)
  testthat::expect_type(p, "list")
})

testthat::test_that("do_LigandReceptorPlot: PASS - from output different keep", {
  p <- SCpubr::do_LigandReceptorPlot(liana_output = liana_output,
                                     top_interactions = 50,
                                     keep_source = "0",
                                     keep_target = "0")
  testthat::expect_type(p, "list")

  p <- SCpubr::do_LigandReceptorPlot(liana_output = liana_output,
                                     top_interactions = 50,
                                     keep_source = "0",
                                     keep_target = "0",
                                     flip = TRUE)
  testthat::expect_type(p, "list")
})

testthat::test_that("do_LigandReceptorPlot: PASS - from output legend.position", {
  p <- SCpubr::do_LigandReceptorPlot(liana_output = liana_output,
                                     top_interactions = 50,
                                     flip = TRUE,
                                     legend.position = "bottom")
  testthat::expect_type(p, "list")
  p <- SCpubr::do_LigandReceptorPlot(liana_output = liana_output,
                                     top_interactions = 50,
                                     flip = TRUE,
                                     legend.position = "right")
  testthat::expect_type(p, "list")
})

testthat::test_that("do_LigandReceptorPlot: FAIL - wrong parameters", {
  testthat::expect_error({SCpubr::do_LigandReceptorPlot(liana_output = liana_output,
                                                        font.type = "wrong")})

  testthat::expect_error({SCpubr::do_LigandReceptorPlot(liana_output = liana_output,
                                                        legend.type = "wrong")})

  testthat::expect_error({SCpubr::do_LigandReceptorPlot(liana_output = liana_output,
                                                        x_labels_angle  = 10)})

  testthat::expect_error({SCpubr::do_LigandReceptorPlot(liana_output = liana_output,
                                                        font.type = "wrong")})

  testthat::expect_error({SCpubr::do_LigandReceptorPlot(liana_output = liana_output,
                                                        legend.position = "wrong")})

  testthat::expect_error({SCpubr::do_LigandReceptorPlot(liana_output = liana_output,
                                                        grid.type = "wrong")})

  testthat::expect_error({SCpubr::do_LigandReceptorPlot(liana_output = liana_output,
                                                        split.by = "wrong")})

})
