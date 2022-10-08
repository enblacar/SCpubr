library(SCpubr)
suppressWarnings(rm(test_list))
test_list <- data("test.data", package = "SCpubr")
sample <- test_list$sample
liana_output <- SCpubr:::test_list$liana_output

testthat::test_that("do_LigandReceptorPlot: PASS - from output", {

  p <- SCpubr::do_LigandReceptorPlot(liana_output = liana_output, add_missing_LR_combinations = FALSE)
  testthat::expect_type(p, "list")

  p <- SCpubr::do_LigandReceptorPlot(liana_output = liana_output, add_missing_LR_combinations = TRUE)
  testthat::expect_type(p, "list")

  p <- SCpubr::do_LigandReceptorPlot(liana_output = liana_output,
                                     dot_border = FALSE,
                                     add_missing_LR_combinations = FALSE)
  testthat::expect_type(p, "list")
  p <- SCpubr::do_LigandReceptorPlot(liana_output = liana_output,
                                     flip = TRUE,
                                     add_missing_LR_combinations = FALSE)
  testthat::expect_type(p, "list")

  p <- SCpubr::do_LigandReceptorPlot(liana_output = liana_output,
                                     rotate_strip_text = TRUE,
                                     add_missing_LR_combinations = FALSE)
  testthat::expect_type(p, "list")
  p <- SCpubr::do_LigandReceptorPlot(liana_output = liana_output,
                                     rotate_strip_text = TRUE,
                                     flip = TRUE,
                                     add_missing_LR_combinations = FALSE)
  testthat::expect_type(p, "list")
  p <- SCpubr::do_LigandReceptorPlot(liana_output = liana_output,
                                     rotate_strip_text = TRUE,
                                     flip = TRUE,
                                     dot_border = FALSE,
                                     add_missing_LR_combinations = FALSE)
  testthat::expect_type(p, "list")
})

testthat::test_that("do_LigandReceptorPlot: PASS - from output different n", {

  p <- SCpubr::do_LigandReceptorPlot(liana_output = liana_output,
                                     top_interactions = 50,
                                     add_missing_LR_combinations = FALSE)
  testthat::expect_type(p, "list")
})

testthat::test_that("do_LigandReceptorPlot: PASS - legend.type", {

  p <- SCpubr::do_LigandReceptorPlot(liana_output = liana_output,
                                     legend.type = "normal",
                                     add_missing_LR_combinations = FALSE)
  testthat::expect_type(p, "list")

  p <- SCpubr::do_LigandReceptorPlot(liana_output = liana_output,
                                     legend.type = "colorbar",
                                     add_missing_LR_combinations = FALSE)
  testthat::expect_type(p, "list")

  p <- SCpubr::do_LigandReceptorPlot(liana_output = liana_output,
                                     legend.type = "colorsteps",
                                     add_missing_LR_combinations = FALSE)
  testthat::expect_type(p, "list")

  p <- SCpubr::do_LigandReceptorPlot(liana_output = liana_output,
                                     legend.type = "colorsteps",
                                     dot_border = FALSE,
                                     add_missing_LR_combinations = FALSE)
  testthat::expect_type(p, "list")

  p <- SCpubr::do_LigandReceptorPlot(liana_output = liana_output,
                                     legend.type = "normal",
                                     dot_border = FALSE,
                                     add_missing_LR_combinations = FALSE)
  testthat::expect_type(p, "list")

  p <- SCpubr::do_LigandReceptorPlot(liana_output = liana_output,
                                     legend.type = "colorbar",
                                     dot_border = FALSE,
                                     add_missing_LR_combinations = FALSE)
  testthat::expect_type(p, "list")

  p <- SCpubr::do_LigandReceptorPlot(liana_output = liana_output,
                                     legend.type = "colorsteps",
                                     add_missing_LR_combinations = FALSE)
  testthat::expect_type(p, "list")
})

testthat::test_that("do_LigandReceptorPlot: PASS - split.by", {

  p <- SCpubr::do_LigandReceptorPlot(liana_output = liana_output,
                                     split.by = "ligand.complex",
                                     dot_border = FALSE,
                                     add_missing_LR_combinations = FALSE)
  testthat::expect_type(p, "list")

  p <- SCpubr::do_LigandReceptorPlot(liana_output = liana_output,
                                     split.by = "receptor.complex",
                                     add_missing_LR_combinations = FALSE)
  testthat::expect_type(p, "list")

  p <- SCpubr::do_LigandReceptorPlot(liana_output = liana_output,
                                     split.by = "ligand.complex",
                                     flip = TRUE,
                                     add_missing_LR_combinations = FALSE)
  testthat::expect_type(p, "list")

  p <- SCpubr::do_LigandReceptorPlot(liana_output = liana_output,
                                     split.by = "receptor.complex",
                                     flip = TRUE)
  testthat::expect_type(p, "list")
})

testthat::test_that("do_LigandReceptorPlot: PASS - from output, angle ", {

  p <- SCpubr::do_LigandReceptorPlot(liana_output = liana_output,
                                     x_labels_angle = 0,
                                     add_missing_LR_combinations = FALSE)
  testthat::expect_type(p, "list")

  p <- SCpubr::do_LigandReceptorPlot(liana_output = liana_output,
                                     x_labels_angle = 45,
                                     add_missing_LR_combinations = FALSE)
  testthat::expect_type(p, "list")

  p <- SCpubr::do_LigandReceptorPlot(liana_output = liana_output,
                                     x_labels_angle = 90,
                                     add_missing_LR_combinations = FALSE)
  testthat::expect_type(p, "list")
})

testthat::test_that("do_LigandReceptorPlot: PASS - from output flip", {

  p <- SCpubr::do_LigandReceptorPlot(liana_output = liana_output,
                                     top_interactions = 50,
                                     flip = TRUE,
                                     add_missing_LR_combinations = FALSE)
  testthat::expect_type(p, "list")
})

testthat::test_that("do_LigandReceptorPlot: PASS - from output different keep", {

  p <- SCpubr::do_LigandReceptorPlot(liana_output = liana_output,
                                     top_interactions = 50,
                                     keep_source = "0",
                                     keep_target = "0",
                                     add_missing_LR_combinations = FALSE)
  testthat::expect_type(p, "list")

  p <- SCpubr::do_LigandReceptorPlot(liana_output = liana_output,
                                     top_interactions = 50,
                                     keep_source = "0",
                                     keep_target = "0",
                                     flip = TRUE,
                                     add_missing_LR_combinations = FALSE)
  testthat::expect_type(p, "list")
})

testthat::test_that("do_LigandReceptorPlot: PASS - from output legend.position", {
  sample <- SCpubr:::test_list$sample
  liana_output <- SCpubr:::test_list$liana_output

  p <- SCpubr::do_LigandReceptorPlot(liana_output = liana_output,
                                     top_interactions = 50,
                                     flip = TRUE,
                                     legend.position = "bottom",
                                     add_missing_LR_combinations = FALSE)
  testthat::expect_type(p, "list")
  p <- SCpubr::do_LigandReceptorPlot(liana_output = liana_output,
                                     top_interactions = 50,
                                     flip = TRUE,
                                     legend.position = "right",
                                     add_missing_LR_combinations = FALSE)
  testthat::expect_type(p, "list")
})

testthat::test_that("do_LigandReceptorPlot: FAIL - wrong parameters", {

  testthat::expect_error({SCpubr::do_LigandReceptorPlot(liana_output = liana_output,
                                                        font.type = "wrong",
                                                        add_missing_LR_combinations = FALSE)})

  testthat::expect_error({SCpubr::do_LigandReceptorPlot(liana_output = liana_output,
                                                        legend.type = "wrong",
                                                        add_missing_LR_combinations = FALSE)})

  testthat::expect_error({SCpubr::do_LigandReceptorPlot(liana_output = liana_output,
                                                        x_labels_angle  = 10,
                                                        add_missing_LR_combinations = FALSE)})

  testthat::expect_error({SCpubr::do_LigandReceptorPlot(liana_output = liana_output,
                                                        font.type = "wrong",
                                                        add_missing_LR_combinations = FALSE)})

  testthat::expect_error({SCpubr::do_LigandReceptorPlot(liana_output = liana_output,
                                                        legend.position = "wrong",
                                                        add_missing_LR_combinations = FALSE)})

  testthat::expect_error({SCpubr::do_LigandReceptorPlot(liana_output = liana_output,
                                                        grid.type = "wrong",
                                                        add_missing_LR_combinations = FALSE)})

  testthat::expect_error({SCpubr::do_LigandReceptorPlot(liana_output = liana_output,
                                                        split.by = "wrong",
                                                        add_missing_LR_combinations = FALSE)})

})

testthat::test_that("do_LigandReceptorPlot: PASS - chord diagrams", {

  out <- SCpubr::do_LigandReceptorPlot(liana_output = liana_output, add_missing_LR_combinations = FALSE, compute_ChordDiagrams = TRUE)
  testthat::expect_type(out, "list")
  testthat::expect_length(out, 3)
  testthat::expect_equal(class(out$dotplot), c("gg", "ggplot"))
  testthat::expect_equal(class(out$chord_total_interactions), c("recordedplot"))
  testthat::expect_equal(class(out$chord_ligand_receptor), c("recordedplot"))
})
