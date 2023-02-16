if(isFALSE(dep_check[["do_LigandReceptorPlot"]])){
  testthat::test_that("do_LigandReceptorPlot: CRAN essentials", {
    p <- SCpubr::do_LigandReceptorPlot(liana_output = liana_output, add_missing_LR_combinations = FALSE)
    testthat::expect_type(p, "list")
  })

  testthat::test_that("do_LigandReceptorPlot: PASS - from output", {
    testthat::skip_on_cran()
    p <- SCpubr::do_LigandReceptorPlot(liana_output = liana_output, add_missing_LR_combinations = FALSE)
    testthat::expect_type(p, "list")

    p <- SCpubr::do_LigandReceptorPlot(liana_output = liana_output, add_missing_LR_combinations = FALSE, plot.grid = FALSE)
    testthat::expect_type(p, "list")

    p <- SCpubr::do_LigandReceptorPlot(liana_output = liana_output, add_missing_LR_combinations = FALSE, plot.grid = TRUE)
    testthat::expect_type(p, "list")

    p <- SCpubr::do_LigandReceptorPlot(liana_output = liana_output, add_missing_LR_combinations = TRUE)
    testthat::expect_type(p, "list")

    p <- SCpubr::do_LigandReceptorPlot(liana_output = liana_output, add_missing_LR_combinations = TRUE, keep_source = c("0", "3", "5"), keep_target = c("0", "2", "4"))
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
    testthat::skip_on_cran()
    p <- SCpubr::do_LigandReceptorPlot(liana_output = liana_output,
                                       top_interactions = 50,
                                       add_missing_LR_combinations = FALSE)
    testthat::expect_type(p, "list")
  })

  testthat::test_that("do_LigandReceptorPlot: PASS - legend.type", {
    testthat::skip_on_cran()
    p <- SCpubr::do_LigandReceptorPlot(liana_output = liana_output,
                                       legend.type = "normal",
                                       add_missing_LR_combinations = FALSE)
    testthat::expect_type(p, "list")

    p <- SCpubr::do_LigandReceptorPlot(liana_output = liana_output,
                                       legend.type = "colorbar",
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

  })

  testthat::test_that("do_LigandReceptorPlot: PASS - split.by", {
    testthat::skip_on_cran()
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
    testthat::skip_on_cran()
    p <- SCpubr::do_LigandReceptorPlot(liana_output = liana_output,
                                       rotate_x_axis_labels = 0,
                                       add_missing_LR_combinations = FALSE)
    testthat::expect_type(p, "list")

    p <- SCpubr::do_LigandReceptorPlot(liana_output = liana_output,
                                       rotate_x_axis_labels = 45,
                                       add_missing_LR_combinations = FALSE)
    testthat::expect_type(p, "list")

    p <- SCpubr::do_LigandReceptorPlot(liana_output = liana_output,
                                       rotate_x_axis_labels = 90,
                                       add_missing_LR_combinations = FALSE)
    testthat::expect_type(p, "list")
  })

  testthat::test_that("do_LigandReceptorPlot: PASS - from output flip", {
    testthat::skip_on_cran()
    p <- SCpubr::do_LigandReceptorPlot(liana_output = liana_output,
                                       top_interactions = 50,
                                       flip = TRUE,
                                       add_missing_LR_combinations = FALSE)
    testthat::expect_type(p, "list")
  })

  testthat::test_that("do_LigandReceptorPlot: PASS - from output different keep", {
    testthat::skip_on_cran()
    p <- SCpubr::do_LigandReceptorPlot(liana_output = liana_output,
                                       top_interactions = 50,
                                       keep_source = "5",
                                       keep_target = "9",
                                       add_missing_LR_combinations = FALSE)
    testthat::expect_type(p, "list")

    p <- SCpubr::do_LigandReceptorPlot(liana_output = liana_output,
                                       top_interactions = 50,
                                       keep_source = "5",
                                       keep_target = "9",
                                       flip = TRUE,
                                       add_missing_LR_combinations = FALSE)
    testthat::expect_type(p, "list")
  })

  testthat::test_that("do_LigandReceptorPlot: PASS - from output legend.position", {
    testthat::skip_on_cran()
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


  testthat::test_that("do_LigandReceptorPlot: PASS - arrange interactions", {
    testthat::skip_on_cran()
    p <- SCpubr::do_LigandReceptorPlot(liana_output = liana_output,
                                       top_interactions = 50,
                                       arrange_interactions_by = "aggregate_rank")
    testthat::expect_type(p, "list")

    p <- SCpubr::do_LigandReceptorPlot(liana_output = liana_output,
                                       top_interactions = 50,
                                       arrange_interactions_by = "specificity")
    testthat::expect_type(p, "list")

    p <- SCpubr::do_LigandReceptorPlot(liana_output = liana_output,
                                       top_interactions = 50,
                                       arrange_interactions_by = "magnitude")
    testthat::expect_type(p, "list")

    p <- SCpubr::do_LigandReceptorPlot(liana_output = liana_output,
                                       top_interactions = 50,
                                       arrange_interactions_by = "both")
    testthat::expect_type(p, "list")
  })

  testthat::test_that("do_LigandReceptorPlot: PASS - sort interactions", {
    testthat::skip_on_cran()
    p <- SCpubr::do_LigandReceptorPlot(liana_output = liana_output,
                                       top_interactions = 50,
                                       sort_interactions_alphabetically =  TRUE)
    testthat::expect_type(p, "list")

    p <- SCpubr::do_LigandReceptorPlot(liana_output = liana_output,
                                       top_interactions = 50,
                                       sort_interactions_alphabetically =  FALSE)
    testthat::expect_type(p, "list")
  })

  testthat::test_that("do_LigandReceptorPlot: FAIL - wrong parameters", {
    testthat::skip_on_cran()
    testthat::expect_error({SCpubr::do_LigandReceptorPlot(liana_output = liana_output,
                                                          font.type = "wrong",
                                                          add_missing_LR_combinations = FALSE)})

    testthat::expect_error({SCpubr::do_LigandReceptorPlot(liana_output = liana_output,
                                                          legend.type = "wrong",
                                                          add_missing_LR_combinations = FALSE)})

    testthat::expect_error({SCpubr::do_LigandReceptorPlot(liana_output = liana_output,
                                                          rotate_x_axis_labels  = 10,
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
    testthat::skip_on_cran()
    out <- SCpubr::do_LigandReceptorPlot(liana_output = liana_output, add_missing_LR_combinations = FALSE, compute_ChordDiagrams = TRUE)
    testthat::expect_type(out, "list")
    testthat::expect_length(out, 3)
    testthat::expect_s3_class(out$dotplot, c("gg", "ggplot"))
    testthat::expect_s3_class(out$chord_total_interactions, c("recordedplot"))
    testthat::expect_s3_class(out$chord_ligand_receptor, c("recordedplot"))
  })
}

