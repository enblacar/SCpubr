if (isFALSE(dep_check[["do_GroupedGOTermPlot"]])){
  testthat::test_that("do_GroupedGOTermPlot: CRAN essential tests", {

    p <- SCpubr::do_GroupedGOTermPlot(genes = c("MBP"),
                                      org.db = org.Hs.eg.db,
                                      GO_ontology = "BP",
                                      levels.use = c(1, 2),
                                      verbose = FALSE,
                                      min.overlap = 1)

    testthat::expect_type(p, "list")
  })

}
