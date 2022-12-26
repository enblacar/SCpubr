if (isFALSE(dep_check[["do_FunctionalAnnotationPlot"]])){
  testthat::test_that("do_FunctionalAnnotationPlot: CRAN essential tests", {

    p <- SCpubr::do_FunctionalAnnotationPlot(genes = c("MBP"),
                                      org.db = org.db,
                                      database = "GO",
                                      GO_ontology = "BP",
                                      min.overlap = 1)

    testthat::expect_type(p, "list")
  })

}
