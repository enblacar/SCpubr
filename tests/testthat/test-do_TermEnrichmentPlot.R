if(base::isFALSE(dep_check[["do_TermEnrichmentPlot"]])){

  testthat::test_that("do_TermEnrichmentPlot: CRAN essentials", {

    p <- SCpubr::do_TermEnrichmentPlot(mat = enriched_terms,
                                       n.terms = 2)
    testthat::expect_true(ggplot2::is_ggplot(p))
  })


  testthat::test_that("do_TermEnrichmentPlot: PASS - legend position = right", {
    testthat::skip_on_cran()
    p <- SCpubr::do_TermEnrichmentPlot(mat = enriched_terms,
                                       legend.position = "right")
    testthat::expect_true(ggplot2::is_ggplot(p))
  })

  testthat::test_that("do_TermEnrichmentPlot: PASS - number of terms", {
    testthat::skip_on_cran()
    p <- SCpubr::do_TermEnrichmentPlot(mat = enriched_terms,
                                       n.terms = 20)
    testthat::expect_true(ggplot2::is_ggplot(p))
  })


  testthat::test_that("do_TermEnrichmentPlot: PASS - length of terms", {
    testthat::skip_on_cran()
    p <- SCpubr::do_TermEnrichmentPlot(mat = enriched_terms,
                                       n.terms = 2,
                                       n.chars = 20)
    testthat::expect_true(ggplot2::is_ggplot(p))
  })

  testthat::test_that("do_TermEnrichmentPlot: PASS - modify colors", {
    testthat::skip_on_cran()
    p <- SCpubr::do_TermEnrichmentPlot(mat = enriched_terms,
                                       n.terms = 2,
                                       sequential.palette = "YlOrRd")
    testthat::expect_true(ggplot2::is_ggplot(p))
  })


  testthat::test_that("do_TermEnrichmentPlot: PASS - single database legend types", {
    testthat::skip_on_cran()
    p <- SCpubr::do_TermEnrichmentPlot(mat = enriched_terms,
                                       n.terms = 2,
                                       legend.type = "normal")
    testthat::expect_true(ggplot2::is_ggplot(p))

    p <- SCpubr::do_TermEnrichmentPlot(mat = enriched_terms,
                                       n.terms = 2,
                                       legend.type = "colorbar")
    testthat::expect_true(ggplot2::is_ggplot(p))
   })


  testthat::test_that("do_TermEnrichmentPlot: FAIL - wrong legend.type", {
    testthat::skip_on_cran()
    testthat::expect_error(SCpubr::do_TermEnrichmentPlot(enriched_terms = enriched_terms,
                                                         n.terms = 2,
                                                         legend.type = "wrong"))
  })

  testthat::test_that("do_TermEnrichmentPlot: FAIL - wrong legend.position", {
    testthat::skip_on_cran()
    testthat::expect_error(SCpubr::do_TermEnrichmentPlot(enriched_terms = enriched_terms,
                                                         n.terms = 2,
                                                         legend.position = "wrong"))
  })

  testthat::test_that("do_TermEnrichmentPlot: FAIL - wrong font.type", {
    testthat::skip_on_cran()
    testthat::expect_error(SCpubr::do_TermEnrichmentPlot(enriched_terms = enriched_terms,
                                                         n.terms = 2,
                                                         font.type = "wrong"))
  })

}

