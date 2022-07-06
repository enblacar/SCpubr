
genes <- c("ABCB1", "ABCG2", "AHR", "AKT1", "AR")


testthat::test_that("do_TermEnrichmentPlot: PASS - A", {
  p <- SCpubr::do_TermEnrichmentPlot(genes = genes,
                                     dbs_use = "A",
                                     nterms = 2)
  testthat::expect_type(p, "list")
  testthat::expect_equal(length(names(p)), 8)
})

testthat::test_that("do_TermEnrichmentPlot: PASS - B", {
  p <- SCpubr::do_TermEnrichmentPlot(genes = genes,
                                     dbs_use = "B",
                                     nterms = 2)
  testthat::expect_type(p, "list")
  testthat::expect_equal(length(names(p)), 4)
})


testthat::test_that("do_TermEnrichmentPlot: PASS - C", {
  p <- SCpubr::do_TermEnrichmentPlot(genes = genes,
                                     dbs_use = "C",
                                     nterms = 2)
  testthat::expect_type(p, "list")
  testthat::expect_equal(length(names(p)), 4)
})

testthat::test_that("do_TermEnrichmentPlot: PASS - legend position = right", {
  p <- SCpubr::do_TermEnrichmentPlot(genes = genes,
                                     nterms = 2,
                                     dbs_use = "GO_Biological_Process_2021",
                                     legend.position = "right")
  testthat::expect_type(p, "list")
})

testthat::test_that("do_TermEnrichmentPlot: PASS - number of terms", {
  p <- SCpubr::do_TermEnrichmentPlot(genes = genes,
                                     dbs_use = "C",
                                     nterms = 2)
  testthat::expect_type(p, "list")
})


testthat::test_that("do_TermEnrichmentPlot: PASS - length of terms", {
  p <- SCpubr::do_TermEnrichmentPlot(genes = genes,
                                     dbs_use = "C",
                                     nterms = 2,
                                     nchar_wrap = 20)
  testthat::expect_type(p, "list")
})

testthat::test_that("do_TermEnrichmentPlot: PASS - modify colors", {
  p <- SCpubr::do_TermEnrichmentPlot(genes = genes,
                                     dbs_use = "C",
                                     nterms = 2,
                                     colors.use = c("#e9d8a6", "#9b2226"))
  testthat::expect_type(p, "list")
})

testthat::test_that("do_TermEnrichmentPlot: PASS - single database", {
  p <- SCpubr::do_TermEnrichmentPlot(genes = genes,
                                     dbs_use = "GO_Biological_Process_2021",
                                     nterms = 2,
                                     colors.use = c("#e9d8a6", "#9b2226"))
  testthat::expect_type(p, "list")
})

testthat::test_that("do_TermEnrichmentPlot: PASS - single database legend normal", {
  p <- SCpubr::do_TermEnrichmentPlot(genes = genes,
                                     dbs_use = "GO_Biological_Process_2021",
                                     nterms = 2,
                                     colors.use = c("#e9d8a6", "#9b2226"),
                                     legend.type = "normal")
  testthat::expect_type(p, "list")
})

testthat::test_that("do_TermEnrichmentPlot: PASS - single database legend colorbar", {
  p <- SCpubr::do_TermEnrichmentPlot(genes = genes,
                                     dbs_use = "GO_Biological_Process_2021",
                                     nterms = 2,
                                     colors.use = c("#e9d8a6", "#9b2226"),
                                     legend.type = "colorbar")
  testthat::expect_type(p, "list")
})


testthat::test_that("do_TermEnrichmentPlot: PASS - single database legend colorsteps", {
  p <- SCpubr::do_TermEnrichmentPlot(genes = genes,
                                     dbs_use = "GO_Biological_Process_2021",
                                     nterms = 2,
                                     colors.use = c("#e9d8a6", "#9b2226"),
                                     legend.type = "colorsteps")
  testthat::expect_type(p, "list")
})

testthat::test_that("do_TermEnrichmentPlot: FAIL - wrong legend.type", {
  testthat::expect_error(SCpubr::do_TermEnrichmentPlot(genes = genes,
                                                       dbs_use = "GO_Biological_Process_2021",
                                                       nterms = 2,
                                                       legend.type = "wrong"))
})

testthat::test_that("do_TermEnrichmentPlot: FAIL - wrong legend.position", {
  testthat::expect_error(SCpubr::do_TermEnrichmentPlot(genes = genes,
                                                       dbs_use = "GO_Biological_Process_2021",
                                                       nterms = 2,
                                                       legend.position = "wrong"))
})

testthat::test_that("do_TermEnrichmentPlot: PASS - single database legend normal - one pvalue", {
  p <- SCpubr::do_TermEnrichmentPlot(genes = genes,
                                     dbs_use = "Azimuth_Cell_Types_2021",
                                     nterms = 2,
                                     colors.use = c("#e9d8a6", "#9b2226"),
                                     legend.type = "normal")
  testthat::expect_type(p, "list")
})

testthat::test_that("do_TermEnrichmentPlot: PASS - single database legend colorbar - one pvalue", {
  p <- SCpubr::do_TermEnrichmentPlot(genes = genes,
                                     dbs_use = "Azimuth_Cell_Types_2021",
                                     nterms = 2,
                                     colors.use = c("#e9d8a6", "#9b2226"),
                                     legend.type = "colorbar")
  testthat::expect_type(p, "list")
})


testthat::test_that("do_TermEnrichmentPlot: PASS - single database legend colorsteps - one pvalue", {
  p <- SCpubr::do_TermEnrichmentPlot(genes = genes,
                                     dbs_use = "Azimuth_Cell_Types_2021",
                                     nterms = 2,
                                     colors.use = c("#e9d8a6", "#9b2226"),
                                     legend.type = "colorsteps")
  testthat::expect_type(p, "list")
})

testthat::test_that("do_TermEnrichmentPlot: FAIL - provide more colors than needed", {
  testthat::expect_error(SCpubr::do_TermEnrichmentPlot(genes = genes,
                                                       dbs_use = "C",
                                                       nterms = 2,
                                                       colors.use = c("#e9d8a6", "#9b2226", "red")))
})

testthat::test_that("do_TermEnrichmentPlot: FAIL - database = NULL", {
  testthat::expect_error(SCpubr::do_TermEnrichmentPlot(genes = genes,
                                                       nterms = 2,
                                                       colors.use = c("#e9d8a6", "#9b2226")))
})

testthat::test_that("do_TermEnrichmentPlot: FAIL - wrong database", {
  testthat::expect_error(SCpubr::do_TermEnrichmentPlot(genes = genes,
                                                       dbs_use = "wrong_database",
                                                       nterms = 2,
                                                       colors.use = c("#e9d8a6", "#9b2226")))
})

testthat::test_that("do_TermEnrichmentPlot: FAIL - more than one database", {
  testthat::expect_error(SCpubr::do_TermEnrichmentPlot(genes = genes,
                                                       dbs_use = c("A", "B"),
                                                       nterms = 2,
                                                       colors.use = c("#e9d8a6", "#9b2226")))
})
