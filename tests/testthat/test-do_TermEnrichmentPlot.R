
testthat::test_that("do_TermEnrichmentPlot: PASS - A", {
  genes <- c("ABCB1", "ABCG2", "AHR", "AKT1", "AR")

  p <- SCpubr::do_TermEnrichmentPlot(genes = genes,
                                     dbs_use = "A",
                                     nterms = 2)
  testthat::expect_type(p, "list")
  testthat::expect_equal(length(names(p)), 8)
})

testthat::test_that("do_TermEnrichmentPlot: PASS - B", {
  genes <- c("ABCB1", "ABCG2", "AHR", "AKT1", "AR")

  p <- SCpubr::do_TermEnrichmentPlot(genes = genes,
                                     dbs_use = "B",
                                     nterms = 2)
  testthat::expect_type(p, "list")
  testthat::expect_equal(length(names(p)), 4)
})


testthat::test_that("do_TermEnrichmentPlot: PASS - C", {
  genes <- c("ABCB1", "ABCG2", "AHR", "AKT1", "AR")

  p <- SCpubr::do_TermEnrichmentPlot(genes = genes,
                                     dbs_use = "C",
                                     nterms = 2)
  testthat::expect_type(p, "list")
  testthat::expect_equal(length(names(p)), 4)
})

testthat::test_that("do_TermEnrichmentPlot: PASS - legend position = right", {
  genes <- c("ABCB1", "ABCG2", "AHR", "AKT1", "AR")

  p <- SCpubr::do_TermEnrichmentPlot(genes = genes,
                                     nterms = 2,
                                     dbs_use = "GO_Biological_Process_2021",
                                     legend.position = "right")
  testthat::expect_type(p, "list")
})

testthat::test_that("do_TermEnrichmentPlot: PASS - number of terms", {
  genes <- c("ABCB1", "ABCG2", "AHR", "AKT1", "AR")

  p <- SCpubr::do_TermEnrichmentPlot(genes = genes,
                                     dbs_use = "C",
                                     nterms = 2)
  testthat::expect_type(p, "list")
})


testthat::test_that("do_TermEnrichmentPlot: PASS - length of terms", {
  genes <- c("ABCB1", "ABCG2", "AHR", "AKT1", "AR")

  p <- SCpubr::do_TermEnrichmentPlot(genes = genes,
                                     dbs_use = "C",
                                     nterms = 2,
                                     nchar_wrap = 20)
  testthat::expect_type(p, "list")
})

testthat::test_that("do_TermEnrichmentPlot: PASS - modify colors", {
  genes <- c("ABCB1", "ABCG2", "AHR", "AKT1", "AR")

  p <- SCpubr::do_TermEnrichmentPlot(genes = genes,
                                     dbs_use = "C",
                                     nterms = 2,
                                     colors.use = c("#e9d8a6", "#9b2226"))
  testthat::expect_type(p, "list")
})

testthat::test_that("do_TermEnrichmentPlot: PASS - single database", {
  genes <- c("ABCB1", "ABCG2", "AHR", "AKT1", "AR")

  p <- SCpubr::do_TermEnrichmentPlot(genes = genes,
                                     dbs_use = "GO_Biological_Process_2021",
                                     nterms = 2,
                                     colors.use = c("#e9d8a6", "#9b2226"))
  testthat::expect_type(p, "list")
})

testthat::test_that("do_TermEnrichmentPlot: PASS - single database legend normal", {
  genes <- c("ABCB1", "ABCG2", "AHR", "AKT1", "AR")

  p <- SCpubr::do_TermEnrichmentPlot(genes = genes,
                                     dbs_use = "GO_Biological_Process_2021",
                                     nterms = 2,
                                     colors.use = c("#e9d8a6", "#9b2226"),
                                     legend.type = "normal")
  testthat::expect_type(p, "list")
})

testthat::test_that("do_TermEnrichmentPlot: PASS - single database legend colorbar", {
  genes <- c("ABCB1", "ABCG2", "AHR", "AKT1", "AR")

  p <- SCpubr::do_TermEnrichmentPlot(genes = genes,
                                     dbs_use = "GO_Biological_Process_2021",
                                     nterms = 2,
                                     colors.use = c("#e9d8a6", "#9b2226"),
                                     legend.type = "colorbar")
  testthat::expect_type(p, "list")
})


testthat::test_that("do_TermEnrichmentPlot: PASS - single database legend colorsteps", {
  genes <- c("ABCB1", "ABCG2", "AHR", "AKT1", "AR")

  p <- SCpubr::do_TermEnrichmentPlot(genes = genes,
                                     dbs_use = "GO_Biological_Process_2021",
                                     nterms = 2,
                                     colors.use = c("#e9d8a6", "#9b2226"),
                                     legend.type = "colorsteps")
  testthat::expect_type(p, "list")
})

testthat::test_that("do_TermEnrichmentPlot: FAIL - wrong legend.type", {
  genes <- c("ABCB1", "ABCG2", "AHR", "AKT1", "AR")

  testthat::expect_error(SCpubr::do_TermEnrichmentPlot(genes = genes,
                                                       dbs_use = "GO_Biological_Process_2021",
                                                       nterms = 2,
                                                       legend.type = "wrong"))
})

testthat::test_that("do_TermEnrichmentPlot: FAIL - wrong legend.position", {
  genes <- c("ABCB1", "ABCG2", "AHR", "AKT1", "AR")

  testthat::expect_error(SCpubr::do_TermEnrichmentPlot(genes = genes,
                                                       dbs_use = "GO_Biological_Process_2021",
                                                       nterms = 2,
                                                       legend.position = "wrong"))
})

testthat::test_that("do_TermEnrichmentPlot: FAIL - wrong font.type", {
  genes <- c("ABCB1", "ABCG2", "AHR", "AKT1", "AR")

  testthat::expect_error(SCpubr::do_TermEnrichmentPlot(genes = genes,
                                                       dbs_use = "GO_Biological_Process_2021",
                                                       nterms = 2,
                                                       font.type = "wrong"))
})

testthat::test_that("do_TermEnrichmentPlot: PASS - single database legend normal - one pvalue", {
  genes <- c("ABCB1", "ABCG2", "AHR", "AKT1", "AR")

  p <- SCpubr::do_TermEnrichmentPlot(genes = genes,
                                     dbs_use = "Azimuth_Cell_Types_2021",
                                     nterms = 2,
                                     colors.use = c("#e9d8a6", "#9b2226"),
                                     legend.type = "normal")
  testthat::expect_type(p, "list")
})

testthat::test_that("do_TermEnrichmentPlot: PASS - single database legend colorbar - one pvalue", {
  genes <- c("ABCB1", "ABCG2", "AHR", "AKT1", "AR")

  p <- SCpubr::do_TermEnrichmentPlot(genes = genes,
                                     dbs_use = "Azimuth_Cell_Types_2021",
                                     nterms = 2,
                                     colors.use = c("#e9d8a6", "#9b2226"),
                                     legend.type = "colorbar")
  testthat::expect_type(p, "list")
})


testthat::test_that("do_TermEnrichmentPlot: PASS - single database legend colorsteps - one pvalue", {
  genes <- c("ABCB1", "ABCG2", "AHR", "AKT1", "AR")

  p <- SCpubr::do_TermEnrichmentPlot(genes = genes,
                                     dbs_use = "Azimuth_Cell_Types_2021",
                                     nterms = 2,
                                     colors.use = c("#e9d8a6", "#9b2226"),
                                     legend.type = "colorsteps")
  testthat::expect_type(p, "list")
})

testthat::test_that("do_TermEnrichmentPlot: FAIL - provide more colors than needed", {
  genes <- c("ABCB1", "ABCG2", "AHR", "AKT1", "AR")

  testthat::expect_error(SCpubr::do_TermEnrichmentPlot(genes = genes,
                                                       dbs_use = "C",
                                                       nterms = 2,
                                                       colors.use = c("#e9d8a6", "#9b2226", "red")))
})

testthat::test_that("do_TermEnrichmentPlot: FAIL - database = NULL", {
  genes <- c("ABCB1", "ABCG2", "AHR", "AKT1", "AR")

  testthat::expect_error(SCpubr::do_TermEnrichmentPlot(genes = genes,
                                                       nterms = 2,
                                                       colors.use = c("#e9d8a6", "#9b2226")))
})

testthat::test_that("do_TermEnrichmentPlot: FAIL - wrong database", {
  genes <- c("ABCB1", "ABCG2", "AHR", "AKT1", "AR")

  testthat::expect_error(SCpubr::do_TermEnrichmentPlot(genes = genes,
                                                       dbs_use = "wrong_database",
                                                       nterms = 2,
                                                       colors.use = c("#e9d8a6", "#9b2226")))
})

testthat::test_that("do_TermEnrichmentPlot: FAIL - more than one database", {
  genes <- c("ABCB1", "ABCG2", "AHR", "AKT1", "AR")

  testthat::expect_error(SCpubr::do_TermEnrichmentPlot(genes = genes,
                                                       dbs_use = c("A", "B"),
                                                       nterms = 2,
                                                       colors.use = c("#e9d8a6", "#9b2226")))
})
