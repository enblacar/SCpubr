utils::data("test.data", package = "SCpubr")

sample <- test.data$monocle_sample
cds <- test.data$monocle_cds

testthat::test_that("do_PseudotimePlot: PASS - default", {

  pseudotime_genes <- Seurat::VariableFeatures(sample)[1:5]

 out <- SCpubr::do_PseudotimePlot(sample = sample,
                                  cds = cds,
                                  pseudotime_genes = pseudotime_genes,
                                  nbin = 10,
                                  ctrl = 10)
  testthat::expect_type(out, "list")
  testthat::expect_length(out, 4)

  out <- SCpubr::do_PseudotimePlot(sample = sample,
                                   cds = cds,
                                   pseudotime_genes = pseudotime_genes,
                                   nbin = 10,
                                   ctrl = 10,
                                   is_max_score_the_start = TRUE)
  testthat::expect_type(out, "list")
  testthat::expect_length(out, 4)

  out <- SCpubr::do_PseudotimePlot(sample = sample,
                                   cds = cds,
                                   pseudotime_genes = pseudotime_genes,
                                   nbin = 10,
                                   ctrl = 10,
                                   is_max_score_the_start = FALSE)
  testthat::expect_type(out, "list")
  testthat::expect_length(out, 4)

  out <- SCpubr::do_PseudotimePlot(sample = sample,
                                   cds = cds,
                                   pseudotime_genes = pseudotime_genes,
                                   nbin = 10,
                                   ctrl = 10,
                                   group.by = "orig.ident")
  testthat::expect_type(out, "list")
  testthat::expect_length(out, 4)

  out <- SCpubr::do_PseudotimePlot(sample = sample,
                                   cds = cds,
                                   pseudotime_genes = pseudotime_genes,
                                   nbin = 10,
                                   ctrl = 10,
                                   compute_monocle_partitions = TRUE,
                                   compute_monocle_clusters = TRUE)
  testthat::expect_type(out, "list")
  testthat::expect_length(out, 4)

  out <- SCpubr::do_PseudotimePlot(sample = sample,
                                   cds = cds,
                                   pseudotime_genes = pseudotime_genes,
                                   nbin = 10,
                                   ctrl = 10,
                                   compute_monocle_partitions = TRUE,
                                   compute_monocle_clusters = FALSE)
  testthat::expect_type(out, "list")
  testthat::expect_length(out, 4)

  out <- SCpubr::do_PseudotimePlot(sample = sample,
                                   cds = cds,
                                   pseudotime_genes = pseudotime_genes,
                                   nbin = 10,
                                   ctrl = 10,
                                   compute_monocle_partitions = FALSE,
                                   compute_monocle_clusters = TRUE)
  testthat::expect_type(out, "list")
  testthat::expect_length(out, 4)

  out <- SCpubr::do_PseudotimePlot(sample = sample,
                                   cds = cds,
                                   pseudotime_genes = pseudotime_genes,
                                   nbin = 10,
                                   ctrl = 10,
                                   compute_monocle_partitions = FALSE,
                                   compute_monocle_clusters = FALSE)
  testthat::expect_type(out, "list")
  testthat::expect_length(out, 4)

  testthat::expect_error(SCpubr::do_PseudotimePlot(sample = sample,
                                                   cds = cds,
                                                   pseudotime_genes = pseudotime_genes,
                                                   nbin = 10,
                                                   ctrl = 10,
                                                   group.by = "wrong"))
})
