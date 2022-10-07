sample <- SCpubr:::test_list$sample
de_genes <- SCpubr:::test_list$de_genes
`%>%` <- magrittr::`%>%`

testthat::test_that("do_VolcanoPlot: PASS - default", {
  p <- SCpubr::do_VolcanoPlot(sample = sample,
                              de_genes = de_genes)
  testthat::expect_type(p, "list")

  de_genes[1, "p_val_adj"] <- 1
  de_genes[2, "avg_log2FC"] <- 0.001

  p <- SCpubr::do_VolcanoPlot(sample = sample,
                              de_genes = de_genes)
  testthat::expect_type(p, "list")

  de_genes <- de_genes %>%
              dplyr::distinct(.data$gene, .keep_all = TRUE) %>%
              tibble::column_to_rownames(var = "gene")

  p <- SCpubr::do_VolcanoPlot(sample = sample,
                              de_genes = de_genes)
  testthat::expect_type(p, "list")
})

testthat::test_that("do_VolcanoPlot: PASS - n_genes", {
  p <- SCpubr::do_VolcanoPlot(sample = sample,
                              de_genes = de_genes,
                              n_genes = 15)
  testthat::expect_type(p, "list")
})

testthat::test_that("do_VolcanoPlot: PASS - use_labels", {
  p <- SCpubr::do_VolcanoPlot(sample = sample,
                              de_genes = de_genes,
                              use_labels = TRUE)
  testthat::expect_type(p, "list")

  p <- SCpubr::do_VolcanoPlot(sample = sample,
                              de_genes = de_genes,
                              use_labels = FALSE)
  testthat::expect_type(p, "list")
})

testthat::test_that("do_VolcanoPlot: PASS - gene tags", {
  p <- SCpubr::do_VolcanoPlot(sample = sample,
                              de_genes = de_genes,
                              add_gene_tags = TRUE)
  testthat::expect_type(p, "list")

  p <- SCpubr::do_VolcanoPlot(sample = sample,
                              de_genes = de_genes,
                              add_gene_tags = FALSE)
  testthat::expect_type(p, "list")
})

testthat::test_that("do_VolcanoPlot: PASS - gene tags order by", {
  p <- SCpubr::do_VolcanoPlot(sample = sample,
                              de_genes = de_genes,
                              add_gene_tags = TRUE,
                              order_tags_by = "both")
  testthat::expect_type(p, "list")

  p <- SCpubr::do_VolcanoPlot(sample = sample,
                              de_genes = de_genes,
                              add_gene_tags = TRUE,
                              order_tags_by = "p_value")
  testthat::expect_type(p, "list")

  p <- SCpubr::do_VolcanoPlot(sample = sample,
                              de_genes = de_genes,
                              add_gene_tags = TRUE,
                              order_tags_by = "logfc")
  testthat::expect_type(p, "list")
})

testthat::test_that("do_VolcanoPlot: FAIL - wrong parameters", {
  testthat::expect_error(SCpubr::do_VolcanoPlot(sample = sample,
                                                de_genes = de_genes,
                                                add_gene_tags = TRUE,
                                                order_tags_by = "wrong"))
})


