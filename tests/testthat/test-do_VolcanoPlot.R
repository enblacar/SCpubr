if(isFALSE(dep_check[["do_VolcanoPlot"]])){
  testthat::test_that("do_VolcanoPlot: PASS - default", {

    `%>%` <- magrittr::`%>%`

    p <- SCpubr::do_VolcanoPlot(sample = sample,
                                de_genes = de_genes)
    testthat::expect_type(p, "list")

    de_genes[1, "p_val_adj"] <- 1
    de_genes[2, "avg_log2FC"] <- 0.001
    de_genes[3, "avg_log2FC"] <- 3
    de_genes[3, "p_val_adj"] <- 0.003

    p <- SCpubr::do_VolcanoPlot(sample = sample,
                                de_genes = de_genes)
    testthat::expect_type(p, "list")

    de_genes <- de_genes %>%
      tibble::as_tibble() %>%
      dplyr::distinct(.data$gene, .keep_all = TRUE) %>%
      tibble::column_to_rownames(var = "gene")

    p <- SCpubr::do_VolcanoPlot(sample = sample,
                                de_genes = de_genes)
    testthat::expect_type(p, "list")
  })

  testthat::test_that("do_VolcanoPlot: PASS - n_genes", {

    `%>%` <- magrittr::`%>%`

    p <- SCpubr::do_VolcanoPlot(sample = sample,
                                de_genes = de_genes,
                                n_genes = 15)
    testthat::expect_type(p, "list")
  })

  testthat::test_that("do_VolcanoPlot: PASS - use_labels", {

    `%>%` <- magrittr::`%>%`

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

    `%>%` <- magrittr::`%>%`

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

    `%>%` <- magrittr::`%>%`

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

    `%>%` <- magrittr::`%>%`

    testthat::expect_error(SCpubr::do_VolcanoPlot(sample = sample,
                                                  de_genes = de_genes,
                                                  add_gene_tags = TRUE,
                                                  order_tags_by = "wrong"))
  })
}



