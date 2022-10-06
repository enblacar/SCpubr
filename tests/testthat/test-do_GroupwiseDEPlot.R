sample <- SCpubr:::test_list$sample

de_genes <- SCpubr:::test_list$de_genes
de_genes_scaled <- SCpubr:::test_list$de_genes_scaled

sample@assays$SCT@scale.data <- as.matrix(sample@assays$SCT@data)

testthat::test_that("do_GroupwiseDEPlot: PASS - default", {
  p <- SCpubr::do_GroupwiseDEPlot(sample = sample,
                                  de_genes = de_genes,
                                  assay = "SCT",
                                  slot = "data")
  testthat::expect_type(p, "S4")

  p <- SCpubr::do_GroupwiseDEPlot(sample = sample,
                                  de_genes = de_genes_scaled,
                                  assay = "SCT",
                                  slot = "scale.data")
  testthat::expect_type(p, "S4")
})

testthat::test_that("do_GroupwiseDEPlot: PASS - heatmap legend side", {
  p <- SCpubr::do_GroupwiseDEPlot(sample = sample,
                                  de_genes = de_genes,
                                  assay = "SCT",
                                  slot = "data",
                                  legend.position = "right")
  testthat::expect_type(p, "S4")

  p <- SCpubr::do_GroupwiseDEPlot(sample = sample,
                                  de_genes = de_genes_scaled,
                                  assay = "SCT",
                                  slot = "scale.data",
                                  legend.position = "right")
  testthat::expect_type(p, "S4")
})


testthat::test_that("do_GroupwiseDEPlot: PASS - multiple grouping", {
  p <- SCpubr::do_GroupwiseDEPlot(sample = sample,
                                  de_genes = de_genes,
                                  assay = "SCT",
                                  slot = "data",
                                  group.by = c("seurat_clusters", "orig.ident"),
                                  row_title_expression = c("", ""))
  testthat::expect_type(p, "S4")

  p <- SCpubr::do_GroupwiseDEPlot(sample = sample,
                                  de_genes = de_genes_scaled,
                                  assay = "SCT",
                                  slot = "scale.data",
                                  group.by = c("seurat_clusters", "orig.ident"),
                                  row_title_expression = c("", ""))
  testthat::expect_type(p, "S4")
})

testthat::test_that("do_GroupwiseDEPlot: FAIL - wrong number of titles", {
  testthat::expect_error({SCpubr::do_GroupwiseDEPlot(sample = sample,
                                                     de_genes = de_genes,
                                                     assay = "SCT",
                                                     slot = "data",
                                                     group.by = c("seurat_clusters", "orig.ident"))})
  testthat::expect_error({SCpubr::do_GroupwiseDEPlot(sample = sample,
                                                     de_genes = de_genes_scaled,
                                                     assay = "SCT",
                                                     slot = "scale.data",
                                                     group.by = c("seurat_clusters", "orig.ident"))})
})

testthat::test_that("do_GroupwiseDEPlot: FAIL - wrong direction", {
  testthat::expect_error({SCpubr::do_GroupwiseDEPlot(sample = sample,
                                                     de_genes = de_genes,
                                                     assay = "SCT",
                                                     slot = "data",
                                                     scale_direction = 0)})
  testthat::expect_error({SCpubr::do_GroupwiseDEPlot(sample = sample,
                                                     de_genes = de_genes_scaled,
                                                     assay = "SCT",
                                                     slot = "scale.data",
                                                     scale_direction = 0)})
})
