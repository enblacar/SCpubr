sample <- use_dataset()
sample$orig.ident <- ifelse(sample$seurat_clusters %in% c("1", "2"), "A", "B")

genes <- list("A" = Seurat::VariableFeatures(sample)[1:5],
              "B" = Seurat::VariableFeatures(sample)[6:10],
              "C" = Seurat::VariableFeatures(sample)[11:15],
              "D" = Seurat::VariableFeatures(sample)[16:20])

testthat::test_that("do_CellularStatesPlot: PASS - 2 variables", {
  p <- SCpubr:::do_CellularStatesPlot(sample = sample,
                                     gene_list = genes,
                                     x1 = "A",
                                     y1 = "B")
  testthat::expect_type(p, "list")
})


testthat::test_that("do_CellularStatesPlot: PASS - 3 variables", {
  p <- SCpubr:::do_CellularStatesPlot(sample = sample,
                                      gene_list = genes,
                                      x1 = "A",
                                      y1 = "B",
                                      x2 = "C")
  testthat::expect_type(p, "list")
})

testthat::test_that("do_CellularStatesPlot: PASS - 4 variables", {
  p <- SCpubr:::do_CellularStatesPlot(sample = sample,
                                      gene_list = genes,
                                      x1 = "A",
                                      y1 = "B",
                                      x2 = "C",
                                      y2 = "D")
  testthat::expect_type(p, "list")
})
