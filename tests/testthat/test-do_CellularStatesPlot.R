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

testthat::test_that("do_CellularStatesPlot: PASS - 2 variables, colors.use", {
  Seurat::Idents(sample) <- sample$orig.ident
  p <- SCpubr:::do_CellularStatesPlot(sample = sample,
                                      gene_list = genes,
                                      colors.use = c("A" = "black", "B" = "red"),
                                      x1 = "A",
                                      y1 = "B")
  testthat::expect_type(p, "list")
})

testthat::test_that("do_CellularStatesPlot: PASS - 2 variables, group.by", {
  p <- SCpubr:::do_CellularStatesPlot(sample = sample,
                                      gene_list = genes,
                                      group.by = "orig.ident",
                                      x1 = "A",
                                      y1 = "B")
  testthat::expect_type(p, "list")
})

testthat::test_that("do_CellularStatesPlot: PASS - 2 variables remove axis ticks", {
  p <- SCpubr:::do_CellularStatesPlot(sample = sample,
                                      gene_list = genes,
                                      x1 = "A",
                                      y1 = "B",
                                      axis.ticks = F)
  testthat::expect_type(p, "list")
})

testthat::test_that("do_CellularStatesPlot: PASS - 2 variables remove axis text", {
  p <- SCpubr:::do_CellularStatesPlot(sample = sample,
                                      gene_list = genes,
                                      x1 = "A",
                                      y1 = "B",
                                      axis.text = F)
  testthat::expect_type(p, "list")
})

testthat::test_that("do_CellularStatesPlot: PASS - 2 variables, group.by and colors.use", {
  p <- SCpubr:::do_CellularStatesPlot(sample = sample,
                                      gene_list = genes,
                                      group.by = "orig.ident",
                                      colors.use = c("A" = "black", "B" = "red"),
                                      x1 = "A",
                                      y1 = "B")
  testthat::expect_type(p, "list")
})

testthat::test_that("do_CellularStatesPlot: FAIL - 2 variables same parameters", {
  testthat::expect_error(SCpubr:::do_CellularStatesPlot(sample = sample,
                                                        gene_list = genes,
                                                        x1 = "A",
                                                        y1 = "A"))
})

testthat::test_that("do_CellularStatesPlot: FAIL - 2 variables x1 not in the list", {
  testthat::expect_error(SCpubr:::do_CellularStatesPlot(sample = sample,
                                                        gene_list = genes,
                                                        x1 = "Not in list",
                                                        y1 = "A"))
})

testthat::test_that("do_CellularStatesPlot: FAIL - 2 variables y1 not in the list", {
  testthat::expect_error(SCpubr:::do_CellularStatesPlot(sample = sample,
                                                        gene_list = genes,
                                                        x1 = "A",
                                                        y1 = "Not in list"))
})

testthat::test_that("do_CellularStatesPlot: PASS - 3 variables", {
  p <- SCpubr:::do_CellularStatesPlot(sample = sample,
                                      gene_list = genes,
                                      x1 = "A",
                                      y1 = "B",
                                      x2 = "C")
  testthat::expect_type(p, "list")
})

testthat::test_that("do_CellularStatesPlot: PASS - 3 variables, verbose", {
  p <- SCpubr:::do_CellularStatesPlot(sample = sample,
                                      gene_list = genes,
                                      x1 = "A",
                                      y1 = "B",
                                      x2 = "C",
                                      verbose = T)
  testthat::expect_type(p, "list")
})

testthat::test_that("do_CellularStatesPlot: FAIL - 3 variables duplicated parameters", {
  testthat::expect_error(SCpubr:::do_CellularStatesPlot(sample = sample,
                                                        gene_list = genes,
                                                        x1 = "A",
                                                        y1 = "A",
                                                        x2 = "B"))
})

testthat::test_that("do_CellularStatesPlot: FAIL - 3 variables x1 not in list", {
  testthat::expect_error(SCpubr:::do_CellularStatesPlot(sample = sample,
                                                        gene_list = genes,
                                                        x1 = "Not in list",
                                                        y1 = "A",
                                                        x2 = "B"))
})

testthat::test_that("do_CellularStatesPlot: FAIL - 3 variables x2 not in list", {
  testthat::expect_error(SCpubr:::do_CellularStatesPlot(sample = sample,
                                                        gene_list = genes,
                                                        x1 = "B",
                                                        y1 = "A",
                                                        x2 = "Not in list"))
})

testthat::test_that("do_CellularStatesPlot: FAIL - 3 variables y1 not in list", {
  testthat::expect_error(SCpubr:::do_CellularStatesPlot(sample = sample,
                                                        gene_list = genes,
                                                        x1 = "A",
                                                        y1 = "Not in list",
                                                        x2 = "B"))
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

testthat::test_that("do_CellularStatesPlot: PASS - 4 variables, verbose", {
  p <- SCpubr:::do_CellularStatesPlot(sample = sample,
                                      gene_list = genes,
                                      x1 = "A",
                                      y1 = "B",
                                      x2 = "C",
                                      y2 = "D",
                                      verbose = T)
  testthat::expect_type(p, "list")
})

testthat::test_that("do_CellularStatesPlot: FAIL - 4 variables repeated parameters", {
  testthat::expect_error(SCpubr:::do_CellularStatesPlot(sample = sample,
                                                        gene_list = genes,
                                                        x1 = "A",
                                                        y1 = "A",
                                                        x2 = "A",
                                                        y2 = "A"))
})

testthat::test_that("do_CellularStatesPlot: FAIL - 4 variables x1 not in list", {
  testthat::expect_error(SCpubr:::do_CellularStatesPlot(sample = sample,
                                                        gene_list = genes,
                                                        x1 = "Not in list",
                                                        y1 = "B",
                                                        x2 = "C",
                                                        y2 = "D"))
})

testthat::test_that("do_CellularStatesPlot: FAIL - 4 variables y1 not in list", {
  testthat::expect_error(SCpubr:::do_CellularStatesPlot(sample = sample,
                                                        gene_list = genes,
                                                        x1 = "A",
                                                        y1 = "Not in list",
                                                        x2 = "C",
                                                        y2 = "D"))
})

testthat::test_that("do_CellularStatesPlot: FAIL - 4 variables x2 not in list", {
  testthat::expect_error(SCpubr:::do_CellularStatesPlot(sample = sample,
                                                        gene_list = genes,
                                                        x1 = "A",
                                                        y1 = "B",
                                                        x2 = "Not in list",
                                                        y2 = "D"))
})

testthat::test_that("do_CellularStatesPlot: FAIL - 4 variables y2 not in list", {
  testthat::expect_error(SCpubr:::do_CellularStatesPlot(sample = sample,
                                                        gene_list = genes,
                                                        x1 = "A",
                                                        y1 = "B",
                                                        x2 = "C",
                                                        y2 = "Not in list"))
})
