
testthat::test_that("do_EnrichmentHeatmap: PASS - normal", {

  sample$orig.ident <- ifelse(sample$seurat_clusters %in% c("1", "2"), "A", "B")

  genes <- list("A" = Seurat::VariableFeatures(sample)[1:5],
                "B" = Seurat::VariableFeatures(sample)[6:10],
                "C" = Seurat::VariableFeatures(sample)[11:15])


  p <- SCpubr::do_EnrichmentHeatmap(sample = sample,
                                    input_gene_list = genes,
                                    nbin = 1,
                                    ctrl = 10)
  testthat::expect_true("HeatmapList" %in% class(p))

  p <- SCpubr::do_EnrichmentHeatmap(sample = sample,
                                    input_gene_list = genes,
                                    legend.position = "top",
                                    nbin = 1,
                                    ctrl = 10)
  testthat::expect_true("HeatmapList" %in% class(p))

  p <- SCpubr::do_EnrichmentHeatmap(sample = sample,
                                    input_gene_list = genes,
                                    legend.position = "right",
                                    nbin = 1,
                                    ctrl = 10)
  testthat::expect_true("HeatmapList" %in% class(p))
})

testthat::test_that("do_EnrichmentHeatmap: PASS - group.by", {

  sample$orig.ident <- ifelse(sample$seurat_clusters %in% c("1", "2"), "A", "B")

  genes <- list("A" = Seurat::VariableFeatures(sample)[1:5],
                "B" = Seurat::VariableFeatures(sample)[6:10],
                "C" = Seurat::VariableFeatures(sample)[11:15])

  p <- SCpubr::do_EnrichmentHeatmap(sample = sample,
                                    input_gene_list = genes,
                                    group.by = "orig.ident",
                                    nbin = 1,
                                    ctrl = 10)
  testthat::expect_true("HeatmapList" %in% class(p))
})

testthat::test_that("do_EnrichmentHeatmap: PASS - group.by and flip", {

  sample$orig.ident <- ifelse(sample$seurat_clusters %in% c("1", "2"), "A", "B")

  genes <- list("A" = Seurat::VariableFeatures(sample)[1:5],
                "B" = Seurat::VariableFeatures(sample)[6:10],
                "C" = Seurat::VariableFeatures(sample)[11:15])

  p <- SCpubr::do_EnrichmentHeatmap(sample = sample,
                                    input_gene_list = genes,
                                    group.by = "orig.ident",
                                    flip = TRUE,
                                    nbin = 1,
                                    ctrl = 10)
  testthat::expect_true("HeatmapList" %in% class(p))
})

testthat::test_that("do_EnrichmentHeatmap: PASS - group.by and flip and column_names_rot", {

  sample$orig.ident <- ifelse(sample$seurat_clusters %in% c("1", "2"), "A", "B")

  genes <- list("A" = Seurat::VariableFeatures(sample)[1:5],
                "B" = Seurat::VariableFeatures(sample)[6:10],
                "C" = Seurat::VariableFeatures(sample)[11:15])

  p <- SCpubr::do_EnrichmentHeatmap(sample = sample,
                                    input_gene_list = genes,
                                    group.by = "orig.ident",
                                    flip = TRUE,
                                    column_names_rot = 0,
                                    nbin = 1,
                                    ctrl = 10)
  testthat::expect_true("HeatmapList" %in% class(p))
})

testthat::test_that("do_EnrichmentHeatmap: PASS - group.by and flip and row_names_rot", {

  sample$orig.ident <- ifelse(sample$seurat_clusters %in% c("1", "2"), "A", "B")

  genes <- list("A" = Seurat::VariableFeatures(sample)[1:5],
                "B" = Seurat::VariableFeatures(sample)[6:10],
                "C" = Seurat::VariableFeatures(sample)[11:15])

  p <- SCpubr::do_EnrichmentHeatmap(sample = sample,
                                    input_gene_list = genes,
                                    group.by = "orig.ident",
                                    flip = TRUE,
                                    row_names_rot = 90,
                                    nbin = 1,
                                    ctrl = 10)
  testthat::expect_true("HeatmapList" %in% class(p))
})



testthat::test_that("do_EnrichmentHeatmap: PASS - multiple variables changing cell size", {

  sample$orig.ident <- ifelse(sample$seurat_clusters %in% c("1", "2"), "A", "B")

  genes <- list("A" = Seurat::VariableFeatures(sample)[1:5],
                "B" = Seurat::VariableFeatures(sample)[6:10],
                "C" = Seurat::VariableFeatures(sample)[11:15])

  p <- SCpubr::do_EnrichmentHeatmap(sample = sample,
                                    input_gene_list = genes,
                                    group.by = "orig.ident",
                                    flip = TRUE,
                                    column_names_rot = 0,
                                    cell_size = 7,
                                    nbin = 1,
                                    ctrl = 10)
  testthat::expect_true("HeatmapList" %in% class(p))
})



testthat::test_that("do_EnrichmentHeatmap: PASS - character list of genes + group by only has 1 entity", {

  sample$orig.ident <- ifelse(sample$seurat_clusters %in% c("1", "2"), "A", "B")

  genes <- list("A" = Seurat::VariableFeatures(sample)[1:5],
                "B" = Seurat::VariableFeatures(sample)[6:10],
                "C" = Seurat::VariableFeatures(sample)[11:15])

  p <- SCpubr::do_EnrichmentHeatmap(sample = sample,
                                    input_gene_list = c("EPC1"),
                                    group.by = "orig.ident",
                                    nbin = 1,
                                    ctrl = 10)
  testthat::expect_true("HeatmapList" %in% class(p))
})

testthat::test_that("do_EnrichmentHeatmap: FAIL - list of genes without name", {

  sample$orig.ident <- ifelse(sample$seurat_clusters %in% c("1", "2"), "A", "B")

  genes <- list("A" = Seurat::VariableFeatures(sample)[1:5],
                "B" = Seurat::VariableFeatures(sample)[6:10],
                "C" = Seurat::VariableFeatures(sample)[11:15])

  testthat::expect_error(SCpubr::do_EnrichmentHeatmap(sample = sample,
                                                      input_gene_list = list("EPC1"),
                                                      group.by = "orig.ident",
                                                      nbin = 1,
                                                      ctrl = 10))

  testthat::expect_error(SCpubr::do_EnrichmentHeatmap(sample = sample,
                                                      input_gene_list = list("A" = "EPC1"),
                                                      group.by = "wrong",
                                                      nbin = 1,
                                                      ctrl = 10))

})

testthat::test_that("do_EnrichmentHeatmap: PASS - group by factor", {

  sample$orig.ident <- ifelse(sample$seurat_clusters %in% c("1", "2"), "A", "B")

  genes <- list("A" = Seurat::VariableFeatures(sample)[1:5],
                "B" = Seurat::VariableFeatures(sample)[6:10],
                "C" = Seurat::VariableFeatures(sample)[11:15])

  sample$seurat_clusters.factor <- factor(sample$seurat_clusters)
  p <- SCpubr::do_EnrichmentHeatmap(sample = sample,
                                    input_gene_list = c("EPC1"),
                                    group.by = "seurat_clusters.factor",
                                    nbin = 1,
                                    ctrl = 10)
  testthat::expect_true("HeatmapList" %in% class(p))
})

testthat::test_that("do_EnrichmentHeatmap: PASS - row title and column title", {

  sample$orig.ident <- ifelse(sample$seurat_clusters %in% c("1", "2"), "A", "B")

  genes <- list("A" = Seurat::VariableFeatures(sample)[1:5],
                "B" = Seurat::VariableFeatures(sample)[6:10],
                "C" = Seurat::VariableFeatures(sample)[11:15])

  sample$seurat_clusters.factor <- factor(sample$seurat_clusters)
  p <- SCpubr::do_EnrichmentHeatmap(sample = sample,
                                    input_gene_list = c("EPC1"),
                                    group.by = "seurat_clusters.factor",
                                    row_title = "A",
                                    column_title = "B",
                                    flip = TRUE,
                                    nbin = 1,
                                    ctrl = 10)
  testthat::expect_true("HeatmapList" %in% class(p))
})

