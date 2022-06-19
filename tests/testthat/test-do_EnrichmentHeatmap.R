sample <- use_dataset()
sample$orig.ident <- ifelse(sample$seurat_clusters %in% c("1", "2"), "A", "B")

genes <- list("A" = Seurat::VariableFeatures(sample)[1:5],
              "B" = Seurat::VariableFeatures(sample)[6:10],
              "C" = Seurat::VariableFeatures(sample)[11:15])

testthat::test_that("do_EnrichmentHeatmap: PASS - normal", {
  p <- SCpubr::do_EnrichmentHeatmap(sample = sample,
                                    list_genes = genes)
  testthat::expect_true("HeatmapList" %in% class(p))
})

testthat::test_that("do_EnrichmentHeatmap: PASS - group.by", {
  p <- SCpubr::do_EnrichmentHeatmap(sample = sample,
                                    list_genes = genes,
                                    group.by = "orig.ident")
  testthat::expect_true("HeatmapList" %in% class(p))
})

testthat::test_that("do_EnrichmentHeatmap: PASS - group.by and transpose", {
  p <- SCpubr::do_EnrichmentHeatmap(sample = sample,
                                    list_genes = genes,
                                    group.by = "orig.ident",
                                    transpose = T)
  testthat::expect_true("HeatmapList" %in% class(p))
})

testthat::test_that("do_EnrichmentHeatmap: PASS - group.by and transpose and column_names_rot", {
  p <- SCpubr::do_EnrichmentHeatmap(sample = sample,
                                    list_genes = genes,
                                    group.by = "orig.ident",
                                    transpose = T,
                                    column_names_rot = 0)
  testthat::expect_true("HeatmapList" %in% class(p))
})

testthat::test_that("do_EnrichmentHeatmap: PASS - group.by and transpose and row_names_rot", {
  p <- SCpubr::do_EnrichmentHeatmap(sample = sample,
                                    list_genes = genes,
                                    group.by = "orig.ident",
                                    transpose = T,
                                    row_names_rot = 90)
  testthat::expect_true("HeatmapList" %in% class(p))
})


testthat::test_that("do_EnrichmentHeatmap: PASS - multiple variables", {
  p <- SCpubr::do_EnrichmentHeatmap(sample = sample,
                                    list_genes = genes,
                                    group.by = "orig.ident",
                                    transpose = T,
                                    column_names_rot = 0,
                                    split.by = "seurat_clusters")
  testthat::expect_true("HeatmapList" %in% class(p))
})

testthat::test_that("do_EnrichmentHeatmap: PASS - multiple variables vertical", {
  p <- SCpubr::do_EnrichmentHeatmap(sample = sample,
                                    list_genes = genes,
                                    group.by = "orig.ident",
                                    transpose = T,
                                    column_names_rot = 0,
                                    split.by = "seurat_clusters",
                                    split.horizontal = F)
  testthat::expect_true("HeatmapList" %in% class(p))
})

testthat::test_that("do_EnrichmentHeatmap: PASS - multiple variables changing cell size", {
  p <- SCpubr::do_EnrichmentHeatmap(sample = sample,
                                    list_genes = genes,
                                    group.by = "orig.ident",
                                    transpose = T,
                                    column_names_rot = 0,
                                    split.by = "seurat_clusters",
                                    cell_size = 7)
  testthat::expect_true("HeatmapList" %in% class(p))
})

testthat::test_that("do_EnrichmentHeatmap: PASS - multiple variables changing color scale", {
  p <- SCpubr::do_EnrichmentHeatmap(sample = sample,
                                    list_genes = genes,
                                    group.by = "orig.ident",
                                    transpose = T,
                                    column_names_rot = 0,
                                    split.by = "seurat_clusters",
                                    colors.use = c("#e9d8a6", "#9b2226"))
  testthat::expect_true("HeatmapList" %in% class(p))
})
