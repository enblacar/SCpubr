library(magrittr)
small_a549_colData_df <- readRDS(system.file("extdata",
                                             "small_a549_dex_pdata.rda",
                                             package = "monocle3"))
small_a549_rowData_df <- readRDS(system.file("extdata",
                                             "small_a549_dex_fdata.rda",
                                             package = "monocle3"))
small_a549_exprs <- readRDS(system.file("extdata",
                                        "small_a549_dex_exprs.rda",
                                        package = "monocle3"))
small_a549_exprs <- small_a549_exprs[,row.names(small_a549_colData_df)]

cds <- monocle3::new_cell_data_set(expression_data = small_a549_exprs,
                                   cell_metadata = small_a549_colData_df,
                                   gene_metadata = small_a549_rowData_df)
suppressMessages({
  cds <- monocle3::preprocess_cds(cds, num_dim = 50)
  cds <- monocle3::reduce_dimension(cds, reduction_method = "PCA")
  cds <- monocle3::reduce_dimension(cds, reduction_method = "UMAP")
  cds <- monocle3::cluster_cells(cds)
})
sample <- Seurat::CreateSeuratObject(counts = small_a549_exprs)
sample <- Seurat::PercentageFeatureSet(sample, pattern = "^MT-", col.name = "percent.mt")
# Normalize.
sample <- suppressWarnings({Seurat::SCTransform(sample, verbose = FALSE)})
sample <- Seurat::RunPCA(sample, verbose = FALSE)
sample <- Seurat::RunUMAP(sample, dims = 1:30, verbose = FALSE)
# Find clusters.
sample <- Seurat::FindNeighbors(sample, dims = 1:30, verbose = FALSE)
sample <- Seurat::FindClusters(sample, resolution = 0.5, verbose = FALSE)



pseudotime_genes <- Seurat::VariableFeatures(sample)[1:5]
testthat::test_that("do_PseudotimePlot: PASS - default", {
 out <- SCpubr::do_PseudotimePlot(sample = sample,
                                  cds = cds,
                                  pseudotime_genes = pseudotime_genes,
                                  nbin = 10,
                                  ctrl = 10)
  testthat::expect_type(out, "list")
  testthat::expect_length(out, 4)
})

testthat::test_that("do_PseudotimePlot: PASS - max score", {
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
})

testthat::test_that("do_PseudotimePlot: PASS - group.by", {
  out <- SCpubr::do_PseudotimePlot(sample = sample,
                                   cds = cds,
                                   pseudotime_genes = pseudotime_genes,
                                   nbin = 10,
                                   ctrl = 10,
                                   group.by = "orig.ident")
  testthat::expect_type(out, "list")
  testthat::expect_length(out, 4)
})

testthat::test_that("do_PseudotimePlot: PASS - combinations", {
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
})


testthat::test_that("do_PseudotimePlot: FAIL - wrong parameters", {
  testthat::expect_error(SCpubr::do_PseudotimePlot(sample = sample,
                                                   cds = cds,
                                                   pseudotime_genes = pseudotime_genes,
                                                   nbin = 10,
                                                   ctrl = 10,
                                                   group.by = "wrong"))
})
