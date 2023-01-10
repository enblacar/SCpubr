if (isFALSE(dep_check[["do_EnrichmentHeatmap"]])){

  testthat::test_that("do_EnrichmentHeatmap: CRAN essential", {
    sample$orig.ident <- ifelse(sample$seurat_clusters %in% c("1", "2"), "A", "B")

    genes <- list("A" = Seurat::VariableFeatures(sample)[1:5],
                  "B" = Seurat::VariableFeatures(sample)[6:10],
                  "C" = Seurat::VariableFeatures(sample)[11:15])


    p <- SCpubr::do_EnrichmentHeatmap(sample = sample,
                                      input_gene_list = genes,
                                      nbin = 1,
                                      ctrl = 10)
    testthat::expect_true("HeatmapList" %in% class(p))
  })

  testthat::test_that("do_EnrichmentHeatmap: PASS -flavors", {
    p <- SCpubr::do_EnrichmentHeatmap(sample = sample,
                                 input_gene_list = c("EPC1"),
                                 flavor = "Seurat",
                                 assay = "SCT",
                                 nbin = 1,
                                 ctrl = 10,
                                 viridis_direction = 1)
    testthat::expect_true("HeatmapList" %in% class(p))

    p <- SCpubr::do_EnrichmentHeatmap(sample = sample,
                                      input_gene_list = c("EPC1"),
                                      flavor = "Seurat",
                                      assay = "SCT",
                                      nbin = 1,
                                      ctrl = 10,
                                      viridis_direction = -1)
    testthat::expect_true("HeatmapList" %in% class(p))

    p <- SCpubr::do_EnrichmentHeatmap(sample = sample,
                                 input_gene_list = c("EPC1"),
                                 flavor = "UCell",
                                 slot = "data",
                                 nbin = 1,
                                 ctrl = 10,
                                 viridis_direction = 1)
    testthat::expect_true("HeatmapList" %in% class(p))

    p <- SCpubr::do_EnrichmentHeatmap(sample = sample,
                                      input_gene_list = c("EPC1"),
                                      flavor = "UCell",
                                      slot = "data",
                                      nbin = 1,
                                      ctrl = 10,
                                      viridis_direction = -1)
    testthat::expect_true("HeatmapList" %in% class(p))

    p <- SCpubr::do_EnrichmentHeatmap(sample = sample,
                                      input_gene_list = c("EPC1"),
                                      flavor = "AUCell",
                                      assay = "SCT",
                                      nbin = 1,
                                      ctrl = 10,
                                      viridis_direction = 1)
    testthat::expect_true("HeatmapList" %in% class(p))

    p <- SCpubr::do_EnrichmentHeatmap(sample = sample,
                                      input_gene_list = c("EPC1"),
                                      flavor = "AUCell",
                                      assay = "SCT",
                                      nbin = 1,
                                      ctrl = 10,
                                      viridis_direction = -1)
    testthat::expect_true("HeatmapList" %in% class(p))
  })

  testthat::test_that("do_EnrichmentHeatmap: PASS - normal", {
    testthat::skip_on_cran()

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
    testthat::skip_on_cran()

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
    testthat::skip_on_cran()

    sample$orig.ident <- ifelse(sample$seurat_clusters %in% c("1", "2"), "A", "B")

    genes <- list("A" = Seurat::VariableFeatures(sample)[1:5],
                  "B" = Seurat::VariableFeatures(sample)[6:10],
                  "C" = Seurat::VariableFeatures(sample)[11:15])

    p <- SCpubr::do_EnrichmentHeatmap(sample = sample,
                                      input_gene_list = genes,
                                      group.by = c("orig.ident", "seurat_clusters"),
                                      flip = TRUE,
                                      nbin = 1,
                                      ctrl = 10)
    testthat::expect_true("HeatmapList" %in% class(p))

    p <- SCpubr::do_EnrichmentHeatmap(sample = sample,
                                      input_gene_list = genes,
                                      group.by = c("orig.ident"),
                                      flip = TRUE,
                                      nbin = 1,
                                      ctrl = 10)
    testthat::expect_true("HeatmapList" %in% class(p))
  })

  testthat::test_that("do_EnrichmentHeatmap: PASS - group.by and flip and column_names_rot", {
    testthat::skip_on_cran()

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
    testthat::skip_on_cran()

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
    testthat::skip_on_cran()

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
    testthat::skip_on_cran()

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

    p <- SCpubr::do_EnrichmentHeatmap(sample = sample,
                                      input_gene_list = genes,
                                      group.by = c("seurat_clusters", "orig.ident"),
                                      nbin = 1,
                                      ctrl = 10)
    testthat::expect_true("HeatmapList" %in% class(p))

    p <- SCpubr::do_EnrichmentHeatmap(sample = sample,
                                      input_gene_list = genes,
                                      group.by = c("seurat_clusters", "orig.ident"),
                                      nbin = 1,
                                      ctrl = 10,
                                      return_object = TRUE,
                                      return_matrix = TRUE)
    testthat::expect_true("list" %in% class(p))

    p <- SCpubr::do_EnrichmentHeatmap(sample = sample,
                                      input_gene_list = genes,
                                      nbin = 1,
                                      ctrl = 10,
                                      return_object = TRUE,
                                      return_matrix = TRUE)
    testthat::expect_true("list" %in% class(p))
  })

  testthat::test_that("do_EnrichmentHeatmap: FAIL - list of genes without name", {
    testthat::skip_on_cran()

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
    testthat::skip_on_cran()

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
    testthat::skip_on_cran()

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

    p <- SCpubr::do_EnrichmentHeatmap(sample = sample,
                                      input_gene_list = c("EPC1"),
                                      group.by = "seurat_clusters.factor",
                                      row_title = "A",
                                      column_title = "B",
                                      flip = FALSE,
                                      nbin = 1,
                                      ctrl = 10)
    testthat::expect_true("HeatmapList" %in% class(p))



    p <- SCpubr::do_EnrichmentHeatmap(sample = sample,
                                      input_gene_list = c("EPC1"),
                                      group.by = c("seurat_clusters", "annotation"),
                                      row_title = c("A", "B"),
                                      column_title = c("C", "D"),
                                      flip = TRUE,
                                      nbin = 1,
                                      ctrl = 10)
    testthat::expect_true("HeatmapList" %in% class(p))

    p <- SCpubr::do_EnrichmentHeatmap(sample = sample,
                                      input_gene_list = c("EPC1"),
                                      group.by = c("seurat_clusters", "annotation"),
                                      row_title = c("A", "B"),
                                      column_title = c("C", "D"),
                                      flip = FALSE,
                                      nbin = 1,
                                      ctrl = 10)
    testthat::expect_true("HeatmapList" %in% class(p))


  })

  testthat::test_that("do_EnrichmentHeatmap: PASS - row title and column title", {
    testthat::skip_on_cran()

    p <- SCpubr::do_EnrichmentHeatmap(sample = sample,
                                      input_gene_list = c("EPC1"),
                                      group.by = "seurat_clusters",
                                      row_title = "A",
                                      column_title = "B",
                                      flip = TRUE,
                                      nbin = 1,
                                      ctrl = 10,
                                      plot_FeaturePlots = TRUE,
                                      plot_GeyserPlots = TRUE,
                                      plot_BeeSwarmPlots = TRUE,
                                      plot_BoxPlots = TRUE,
                                      plot_ViolinPlots = TRUE)
    testthat::expect_type(p, "list")

    p <- SCpubr::do_EnrichmentHeatmap(sample = sample,
                                      input_gene_list = c("EPC1"),
                                      group.by = "seurat_clusters",
                                      row_title = "A",
                                      column_title = "B",
                                      flip = TRUE,
                                      nbin = 1,
                                      ctrl = 10,
                                      plot_FeaturePlots = TRUE,
                                      plot_GeyserPlots = TRUE,
                                      plot_BeeSwarmPlots = TRUE,
                                      plot_BoxPlots = TRUE,
                                      plot_ViolinPlots = TRUE,
                                      assay = "SCT",
                                      slot = "data",
                                      reduction = "umap",
                                      flavor = "AUCell")
    testthat::expect_type(p, "list")

    p <- SCpubr::do_EnrichmentHeatmap(sample = sample,
                                      input_gene_list = c("EPC1"),
                                      group.by = "seurat_clusters",
                                      row_title = "A",
                                      column_title = "B",
                                      flip = TRUE,
                                      nbin = 1,
                                      ctrl = 10,
                                      plot_FeaturePlots = TRUE,
                                      plot_GeyserPlots = TRUE,
                                      min.cutoff = 0.1,
                                      max.cutoff = 0.5,
                                      plot_BeeSwarmPlots = TRUE,
                                      plot_BoxPlots = TRUE,
                                      plot_ViolinPlots = TRUE)
    testthat::expect_type(p, "list")


  })

  testthat::test_that("do_EnrichmentHeatmap: ERROR - wrong arguments", {
    testthat::skip_on_cran()

    sample$orig.ident <- ifelse(sample$seurat_clusters %in% c("1", "2"), "A", "B")

    genes <- list("A" = Seurat::VariableFeatures(sample)[1:5],
                  "B" = Seurat::VariableFeatures(sample)[6:10],
                  "C" = Seurat::VariableFeatures(sample)[11:15])

    testthat::expect_error({SCpubr::do_EnrichmentHeatmap(sample = sample,
                                      input_gene_list = c("EPC1"),
                                      flavor = "Seurat",
                                      slot = "data",
                                      nbin = 1,
                                      ctrl = 10)})

    testthat::expect_error({SCpubr::do_EnrichmentHeatmap(sample = sample,
                                                        input_gene_list = c("EPC1"),
                                                        flavor = "UCell",
                                                        assay = "SCT",
                                                        nbin = 1,
                                                        ctrl = 10)})

    testthat::expect_error({SCpubr::do_EnrichmentHeatmap(sample = sample,
                                                         input_gene_list = c("EPC1"),
                                                         min.cutoff = -10,
                                                         nbin = 1,
                                                         ctrl = 10)})

    testthat::expect_error({SCpubr::do_EnrichmentHeatmap(sample = sample,
                                                        input_gene_list = c("EPC1"),
                                                        max.cutoff = 200,
                                                        nbin = 1,
                                                        ctrl = 10)})

    testthat::expect_error({SCpubr::do_EnrichmentHeatmap(sample = sample,
                                                         input_gene_list = c("EPC1"),
                                                         max.cutoff = 1,
                                                         min.cutoff = 2,
                                                         nbin = 1,
                                                         ctrl = 10)})

    testthat::expect_error({SCpubr::do_EnrichmentHeatmap(sample = sample,
                                                         input_gene_list = c("EPC1"),
                                                         group.by = c("seurat_clusters", "annotation"),
                                                         column_title = c("A"),
                                                         nbin = 1,
                                                         ctrl = 10,
                                                         flip = TRUE)})

    testthat::expect_error({SCpubr::do_EnrichmentHeatmap(sample = sample,
                                                         input_gene_list = c("EPC1"),
                                                         group.by = c("seurat_clusters", "annotation"),
                                                         row_title = c("A"),
                                                         nbin = 1,
                                                         ctrl = 10,
                                                         flip = FALSE)})
  })
}


