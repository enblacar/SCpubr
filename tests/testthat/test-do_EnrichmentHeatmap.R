if (isFALSE(dep_check[["do_EnrichmentHeatmap"]])){

  testthat::test_that("do_EnrichmentHeatmap: CRAN essential", {
    sample$orig.ident <- ifelse(sample$seurat_clusters %in% c("1", "2"), "A", "B")

    genes <- list("A" = rownames(sample)[1:5],
                  "B" = rownames(sample)[6:10],
                  "C" = rownames(sample)[11:15])


    p <- SCpubr::do_EnrichmentHeatmap(sample = sample,
                                      input_gene_list = genes,
                                      nbin = 1,
                                      ctrl = 10)
    testthat::expect_true("ggplot" %in% class(p))
  })

  testthat::test_that("do_EnrichmentHeatmap: PASS -flavors", {
    testthat::skip_on_cran()
    genes <- list("A" = rownames(sample)[1:5],
                  "B" = rownames(sample)[6:10],
                  "C" = rownames(sample)[11:15])
    
    p <- SCpubr::do_EnrichmentHeatmap(sample = sample,
                                      input_gene_list = genes,
                                      group.by = c("seurat_clusters", "orig.ident"),
                                      nbin = 1,
                                      ctrl = 10,
                                      enforce_symmetry = TRUE)
    testthat::expect_true("ggplot" %in% class(p))
    
    p <- SCpubr::do_EnrichmentHeatmap(sample = sample,
                                      input_gene_list = genes,
                                      group.by = c("seurat_clusters", "orig.ident"),
                                      nbin = 1,
                                      ctrl = 10,
                                      enforce_symmetry = FALSE,
                                      use_viridis = TRUE,
                                      viridis.direction = 1)
    testthat::expect_true("ggplot" %in% class(p))
    
    p <- SCpubr::do_EnrichmentHeatmap(sample = sample,
                                      input_gene_list = genes,
                                      group.by = c("seurat_clusters", "orig.ident"),
                                      nbin = 1,
                                      ctrl = 10,
                                      enforce_symmetry = FALSE,
                                      use_viridis = TRUE,
                                      viridis.direction = -1)
    testthat::expect_true("ggplot" %in% class(p))
    
    p <- SCpubr::do_EnrichmentHeatmap(sample = sample,
                                      input_gene_list = genes,
                                      group.by = c("seurat_clusters", "orig.ident"),
                                      nbin = 1,
                                      ctrl = 10,
                                      enforce_symmetry = FALSE,
                                      use_viridis = FALSE,
                                      sequential.direction = 1)
    testthat::expect_true("ggplot" %in% class(p))
    
    p <- SCpubr::do_EnrichmentHeatmap(sample = sample,
                                      input_gene_list = genes,
                                      group.by = c("seurat_clusters", "orig.ident"),
                                      nbin = 1,
                                      ctrl = 10,
                                      enforce_symmetry = FALSE,
                                      use_viridis = FALSE,
                                      sequential.direction = -1)
    testthat::expect_true("ggplot" %in% class(p))
    
    genes <- list("A" = rownames(sample)[1:5])
    
    p <- SCpubr::do_EnrichmentHeatmap(sample = sample,
                                      input_gene_list = genes,
                                      nbin = 1,
                                      ctrl = 10,
                                      flip = FALSE)
    testthat::expect_true("ggplot" %in% class(p))
    
    p <- SCpubr::do_EnrichmentHeatmap(sample = sample,
                                      input_gene_list = genes,
                                      nbin = 1,
                                      ctrl = 10,
                                      flip = TRUE)
    testthat::expect_true("ggplot" %in% class(p))
    
    genes <- list("A" = rownames(sample)[1:5],
                  "B" = rownames(sample)[6:10],
                  "C" = rownames(sample)[11:15])
    
    p <- SCpubr::do_EnrichmentHeatmap(sample = sample,
                                      input_gene_list = genes,
                                      nbin = 1,
                                      ctrl = 10,
                                      geneset.order = c("A", "B", "C"))
    testthat::expect_true("ggplot" %in% class(p))
    
    
    p <- SCpubr::do_EnrichmentHeatmap(sample = sample,
                                 input_gene_list = genes,
                                 flavor = "Seurat",
                                 assay = "SCT",
                                 nbin = 1,
                                 ctrl = 10,
                                 viridis.direction = 1)
    testthat::expect_true("ggplot" %in% class(p))

    p <- SCpubr::do_EnrichmentHeatmap(sample = sample,
                                      input_gene_list = genes,
                                      flavor = "Seurat",
                                      assay = "SCT",
                                      nbin = 1,
                                      ctrl = 10,
                                      viridis.direction = -1)
    testthat::expect_true("ggplot" %in% class(p))

    p <- SCpubr::do_EnrichmentHeatmap(sample = sample,
                                 input_gene_list = genes,
                                 flavor = "UCell",
                                 slot = "data",
                                 nbin = 1,
                                 ctrl = 10,
                                 viridis.direction = 1)
    testthat::expect_true("ggplot" %in% class(p))

    p <- SCpubr::do_EnrichmentHeatmap(sample = sample,
                                      input_gene_list = genes,
                                      flavor = "UCell",
                                      slot = "data",
                                      nbin = 1,
                                      ctrl = 10,
                                      viridis.direction = -1)
    testthat::expect_true("ggplot" %in% class(p))

    p <- SCpubr::do_EnrichmentHeatmap(sample = sample,
                                      input_gene_list = genes,
                                      flavor = "AUCell",
                                      assay = "SCT",
                                      nbin = 1,
                                      ctrl = 10,
                                      viridis.direction = 1)
    testthat::expect_true("ggplot" %in% class(p))

    p <- SCpubr::do_EnrichmentHeatmap(sample = sample,
                                      input_gene_list = genes,
                                      flavor = "AUCell",
                                      assay = "SCT",
                                      nbin = 1,
                                      ctrl = 10,
                                      viridis.direction = -1)
    testthat::expect_true("ggplot" %in% class(p))
  })

  testthat::test_that("do_EnrichmentHeatmap: PASS - normal", {
    testthat::skip_on_cran()

    sample$orig.ident <- ifelse(sample$seurat_clusters %in% c("1", "2"), "A", "B")

    genes <- list("A" = rownames(sample)[1:5],
                  "B" = rownames(sample)[6:10],
                  "C" = rownames(sample)[11:15])


    p <- SCpubr::do_EnrichmentHeatmap(sample = sample,
                                      input_gene_list = genes,
                                      nbin = 1,
                                      ctrl = 10)
    testthat::expect_true("ggplot" %in% class(p))

    p <- SCpubr::do_EnrichmentHeatmap(sample = sample,
                                      input_gene_list = genes,
                                      legend.position = "top",
                                      nbin = 1,
                                      ctrl = 10)
    testthat::expect_true("ggplot" %in% class(p))

    p <- SCpubr::do_EnrichmentHeatmap(sample = sample,
                                      input_gene_list = genes,
                                      legend.position = "right",
                                      nbin = 1,
                                      ctrl = 10)
    testthat::expect_true("ggplot" %in% class(p))
  })

  testthat::test_that("do_EnrichmentHeatmap: PASS - group.by", {
    testthat::skip_on_cran()

    sample$orig.ident <- ifelse(sample$seurat_clusters %in% c("1", "2"), "A", "B")

    genes <- list("A" = rownames(sample)[1:5],
                  "B" = rownames(sample)[6:10],
                  "C" = rownames(sample)[11:15])

    p <- SCpubr::do_EnrichmentHeatmap(sample = sample,
                                      input_gene_list = genes,
                                      group.by = "orig.ident",
                                      nbin = 1,
                                      ctrl = 10)
    testthat::expect_true("ggplot" %in% class(p))
  })

  testthat::test_that("do_EnrichmentHeatmap: PASS - group.by and flip", {
    testthat::skip_on_cran()

    sample$orig.ident <- ifelse(sample$seurat_clusters %in% c("1", "2"), "A", "B")

    genes <- list("A" = rownames(sample)[1:5],
                  "B" = rownames(sample)[6:10],
                  "C" = rownames(sample)[11:15])

    p <- SCpubr::do_EnrichmentHeatmap(sample = sample,
                                      input_gene_list = genes,
                                      group.by = c("orig.ident", "seurat_clusters"),
                                      flip = TRUE,
                                      nbin = 1,
                                      ctrl = 10)
    testthat::expect_true("ggplot" %in% class(p))

    p <- SCpubr::do_EnrichmentHeatmap(sample = sample,
                                      input_gene_list = genes,
                                      group.by = c("orig.ident"),
                                      flip = TRUE,
                                      nbin = 1,
                                      ctrl = 10)
    testthat::expect_true("ggplot" %in% class(p))
  })









  testthat::test_that("do_EnrichmentHeatmap: PASS - character list of genes + group by only has 1 entity", {
    testthat::skip_on_cran()

    sample$orig.ident <- ifelse(sample$seurat_clusters %in% c("1", "2"), "A", "B")

    genes <- list("A" = rownames(sample)[1:5],
                  "B" = rownames(sample)[6:10],
                  "C" = rownames(sample)[11:15])

    p <- SCpubr::do_EnrichmentHeatmap(sample = sample,
                                      input_gene_list = genes,
                                      group.by = "orig.ident",
                                      nbin = 1,
                                      ctrl = 10)
    testthat::expect_true("ggplot" %in% class(p))

    p <- SCpubr::do_EnrichmentHeatmap(sample = sample,
                                      input_gene_list = genes,
                                      group.by = c("seurat_clusters", "orig.ident"),
                                      nbin = 1,
                                      ctrl = 10)
    testthat::expect_true("ggplot" %in% class(p))

    p <- SCpubr::do_EnrichmentHeatmap(sample = sample,
                                      input_gene_list = genes,
                                      group.by = c("seurat_clusters", "orig.ident"),
                                      nbin = 1,
                                      ctrl = 10,
                                      return_object = TRUE)
    testthat::expect_true("list" %in% class(p))

    p <- SCpubr::do_EnrichmentHeatmap(sample = sample,
                                      input_gene_list = genes,
                                      nbin = 1,
                                      ctrl = 10,
                                      return_object = TRUE)
    testthat::expect_true("list" %in% class(p))
  })

  testthat::test_that("do_EnrichmentHeatmap: FAIL - list of genes without name", {
    testthat::skip_on_cran()

    sample$orig.ident <- ifelse(sample$seurat_clusters %in% c("1", "2"), "A", "B")

    genes <- list("A" = rownames(sample)[1:5],
                  "B" = rownames(sample)[6:10],
                  "C" = rownames(sample)[11:15])

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

    genes <- list("A" = rownames(sample)[1:5],
                  "B" = rownames(sample)[6:10],
                  "C" = rownames(sample)[11:15])

    sample$seurat_clusters.factor <- factor(sample$seurat_clusters)
    p <- SCpubr::do_EnrichmentHeatmap(sample = sample,
                                      input_gene_list = genes,
                                      group.by = "seurat_clusters.factor",
                                      nbin = 1,
                                      ctrl = 10)
    testthat::expect_true("ggplot" %in% class(p))
  })


  testthat::test_that("do_EnrichmentHeatmap: ERROR - wrong arguments", {
    testthat::skip_on_cran()

    sample$orig.ident <- ifelse(sample$seurat_clusters %in% c("1", "2"), "A", "B")

    genes <- list("A" = rownames(sample)[1:5],
                  "B" = rownames(sample)[6:10],
                  "C" = rownames(sample)[11:15])

    testthat::expect_warning({SCpubr::do_EnrichmentHeatmap(sample = sample,
                                      input_gene_list = genes,
                                      flavor = "Seurat",
                                      slot = "data",
                                      nbin = 1,
                                      ctrl = 10)})

    testthat::expect_warning({SCpubr::do_EnrichmentHeatmap(sample = sample,
                                                        input_gene_list = genes,
                                                        flavor = "UCell",
                                                        assay = "SCT",
                                                        nbin = 1,
                                                        ctrl = 10)})

    testthat::expect_error({SCpubr::do_EnrichmentHeatmap(sample = sample,
                                                         input_gene_list = genes,
                                                         min.cutoff = -10,
                                                         nbin = 1,
                                                         ctrl = 10)})

    testthat::expect_error({SCpubr::do_EnrichmentHeatmap(sample = sample,
                                                        input_gene_list = genes,
                                                        max.cutoff = 200,
                                                        nbin = 1,
                                                        ctrl = 10)})

    testthat::expect_error({SCpubr::do_EnrichmentHeatmap(sample = sample,
                                                         input_gene_list = genes,
                                                         max.cutoff = 1,
                                                         min.cutoff = 2,
                                                         nbin = 1,
                                                         ctrl = 10)})



    testthat::expect_error({SCpubr::do_EnrichmentHeatmap(sample = sample,
                                                         input_gene_list = list("A" = c("EPC1")),
                                                         group.by = c("seurat_clusters", "annotation"),
                                                         nbin = 1,
                                                         ctrl = 10,
                                                         flip = FALSE,
                                                         geneset.order = "wrong")})
    
    testthat::expect_warning({SCpubr::do_EnrichmentHeatmap(sample = sample,
                                                         input_gene_list = list("A_A" = c("EPC1")),
                                                         group.by = c("seurat_clusters", "annotation"),
                                                         nbin = 1,
                                                         ctrl = 10,
                                                         flip = FALSE)})
    
    testthat::expect_error({SCpubr::do_EnrichmentHeatmap(sample = sample,
                                                         input_gene_list = c("EPC1"),
                                                         nbin = 1,
                                                         ctrl = 10)})
  })
}


