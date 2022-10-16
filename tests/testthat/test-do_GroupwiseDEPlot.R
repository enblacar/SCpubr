if(isFALSE(dep_check[["do_GroupwiseDEPlot"]])){
  testthat::test_that("do_GroupwiseDEPlot: PASS - default", {



    sample@assays$SCT@scale.data <- as.matrix(sample@assays$SCT@data)

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


    sample@assays$SCT@scale.data <- as.matrix(sample@assays$SCT@data)

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


    sample@assays$SCT@scale.data <- as.matrix(sample@assays$SCT@data)

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


    sample@assays$SCT@scale.data <- as.matrix(sample@assays$SCT@data)

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


    sample@assays$SCT@scale.data <- as.matrix(sample@assays$SCT@data)

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

}
