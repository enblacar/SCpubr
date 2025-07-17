if(base::isFALSE(dep_check[["do_GroupwiseDEHeatmap"]])){
  testthat::test_that("do_GroupwiseDEHeatmap: CRAN essentials", {
    suppressWarnings({
    sample <- SeuratObject::SetAssayData(object = sample,
                                         assay = "SCT",
                                         slot = "scale.data",
                                         new.data = as.matrix(SeuratObject::GetAssayData(object = sample,
                                                                                         assay = "SCT",
                                                                                         slot = "data")))
    })

    p <- SCpubr::do_GroupwiseDEHeatmap(sample = sample,
                                    de_genes = de_genes,
                                    assay = "SCT",
                                    slot = "data")
    testthat::expect_true(ggplot2::is_ggplot(p))

    p <- SCpubr::do_GroupwiseDEHeatmap(sample = sample,
                                    de_genes = de_genes_scaled,
                                    assay = "SCT",
                                    slot = "scale.data")
    testthat::expect_true(ggplot2::is_ggplot(p))
  })

  testthat::test_that("do_GroupwiseDEHeatmap: PASS - default", {
    testthat::skip_on_cran()


    suppressWarnings({
      sample <- SeuratObject::SetAssayData(object = sample,
                                           assay = "SCT",
                                           slot = "scale.data",
                                           new.data = as.matrix(SeuratObject::GetAssayData(object = sample,
                                                                                       assay = "SCT",
                                                                                       slot = "data")))
    })

    p <- SCpubr::do_GroupwiseDEHeatmap(sample = sample,
                                    de_genes = de_genes,
                                    assay = "SCT",
                                    slot = "data",
                                    use_viridis = FALSE,
                                    sequential.direction = 1)
    testthat::expect_true(ggplot2::is_ggplot(p))
    
    p <- SCpubr::do_GroupwiseDEHeatmap(sample = sample,
                                    de_genes = de_genes,
                                    assay = "SCT",
                                    slot = "data",
                                    use_viridis = FALSE,
                                    sequential.direction = -1)
    testthat::expect_true(ggplot2::is_ggplot(p))
    
    p <- SCpubr::do_GroupwiseDEHeatmap(sample = sample,
                                    de_genes = de_genes,
                                    assay = "SCT",
                                    slot = "data",
                                    use_viridis = TRUE,
                                    viridis.direction = 1)
    testthat::expect_true(ggplot2::is_ggplot(p))
    
    p <- SCpubr::do_GroupwiseDEHeatmap(sample = sample,
                                    de_genes = de_genes,
                                    assay = "SCT",
                                    slot = "data",
                                    use_viridis = TRUE,
                                    viridis.direction = -1)
    testthat::expect_true(ggplot2::is_ggplot(p))

    p <- SCpubr::do_GroupwiseDEHeatmap(sample = sample,
                                    de_genes = de_genes_scaled,
                                    assay = "SCT",
                                    slot = "scale.data")
    testthat::expect_true(ggplot2::is_ggplot(p))


    p <- SCpubr::do_GroupwiseDEHeatmap(sample = sample,
                                    de_genes = de_genes,
                                    assay = "SCT",
                                    slot = "data",
                                    viridis.direction = 1,
                                    max.cutoff = 1.2,
                                    min.cutoff = 1)
    testthat::expect_true(ggplot2::is_ggplot(p))

    p <- SCpubr::do_GroupwiseDEHeatmap(sample = sample,
                                    de_genes = de_genes,
                                    assay = "SCT",
                                    slot = "data",
                                    viridis.direction = 1,
                                    min.cutoff = 1)
    testthat::expect_true(ggplot2::is_ggplot(p))

    p <- SCpubr::do_GroupwiseDEHeatmap(sample = sample,
                                    de_genes = de_genes,
                                    assay = "SCT",
                                    slot = "data",
                                    viridis.direction = 1,
                                    max.cutoff = 1.2)
    testthat::expect_true(ggplot2::is_ggplot(p))

  })

  testthat::test_that("do_GroupwiseDEHeatmap: PASS - heatmap legend side", {
    testthat::skip_on_cran()

    suppressWarnings({
      sample <- SeuratObject::SetAssayData(object = sample,
                                           assay = "SCT",
                                           slot = "scale.data",
                                           new.data = as.matrix(SeuratObject::GetAssayData(object = sample,
                                                                                       assay = "SCT",
                                                                                       slot = "data")))
    })

    p <- SCpubr::do_GroupwiseDEHeatmap(sample = sample,
                                    de_genes = de_genes,
                                    assay = "SCT",
                                    slot = "data",
                                    legend.position = "right")
    testthat::expect_true(ggplot2::is_ggplot(p))

    p <- SCpubr::do_GroupwiseDEHeatmap(sample = sample,
                                    de_genes = de_genes_scaled,
                                    assay = "SCT",
                                    slot = "scale.data",
                                    legend.position = "right")
    testthat::expect_true(ggplot2::is_ggplot(p))
  })


  testthat::test_that("do_GroupwiseDEHeatmap: FAIL - wrong direction", {
    testthat::skip_on_cran()

    testthat::expect_error({SCpubr::do_GroupwiseDEHeatmap(sample = sample,
                                                       de_genes = de_genes,
                                                       assay = "SCT",
                                                       slot = "data",
                                                       viridis.direction = 0)})
  })

  testthat::test_that("do_ExpressionHeatmap: FAIL", {
    testthat::skip_on_cran()
    testthat::expect_error({SCpubr::do_GroupwiseDEHeatmap(sample = sample,
                                                       de_genes = de_genes,
                                                       assay = "SCT",
                                                       slot = "data",
                                                       viridis.direction = 1,
                                                       min.cutoff = -10)})

    testthat::expect_error({SCpubr::do_GroupwiseDEHeatmap(sample = sample,
                                                       de_genes = de_genes,
                                                       assay = "SCT",
                                                       slot = "data",
                                                       viridis.direction = 1,
                                                       max.cutoff = 200)})

    testthat::expect_error({SCpubr::do_GroupwiseDEHeatmap(sample = sample,
                                                       de_genes = de_genes,
                                                       assay = "SCT",
                                                       slot = "data",
                                                       viridis.direction = 1,
                                                       max.cutoff = 1,
                                                       min.cutoff = 2)})

  })

}

