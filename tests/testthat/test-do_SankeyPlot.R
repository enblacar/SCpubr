if(isFALSE(dep_check[["do_SankeyPlot"]])){
  testthat::test_that("do_SankeyPlot: CRAN essentials", {

    p <- SCpubr::do_SankeyPlot(sample = sample,
                               first_group = "seurat_clusters",
                               last_group = "orig.ident")
    testthat::expect_type(p, "list")
  })

  testthat::test_that("do_SankeyPlot: PASS - default", {
    testthat::skip_on_cran()


    p <- SCpubr::do_SankeyPlot(sample = sample,
                               first_group = "seurat_clusters",
                               last_group = "orig.ident")
    testthat::expect_type(p, "list")
  })

  testthat::test_that("do_SankeyPlot: PASS - test colors when colors.use is null", {
    testthat::skip_on_cran()
    sample$first_group <- as.character(sample$orig.ident)
    sample$middle_group <- as.character(sample$orig.ident)
    sample$middle_group2 <- as.character(ifelse(sample$seurat_clusters %in% c("0"), "A", "B"))
    sample$last_group <- as.character(sample$orig.ident)
    sample$first_group_factor <- factor(sample$orig.ident)
    sample$middle_group_factor <- factor(sample$orig.ident)
    sample$middle_group2_factor <- factor(ifelse(sample$seurat_clusters %in% c("0"), "A", "B"))
    sample$last_group_factor <- factor(sample$orig.ident)

    p <- SCpubr::do_SankeyPlot(sample = sample,
                               first_group = "first_group",
                               middle_groups = "middle_group",
                               last_group = "last_group")
    testthat::expect_type(p, "list")

    p <- SCpubr::do_SankeyPlot(sample = sample,
                               first_group = "first_group",
                               middle_groups = "middle_group",
                               last_group = "last_group",
                               colors.middle = c("Cell" = "red"))
    testthat::expect_type(p, "list")


    testthat::expect_error({SCpubr::do_SankeyPlot(sample = sample,
                                                  first_group = "first_group",
                                                  middle_groups = "middle_group",
                                                  last_group = "last_group",
                                                  colors.middle = c("Cell?" = "red"))})

    testthat::expect_error({SCpubr::do_SankeyPlot(sample = sample,
                                                  first_group = "first_group",
                                                  middle_groups = "middle_group2",
                                                  last_group = "last_group",
                                                  colors.middle = c("A" = "red"))})

    p <- SCpubr::do_SankeyPlot(sample = sample,
                               first_group = "first_group_factor",
                               middle_groups = "middle_group_factor",
                               last_group = "last_group_factor")
    testthat::expect_type(p, "list")

    p <- SCpubr::do_SankeyPlot(sample = sample,
                               first_group = "first_group_factor",
                               middle_groups = "middle_group_factor",
                               last_group = "last_group_factor",
                               colors.middle = c("Cell" = "red"))
    testthat::expect_type(p, "list")

    testthat::expect_error({SCpubr::do_SankeyPlot(sample = sample,
                                                  first_group = "first_group_factor",
                                                  middle_groups = "middle_group_factor",
                                                  last_group = "last_group_factor",
                                                  colors.middle = c("Cell?" = "red"))})

    testthat::expect_error({SCpubr::do_SankeyPlot(sample = sample,
                                                  first_group = "first_group_factor",
                                                  middle_groups = "middle_group2_factor",
                                                  last_group = "last_group_factor",
                                                  colors.middle = c("A" = "red"))})
  })

  testthat::test_that("do_SankeyPlot: PASS - middle_groups", {
    testthat::skip_on_cran()


    sample$middle_group <- sample$seurat_clusters
    p <- SCpubr::do_SankeyPlot(sample = sample,
                               first_group = "seurat_clusters",
                               last_group = "orig.ident",
                               middle_group = "middle_group")
    testthat::expect_type(p, "list")
  })

  testthat::test_that("do_SankeyPlot: PASS - hjust", {

    testthat::skip_on_cran()

    sample$middle_group <- sample$seurat_clusters
    p <- SCpubr::do_SankeyPlot(sample = sample,
                               first_group = "seurat_clusters",
                               last_group = "orig.ident",
                               hjust = 0.5)
    testthat::expect_type(p, "list")
  })

  testthat::test_that("do_SankeyPlot: FAILS", {
    testthat::skip_on_cran()


    testthat::expect_error({SCpubr::do_SankeyPlot(sample = sample,
                                                  first_group = "seurat_clusters",
                                                  last_group = "orig.ident",
                                                  type = "wrong")})

    testthat::expect_error({SCpubr::do_SankeyPlot(sample = sample,
                                                  first_group = "seurat_clusters",
                                                  last_group = "orig.ident",
                                                  position = "wrong")})

    testthat::expect_error({SCpubr::do_SankeyPlot(sample = sample,
                                                  first_group = "wrong",
                                                  last_group = "orig.ident")})

    testthat::expect_error({SCpubr::do_SankeyPlot(sample = sample,
                                                  first_group = "nCount_RNA",
                                                  last_group = "orig.ident")})

    testthat::expect_error({SCpubr::do_SankeyPlot(sample = sample,
                                                  last_group = "wrong",
                                                  first_group = "orig.ident")})

    testthat::expect_error({SCpubr::do_SankeyPlot(sample = sample,
                                                  last_group = "nCount_RNA",
                                                  first_group = "orig.ident")})
    testthat::expect_error({SCpubr::do_SankeyPlot(sample = sample,
                                                  last_group = "seurat_clusters",
                                                  middle_groups = "wrong",
                                                  first_group = "orig.ident")})

    testthat::expect_error({SCpubr::do_SankeyPlot(sample = sample,
                                                  last_group = "seurat_clusters",
                                                  middle_groups = "nCount_RNA",
                                                  first_group = "orig.ident")})

    testthat::expect_error({SCpubr::do_SankeyPlot(sample = sample,
                                                  first_group = "seurat_clusters",
                                                  last_group = "orig.ident",
                                                  colors.first = c("A" = "red"))})
    testthat::expect_error({SCpubr::do_SankeyPlot(sample = sample,
                                                  first_group = "seurat_clusters",
                                                  last_group = "orig.ident",
                                                  colors.first = c("0" = "red"))})

    testthat::expect_error({SCpubr::do_SankeyPlot(sample = sample,
                                                  first_group = "seurat_clusters",
                                                  middle_groups = "middle_group",
                                                  last_group = "orig.ident",
                                                  colors.middle = c("A" = "red"))})
    testthat::expect_error({SCpubr::do_SankeyPlot(sample = sample,
                                                  first_group = "seurat_clusters",
                                                  middle_groups = "middle_group",
                                                  last_group = "orig.ident",
                                                  colors.middle = c("0" = "red"))})
    sample$orig.ident <- ifelse(sample$seurat_clusters == "0", "C", "B")
    testthat::expect_error({SCpubr::do_SankeyPlot(sample = sample,
                                                  first_group = "seurat_clusters",
                                                  last_group = "orig.ident",
                                                  colors.last = c("A" = "red"))})
    testthat::expect_error({SCpubr::do_SankeyPlot(sample = sample,
                                                  first_group = "seurat_clusters",
                                                  last_group = "orig.ident",
                                                  colors.last = c("B" = "red"))})
  })
}


