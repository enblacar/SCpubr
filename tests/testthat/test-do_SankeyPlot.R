sample <- SCpubr:::test_list$sample

testthat::test_that("do_SankeyPlot: FAIL - no dplyr", {
  testthat::expect_error({SCpubr::do_SankeyPlot(sample = sample,
                                                first_group = "seurat_clusters",
                                                last_group = "orig.ident")})
})
suppressMessages(library(dplyr))
testthat::test_that("do_SankeyPlot: PASS - default", {
  p <- SCpubr::do_SankeyPlot(sample = sample,
                             first_group = "seurat_clusters",
                             last_group = "orig.ident")
  testthat::expect_type(p, "list")
})


testthat::test_that("do_SankeyPlot: PASS - middle_groups", {
  sample$middle_group <- sample$seurat_clusters
  p <- SCpubr::do_SankeyPlot(sample = sample,
                             first_group = "seurat_clusters",
                             last_group = "orig.ident",
                             middle_group = "middle_group")
  testthat::expect_type(p, "list")
})

testthat::test_that("do_SankeyPlot: PASS - hjust", {
  sample$middle_group <- sample$seurat_clusters
  p <- SCpubr::do_SankeyPlot(sample = sample,
                             first_group = "seurat_clusters",
                             last_group = "orig.ident",
                             hjust = 0.5)
  testthat::expect_type(p, "list")
})

testthat::test_that("do_SankeyPlot: FAILS", {
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
