if(isFALSE(dep_check[["do_GeyserPlot"]])){
  testthat::test_that("do_GeyserPlot: PASS - default parameters", {



    p <- SCpubr::do_GeyserPlot(sample = sample,
                               features = "EPC1")
    testthat::expect_type(p, "list")
  })

  testthat::test_that("do_GeyserPlot: PASS - default parameters = symmetrical scale", {



    p <- SCpubr::do_GeyserPlot(sample = sample,
                               features = "EPC1",
                               enforce_symmetry = TRUE)
    testthat::expect_type(p, "list")

    p <- SCpubr::do_GeyserPlot(sample = sample,
                               features = "EPC1",
                               enforce_symmetry = FALSE)
    testthat::expect_type(p, "list")
  })


  testthat::test_that("do_GeyserPlot: PASS - default parameters = categorical scale", {



    p <- SCpubr::do_GeyserPlot(sample = sample,
                               features = "EPC1",
                               scale_type = "categorical")
    testthat::expect_type(p, "list")
  })

  testthat::test_that("do_GeyserPlot: PASS - default parameters = color.by", {



    p <- SCpubr::do_GeyserPlot(sample = sample,
                               features = "EPC1",
                               scale_type = "categorical",
                               color.by = "orig.ident")
    testthat::expect_type(p, "list")

    p <- SCpubr::do_GeyserPlot(sample = sample,
                               features = "EPC1",
                               scale_type = "continuous",
                               color.by = "PC_1")
    testthat::expect_type(p, "list")

    p <- SCpubr::do_GeyserPlot(sample = sample,
                               features = "PC_2",
                               scale_type = "categorical",
                               color.by = "orig.ident")
    testthat::expect_type(p, "list")

    p <- SCpubr::do_GeyserPlot(sample = sample,
                               features = "PC_2",
                               scale_type = "continuous",
                               color.by = "nCount_RNA")
    testthat::expect_type(p, "list")

    p <- SCpubr::do_GeyserPlot(sample = sample,
                               features = "nCount_RNA",
                               scale_type = "categorical",
                               color.by = "orig.ident")
    testthat::expect_type(p, "list")

    p <- SCpubr::do_GeyserPlot(sample = sample,
                               features = "nCount_RNA",
                               scale_type = "continuous",
                               color.by = "EPC1")
    testthat::expect_type(p, "list")
  })


  testthat::test_that("do_GeyserPlot: PASS - categorical colors.use", {



    p <- SCpubr::do_GeyserPlot(sample = sample,
                               features = "EPC1",
                               scale_type = "categorical",
                               group.by = "orig.ident",
                               colors.use = c("Cell" = "green"))
    testthat::expect_type(p, "list")

    p <- SCpubr::do_GeyserPlot(sample = sample,
                               features = "EPC1",
                               scale_type = "categorical",
                               group.by = "orig.ident",
                               colors.use = NULL)
    testthat::expect_type(p, "list")
  })

  testthat::test_that("do_GeyserPlot: PASS - order by mean", {



    p <- SCpubr::do_GeyserPlot(sample = sample,
                               features = "EPC1",
                               scale_type = "continuous",
                               order_by_mean = TRUE)
    testthat::expect_type(p, "list")

    p <- SCpubr::do_GeyserPlot(sample = sample,
                               features = "EPC1",
                               scale_type = "continuous",
                               order_by_mean = FALSE)
    testthat::expect_type(p, "list")
  })

  testthat::test_that("do_GeyserPlot: PASS - plot cell borders", {



    p <- SCpubr::do_GeyserPlot(sample = sample,
                               features = "EPC1",
                               scale_type = "continuous",
                               plot_cell_borders = TRUE)
    testthat::expect_type(p, "list")

    p <- SCpubr::do_GeyserPlot(sample = sample,
                               features = "EPC1",
                               scale_type = "continuous",
                               plot_cell_borders = FALSE)
    testthat::expect_type(p, "list")
  })

  testthat::test_that("do_GeyserPlot: PASS - split.by", {



    p <- SCpubr::do_GeyserPlot(sample = sample,
                               features = "EPC1",
                               scale_type = "continuous",
                               split.by = "seurat_clusters",
                               plot_cell_borders = TRUE)
    testthat::expect_type(p, "list")

    p <- SCpubr::do_GeyserPlot(sample = sample,
                               features = "EPC1",
                               scale_type = "continuous",
                               split.by = "seurat_clusters",
                               plot_cell_borders = FALSE)
    testthat::expect_type(p, "list")
  })


  testthat::test_that("do_BarPlot: FAIL - wrong paramters", {



    testthat::expect_error(SCpubr::do_GeyserPlot(sample = sample,
                                                 features = "EPC1",
                                                 scale_type = "wrong"))
    testthat::expect_error(SCpubr::do_GeyserPlot(sample = sample,
                                                 features = "EPC1",
                                                 split.by = "wrong"))
    testthat::expect_error(SCpubr::do_GeyserPlot(sample = sample,
                                                 features = "EPC1",
                                                 group.by = "wrong"))
    testthat::expect_error(SCpubr::do_GeyserPlot(sample = sample,
                                                 features = "EPC1",
                                                 jitter = 1))
    testthat::expect_error(SCpubr::do_GeyserPlot(sample = sample,
                                                 features = "EPC1",
                                                 scale_type = "categorical",
                                                 color.by = "EPC1"))
  })

  testthat::test_that("do_GeyserPlot: PASS - show legend", {



    p <- SCpubr::do_GeyserPlot(sample = sample,
                               features = "EPC1")
    testthat::expect_type(p, "list")

    p <- SCpubr::do_GeyserPlot(sample = sample,
                               features = "EPC1")
    testthat::expect_type(p, "list")
  })

  testthat::test_that("do_GeyserPlot: PASS - several features", {



    p <- SCpubr::do_GeyserPlot(sample = sample,
                               features = c("EPC1", "PC_1"))
    testthat::expect_type(p, "list")
    testthat::expect_length(p, 2)
  })

  testthat::test_that("do_GeyserPlot: PASS - color.by factor", {



    sample$seurat_clusters_factor <- as.factor(sample$seurat_clusters)
    sample$seurat_clusters_character <- as.character(sample$seurat_clusters)
    p <- SCpubr::do_GeyserPlot(sample = sample,
                               features = "EPC1",
                               color.by = "seurat_clusters_factor")
    testthat::expect_type(p, "list")

    p <- SCpubr::do_GeyserPlot(sample = sample,
                               features = "EPC1",
                               color.by = "seurat_clusters_character")
    testthat::expect_type(p, "list")

    p <- SCpubr::do_GeyserPlot(sample = sample,
                               features = "EPC1",
                               scale_type = "categorical",
                               color.by = "seurat_clusters_character")
    testthat::expect_type(p, "list")

    p <- SCpubr::do_GeyserPlot(sample = sample,
                               features = "EPC1",
                               scale_type = "categorical",
                               color.by = "seurat_clusters_factor")
    testthat::expect_type(p, "list")
  })

}



