if(isFALSE(dep_check[["do_NebulosaPlot"]])){
  testthat::test_that("do_NebulosaPlot: PASS - single feature", {



    p <- SCpubr::do_NebulosaPlot(sample = sample,
                                 features = c("EPC1"))
    testthat::expect_type(p, "list")
  })

  testthat::test_that("do_NebulosaPlot: PASS - cell_borders", {



    p <- SCpubr::do_NebulosaPlot(sample = sample, features = "EPC1", plot_cell_borders = TRUE)
    testthat::expect_type(p, "list")
    p <- suppressWarnings({SCpubr::do_NebulosaPlot(sample = sample, features = c("EPC1", "PC_1"), plot_cell_borders = TRUE)})
    testthat::expect_type(p, "list")
  })

  testthat::test_that("do_NebulosaPlot: PASS - single feature legend normal", {



    p <- SCpubr::do_NebulosaPlot(sample = sample,
                                 features = c("EPC1"),
                                 legend.type = "normal")
    testthat::expect_type(p, "list")
  })

  testthat::test_that("do_NebulosaPlot: PASS - single feature legend colorbar", {



    p <- SCpubr::do_NebulosaPlot(sample = sample,
                                 features = c("EPC1"),
                                 legend.type = "colorbar")
    testthat::expect_type(p, "list")
  })

  testthat::test_that("do_NebulosaPlot: PASS - single feature legend colorsteps", {



    p <- SCpubr::do_NebulosaPlot(sample = sample,
                                 features = c("EPC1"),
                                 legend.type = "colorsteps")
    testthat::expect_type(p, "list")
  })

  testthat::test_that("do_NebulosaPlot: FAIL - wrong legend type ", {



    testthat::expect_error(SCpubr::do_NebulosaPlot(sample = sample,
                                                   features = c("EPC1"),
                                                   legend.type = "wrong"))
  })

  testthat::test_that("do_NebulosaPlot: FAIL - wrong legend position ", {



    testthat::expect_error(SCpubr::do_NebulosaPlot(sample = sample,
                                                   features = c("EPC1"),
                                                   legend.position = "wrong"))
  })

  testthat::test_that("do_NebulosaPlot: FAIL - wrong font.type", {



    testthat::expect_error(SCpubr::do_NebulosaPlot(sample = sample,
                                                   features = c("EPC1"),
                                                   font.type = "wrong"))
  })

  testthat::test_that("do_NebulosaPlot: PASS - single feature distinct dims", {



    p <- SCpubr::do_NebulosaPlot(sample = sample,
                                 features = c("EPC1"),
                                 dims = c(2, 1))
    testthat::expect_type(p, "list")
  })

  testthat::test_that("do_FeaturePlot: PASS - diffusion", {



    test <- sample@reductions$umap[[]]
    colnames(test) <- c("DC_1", "DC_2")
    obj <- Seurat::CreateDimReducObject(test, assay = "SCT", key = "DC_")
    sample@reductions$diffusion <- obj
    p <- suppressWarnings(SCpubr::do_NebulosaPlot(sample,
                                                  features = c("PC_1"),
                                                  reduction = "diffusion"))
    testthat::expect_type(p, "list")
  })

  testthat::test_that("do_NebulosaPlot: PASS - several", {



    p <- SCpubr::do_NebulosaPlot(sample = sample,
                                 features = c("EPC1", "LTV1"))
    testthat::expect_type(p, "list")
  })

  testthat::test_that("do_NebulosaPlot: PASS - several, joint", {



    p <- SCpubr::do_NebulosaPlot(sample = sample,
                                 features = c("EPC1", "LTV1"),
                                 joint = TRUE)
    testthat::expect_type(p, "list")
  })

  testthat::test_that("do_NebulosaPlot: PASS - several, joint only joint", {



    p <- SCpubr::do_NebulosaPlot(sample = sample,
                                 features = c("EPC1", "LTV1"),
                                 joint = TRUE,
                                 return_only_joint = TRUE)
    testthat::expect_type(p, "list")
  })

  testthat::test_that("do_NebulosaPlot: PASS - title", {



    p <- SCpubr::do_NebulosaPlot(sample = sample,
                                 features = c("EPC1", "LTV1"),
                                 joint = TRUE,
                                 return_only_joint = TRUE,
                                 plot.title = "Title")
    testthat::expect_type(p, "list")
  })

  testthat::test_that("do_NebulosaPlot: PASS - subtitle", {



    p <- SCpubr::do_NebulosaPlot(sample = sample,
                                 features = c("EPC1", "LTV1"),
                                 joint = TRUE,
                                 return_only_joint = TRUE,
                                 plot.subtitle = "Subtitle")
    testthat::expect_type(p, "list")
  })

  testthat::test_that("do_NebulosaPlot: PASS - caption", {



    p <- SCpubr::do_NebulosaPlot(sample = sample,
                                 features = c("EPC1", "LTV1"),
                                 joint = TRUE,
                                 return_only_joint = TRUE,
                                 plot.caption = "Caption")
    testthat::expect_type(p, "list")
  })


  testthat::test_that("do_NebulosaPlot: PASS - color map", {



    p <- SCpubr::do_NebulosaPlot(sample = sample,
                                 features = c("EPC1"),
                                 viridis_color_map = "F")
    testthat::expect_type(p, "list")
  })

  testthat::test_that("do_NebulosaPlot: PASS - legend top", {



    p <- SCpubr::do_NebulosaPlot(sample = sample,
                                 features = c("EPC1"),
                                 legend.position = "left")
    testthat::expect_type(p, "list")
  })

  testthat::test_that("do_NebulosaPlot: PASS - legend top", {



    p <- SCpubr::do_NebulosaPlot(sample = sample,
                                 features = c("EPC1"),
                                 legend.position = "top")
    testthat::expect_type(p, "list")
  })

  testthat::test_that("do_NebulosaPlot: WARNING - features as list", {



    testthat::expect_warning(SCpubr::do_NebulosaPlot(sample = sample,
                                                     features = list("EPC1"),
                                                     viridis_color_map = "F"))
  })


  testthat::test_that("do_NebulosaPlot: PASS - no legend", {



    p <- SCpubr::do_NebulosaPlot(sample = sample,
                                 features = c("EPC1"),
                                 legend.position = "none")
    testthat::expect_type(p, "list")
  })

  testthat::test_that("do_NebulosaPlot: PASS - patchwork title, subtitle and caption", {



    p <- SCpubr::do_NebulosaPlot(sample = sample,
                                 features = c("EPC1", "LTV1"),
                                 plot.title = "A",
                                 plot.subtitle = "B",
                                 plot.caption = "C")
    testthat::expect_type(p, "list")
  })

  testthat::test_that("do_NebulosaPlot: PASS - plot axis", {



    p <- SCpubr::do_NebulosaPlot(sample = sample, plot.axes = TRUE, features = "nCount_RNA")
    testthat::expect_type(p, "list")

    p <- SCpubr::do_NebulosaPlot(sample = sample, reduction = "pca", plot.axes = TRUE, features = "nCount_RNA")
    testthat::expect_type(p, "list")

    p <- SCpubr::do_NebulosaPlot(sample = sample, dims = c(2, 1), plot.axes = TRUE, features = "nCount_RNA")
    testthat::expect_type(p, "list")

    sample@reductions$diffusion <- sample@reductions$umap
    p <- SCpubr::do_NebulosaPlot(sample = sample,
                                 reduction = "diffusion",
                                 plot.axes = TRUE,
                                 features = "nCount_RNA")
    testthat::expect_type(p, "list")
  })
}

