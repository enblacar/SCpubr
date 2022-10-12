
testthat::test_that("do_FeaturePlot: PASS - single feature", {


  p <- SCpubr::do_FeaturePlot(sample = sample,
                              features = "nCount_RNA")
  testthat::expect_type(p, "list")
})

testthat::test_that("do_FeaturePlot: PASS - legend.title", {


  p <- SCpubr::do_FeaturePlot(sample = sample,
                              features = "nCount_RNA",
                              legend.title = "pepe")
  testthat::expect_type(p, "list")

  p <- SCpubr::do_FeaturePlot(sample = sample,
                              features = "nCount_RNA",
                              legend.title = "pepe",
                              split.by = "seurat_clusters")
  testthat::expect_type(p, "list")

  p <- SCpubr::do_FeaturePlot(sample = sample,
                              features = "nCount_RNA",
                              split.by = "seurat_clusters")
  testthat::expect_type(p, "list")
})

testthat::test_that("do_FeaturePlot: PASS - cell_borders", {


  p <- SCpubr::do_FeaturePlot(sample = sample, features = "EPC1", plot_cell_borders = TRUE)
  testthat::expect_type(p, "list")
  p <- SCpubr::do_FeaturePlot(sample = sample, features = "EPC1", plot_cell_borders = TRUE, raster = TRUE, pt.size = 1)
  testthat::expect_type(p, "list")
  p <- SCpubr::do_FeaturePlot(sample = sample, features = "EPC1", plot_cell_borders = TRUE, idents.highlight = "1")
  testthat::expect_type(p, "list")
  p <- SCpubr::do_FeaturePlot(sample = sample, features = "EPC1", plot_cell_borders = TRUE, raster = TRUE, idents.highlight = "1", pt.size = 1)
  testthat::expect_type(p, "list")
  p <- SCpubr::do_FeaturePlot(sample = sample, features = "EPC1", plot_cell_borders = TRUE, split.by = "seurat_clusters")
  testthat::expect_type(p, "list")
  p <- SCpubr::do_FeaturePlot(sample = sample, features = "EPC1", plot_cell_borders = TRUE, raster = TRUE, split.by = "seurat_clusters", pt.size = 1)
  testthat::expect_type(p, "list")
})

testthat::test_that("do_FeaturePlot: PASS - enforce_symmetry", {


  p <- SCpubr::do_FeaturePlot(sample = sample, features = "EPC1", enforce_symmetry = TRUE)
  testthat::expect_type(p, "list")

  p <- SCpubr::do_FeaturePlot(sample = sample, features = c("EPC1", "nCount_RNA"), enforce_symmetry = TRUE)
  testthat::expect_type(p, "list")

  p <- SCpubr::do_FeaturePlot(sample = sample, features = "EPC1", enforce_symmetry = TRUE, idents.highlight = c("1", "3"))
  testthat::expect_type(p, "list")

  p <- SCpubr::do_FeaturePlot(sample = sample, features = "EPC1", enforce_symmetry = TRUE, split.by = "seurat_clusters")
  testthat::expect_type(p, "list")
})

testthat::test_that("do_FeaturePlot: PASS - multiple features", {


  p <- SCpubr::do_FeaturePlot(sample = sample,
                              features = c("nCount_RNA", "nFeature_RNA"))
  testthat::expect_type(p, "list")
})

testthat::test_that("do_FeaturePlot: PASS - title", {


  p <- SCpubr::do_FeaturePlot(sample = sample,
                              features = "nCount_RNA",
                              plot.title = "Feature Plot")
  testthat::expect_type(p, "list")
})

testthat::test_that("do_FeaturePlot: PASS - subtitle", {


  p <- SCpubr::do_FeaturePlot(sample = sample,
                              features = "nCount_RNA",
                              plot.subtitle = "Feature Plot")
  testthat::expect_type(p, "list")
})

testthat::test_that("do_FeaturePlot: PASS - caption", {


  p <- SCpubr::do_FeaturePlot(sample = sample,
                              features = "nCount_RNA",
                              plot.caption = "Feature Plot")
  testthat::expect_type(p, "list")
})

testthat::test_that("do_FeaturePlot: PASS - individual titles", {


  p <- SCpubr::do_FeaturePlot(sample = sample,
                              features = c("nCount_RNA", "nFeature_RNA"),
                              individual.titles = c("A", NA))
  testthat::expect_type(p, "list")
})

testthat::test_that("do_FeaturePlot: PASS - individual subtitles ", {


  p <- SCpubr::do_FeaturePlot(sample = sample,
                              features = c("nCount_RNA", "nFeature_RNA"),
                              individual.subtitles = c("A", NA))
  testthat::expect_type(p, "list")
})

testthat::test_that("do_FeaturePlot: PASS - individual captions", {


  p <- SCpubr::do_FeaturePlot(sample = sample,
                              features = c("nCount_RNA", "nFeature_RNA"),
                              individual.captions = c("A", NA))
  testthat::expect_type(p, "list")
})

testthat::test_that("do_FeaturePlot: PASS - dims", {


  p <- SCpubr::do_FeaturePlot(sample = sample,
                              features = "nCount_RNA",
                              dims = c(2, 1))
  testthat::expect_type(p, "list")
})

testthat::test_that("do_FeaturePlot: PASS - subset of cells", {


  cells.plot <- colnames(sample[, !(sample$seurat_clusters %in% c("2", "5", "8"))])
  p <- SCpubr::do_FeaturePlot(sample,
                              features = c("EPC1"),
                              cells.highlight = cells.plot)
  testthat::expect_type(p, "list")
})

testthat::test_that("do_FeaturePlot: PASS - subset of identities", {


  p <- SCpubr::do_FeaturePlot(sample,
                              features = c("EPC1"),
                              idents.highlight = c("1", "2"))
  testthat::expect_type(p, "list")
})

testthat::test_that("do_FeaturePlot: PASS - subset of cells and identities", {


  cells.plot <- colnames(sample[, !(sample$seurat_clusters %in% c("2", "5", "8"))])
  p <- SCpubr::do_FeaturePlot(sample,
                              features = c("EPC1"),
                              cells.highlight = cells.plot,
                              idents.highlight = c("1", "2"))
  testthat::expect_type(p, "list")
})

testthat::test_that("do_FeaturePlot: PASS - split.by", {


  p <- SCpubr::do_FeaturePlot(sample,
                              features = c("EPC1"),
                              split.by = "seurat_clusters")
  testthat::expect_type(p, "list")
})

testthat::test_that("do_FeaturePlot: PASS - split.by and split.by.idents", {


  p <- SCpubr::do_FeaturePlot(sample,
                              features = c("EPC1"),
                              split.by = "seurat_clusters",
                              split.by.idents = c("1", "2"))
  testthat::expect_type(p, "list")
})

testthat::test_that("do_FeaturePlot: PASS - split.by and split.by.idents multiple features", {


  p <- SCpubr::do_FeaturePlot(sample,
                              features = c("EPC1", "nCount_RNA"),
                              split.by = "seurat_clusters",
                              split.by.idents = c("1", "2"))
  testthat::expect_type(p, "list")
})

testthat::test_that("do_FeaturePlot: PASS - modify color maps", {


  p <- SCpubr::do_FeaturePlot(sample,
                              features = c("nCount_RNA"),
                              viridis_color_map = "F")
  testthat::expect_type(p, "list")
})

testthat::test_that("do_FeaturePlot: WARNING - features as a list", {


  testthat::expect_warning(SCpubr::do_FeaturePlot(sample,
                                                  features = list("A" = c("nCount_RNA"))))
})

testthat::test_that("do_FeaturePlot: FAIL - individual titles, subtitles or captions do not match with number of features", {


  testthat::expect_error(SCpubr::do_FeaturePlot(sample,
                                                features = c("nCount_RNA", "EPC1"),
                                                individual.titles = "A"))
  testthat::expect_error(SCpubr::do_FeaturePlot(sample,
                                                features = c("nCount_RNA", "EPC1"),
                                                individual.subtitles = "A"))
  testthat::expect_error(SCpubr::do_FeaturePlot(sample,
                                                features = c("nCount_RNA", "EPC1"),
                                                individual.captions = "A"))
})

testthat::test_that("do_FeaturePlot: PASS - legend position = right", {


  p <- SCpubr::do_FeaturePlot(sample,
                              features = c("nCount_RNA"),
                              legend.position = "right")
  testthat::expect_type(p, "list")
})

testthat::test_that("do_FeaturePlot: WARNING - raster and small point size", {


  testthat::expect_warning(SCpubr::do_FeaturePlot(sample,
                                                  features = c("nCount_RNA"),
                                                  raster = TRUE,
                                                  pt.size = 0.5))
})

testthat::test_that("do_FeaturePlot: PASS - ussing diffusion reduction", {


  test <- sample@reductions$umap[[]]
  colnames(test) <- c("DC_1", "DC_2")
  obj <- Seurat::CreateDimReducObject(test, assay = "SCT", key = "DC_")
  sample@reductions$diffusion <- obj
  p <- SCpubr::do_FeaturePlot(sample,
                              features = c("nCount_RNA"),
                              reduction = "diffusion")
  testthat::expect_type(p, "list")
})

testthat::test_that("do_FeaturePlot: PASS - duplicated split.by.idents", {


  testthat::expect_message(SCpubr::do_FeaturePlot(sample,
                                                  features = c("nCount_RNA"),
                                                  split.by = "seurat_clusters",
                                                  split.by.idents = c("2", "2")))
})


testthat::test_that("do_FeaturePlot: PASS - plotting a Dimensional reduction component", {


  p <- SCpubr::do_FeaturePlot(sample,
                              features = c("PC_1"))
  testthat::expect_type(p, "list")
})

testthat::test_that("do_FeaturePlot: PASS - split.by factor", {


  sample$factor_seurat_clusters <- factor(sample$seurat_clusters, levels = c("2", "0", "1", "3", "4", "5", "6", "7", "8"))
  p <- SCpubr::do_FeaturePlot(sample,
                              features = c("PC_1"),
                              split.by = "factor_seurat_clusters")
  testthat::expect_type(p, "list")
})

testthat::test_that("do_FeaturePlot: PASS - split.by and plot.title", {


  p <- SCpubr::do_FeaturePlot(sample,
                              features = c("PC_1"),
                              split.by = "seurat_clusters",
                              plot.title = "Title",
                              plot.subtitle = "Subtitle",
                              plot.caption = "Caption")
  testthat::expect_type(p, "list")
})

testthat::test_that("do_FeaturePlot: PASS - split.by and pca", {


  p <- SCpubr::do_FeaturePlot(sample,
                              features = c("PC_1"),
                              split.by = "seurat_clusters",
                              reduction = "pca")
  testthat::expect_type(p, "list")
})

testthat::test_that("do_FeaturePlot: PASS - split.by and diffusion", {


  test <- sample@reductions$umap[[]]
  colnames(test) <- c("DC_1", "DC_2")
  obj <- Seurat::CreateDimReducObject(test, assay = "SCT", key = "DC_")
  sample@reductions$diffusion <- obj
  p <- SCpubr::do_FeaturePlot(sample,
                              features = c("PC_1"),
                              split.by = "seurat_clusters",
                              reduction = "diffusion")
  testthat::expect_type(p, "list")
})

testthat::test_that("do_FeaturePlot: PASS - remove legend", {


  p <- SCpubr::do_FeaturePlot(sample,
                              features = c("PC_1"),
                              legend.position = "none")
  testthat::expect_type(p, "list")
})

testthat::test_that("do_FeaturePlot: PASS - normal legend", {


  p <- SCpubr::do_FeaturePlot(sample,
                              features = c("PC_1"),
                              legend.type = "normal")
  testthat::expect_type(p, "list")
})

testthat::test_that("do_FeaturePlot: PASS - colorbar legend", {


  p <- SCpubr::do_FeaturePlot(sample,
                              features = c("PC_1"),
                              legend.type = "colorbar")
  testthat::expect_type(p, "list")
})

testthat::test_that("do_FeaturePlot: PASS - colorsteps legend", {


  p <- SCpubr::do_FeaturePlot(sample,
                              features = c("PC_1"),
                              legend.type = "colorsteps")
  testthat::expect_type(p, "list")
})

testthat::test_that("do_FeaturePlot: PASS - normal legend - split.by", {


  p <- SCpubr::do_FeaturePlot(sample,
                              features = c("PC_1"),
                              legend.type = "normal",
                              split.by = "seurat_clusters")
  testthat::expect_type(p, "list")
})

testthat::test_that("do_FeaturePlot: PASS - colorbar legend - split.by", {


  p <- SCpubr::do_FeaturePlot(sample,
                              features = c("PC_1"),
                              legend.type = "colorbar",
                              split.by = "seurat_clusters")
  testthat::expect_type(p, "list")
})

testthat::test_that("do_FeaturePlot: PASS - colorsteps legend - split.by", {


  p <- SCpubr::do_FeaturePlot(sample,
                              features = c("PC_1"),
                              legend.type = "colorsteps",
                              split.by = "seurat_clusters")
  testthat::expect_type(p, "list")
})

testthat::test_that("do_FeaturePlot: FAIL - wrong legend type", {


  testthat::expect_error(SCpubr::do_FeaturePlot(sample,
                                                features = c("PC_1"),
                                                legend.type = "wrong"))
})

testthat::test_that("do_FeaturePlot: FAIL - wrong legend position", {


  testthat::expect_error(SCpubr::do_FeaturePlot(sample,
                                                features = c("PC_1"),
                                                legend.position = "wrong"))
})
testthat::test_that("do_FeaturePlot: FAIL - wrong font.type", {


  testthat::expect_error(SCpubr::do_FeaturePlot(sample,
                                                features = c("PC_1"),
                                                font.type = "wrong"))
})

testthat::test_that("do_FeaturePlot: PASS - plot axis", {


  p <- SCpubr::do_FeaturePlot(sample = sample, plot.axes = TRUE, features = "nCount_RNA")
  testthat::expect_type(p, "list")

  p <- SCpubr::do_FeaturePlot(sample = sample, reduction = "pca", plot.axes = TRUE, features = "nCount_RNA")
  testthat::expect_type(p, "list")

  p <- SCpubr::do_FeaturePlot(sample = sample, dims = c(2, 1), plot.axes = TRUE, features = "nCount_RNA")
  testthat::expect_type(p, "list")

  sample@reductions$diffusion <- sample@reductions$umap
  p <- SCpubr::do_FeaturePlot(sample = sample,
                          reduction = "diffusion",
                          plot.axes = TRUE,
                          features = "nCount_RNA")
  testthat::expect_type(p, "list")
})
