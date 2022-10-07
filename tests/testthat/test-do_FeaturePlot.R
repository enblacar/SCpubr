testthat::test_that("do_FeaturePlot: PASS - single feature", {
  sample <- SCpubr:::test_list$sample

  p <- SCpubr::do_FeaturePlot(sample = sample,
                              features = "nCount_RNA")
  testthat::expect_type(p, "list")
})

testthat::test_that("do_FeaturePlot: PASS - legend.title", {
  sample <- SCpubr:::test_list$sample

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
  sample <- SCpubr:::test_list$sample

  p <- SCpubr::do_FeaturePlot(sample = sample, features = "CD14", plot_cell_borders = T)
  testthat::expect_type(p, "list")
  p <- SCpubr::do_FeaturePlot(sample = sample, features = "CD14", plot_cell_borders = T, raster = T, pt.size = 1)
  testthat::expect_type(p, "list")
  p <- SCpubr::do_FeaturePlot(sample = sample, features = "CD14", plot_cell_borders = T, idents.highlight = "1")
  testthat::expect_type(p, "list")
  p <- SCpubr::do_FeaturePlot(sample = sample, features = "CD14", plot_cell_borders = T, raster = T, idents.highlight = "1", pt.size = 1)
  testthat::expect_type(p, "list")
  p <- SCpubr::do_FeaturePlot(sample = sample, features = "CD14", plot_cell_borders = T, split.by = "seurat_clusters")
  testthat::expect_type(p, "list")
  p <- SCpubr::do_FeaturePlot(sample = sample, features = "CD14", plot_cell_borders = T, raster = T, split.by = "seurat_clusters", pt.size = 1)
  testthat::expect_type(p, "list")
})

testthat::test_that("do_FeaturePlot: PASS - enforce_symmetry", {
  sample <- SCpubr:::test_list$sample

  p <- SCpubr::do_FeaturePlot(sample = sample, features = "CD14", enforce_symmetry = T)
  testthat::expect_type(p, "list")

  p <- SCpubr::do_FeaturePlot(sample = sample, features = c("CD14", "nCount_RNA"), enforce_symmetry = T)
  testthat::expect_type(p, "list")

  p <- SCpubr::do_FeaturePlot(sample = sample, features = "CD14", enforce_symmetry = T, idents.highlight = c("1", "3"))
  testthat::expect_type(p, "list")

  p <- SCpubr::do_FeaturePlot(sample = sample, features = "CD14", enforce_symmetry = T, split.by = "seurat_clusters")
  testthat::expect_type(p, "list")
})

testthat::test_that("do_FeaturePlot: PASS - multiple features", {
  sample <- SCpubr:::test_list$sample

  p <- SCpubr::do_FeaturePlot(sample = sample,
                              features = c("nCount_RNA", "nFeature_RNA"))
  testthat::expect_type(p, "list")
})

testthat::test_that("do_FeaturePlot: PASS - title", {
  sample <- SCpubr:::test_list$sample

  p <- SCpubr::do_FeaturePlot(sample = sample,
                              features = "nCount_RNA",
                              plot.title = "Feature Plot")
  testthat::expect_type(p, "list")
})

testthat::test_that("do_FeaturePlot: PASS - subtitle", {
  sample <- SCpubr:::test_list$sample

  p <- SCpubr::do_FeaturePlot(sample = sample,
                              features = "nCount_RNA",
                              plot.subtitle = "Feature Plot")
  testthat::expect_type(p, "list")
})

testthat::test_that("do_FeaturePlot: PASS - caption", {
  sample <- SCpubr:::test_list$sample

  p <- SCpubr::do_FeaturePlot(sample = sample,
                              features = "nCount_RNA",
                              plot.caption = "Feature Plot")
  testthat::expect_type(p, "list")
})

testthat::test_that("do_FeaturePlot: PASS - individual titles", {
  sample <- SCpubr:::test_list$sample

  p <- SCpubr::do_FeaturePlot(sample = sample,
                              features = c("nCount_RNA", "nFeature_RNA"),
                              individual.titles = c("A", NA))
  testthat::expect_type(p, "list")
})

testthat::test_that("do_FeaturePlot: PASS - individual subtitles ", {
  sample <- SCpubr:::test_list$sample

  p <- SCpubr::do_FeaturePlot(sample = sample,
                              features = c("nCount_RNA", "nFeature_RNA"),
                              individual.subtitles = c("A", NA))
  testthat::expect_type(p, "list")
})

testthat::test_that("do_FeaturePlot: PASS - individual captions", {
  sample <- SCpubr:::test_list$sample

  p <- SCpubr::do_FeaturePlot(sample = sample,
                              features = c("nCount_RNA", "nFeature_RNA"),
                              individual.captions = c("A", NA))
  testthat::expect_type(p, "list")
})

testthat::test_that("do_FeaturePlot: PASS - dims", {
  sample <- SCpubr:::test_list$sample

  p <- SCpubr::do_FeaturePlot(sample = sample,
                              features = "nCount_RNA",
                              dims = c(2, 1))
  testthat::expect_type(p, "list")
})

testthat::test_that("do_FeaturePlot: PASS - subset of cells", {
  sample <- SCpubr:::test_list$sample

  cells.plot <- colnames(sample[, !(sample$seurat_clusters %in% c("2", "5", "8"))])
  p <- SCpubr::do_FeaturePlot(sample,
                              features = c("CD14"),
                              cells.highlight = cells.plot)
  testthat::expect_type(p, "list")
})

testthat::test_that("do_FeaturePlot: PASS - subset of identities", {
  sample <- SCpubr:::test_list$sample

  p <- SCpubr::do_FeaturePlot(sample,
                              features = c("CD14"),
                              idents.highlight = c("1", "2"))
  testthat::expect_type(p, "list")
})

testthat::test_that("do_FeaturePlot: PASS - subset of cells and identities", {
  sample <- SCpubr:::test_list$sample

  cells.plot <- colnames(sample[, !(sample$seurat_clusters %in% c("2", "5", "8"))])
  p <- SCpubr::do_FeaturePlot(sample,
                              features = c("CD14"),
                              cells.highlight = cells.plot,
                              idents.highlight = c("1", "2"))
  testthat::expect_type(p, "list")
})

testthat::test_that("do_FeaturePlot: PASS - split.by", {
  sample <- SCpubr:::test_list$sample

  p <- SCpubr::do_FeaturePlot(sample,
                              features = c("CD14"),
                              split.by = "seurat_clusters")
  testthat::expect_type(p, "list")
})

testthat::test_that("do_FeaturePlot: PASS - split.by and split.by.idents", {
  sample <- SCpubr:::test_list$sample

  p <- SCpubr::do_FeaturePlot(sample,
                              features = c("CD14"),
                              split.by = "seurat_clusters",
                              split.by.idents = c("1", "2"))
  testthat::expect_type(p, "list")
})

testthat::test_that("do_FeaturePlot: PASS - split.by and split.by.idents multiple features", {
  sample <- SCpubr:::test_list$sample

  p <- SCpubr::do_FeaturePlot(sample,
                              features = c("CD14", "nCount_RNA"),
                              split.by = "seurat_clusters",
                              split.by.idents = c("1", "2"))
  testthat::expect_type(p, "list")
})

testthat::test_that("do_FeaturePlot: PASS - modify color maps", {
  sample <- SCpubr:::test_list$sample

  p <- SCpubr::do_FeaturePlot(sample,
                              features = c("nCount_RNA"),
                              viridis_color_map = "F")
  testthat::expect_type(p, "list")
})

testthat::test_that("do_FeaturePlot: WARNING - features as a list", {
  sample <- SCpubr:::test_list$sample

  testthat::expect_warning(SCpubr::do_FeaturePlot(sample,
                                                  features = list("A" = c("nCount_RNA"))))
})

testthat::test_that("do_FeaturePlot: FAIL - individual titles, subtitles or captions do not match with number of features", {
  sample <- SCpubr:::test_list$sample

  testthat::expect_error(SCpubr::do_FeaturePlot(sample,
                                                features = c("nCount_RNA", "CD14"),
                                                individual.titles = "A"))
  testthat::expect_error(SCpubr::do_FeaturePlot(sample,
                                                features = c("nCount_RNA", "CD14"),
                                                individual.subtitles = "A"))
  testthat::expect_error(SCpubr::do_FeaturePlot(sample,
                                                features = c("nCount_RNA", "CD14"),
                                                individual.captions = "A"))
})

testthat::test_that("do_FeaturePlot: PASS - legend position = right", {
  sample <- SCpubr:::test_list$sample

  p <- SCpubr::do_FeaturePlot(sample,
                              features = c("nCount_RNA"),
                              legend.position = "right")
  testthat::expect_type(p, "list")
})

testthat::test_that("do_FeaturePlot: WARNING - raster and small point size", {
  sample <- SCpubr:::test_list$sample

  testthat::expect_warning(SCpubr::do_FeaturePlot(sample,
                                                  features = c("nCount_RNA"),
                                                  raster = T,
                                                  pt.size = 0.5))
})

testthat::test_that("do_FeaturePlot: PASS - ussing diffusion reduction", {
  sample <- SCpubr:::test_list$sample

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
  sample <- SCpubr:::test_list$sample

  testthat::expect_message(SCpubr::do_FeaturePlot(sample,
                                                  features = c("nCount_RNA"),
                                                  split.by = "seurat_clusters",
                                                  split.by.idents = c("2", "2")))
})


testthat::test_that("do_FeaturePlot: PASS - plotting a Dimensional reduction component", {
  sample <- SCpubr:::test_list$sample

  p <- SCpubr::do_FeaturePlot(sample,
                              features = c("PC_1"))
  testthat::expect_type(p, "list")
})

testthat::test_that("do_FeaturePlot: PASS - split.by factor", {
  sample <- SCpubr:::test_list$sample

  sample$factor_seurat_clusters <- factor(sample$seurat_clusters, levels = c("2", "0", "1", "3","4", "5", "6", "7", "8"))
  p <- SCpubr::do_FeaturePlot(sample,
                              features = c("PC_1"),
                              split.by = "factor_seurat_clusters")
  testthat::expect_type(p, "list")
})

testthat::test_that("do_FeaturePlot: PASS - split.by and plot.title", {
  sample <- SCpubr:::test_list$sample

  p <- SCpubr::do_FeaturePlot(sample,
                              features = c("PC_1"),
                              split.by = "seurat_clusters",
                              plot.title = "Title",
                              plot.subtitle = "Subtitle",
                              plot.caption = "Caption")
  testthat::expect_type(p, "list")
})

testthat::test_that("do_FeaturePlot: PASS - split.by and pca", {
  sample <- SCpubr:::test_list$sample

  p <- SCpubr::do_FeaturePlot(sample,
                              features = c("PC_1"),
                              split.by = "seurat_clusters",
                              reduction = "pca")
  testthat::expect_type(p, "list")
})

testthat::test_that("do_FeaturePlot: PASS - split.by and diffusion", {
  sample <- SCpubr:::test_list$sample

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
  sample <- SCpubr:::test_list$sample

  p <- SCpubr::do_FeaturePlot(sample,
                              features = c("PC_1"),
                              legend.position = "none")
  testthat::expect_type(p, "list")
})

testthat::test_that("do_FeaturePlot: PASS - normal legend", {
  sample <- SCpubr:::test_list$sample

  p <- SCpubr::do_FeaturePlot(sample,
                              features = c("PC_1"),
                              legend.type = "normal")
  testthat::expect_type(p, "list")
})

testthat::test_that("do_FeaturePlot: PASS - colorbar legend", {
  sample <- SCpubr:::test_list$sample

  p <- SCpubr::do_FeaturePlot(sample,
                              features = c("PC_1"),
                              legend.type = "colorbar")
  testthat::expect_type(p, "list")
})

testthat::test_that("do_FeaturePlot: PASS - colorsteps legend", {
  sample <- SCpubr:::test_list$sample

  p <- SCpubr::do_FeaturePlot(sample,
                              features = c("PC_1"),
                              legend.type = "colorsteps")
  testthat::expect_type(p, "list")
})

testthat::test_that("do_FeaturePlot: PASS - normal legend - split.by", {
  sample <- SCpubr:::test_list$sample

  p <- SCpubr::do_FeaturePlot(sample,
                              features = c("PC_1"),
                              legend.type = "normal",
                              split.by = "seurat_clusters")
  testthat::expect_type(p, "list")
})

testthat::test_that("do_FeaturePlot: PASS - colorbar legend - split.by", {
  sample <- SCpubr:::test_list$sample

  p <- SCpubr::do_FeaturePlot(sample,
                              features = c("PC_1"),
                              legend.type = "colorbar",
                              split.by = "seurat_clusters")
  testthat::expect_type(p, "list")
})

testthat::test_that("do_FeaturePlot: PASS - colorsteps legend - split.by", {
  sample <- SCpubr:::test_list$sample

  p <- SCpubr::do_FeaturePlot(sample,
                              features = c("PC_1"),
                              legend.type = "colorsteps",
                              split.by = "seurat_clusters")
  testthat::expect_type(p, "list")
})

testthat::test_that("do_FeaturePlot: FAIL - wrong legend type", {
  sample <- SCpubr:::test_list$sample

  testthat::expect_error(SCpubr::do_FeaturePlot(sample,
                                                features = c("PC_1"),
                                                legend.type = "wrong"))
})

testthat::test_that("do_FeaturePlot: FAIL - wrong legend position", {
  sample <- SCpubr:::test_list$sample

  testthat::expect_error(SCpubr::do_FeaturePlot(sample,
                                                features = c("PC_1"),
                                                legend.position = "wrong"))
})
testthat::test_that("do_FeaturePlot: FAIL - wrong font.type", {
  sample <- SCpubr:::test_list$sample

  testthat::expect_error(SCpubr::do_FeaturePlot(sample,
                                                features = c("PC_1"),
                                                font.type = "wrong"))
})

testthat::test_that("do_FeaturePlot: PASS - plot axis", {
  sample <- SCpubr:::test_list$sample

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
