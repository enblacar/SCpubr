
testthat::test_that("do_DimPlot: PASS - sample", {
  sample <- SCpubr:::test_list$sample

  p <- SCpubr::do_DimPlot(sample = sample)
  testthat::expect_type(p, "list")
})

testthat::test_that("do_DimPlot: PASS - plot axis", {
  sample <- SCpubr:::test_list$sample


  p <- SCpubr::do_DimPlot(sample = sample, plot.axes = TRUE)
  testthat::expect_type(p, "list")

  p <- SCpubr::do_DimPlot(sample = sample, reduction = "pca", plot.axes = TRUE)
  testthat::expect_type(p, "list")

  p <- SCpubr::do_DimPlot(sample = sample, dims = c(2, 1), plot.axes = TRUE)
  testthat::expect_type(p, "list")

  sample@reductions$diffusion <- sample@reductions$umap
  p <- SCpubr::do_DimPlot(sample = sample,
                          reduction = "diffusion",
                          plot.axes = TRUE)
  testthat::expect_type(p, "list")
})

testthat::test_that("do_DimPlot: PASS - sample cell_borders", {
  sample <- SCpubr:::test_list$sample


  p <- SCpubr::do_DimPlot(sample = sample, plot_cell_borders = T)
  testthat::expect_type(p, "list")
  p <- SCpubr::do_DimPlot(sample = sample, plot_cell_borders = T, raster = T, pt.size = 1)
  testthat::expect_type(p, "list")
  p <- SCpubr::do_DimPlot(sample = sample, plot_cell_borders = T, idents.highlight = "1")
  testthat::expect_type(p, "list")
  p <- SCpubr::do_DimPlot(sample = sample, plot_cell_borders = T, raster = T, idents.highlight = "1", pt.size = 1)
  testthat::expect_type(p, "list")
  p <- SCpubr::do_DimPlot(sample = sample, plot_cell_borders = T, split.by = "seurat_clusters")
  testthat::expect_type(p, "list")
  p <- SCpubr::do_DimPlot(sample = sample, plot_cell_borders = T, raster = T, split.by = "seurat_clusters", pt.size = 1)
  testthat::expect_type(p, "list")
})

testthat::test_that("do_DimPlot: PASS - sample marginal", {
  sample <- SCpubr:::test_list$sample


  p <- SCpubr::do_DimPlot(sample = sample,
                          plot_marginal_distributions = TRUE,
                          marginal.type = "density",
                          plot_cell_borders = FALSE)
  testthat::expect_type(p, "list")

  p <- SCpubr::do_DimPlot(sample = sample,
                          plot_marginal_distributions = TRUE,
                          marginal.type = "histogram",
                          plot_cell_borders = FALSE)
  testthat::expect_type(p, "list")

  p <- suppressWarnings({SCpubr::do_DimPlot(sample = sample,
                                            plot_marginal_distributions = TRUE,
                                            plot_cell_borders = FALSE,
                                            marginal.type = "violin")})
  testthat::expect_type(p, "list")

  p <- SCpubr::do_DimPlot(sample = sample,
                          plot_marginal_distributions = TRUE,
                          marginal.type = "boxplot",
                          plot_cell_borders = FALSE)
  testthat::expect_type(p, "list")

  p <- SCpubr::do_DimPlot(sample = sample,
                          plot_marginal_distributions = TRUE,
                          marginal.type = "densigram",
                          plot_cell_borders = FALSE)
  testthat::expect_type(p, "list")
})

testthat::test_that("do_DimPlot: PASS - sample marginal size", {
  sample <- SCpubr:::test_list$sample

  p <- SCpubr::do_DimPlot(sample = sample,
                          plot_marginal_distributions = TRUE,
                          marginal.size = 9,
                          plot_cell_borders = FALSE)
  testthat::expect_type(p, "list")
})

testthat::test_that("do_DimPlot: PASS - sample marginal group", {
  sample <- SCpubr:::test_list$sample

  p <- SCpubr::do_DimPlot(sample = sample,
                          plot_marginal_distributions = TRUE,
                          marginal.group = F,
                          plot_cell_borders = FALSE)
  testthat::expect_type(p, "list")
})

testthat::test_that("do_DimPlot: FAIL - sample marginal wrong marginal.type", {
  sample <- SCpubr:::test_list$sample

  testthat::expect_error({SCpubr::do_DimPlot(sample = sample,
                                             plot_marginal_distributions = TRUE,
                                             plot_cell_borders = FALSE,
                                             marginal.type = "wrong")})
})

testthat::test_that("do_DimPlot: FAIL - sample marginal used alongside split.by or cells.highlight", {
  sample <- SCpubr:::test_list$sample

  testthat::expect_error({SCpubr::do_DimPlot(sample = sample,
                                             plot_marginal_distributions = TRUE,
                                             plot_cell_borders = FALSE,
                                             split.by = "seurat_clusters")})

  testthat::expect_error({SCpubr::do_DimPlot(sample = sample,
                                             plot_marginal_distributions = TRUE,
                                             plot_cell_borders = FALSE,
                                             cells.highlight =  colnames(sample))})

  testthat::expect_error({SCpubr::do_DimPlot(sample = sample,
                                             plot_marginal_distributions = TRUE,
                                             plot_cell_borders = FALSE,
                                             idents.highlight =  c("1"))})
})

testthat::test_that("do_DimPlot: PASS - title", {
  sample <- SCpubr:::test_list$sample

  p <- SCpubr::do_DimPlot(sample = sample,
                          plot.title = "My awesome SC data set")
  testthat::expect_type(p, "list")
})

testthat::test_that("do_DimPlot: PASS - subtitle", {
  sample <- SCpubr:::test_list$sample

  p <- SCpubr::do_DimPlot(sample = sample,
                          plot.subtitle = "My awesome SC data set")
  testthat::expect_type(p, "list")
})

testthat::test_that("do_DimPlot: PASS - caption", {
  sample <- SCpubr:::test_list$sample

  p <- SCpubr::do_DimPlot(sample = sample,
                          plot.caption = "My awesome SC data set")
  testthat::expect_type(p, "list")
})

testthat::test_that("do_DimPlot: PASS - sample + group.by", {
  sample <- SCpubr:::test_list$sample

  p <- SCpubr::do_DimPlot(sample = sample, group.by = "seurat_clusters")
  testthat::expect_type(p, "list")
})

testthat::test_that("do_DimPlot: PASS - sample + split.by", {
  sample <- SCpubr:::test_list$sample

  p <- SCpubr::do_DimPlot(sample = sample, split.by = "seurat_clusters")
  testthat::expect_type(p, "list")
})

testthat::test_that("do_DimPlot: PASS - sample + split.by + idents.keep", {
  sample <- SCpubr:::test_list$sample

  p <- SCpubr::do_DimPlot(sample = sample, split.by = "seurat_clusters", idents.keep = c("1", "3", "5"))
  testthat::expect_type(p, "list")
})

testthat::test_that("do_DimPlot: PASS - dims", {
  sample <- SCpubr:::test_list$sample

  p <- SCpubr::do_DimPlot(sample = sample, dims = c(1, 2))
  testthat::expect_type(p, "list")
})

testthat::test_that("do_DimPlot: PASS - legend.position", {
  sample <- SCpubr:::test_list$sample

  p <- SCpubr::do_DimPlot(sample = sample, legend.position = "top")
  testthat::expect_type(p, "list")
})

testthat::test_that("do_DimPlot: PASS - legend.ncol", {
  sample <- SCpubr:::test_list$sample

  p <- SCpubr::do_DimPlot(sample = sample, legend.ncol = 2)
  testthat::expect_type(p, "list")
})

testthat::test_that("do_DimPlot: PASS - legend.nrow", {
  sample <- SCpubr:::test_list$sample

  p <- SCpubr::do_DimPlot(sample = sample, legend.nrow = 2)
  testthat::expect_type(p, "list")
})

testthat::test_that("do_DimPlot: PASS - label", {
  sample <- SCpubr:::test_list$sample

  p <- SCpubr::do_DimPlot(sample = sample, label = T)
  testthat::expect_type(p, "list")
})

testthat::test_that("do_DimPlot: PASS - order", {
  sample <- SCpubr:::test_list$sample

  p <- SCpubr::do_DimPlot(sample = sample, order = "5", shuffle = F)
  testthat::expect_type(p, "list")
})

testthat::test_that("do_DimPlot: PASS - colors.use", {
  sample <- SCpubr:::test_list$sample

  p <- SCpubr::do_DimPlot(sample = sample, colors.use = c("0" = "#001219",
                                                          "1" = "#005f73",
                                                          "2" = "#0a9396",
                                                          "3" = "#94d2bd",
                                                          "4" = "#e9d8a6",
                                                          "5" = "#ee9b00",
                                                          "6" = "#ca6702",
                                                          "7" = "#bb3e03",
                                                          "8" = "#ae2012"))
  testthat::expect_type(p, "list")
})

testthat::test_that("do_DimPlot: PASS - cells.highlight", {
  sample <- SCpubr:::test_list$sample

  p <- SCpubr::do_DimPlot(sample = sample, cells.highlight = sample(colnames(sample), 50))
  testthat::expect_type(p, "list")
})

testthat::test_that("do_DimPlot: PASS - idents.highlight", {
  sample <- SCpubr:::test_list$sample

  p <- SCpubr::do_DimPlot(sample = sample, idents.highlight = "5")
  testthat::expect_type(p, "list")
})

testthat::test_that("do_DimPlot: PASS - cells.highlight and idents.highlight", {
  sample <- SCpubr:::test_list$sample

  p <- SCpubr::do_DimPlot(sample = sample, cells.highlight = sample(colnames(sample), 50), idents.highlight = "2")
  testthat::expect_type(p, "list")
})

testthat::test_that("do_DimPlot: PASS - idents.keep", {
  sample <- SCpubr:::test_list$sample

  p <- SCpubr::do_DimPlot(sample = sample, idents.keep = "5")
  testthat::expect_type(p, "list")
})


testthat::test_that("do_DimPlot: FAIL - group.by and split.by used", {
  sample <- SCpubr:::test_list$sample

  testthat::expect_error(SCpubr::do_DimPlot(sample = sample, group.by = "seurat_clusters", split.by = "seurat_clusters"))
})

testthat::test_that("do_DimPlot: FAIL - group.by and cells.highlights used", {
  sample <- SCpubr:::test_list$sample

  testthat::expect_error(SCpubr::do_DimPlot(sample = sample, group.by = "seurat_clusters", cells.highlight = colnames(sample)))
})

testthat::test_that("do_DimPlot: FAIL - split.by and cells.highlights used", {
  sample <- SCpubr:::test_list$sample

  testthat::expect_error(SCpubr::do_DimPlot(sample = sample, split.by = "seurat_clusters", cells.highlight = colnames(sample)))
})

testthat::test_that("do_DimPlot: WARNING - order and shuffle used", {
  sample <- SCpubr:::test_list$sample

  testthat::expect_warning(SCpubr::do_DimPlot(sample = sample, order = "4", shuffle = T))
})

testthat::test_that("do_DimPlot: FAIL - more than one NA values", {
  sample <- SCpubr:::test_list$sample

  testthat::expect_error(SCpubr::do_DimPlot(sample = sample, na.value = c("red", "blue")))
})

testthat::test_that("do_DimPlot: WARNING - raster = T and pt.size lower than 1", {
  sample <- SCpubr:::test_list$sample

  testthat::expect_warning(SCpubr::do_DimPlot(sample = sample, raster = T, pt.size = 0.5))
})

colors <- c("0" = "#001219",
            "1" = "#005f73",
            "2" = "#0a9396",
            "3" = "#94d2bd",
            "4" = "#e9d8a6",
            "5" = "#ee9b00",
            "6" = "#ca6702",
            "7" = "#bb3e03",
            "8" = "#ae2012")

testthat::test_that("do_DimPlot: PASS - group.by + colors", {
  sample <- SCpubr:::test_list$sample

  p <- SCpubr::do_DimPlot(sample = sample, group.by = "seurat_clusters", colors.use = colors)
  testthat::expect_type(p, "list")
})

testthat::test_that("do_DimPlot: PASS - split.by + colors", {
  sample <- SCpubr:::test_list$sample

  p <- SCpubr::do_DimPlot(sample = sample, split.by = "seurat_clusters", colors.use = colors)
  testthat::expect_type(p, "list")
})

testthat::test_that("do_DimPlot: FAIL - more than 1 color with cells highlight", {
  sample <- SCpubr:::test_list$sample

  testthat::expect_error(SCpubr::do_DimPlot(sample = sample, colors.use = colors, idents.highlight = "4"))
})

testthat::test_that("do_DimPlot: FAIL - idents.keep not in the levels of the sample", {
  sample <- SCpubr:::test_list$sample

  testthat::expect_error(SCpubr::do_DimPlot(sample = sample, idents.keep = c("4", "Not an ident")))
})

testthat::test_that("do_DimPlot: FAIL - idents.keep not in the unique values of group.by", {
  sample <- SCpubr:::test_list$sample

  testthat::expect_error(SCpubr::do_DimPlot(sample = sample, group.by = "orig.ident", idents.keep = c("4", "Not an ident")))
})

testthat::test_that("do_DimPlot: FAIL - idents.keep not in the unique values of split.by", {
  sample <- SCpubr:::test_list$sample

  testthat::expect_error(SCpubr::do_DimPlot(sample = sample, split.by = "orig.ident", idents.keep = c("4", "Not an ident")))
})


testthat::test_that("do_DimPlot: PASS - split.by + plot.title, subtitle and caption", {
  sample <- SCpubr:::test_list$sample

  p <- SCpubr::do_DimPlot(sample = sample, split.by = "orig.ident",
                          plot.title = "Plot title",
                          plot.subtitle = "Plot subtitle",
                          plot.caption = "Plot caption")
  testthat::expect_type(p, "list")
})

testthat::test_that("do_DimPlot: PASS - legend.position none", {
  sample <- SCpubr:::test_list$sample

  p <- SCpubr::do_DimPlot(sample = sample,
                          legend.position = "none")
  testthat::expect_type(p, "list")
})

testthat::test_that("do_DimPlot: PASS - dims different", {
  sample <- SCpubr:::test_list$sample

  p <- SCpubr::do_DimPlot(sample = sample,
                          dims = c(2, 1))
  testthat::expect_type(p, "list")
})

testthat::test_that("do_DimPlot: PASS - diffusion maps", {
  sample <- SCpubr:::test_list$sample

  sample@reductions$diffusion <- sample@reductions$umap
  p <- SCpubr::do_DimPlot(sample = sample,
                          reduction = "diffusion")
  testthat::expect_type(p, "list")
})

testthat::test_that("do_DimPlot: PASS - group.by + idents.keep", {
  sample <- SCpubr:::test_list$sample

  p <- SCpubr::do_DimPlot(sample = sample,
                          group.by = "seurat_clusters",
                          idents.keep = "4")
  testthat::expect_type(p, "list")
})


testthat::test_that("do_DimPlot: PASS - split.by + factor", {
  sample <- SCpubr:::test_list$sample

  sample$seurat_clusters <- factor(sample$seurat_clusters)
  p <- SCpubr::do_DimPlot(sample = sample, split.by = "seurat_clusters", colors.use = colors)
  testthat::expect_type(p, "list")
})

testthat::test_that("do_DimPlot: PASS - split.by + factor + idents.keep", {
  sample <- SCpubr:::test_list$sample

  sample$seurat_clusters <- factor(sample$seurat_clusters)
  p <- SCpubr::do_DimPlot(sample = sample, split.by = "seurat_clusters", colors.use = colors, idents.keep = "4")
  testthat::expect_type(p, "list")
})

testthat::test_that("do_DimPlot: FAIL - wrong font.type", {
  sample <- SCpubr:::test_list$sample

  testthat::expect_error(SCpubr::do_DimPlot(sample = sample, font.type = "wrong"))
})
