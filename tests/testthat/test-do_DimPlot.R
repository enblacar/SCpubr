sample <- use_dataset()
testthat::test_that("do_DimPlot: PASS - sample", {
  p <- SCpubr::do_DimPlot(sample = sample)
  testthat::expect_type(p, "list")
})

testthat::test_that("do_DimPlot: PASS - title", {
  p <- SCpubr::do_DimPlot(sample = sample,
                          plot.title = "My awesome SC data set")
  testthat::expect_type(p, "list")
})

testthat::test_that("do_DimPlot: PASS - subtitle", {
  p <- SCpubr::do_DimPlot(sample = sample,
                          plot.subtitle = "My awesome SC data set")
  testthat::expect_type(p, "list")
})

testthat::test_that("do_DimPlot: PASS - caption", {
  p <- SCpubr::do_DimPlot(sample = sample,
                          plot.caption = "My awesome SC data set")
  testthat::expect_type(p, "list")
})

testthat::test_that("do_DimPlot: PASS - sample + group.by", {
  p <- SCpubr::do_DimPlot(sample = sample, group.by = "seurat_clusters")
  testthat::expect_type(p, "list")
})

testthat::test_that("do_DimPlot: PASS - sample + split.by", {
  p <- SCpubr::do_DimPlot(sample = sample, split.by = "seurat_clusters")
  testthat::expect_type(p, "list")
})

testthat::test_that("do_DimPlot: PASS - sample + split.by + idents.keep", {
  p <- SCpubr::do_DimPlot(sample = sample, split.by = "seurat_clusters", idents.keep = c("1", "3", "5"))
  testthat::expect_type(p, "list")
})

testthat::test_that("do_DimPlot: PASS - dims", {
  p <- SCpubr::do_DimPlot(sample = sample, dims = c(1, 2))
  testthat::expect_type(p, "list")
})

testthat::test_that("do_DimPlot: PASS - legend.position", {
  p <- SCpubr::do_DimPlot(sample = sample, legend.position = "top")
  testthat::expect_type(p, "list")
})

testthat::test_that("do_DimPlot: PASS - legend.ncol", {
  p <- SCpubr::do_DimPlot(sample = sample, legend.ncol = 2)
  testthat::expect_type(p, "list")
})

testthat::test_that("do_DimPlot: PASS - legend.nrow", {
  p <- SCpubr::do_DimPlot(sample = sample, legend.nrow = 2)
  testthat::expect_type(p, "list")
})

testthat::test_that("do_DimPlot: PASS - label", {
  p <- SCpubr::do_DimPlot(sample = sample, label = T)
  testthat::expect_type(p, "list")
})

testthat::test_that("do_DimPlot: PASS - order", {
  p <- SCpubr::do_DimPlot(sample = sample, order = "5", shuffle = F)
  testthat::expect_type(p, "list")
})

testthat::test_that("do_DimPlot: PASS - colors.use", {
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
  p <- SCpubr::do_DimPlot(sample = sample, cells.highlight = sample(colnames(sample), 50))
  testthat::expect_type(p, "list")
})

testthat::test_that("do_DimPlot: PASS - idents.highlight", {
  p <- SCpubr::do_DimPlot(sample = sample, idents.highlight = "5")
  testthat::expect_type(p, "list")
})

testthat::test_that("do_DimPlot: PASS - cells.highlight and idents.highlight", {
  p <- SCpubr::do_DimPlot(sample = sample, cells.highlight = sample(colnames(sample), 50), idents.highlight = "2")
  testthat::expect_type(p, "list")
})

testthat::test_that("do_DimPlot: PASS - idents.keep", {
  p <- SCpubr::do_DimPlot(sample = sample, idents.keep = "5")
  testthat::expect_type(p, "list")
})


testthat::test_that("do_DimPlot: FAIL - group.by and split.by used", {
  testthat::expect_error(SCpubr::do_DimPlot(sample = sample, group.by = "seurat_clusters", split.by = "seurat_clusters"))
})

testthat::test_that("do_DimPlot: FAIL - group.by and cells.highlights used", {
  testthat::expect_error(SCpubr::do_DimPlot(sample = sample, group.by = "seurat_clusters", cells.highlight = colnames(sample)))
})

testthat::test_that("do_DimPlot: FAIL - split.by and cells.highlights used", {
  testthat::expect_error(SCpubr::do_DimPlot(sample = sample, split.by = "seurat_clusters", cells.highlight = colnames(sample)))
})

testthat::test_that("do_DimPlot: WARNING - order and shuffle used", {
  testthat::expect_warning(SCpubr::do_DimPlot(sample = sample, order = "4", shuffle = T))
})

testthat::test_that("do_DimPlot: FAIL - more than one NA values", {
  testthat::expect_error(SCpubr::do_DimPlot(sample = sample, na.value = c("red", "blue")))
})

testthat::test_that("do_DimPlot: WARNING - raster = T and pt.size lower than 1", {
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
  p <- SCpubr::do_DimPlot(sample = sample, group.by = "seurat_clusters", colors.use = colors)
  testthat::expect_type(p, "list")
})

testthat::test_that("do_DimPlot: PASS - split.by + colors", {
  p <- SCpubr::do_DimPlot(sample = sample, split.by = "seurat_clusters", colors.use = colors)
  testthat::expect_type(p, "list")
})

testthat::test_that("do_DimPlot: FAIL - more than 1 color with cells highlight", {
  testthat::expect_error(SCpubr::do_DimPlot(sample = sample, colors.use = colors, idents.highlight = "4"))
})

testthat::test_that("do_DimPlot: FAIL - idents.keep not in the levels of the sample", {
  testthat::expect_error(SCpubr::do_DimPlot(sample = sample, idents.keep = c("4", "Not an ident")))
})

testthat::test_that("do_DimPlot: FAIL - idents.keep not in the unique values of group.by", {
  testthat::expect_error(SCpubr::do_DimPlot(sample = sample, group.by = "orig.ident", idents.keep = c("4", "Not an ident")))
})

testthat::test_that("do_DimPlot: FAIL - idents.keep not in the unique values of split.by", {
  testthat::expect_error(SCpubr::do_DimPlot(sample = sample, split.by = "orig.ident", idents.keep = c("4", "Not an ident")))
})


testthat::test_that("do_DimPlot: PASS - split.by + plot.title, subtitle and caption", {
  p <- SCpubr::do_DimPlot(sample = sample, split.by = "orig.ident",
                          plot.title = "Plot title",
                          plot.subtitle = "Plot subtitle",
                          plot.caption = "Plot caption")
  testthat::expect_type(p, "list")
})

testthat::test_that("do_DimPlot: PASS - legend = F", {
  p <- SCpubr::do_DimPlot(sample = sample,
                          legend = F)
  testthat::expect_type(p, "list")
})

testthat::test_that("do_DimPlot: PASS - dims different", {
  p <- SCpubr::do_DimPlot(sample = sample,
                          dims = c(2, 1))
  testthat::expect_type(p, "list")
})

testthat::test_that("do_DimPlot: PASS - diffusion maps", {
  sample@reductions$diffusion <- sample@reductions$umap
  p <- SCpubr::do_DimPlot(sample = sample,
                          reduction = "diffusion")
  testthat::expect_type(p, "list")
})


testthat::test_that("do_DimPlot: PASS - split.by + factor", {
  sample$seurat_clusters <- factor(sample$seurat_clusters)
  p <- SCpubr::do_DimPlot(sample = sample, split.by = "seurat_clusters", colors.use = colors)
  testthat::expect_type(p, "list")
})

testthat::test_that("do_DimPlot: PASS - split.by + factor + idents.keep", {
  sample$seurat_clusters <- factor(sample$seurat_clusters)
  p <- SCpubr::do_DimPlot(sample = sample, split.by = "seurat_clusters", colors.use = colors, idents.keep = "4")
  testthat::expect_type(p, "list")
})


