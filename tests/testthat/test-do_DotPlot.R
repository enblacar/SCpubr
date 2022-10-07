testthat::test_that("do_DotPlot: PASS - one variable", {
  sample <- SCpubr:::test_list$sample

  p <- SCpubr::do_DotPlot(sample = sample,
                          features = "CD14")
  testthat::expect_type(p, "list")

  p <- SCpubr::do_DotPlot(sample = sample,
                          features = "CD14",
                          dot_border = FALSE)
  testthat::expect_type(p, "list")
})

testthat::test_that("do_DotPlot: PASS - plot grid", {
  sample <- SCpubr:::test_list$sample

  p <- SCpubr::do_DotPlot(sample = sample,
                          features = "CD14",
                          plot.grid = TRUE)
  testthat::expect_type(p, "list")

  p <- SCpubr::do_DotPlot(sample = sample,
                          features = "CD14",
                          plot.grid = FALSE)
  testthat::expect_type(p, "list")
})

testthat::test_that("do_DotPlot: PASS - use_viridis", {
  sample <- SCpubr:::test_list$sample

  p <- SCpubr::do_DotPlot(sample = sample,
                          features = "CD14",
                          use_viridis = TRUE,
                          dot_border = FALSE)
  testthat::expect_type(p, "list")

  p <- SCpubr::do_DotPlot(sample = sample,
                          features = "CD14",
                          use_viridis = TRUE,
                          dot_border = TRUE)
  testthat::expect_type(p, "list")
})

testthat::test_that("do_DotPlot: PASS - one variable legend normal", {
  sample <- SCpubr:::test_list$sample

  p <- SCpubr::do_DotPlot(sample = sample,
                          features = "CD14",
                          legend.type = "normal")
  testthat::expect_type(p, "list")
})

testthat::test_that("do_DotPlot: PASS - one variable legend colorbar", {
  sample <- SCpubr:::test_list$sample

  p <- SCpubr::do_DotPlot(sample = sample,
                          features = "CD14",
                          legend.type = "colorbar")
  testthat::expect_type(p, "list")
})

testthat::test_that("do_DotPlot: PASS - one variable legend colorsteps", {
  sample <- SCpubr:::test_list$sample

  p <- SCpubr::do_DotPlot(sample = sample,
                          features = "CD14",
                          legend.type = "colorsteps")
  testthat::expect_type(p, "list")
})

testthat::test_that("do_DotPlot: PASS - one variable legend colorsteps legend to the right", {
  sample <- SCpubr:::test_list$sample

  p <- SCpubr::do_DotPlot(sample = sample,
                          features = "CD14",
                          legend.type = "colorsteps",
                          legend.position = "right")
  testthat::expect_type(p, "list")
})

testthat::test_that("do_DotPlot: FAIL - wrong legend type", {
  sample <- SCpubr:::test_list$sample

  testthat::expect_error(suppressWarnings({SCpubr::do_DotPlot(sample = sample,
                                                              features = "CD14",
                                                              flip = T,
                                                              legend.type = "wrong")}))

})

testthat::test_that("do_DotPlot: FAIL - wrong legend position", {
  sample <- SCpubr:::test_list$sample

  testthat::expect_error(suppressWarnings({SCpubr::do_DotPlot(sample = sample,
                                                              features = "CD14",
                                                              flip = T,
                                                              legend.position = "wrong")}))

})

testthat::test_that("do_DotPlot: FAIL - wrong font.type", {
  sample <- SCpubr:::test_list$sample

  testthat::expect_error(suppressWarnings({SCpubr::do_DotPlot(sample = sample,
                                                              features = "CD14",
                                                              flip = T,
                                                              font.type = "wrong")}))

})

testthat::test_that("do_DotPlot: PASS - one variable flip", {
  sample <- SCpubr:::test_list$sample

  p <- SCpubr::do_DotPlot(sample = sample,
                          features = "CD14",
                          flip = T)
  testthat::expect_type(p, "list")
})

testthat::test_that("do_DotPlot: PASS - multiple features", {
  sample <- SCpubr:::test_list$sample

  genes <- c("IL7R", "CCR7", "CD14", "LYZ",
             "S100A4", "MS4A1", "CD8A", "FCGR3A",
             "MS4A7", "GNLY", "NKG7", "FCER1A",
             "CST3", "PPBP")
  p <- suppressWarnings({SCpubr::do_DotPlot(sample = sample,
                          features = genes)})
  testthat::expect_type(p, "list")
})

testthat::test_that("do_DotPlot: PASS - multiple features flip", {
  sample <- SCpubr:::test_list$sample

  genes <- c("IL7R", "CCR7", "CD14", "LYZ",
             "S100A4", "MS4A1", "CD8A", "FCGR3A",
             "MS4A7", "GNLY", "NKG7", "FCER1A",
             "CST3", "PPBP")
  p <- suppressWarnings({SCpubr::do_DotPlot(sample = sample,
                                            features = genes,
                                            flip = T)})
  testthat::expect_type(p, "list")
})

testthat::test_that("do_DotPlot: PASS - multiple features flip rotate x labels", {
  sample <- SCpubr:::test_list$sample

  genes <- c("IL7R", "CCR7", "CD14", "LYZ",
             "S100A4", "MS4A1", "CD8A", "FCGR3A",
             "MS4A7", "GNLY", "NKG7", "FCER1A",
             "CST3", "PPBP")
  p <- suppressWarnings({SCpubr::do_DotPlot(sample = sample,
                                            features = genes,
                                            flip = T,
                                            rotate_x_axis_labels = T)})
  testthat::expect_type(p, "list")
})

testthat::test_that("do_DotPlot: PASS - list of features", {
  sample <- SCpubr:::test_list$sample

  genes <- list("Naive CD4+ T" = c("IL7R", "CCR7"),
                "CD14+ Mono" = c("CD14", "LYZ"),
                "Memory CD4+" = c("S100A4"),
                "B" = c("MS4A1"),
                "CD8+ T" = c("CD8A"),
                "FCGR3A+ Mono" = c("FCGR3A", "MS4A7"),
                "NK" = c("GNLY", "NKG7"),
                "DC" = c("FCER1A", "CST3"),
                "Platelet" = c("PPBP"))
  p <- suppressWarnings({SCpubr::do_DotPlot(sample = sample,
                          features = genes)})
  testthat::expect_type(p, "list")
})

testthat::test_that("do_DotPlot: FAIL - list of features flip", {
  sample <- SCpubr:::test_list$sample

  genes <- list("Naive CD4+ T" = c("IL7R", "CCR7"),
                "CD14+ Mono" = c("CD14", "LYZ"),
                "Memory CD4+" = c("S100A4"),
                "B" = c("MS4A1"),
                "CD8+ T" = c("CD8A"),
                "FCGR3A+ Mono" = c("FCGR3A", "MS4A7"),
                "NK" = c("GNLY", "NKG7"),
                "DC" = c("FCER1A", "CST3"),
                "Platelet" = c("PPBP"))
  testthat::expect_error(suppressWarnings({SCpubr::do_DotPlot(sample = sample,
                                            features = genes,
                                            flip = T)}))

})

testthat::test_that("do_DotPlot: PASS - list of features cluster idents", {
  sample <- SCpubr:::test_list$sample

  genes <- list("Naive CD4+ T" = c("IL7R", "CCR7"),
                "CD14+ Mono" = c("CD14", "LYZ"),
                "Memory CD4+" = c("S100A4"),
                "B" = c("MS4A1"),
                "CD8+ T" = c("CD8A"),
                "FCGR3A+ Mono" = c("FCGR3A", "MS4A7"),
                "NK" = c("GNLY", "NKG7"),
                "DC" = c("FCER1A", "CST3"),
                "Platelet" = c("PPBP"))
  p <- suppressWarnings({SCpubr::do_DotPlot(sample = sample,
                                            features = genes,
                                            cluster.idents = T)})
  testthat::expect_type(p, "list")
})

testthat::test_that("do_DotPlot: PASS - list of features cluster idents modify colors", {
  sample <- SCpubr:::test_list$sample

  genes <- list("Naive CD4+ T" = c("IL7R", "CCR7"),
                "CD14+ Mono" = c("CD14", "LYZ"),
                "Memory CD4+" = c("S100A4"),
                "B" = c("MS4A1"),
                "CD8+ T" = c("CD8A"),
                "FCGR3A+ Mono" = c("FCGR3A", "MS4A7"),
                "NK" = c("GNLY", "NKG7"),
                "DC" = c("FCER1A", "CST3"),
                "Platelet" = c("PPBP"))
  p <- suppressWarnings({SCpubr::do_DotPlot(sample = sample,
                                            features = genes,
                                            cluster.idents = T,
                                            colors.use = c("#001219", "#e9d8a6"))})
  testthat::expect_type(p, "list")
})

testthat::test_that("do_DotPlot: PASS - one variable split.by", {
  sample <- SCpubr:::test_list$sample

  p <- SCpubr::do_DotPlot(sample = sample,
                          features = "CD14",
                          split.by = "seurat_clusters")
  testthat::expect_type(p, "list")
})

testthat::test_that("do_DotPlot: PASS - one variable split.by factor", {
  sample <- SCpubr:::test_list$sample

  sample$seurat_clusters.factor <- factor(sample$seurat_clusters)
  p <- SCpubr::do_DotPlot(sample = sample,
                          features = "CD14",
                          split.by = "seurat_clusters.factor")
  testthat::expect_type(p, "list")
})


testthat::test_that("do_DotPlot: PASS - one variable xlab, ylab, title, subtitle, caption", {
  sample <- SCpubr:::test_list$sample

  p <- SCpubr::do_DotPlot(sample = sample,
                          features = "CD14",
                          xlab = "A",
                          ylab = "B",
                          plot.title = "C",
                          plot.subtitle = "D",
                          plot.caption = "E")
  testthat::expect_type(p, "list")
})
