library(SCpubr)
suppressWarnings(rm(test_list))
test_list <- SCpubr::test.data
sample <- test_list$sample

testthat::test_that("do_DotPlot: PASS - one variable", {


  p <- SCpubr::do_DotPlot(sample = sample,
                          features = "CD14")
  testthat::expect_type(p, "list")

  p <- SCpubr::do_DotPlot(sample = sample,
                          features = "CD14",
                          dot_border = FALSE)
  testthat::expect_type(p, "list")
})

testthat::test_that("do_DotPlot: PASS - plot grid", {


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


  p <- SCpubr::do_DotPlot(sample = sample,
                          features = "CD14",
                          legend.type = "normal")
  testthat::expect_type(p, "list")
})

testthat::test_that("do_DotPlot: PASS - one variable legend colorbar", {


  p <- SCpubr::do_DotPlot(sample = sample,
                          features = "CD14",
                          legend.type = "colorbar")
  testthat::expect_type(p, "list")
})

testthat::test_that("do_DotPlot: PASS - one variable legend colorsteps", {


  p <- SCpubr::do_DotPlot(sample = sample,
                          features = "CD14",
                          legend.type = "colorsteps")
  testthat::expect_type(p, "list")
})

testthat::test_that("do_DotPlot: PASS - one variable legend colorsteps legend to the right", {


  p <- SCpubr::do_DotPlot(sample = sample,
                          features = "CD14",
                          legend.type = "colorsteps",
                          legend.position = "right")
  testthat::expect_type(p, "list")
})

testthat::test_that("do_DotPlot: FAIL - wrong legend type", {


  testthat::expect_error(suppressWarnings({SCpubr::do_DotPlot(sample = sample,
                                                              features = "CD14",
                                                              flip = T,
                                                              legend.type = "wrong")}))

})

testthat::test_that("do_DotPlot: FAIL - wrong legend position", {


  testthat::expect_error(suppressWarnings({SCpubr::do_DotPlot(sample = sample,
                                                              features = "CD14",
                                                              flip = T,
                                                              legend.position = "wrong")}))

})

testthat::test_that("do_DotPlot: FAIL - wrong font.type", {


  testthat::expect_error(suppressWarnings({SCpubr::do_DotPlot(sample = sample,
                                                              features = "CD14",
                                                              flip = T,
                                                              font.type = "wrong")}))

})

testthat::test_that("do_DotPlot: PASS - one variable flip", {


  p <- SCpubr::do_DotPlot(sample = sample,
                          features = "CD14",
                          flip = T)
  testthat::expect_type(p, "list")
})

testthat::test_that("do_DotPlot: PASS - multiple features", {


  genes <- c("IL7R", "CCR7", "CD14", "LYZ",
             "S100A4", "MS4A1", "CD8A", "FCGR3A",
             "MS4A7", "GNLY", "NKG7", "FCER1A",
             "CST3", "PPBP")
  p <- suppressWarnings({SCpubr::do_DotPlot(sample = sample,
                          features = genes)})
  testthat::expect_type(p, "list")
})

testthat::test_that("do_DotPlot: PASS - multiple features flip", {


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


  p <- SCpubr::do_DotPlot(sample = sample,
                          features = "CD14",
                          split.by = "seurat_clusters")
  testthat::expect_type(p, "list")
})

testthat::test_that("do_DotPlot: PASS - one variable split.by factor", {


  sample$seurat_clusters.factor <- factor(sample$seurat_clusters)
  p <- SCpubr::do_DotPlot(sample = sample,
                          features = "CD14",
                          split.by = "seurat_clusters.factor")
  testthat::expect_type(p, "list")
})


testthat::test_that("do_DotPlot: PASS - one variable xlab, ylab, title, subtitle, caption", {


  p <- SCpubr::do_DotPlot(sample = sample,
                          features = "CD14",
                          xlab = "A",
                          ylab = "B",
                          plot.title = "C",
                          plot.subtitle = "D",
                          plot.caption = "E")
  testthat::expect_type(p, "list")
})
