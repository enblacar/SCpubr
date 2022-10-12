
testthat::test_that("do_DotPlot: PASS - one variable", {


  p <- SCpubr::do_DotPlot(sample = sample,
                          features = "EPC1")
  testthat::expect_type(p, "list")

  p <- SCpubr::do_DotPlot(sample = sample,
                          features = "EPC1",
                          dot_border = FALSE)
  testthat::expect_type(p, "list")
})

testthat::test_that("do_DotPlot: PASS - plot grid", {


  p <- SCpubr::do_DotPlot(sample = sample,
                          features = "EPC1",
                          plot.grid = TRUE)
  testthat::expect_type(p, "list")

  p <- SCpubr::do_DotPlot(sample = sample,
                          features = "EPC1",
                          plot.grid = FALSE)
  testthat::expect_type(p, "list")
})

testthat::test_that("do_DotPlot: PASS - use_viridis", {


  p <- SCpubr::do_DotPlot(sample = sample,
                          features = "EPC1",
                          use_viridis = TRUE,
                          dot_border = FALSE)
  testthat::expect_type(p, "list")

  p <- SCpubr::do_DotPlot(sample = sample,
                          features = "EPC1",
                          use_viridis = TRUE,
                          dot_border = TRUE)
  testthat::expect_type(p, "list")
})

testthat::test_that("do_DotPlot: PASS - one variable legend normal", {


  p <- SCpubr::do_DotPlot(sample = sample,
                          features = "EPC1",
                          legend.type = "normal")
  testthat::expect_type(p, "list")
})

testthat::test_that("do_DotPlot: PASS - one variable legend colorbar", {


  p <- SCpubr::do_DotPlot(sample = sample,
                          features = "EPC1",
                          legend.type = "colorbar")
  testthat::expect_type(p, "list")
})

testthat::test_that("do_DotPlot: PASS - one variable legend colorsteps", {


  p <- SCpubr::do_DotPlot(sample = sample,
                          features = "EPC1",
                          legend.type = "colorsteps")
  testthat::expect_type(p, "list")
})

testthat::test_that("do_DotPlot: PASS - one variable legend colorsteps legend to the right", {


  p <- SCpubr::do_DotPlot(sample = sample,
                          features = "EPC1",
                          legend.type = "colorsteps",
                          legend.position = "right")
  testthat::expect_type(p, "list")
})

testthat::test_that("do_DotPlot: FAIL - wrong legend type", {


  testthat::expect_error(suppressWarnings({SCpubr::do_DotPlot(sample = sample,
                                                              features = "EPC1",
                                                              flip = TRUE,
                                                              legend.type = "wrong")}))

})

testthat::test_that("do_DotPlot: FAIL - wrong legend position", {


  testthat::expect_error(suppressWarnings({SCpubr::do_DotPlot(sample = sample,
                                                              features = "EPC1",
                                                              flip = TRUE,
                                                              legend.position = "wrong")}))

})

testthat::test_that("do_DotPlot: FAIL - wrong font.type", {


  testthat::expect_error(suppressWarnings({SCpubr::do_DotPlot(sample = sample,
                                                              features = "EPC1",
                                                              flip = TRUE,
                                                              font.type = "wrong")}))

})

testthat::test_that("do_DotPlot: PASS - one variable flip", {


  p <- SCpubr::do_DotPlot(sample = sample,
                          features = "EPC1",
                          flip = TRUE)
  testthat::expect_type(p, "list")
})

testthat::test_that("do_DotPlot: PASS - multiple features", {


  genes <- list("Naive CD4+ T" = Seurat::VariableFeatures(sample)[1:2],
                "EPC1+ Mono" = Seurat::VariableFeatures(sample)[3:4],
                "Memory CD4+" = Seurat::VariableFeatures(sample)[5],
                "B" = Seurat::VariableFeatures(sample)[6],
                "CD8+ T" = Seurat::VariableFeatures(sample)[7],
                "FCGR3A+ Mono" = Seurat::VariableFeatures(sample)[8:9],
                "NK" = Seurat::VariableFeatures(sample)[10:11],
                "DC" = Seurat::VariableFeatures(sample)[12:13],
                "Platelet" = Seurat::VariableFeatures(sample)[14])
  genes <- unname(unlist(genes))
  p <- suppressWarnings({SCpubr::do_DotPlot(sample = sample,
                          features = genes)})
  testthat::expect_type(p, "list")
})

testthat::test_that("do_DotPlot: PASS - multiple features flip", {


  genes <- list("Naive CD4+ T" = Seurat::VariableFeatures(sample)[1:2],
                "EPC1+ Mono" = Seurat::VariableFeatures(sample)[3:4],
                "Memory CD4+" = Seurat::VariableFeatures(sample)[5],
                "B" = Seurat::VariableFeatures(sample)[6],
                "CD8+ T" = Seurat::VariableFeatures(sample)[7],
                "FCGR3A+ Mono" = Seurat::VariableFeatures(sample)[8:9],
                "NK" = Seurat::VariableFeatures(sample)[10:11],
                "DC" = Seurat::VariableFeatures(sample)[12:13],
                "Platelet" = Seurat::VariableFeatures(sample)[14])
  genes <- unname(unlist(genes))
  p <- suppressWarnings({SCpubr::do_DotPlot(sample = sample,
                                            features = genes,
                                            flip = TRUE)})
  testthat::expect_type(p, "list")
})

testthat::test_that("do_DotPlot: PASS - multiple features flip rotate x labels", {


  genes <- list("Naive CD4+ T" = Seurat::VariableFeatures(sample)[1:2],
                "EPC1+ Mono" = Seurat::VariableFeatures(sample)[3:4],
                "Memory CD4+" = Seurat::VariableFeatures(sample)[5],
                "B" = Seurat::VariableFeatures(sample)[6],
                "CD8+ T" = Seurat::VariableFeatures(sample)[7],
                "FCGR3A+ Mono" = Seurat::VariableFeatures(sample)[8:9],
                "NK" = Seurat::VariableFeatures(sample)[10:11],
                "DC" = Seurat::VariableFeatures(sample)[12:13],
                "Platelet" = Seurat::VariableFeatures(sample)[14])
  genes <- unname(unlist(genes))
  p <- suppressWarnings({SCpubr::do_DotPlot(sample = sample,
                                            features = genes,
                                            flip = TRUE,
                                            rotate_x_axis_labels = TRUE)})
  testthat::expect_type(p, "list")
})

testthat::test_that("do_DotPlot: PASS - list of features", {


  genes <- list("Naive CD4+ T" = Seurat::VariableFeatures(sample)[1:2],
                "EPC1+ Mono" = Seurat::VariableFeatures(sample)[3:4],
                "Memory CD4+" = Seurat::VariableFeatures(sample)[5],
                "B" = Seurat::VariableFeatures(sample)[6],
                "CD8+ T" = Seurat::VariableFeatures(sample)[7],
                "FCGR3A+ Mono" = Seurat::VariableFeatures(sample)[8:9],
                "NK" = Seurat::VariableFeatures(sample)[10:11],
                "DC" = Seurat::VariableFeatures(sample)[12:13],
                "Platelet" = Seurat::VariableFeatures(sample)[14])
  p <- suppressWarnings({SCpubr::do_DotPlot(sample = sample,
                          features = genes)})
  testthat::expect_type(p, "list")
})

testthat::test_that("do_DotPlot: FAIL - list of features flip", {


  genes <- list("Naive CD4+ T" = Seurat::VariableFeatures(sample)[1:2],
                "EPC1+ Mono" = Seurat::VariableFeatures(sample)[3:4],
                "Memory CD4+" = Seurat::VariableFeatures(sample)[5],
                "B" = Seurat::VariableFeatures(sample)[6],
                "CD8+ T" = Seurat::VariableFeatures(sample)[7],
                "FCGR3A+ Mono" = Seurat::VariableFeatures(sample)[8:9],
                "NK" = Seurat::VariableFeatures(sample)[10:11],
                "DC" = Seurat::VariableFeatures(sample)[12:13],
                "Platelet" = Seurat::VariableFeatures(sample)[14])
  testthat::expect_error(suppressWarnings({SCpubr::do_DotPlot(sample = sample,
                                            features = genes,
                                            flip = TRUE)}))

})

testthat::test_that("do_DotPlot: PASS - list of features cluster idents", {


  genes <- list("Naive CD4+ T" = Seurat::VariableFeatures(sample)[1:2],
                "EPC1+ Mono" = Seurat::VariableFeatures(sample)[3:4],
                "Memory CD4+" = Seurat::VariableFeatures(sample)[5],
                "B" = Seurat::VariableFeatures(sample)[6],
                "CD8+ T" = Seurat::VariableFeatures(sample)[7],
                "FCGR3A+ Mono" = Seurat::VariableFeatures(sample)[8:9],
                "NK" = Seurat::VariableFeatures(sample)[10:11],
                "DC" = Seurat::VariableFeatures(sample)[12:13],
                "Platelet" = Seurat::VariableFeatures(sample)[14])
  p <- suppressWarnings({SCpubr::do_DotPlot(sample = sample,
                                            features = genes,
                                            cluster.idents = TRUE)})
  testthat::expect_type(p, "list")
})

testthat::test_that("do_DotPlot: PASS - list of features cluster idents modify colors", {


  genes <- list("Naive CD4+ T" = Seurat::VariableFeatures(sample)[1:2],
                "EPC1+ Mono" = Seurat::VariableFeatures(sample)[3:4],
                "Memory CD4+" = Seurat::VariableFeatures(sample)[5],
                "B" = Seurat::VariableFeatures(sample)[6],
                "CD8+ T" = Seurat::VariableFeatures(sample)[7],
                "FCGR3A+ Mono" = Seurat::VariableFeatures(sample)[8:9],
                "NK" = Seurat::VariableFeatures(sample)[10:11],
                "DC" = Seurat::VariableFeatures(sample)[12:13],
                "Platelet" = Seurat::VariableFeatures(sample)[14])
  p <- suppressWarnings({SCpubr::do_DotPlot(sample = sample,
                                            features = genes,
                                            cluster.idents = TRUE,
                                            colors.use = c("#001219", "#e9d8a6"))})
  testthat::expect_type(p, "list")
})

testthat::test_that("do_DotPlot: PASS - one variable split.by", {


  p <- SCpubr::do_DotPlot(sample = sample,
                          features = "EPC1",
                          split.by = "seurat_clusters")
  testthat::expect_type(p, "list")
})

testthat::test_that("do_DotPlot: PASS - one variable split.by factor", {


  sample$seurat_clusters.factor <- factor(sample$seurat_clusters)
  p <- SCpubr::do_DotPlot(sample = sample,
                          features = "EPC1",
                          split.by = "seurat_clusters.factor")
  testthat::expect_type(p, "list")
})


testthat::test_that("do_DotPlot: PASS - one variable xlab, ylab, title, subtitle, caption", {


  p <- SCpubr::do_DotPlot(sample = sample,
                          features = "EPC1",
                          xlab = "A",
                          ylab = "B",
                          plot.title = "C",
                          plot.subtitle = "D",
                          plot.caption = "E")
  testthat::expect_type(p, "list")
})
