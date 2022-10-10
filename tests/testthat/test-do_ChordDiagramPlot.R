
testthat::test_that("do_ChordDiagramPlot: PASS - default", {

  sample$seurat_clusters_char <- as.character(sample$seurat_clusters)
  sample$orig.ident_char <- sample$orig.ident
  sample$orig.ident <- factor(sample$orig.ident)
  p <- SCpubr::do_ChordDiagramPlot(sample = sample,
                                   from = "seurat_clusters",
                                   to = "orig.ident")
  testthat::expect_type(p, "list")

  p <- SCpubr::do_ChordDiagramPlot(sample = sample,
                                   from = "seurat_clusters_char",
                                   to = "orig.ident_char")
  testthat::expect_type(p, "list")
})


testthat::test_that("do_ChordDiagramPlot: PASS - colors", {
  sample$seurat_clusters_char <- as.character(sample$seurat_clusters)
  sample$orig.ident_char <- sample$orig.ident
  sample$orig.ident <- factor(sample$orig.ident)

  p <- SCpubr::do_ChordDiagramPlot(sample = sample,
                                   from = "orig.ident",
                                   to = "seurat_clusters",
                                   colors.from = c("Cell" = "blue"))
  testthat::expect_type(p, "list")

  p <- SCpubr::do_ChordDiagramPlot(sample = sample,
                                   from = "orig.ident_char",
                                   to = "seurat_clusters",
                                   colors.from = c("Cell" = "blue"))
  testthat::expect_type(p, "list")

  p <- SCpubr::do_ChordDiagramPlot(sample = sample,
                                   from = "orig.ident",
                                   to = "seurat_clusters")
  testthat::expect_type(p, "list")

  p <- SCpubr::do_ChordDiagramPlot(sample = sample,
                                   from = "orig.ident_char",
                                   to = "seurat_clusters")
  testthat::expect_type(p, "list")

  p <- SCpubr::do_ChordDiagramPlot(sample = sample,
                                   from = "seurat_clusters",
                                   to = "orig.ident",
                                   colors.to = c("Cell" = "blue"))
  testthat::expect_type(p, "list")

  p <- SCpubr::do_ChordDiagramPlot(sample = sample,
                                   from = "seurat_clusters",
                                   to = "orig.ident_char",
                                   colors.to = c("Cell" = "blue"))
  testthat::expect_type(p, "list")

  p <- SCpubr::do_ChordDiagramPlot(sample = sample,
                                   from = "seurat_clusters",
                                   to = "orig.ident")
  testthat::expect_type(p, "list")

  p <- SCpubr::do_ChordDiagramPlot(sample = sample,
                                   from = "seurat_clusters",
                                   to = "orig.ident_char")
  testthat::expect_type(p, "list")

  p <- SCpubr::do_ChordDiagramPlot(sample = sample,
                                   from = "orig.ident",
                                   to = "seurat_clusters",
                                   colors.from = c("Cell" = "#345211"))
  testthat::expect_type(p, "list")

  p <- SCpubr::do_ChordDiagramPlot(sample = sample,
                                   from = "orig.ident",
                                   to = "seurat_clusters",
                                   colors.from = c("Cell" = "#345211FF"))
  testthat::expect_type(p, "list")

  sample$orig.ident <- ifelse(sample$seurat_clusters %in% c("0"), "A", "B")
  p <- SCpubr::do_ChordDiagramPlot(sample = sample,
                                   from = "orig.ident",
                                   to = "seurat_clusters",
                                   colors.from = c("A" = "#345211", "B" = "#345222"),
                                   highlight_group = "A")
  testthat::expect_type(p, "list")

  p <- SCpubr::do_ChordDiagramPlot(sample = sample,
                                   from = "orig.ident",
                                   to = "seurat_clusters",
                                   colors.from = c("Cell" = "#345211FF"),
                                   highlight_group = "Cell")
  testthat::expect_type(p, "list")

})

testthat::test_that("do_ChordDiagramPlot: PASS - link border color", {


  p <- SCpubr::do_ChordDiagramPlot(sample = sample,
                                   from = "seurat_clusters",
                                   to = "orig.ident",
                                   link.border.color = "black")
  testthat::expect_type(p, "list")
})

testthat::test_that("do_ChordDiagramPlot: PASS - alignment", {


  p <- SCpubr::do_ChordDiagramPlot(sample = sample,
                                   from = "seurat_clusters",
                                   to = "orig.ident",
                                   alignment = "vertical")
  testthat::expect_type(p, "list")

  p <- SCpubr::do_ChordDiagramPlot(sample = sample,
                                   from = "seurat_clusters",
                                   to = "orig.ident",
                                   alignment = "horizontal")
  testthat::expect_type(p, "list")
})

testthat::test_that("do_ChordDiagramPlot: PASS - highlight group", {


  p <- SCpubr::do_ChordDiagramPlot(sample = sample,
                                   from = "seurat_clusters",
                                   to = "orig.ident",
                                   highlight_group = "0")
  testthat::expect_type(p, "list")
})

testthat::test_that("do_ChordDiagramPlot: FAILS", {


  testthat::expect_error({SCpubr::do_ChordDiagramPlot(sample = sample,
                                                      from = "seurat_clusters",
                                                      to = "orig.ident",
                                                      alignment = "wrong")})

  testthat::expect_error({SCpubr::do_ChordDiagramPlot(sample = sample,
                                                      from = "seurat_clusters",
                                                      to = "orig.ident",
                                                      link.arr.type = "wrong")})

  testthat::expect_error({SCpubr::do_ChordDiagramPlot(sample = sample,
                                                      from = "seurat_clusters",
                                                      to = "orig.ident",
                                                      highlight_group = "0",
                                                      alpha.highlight = 120)})

  testthat::expect_error({SCpubr::do_ChordDiagramPlot(sample = NULL,
                                                      from = "seurat_clusters",
                                                      to = "orig.ident")})

  testthat::expect_error({SCpubr::do_ChordDiagramPlot(sample = sample,
                                                      from = NULL,
                                                      to = "orig.ident")})

  testthat::expect_error({SCpubr::do_ChordDiagramPlot(sample = sample,
                                                      from = "seurat_clusters",
                                                      to = NULL)})

  testthat::expect_error({SCpubr::do_ChordDiagramPlot(sample = sample,
                                                      from = "seurat_clusters",
                                                      to = "orig.ident",
                                                      directional = 4)})

  testthat::expect_error({SCpubr::do_ChordDiagramPlot(sample = sample,
                                                      from = "seurat_clusters",
                                                      to = "orig.ident",
                                                      direction.type = "wrong")})

  testthat::expect_error({SCpubr::do_ChordDiagramPlot(sample = sample,
                                                      from = "seurat_clusters",
                                                      to = "orig.ident",
                                                      self.link = 3 )})

  testthat::expect_error({SCpubr::do_ChordDiagramPlot(sample = sample,
                                                      from = "CD14",
                                                      to = "orig.ident",
                                                      direction.type = "wrong")})

  testthat::expect_error({SCpubr::do_ChordDiagramPlot(sample = sample,
                                                      from = "seurat_clusters",
                                                      to = "CD14",
                                                      direction.type = "wrong")})

  testthat::expect_error({SCpubr::do_ChordDiagramPlot(sample = sample,
                                                      from = "nCount_RNA",
                                                      to = "orig.ident",
                                                      direction.type = "wrong")})

  testthat::expect_error({SCpubr::do_ChordDiagramPlot(sample = sample,
                                                      from = "seurat_clusters",
                                                      to = "nCount_RNA",
                                                      direction.type = "wrong")})
})


testthat::test_that("do_ChordDiagramPlot: internal_use", {
  liana_output <- liana_output %>%
                  liana::liana_aggregate(verbose = FALSE) %>%
                  dplyr::mutate(magnitude = .data$sca.LRscore) %>%
                  dplyr::mutate(specificity = .data$natmi.edge_specificity) %>%
                  dplyr::arrange(dplyr::desc(.data$specificity), dplyr::desc(.data$magnitude)) %>%
                  tidyr::unite(c("ligand.complex", "receptor.complex"),
                               col = "interaction",
                               sep = " | ",
                               remove = FALSE) %>%
                  # Merge source and target column into one, for future filtering.
                  tidyr::unite(c("source", "target"),
                               col = "interacting_clusters",
                               remove = FALSE)
  output_copy <- liana_output %>% dplyr::filter(.data$aggregate_rank <= 0.05)
  liana_output <- liana_output %>%
                  # Filter based on the top X interactions of ascending sensibilities.
                  dplyr::inner_join(y = {liana_output %>%
                                         dplyr::distinct_at(dplyr::all_of(c("ligand.complex", "receptor.complex"))) %>%
                                         dplyr::slice_head(n = 25)},
                                         by = c("ligand.complex", "receptor.complex"))

  testthat::expect_error({SCpubr::do_ChordDiagramPlot(sample = sample,
                                                      from = "seurat_clusters",
                                                      to = "orig.ident",
                                                      from_df = TRUE,
                                                      df = 3)})

  testthat::expect_error({SCpubr::do_ChordDiagramPlot(sample = sample,
                                                      from = "seurat_clusters",
                                                      to = "orig.ident",
                                                      from_df = TRUE,
                                                      df = as.data.frame(liana_output))})

  data <- output_copy %>%
          dplyr::select(dplyr::all_of(c("source", "target"))) %>%
          dplyr::group_by(.data$target, .data$source) %>%
          dplyr::summarise(value = dplyr::n()) %>%
          dplyr::rename("from" = .data[["source"]],
                        "to" = .data[["target"]]) %>%
          dplyr::select(dplyr::all_of(c("from", "to", "value")))
  data.test <- data
  colnames(data.test) <- c("A", "B", "C")
  testthat::expect_error({SCpubr::do_ChordDiagramPlot(sample = sample,
                                                      from = "seurat_clusters",
                                                      to = "orig.ident",
                                                      from_df = TRUE,
                                                      df = data.test)})

  data.test <- data
  data.test$from <- 3
  testthat::expect_error({SCpubr::do_ChordDiagramPlot(sample = sample,
                                                      from = "seurat_clusters",
                                                      to = "orig.ident",
                                                      from_df = TRUE,
                                                      df = data.test)})

  data.test <- data
  data.test$to <- 3
  testthat::expect_error({SCpubr::do_ChordDiagramPlot(sample = sample,
                                                      from = "seurat_clusters",
                                                      to = "orig.ident",
                                                      from_df = TRUE,
                                                      df = data.test)})

  data.test <- data
  data.test$value <- "B"
  testthat::expect_error({SCpubr::do_ChordDiagramPlot(sample = sample,
                                                      from = "seurat_clusters",
                                                      to = "orig.ident",
                                                      from_df = TRUE,
                                                      df = data.test)})


  testthat::expect_error({SCpubr::do_ChordDiagramPlot(sample = sample,
                                                      from = "seurat_clusters",
                                                      to = "orig.ident",
                                                      link.arr.type = "wrong")})

  testthat::expect_error({SCpubr::do_ChordDiagramPlot(sample = sample,
                                                      from = "seurat_clusters",
                                                      to = "orig.ident",
                                                      highlight_group = "0",
                                                      alpha.highlight = 120)})
})
