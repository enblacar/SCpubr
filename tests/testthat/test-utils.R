# CHECK SUGGESTS
testthat::test_that("utils: check_suggests - FAIL - Wrong function", {
  testthat::expect_error(check_suggests("wrong_name"))
})

testthat::test_that("utils: check_suggests - PASS - Correct function", {
  testthat::expect_silent(check_suggests("do_DimPlot"))
})


# STATE DEPENDENCIES

testthat::test_that("utils: state_dependencies - FAIL - Wrong function", {
  testthat::expect_error(state_dependencies("wrong_name"))
})

testthat::test_that("utils: state_dependencies - PASS - Correct function, one name", {
  suppressMessages({testthat::expect_message(state_dependencies("do_DimPlot"))})
})

testthat::test_that("utils: state_dependencies - PASS - Correct function, several names", {
  suppressMessages({testthat::expect_message(state_dependencies(c("do_DimPlot", "do_FeaturePlot")))})
})

testthat::test_that("utils: state_dependencies - PASS - Correct function, no parameters provided", {
  suppressMessages({testthat::expect_message(state_dependencies())})
})

# CHECK SEURAT

# CHECK SUGGESTS
testthat::test_that("utils: check_Seurat - FAIL - Not Seurat object", {
  testthat::expect_error(check_Seurat("not a Seurat object"))
})

testthat::test_that("utils: check_suggests - PASS - Seurat object", {
  testthat::expect_silent(check_Seurat(sample))
})


# CHECK COLORS
testthat::test_that("utils: check_colors - FAIL - wrong color", {
  testthat::expect_error(check_colors("not_a_color"))

})

testthat::test_that("utils: check_colors - FAIL - wrong color in a vector of colors", {
  testthat::expect_error(check_colors(c("not_a_color", "red", "blue")))
})

testthat::test_that("utils: check_colors - PASS - One color", {
  testthat::expect_silent(check_colors("red"))
})

testthat::test_that("utils: check_colors - PASS - Several colors", {
  testthat::expect_silent(check_colors(c("red", "blue")))
})




# CHECK CONSISTENCY COLORS AND NAMES

testthat::test_that("utils: check_consistency_colors_and_names - FAIL - more colors provided", {
  testthat::expect_error(check_consistency_colors_and_names(sample = sample,
                                                            colors = c("a" = "red", "b" = "blue"),
                                                            grouping_variable = "orig.ident"))
})

testthat::test_that("utils: check_consistency_colors_and_names - FAIL - names of colors not matching", {
  testthat::expect_error(check_consistency_colors_and_names(sample = sample,
                                                            colors = c("a" = "red"),
                                                            grouping_variable = "orig.ident"))
})

testthat::test_that("utils: check_consistency_colors_and_names - FAIL - less colors provided", {
  testthat::expect_error(check_consistency_colors_and_names(sample = sample,
                                                            colors = c("a" = "red"),
                                                            grouping_variable = "seurat_clusters"))
})

testthat::test_that("utils: check_consistency_colors_and_names - PASS - Colors matching", {
  testthat::expect_silent(check_consistency_colors_and_names(sample = sample,
                                                             colors = c("0" = "red",
                                                                        "1" = "red",
                                                                        "2" = "red",
                                                                        "3" = "red",
                                                                        "4" = "red",
                                                                        "5" = "red",
                                                                        "6" = "red",
                                                                        "7" = "red",
                                                                        "8" = "red")))
})

testthat::test_that("utils: check_consistency_colors_and_names - PASS - Colors matching, grouping variable", {
  testthat::expect_silent(check_consistency_colors_and_names(sample = sample,
                                                             colors = c("Cell" = "red"),
                                                             grouping_variable = "orig.ident"))
})


# GENERATE COLOR SCALE
testthat::test_that("utils: generate_color_scale - PASS - equal length of output", {
  names_use <- c("a", "b", "c")
  colors <- colortools::setColors("#457b9d", length(names_use))
  testthat::expect_length(colors, length(names_use))
})

# COMPUTE SCALES LIMITS

testthat::test_that("utils: compute_scale_limits - PASS - using a gene", {
  output <- compute_scale_limits(sample = sample,
                                 feature = "CD14")
  testthat::expect_length(output, 2)
})

testthat::test_that("utils: compute_scale_limits - PASS - using a metadata variable", {
  output <- compute_scale_limits(sample = sample,
                                 feature = "orig.ident")
  testthat::expect_length(output, 2)
})

testthat::test_that("utils: compute_scale_limits - PASS - using dimensional reduction variable", {
  output <- compute_scale_limits(sample = sample,
                                 feature = "PC_1")
  testthat::expect_length(output, 2)
})

# CHECK FEATURE

testthat::test_that("utils: check_feature - FAIL - using the wrong gene", {
  testthat::expect_error(check_feature(sample = sample,
                                       features = "NOTCD14"))
})

testthat::test_that("utils: check_feature - FAIL - using the wrong metadata", {
  testthat::expect_error(check_feature(sample = sample,
                                       features = "oris.ident"))
})

testthat::test_that("utils: check_feature - FAIL - using the wrong dimensional reduction variable", {
  testthat::expect_error(check_feature(sample = sample,
                                       features = "UMAP_38"))
})

testthat::test_that("utils: check_feature - FAIL - all features failing while in permissive mode", {
  testthat::expect_error(check_feature(sample = sample,
                                       features = c("NOTCD14", "UMAP_38"),
                                       permissive = TRUE))
})

testthat::test_that("utils: check_feature - WARNING - using one wrong gene and one good", {
  testthat::expect_warning(check_feature(sample = sample,
                                         features = c("NOTCD14", "CD14"),
                                         permissive = TRUE))
})

testthat::test_that("utils: check_feature - WARNING - using one wrong metadata variable and one good", {
  testthat::expect_warning(check_feature(sample = sample,
                                         features = c("oris.ident", "orig.ident"),
                                         permissive = TRUE))
})

testthat::test_that("utils: check_feature - WARNING - using one wrong dimensional reduction variable and one good", {
  testthat::expect_warning(check_feature(sample = sample,
                                         features = c("UMAP_38", "PC_1"),
                                         permissive = TRUE))
})

testthat::test_that("utils: check_feature - PASS - dump reduction names", {
  dim_names <- check_feature(sample = sample,
                             features = c("PC_1"),
                             dump_reduction_names = TRUE)
  expected_output <- 0
  for (dim_red in names(sample@reductions)){
    expected_output <- expected_output + length(colnames(sample@reductions[[dim_red]][[]]))
  }
  testthat::expect_length(dim_names, expected_output)
})

testthat::test_that("utils: check_feature - PASS - permissive check length of output", {
  testthat::expect_warning({
    features <- check_feature(sample = sample,
                              features = c("PC_1", "PC_99"),
                              permissive = TRUE)
    testthat::expect_length(features, 1)
  })
})

testthat::test_that("utils: check_feature - PASS - permissive check length of output when both permissive and dump_reduction_names are present.", {
  output <- check_feature(sample = sample,
                          features = c("PC_1"),
                          dump_reduction_names = TRUE,
                          permissive = TRUE)
  testthat::expect_length(output, 2)
})

testthat::test_that("utils: check_feature - ERROR - using the wrong enforcer", {
  testthat::expect_error(check_feature(sample = sample,
                                       features = c("CD14"),
                                       enforce_check = TRUE,
                                       enforce_parameter = "Gene"))
})

testthat::test_that("utils: check_feature - ERROR - using the wrong feature for the selected enforcer", {
  testthat::expect_error(check_feature(sample = sample,
                                       features = c("CD14"),
                                       enforce_check = TRUE,
                                       enforce_parameter = "metadata"))
})


# REMOVE NOT FOUND FEATURES
testthat::test_that("utils: remove_not_found_features - PASS - 0 features removed - character", {
  features <- c("a", "b")
  not_found_features <- ""
  output <- remove_not_found_features(features = features, not_found_features = not_found_features)
  testthat::expect_length(output, 2)
  testthat::expect_type(output, "character")
})

testthat::test_that("utils: remove_not_found_features - PASS - 1 features removed - character", {
  features <- c("a", "b")
  not_found_features <- "a"
  output <- remove_not_found_features(features = features, not_found_features = not_found_features)
  testthat::expect_length(output, 1)
  testthat::expect_type(output, "character")
})

testthat::test_that("utils: remove_not_found_features - PASS - 2 features removed - character", {
  features <- c("a", "b")
  not_found_features <- c("a", "b")
  output <- remove_not_found_features(features = features, not_found_features = not_found_features)
  testthat::expect_length(output, 0)
  testthat::expect_type(output, "character")
})

testthat::test_that("utils: remove_not_found_features - PASS - 0 features removed - list", {
  features <- list("A" = c("a"),
                   "B" = c("b"))
  not_found_features <- ""
  output <- remove_not_found_features(features = features, not_found_features = not_found_features)
  testthat::expect_length(output$A, 1)
  testthat::expect_length(output$B, 1)
  testthat::expect_type(output, "list")
})

testthat::test_that("utils: remove_not_found_features - PASS - 1 features removed - list", {
  features <- list("A" = c("a"),
                   "B" = c("b"))
  not_found_features <- "a"
  output <- remove_not_found_features(features = features, not_found_features = not_found_features)
  testthat::expect_length(output$A, 0)
  testthat::expect_length(output$B, 1)
  testthat::expect_type(output, "list")
})

testthat::test_that("utils: remove_not_found_features - PASS - 2 features removed - list", {
  features <- list("A" = c("a"),
                   "B" = c("b"))
  not_found_features <- c("a", "b")
  output <- remove_not_found_features(features = features, not_found_features = not_found_features)
  testthat::expect_length(output$A, 0)
  testthat::expect_length(output$B, 0)
  testthat::expect_type(output, "list")
})


# REMOVE DUPLICATED FEATURES

testthat::test_that("utils: remove_duplicated_features - WARNING - having duplicated features - character", {
  features <- c("a", "a")
  testthat::expect_warning(remove_duplicated_features(features))
  output <- suppressWarnings({remove_duplicated_features(features)})
  testthat::expect_type(output, "character")
})

testthat::test_that("utils: remove_duplicated_features - WARNING - having duplicated features across lists - list", {
  features <- list("A" = c("a"),
                   "B" = c("a"))
  testthat::expect_warning(remove_duplicated_features(features))
  output <- suppressWarnings({remove_duplicated_features(features)})
  testthat::expect_type(output, "list")
})

testthat::test_that("utils: remove_duplicated_features - WARNING - having duplicated features within lists - list", {
  features <- list("A" = c("a", "a"),
                   "B" = c("b"))
  testthat::expect_warning(remove_duplicated_features(features))
  output <- suppressWarnings({remove_duplicated_features(features)})
  testthat::expect_type(output, "list")
})

testthat::test_that("utils: remove_duplicated_features - WARNING - having duplicated features across and between lists - list", {
  features <- list("A" = c("a", "a"),
                   "B" = c("a"))
  suppressWarnings({testthat::expect_warning(remove_duplicated_features(features))})
  output <- suppressWarnings({remove_duplicated_features(features)})
  testthat::expect_type(output, "list")

})
