if (isFALSE(dep_check[["utils"]])){
  # CHECK SUGGESTS
  testthat::test_that("utils: check_suggests - FAIL - Wrong function", {
    testthat::expect_error(SCpubr:::check_suggests("wrong_name"))
  })

  testthat::test_that("utils: check_suggests - FAIL - Package not installed", {
    testthat::expect_error(SCpubr:::check_suggests("testing"))
  })

  testthat::test_that("utils: check_suggests - PASS - Correct function", {
    testthat::expect_silent(SCpubr:::check_suggests("do_DimPlot"))
  })


  # STATE DEPENDENCIES

  testthat::test_that("utils: state_dependencies - FAIL - Wrong function", {
    testthat::expect_error(SCpubr::state_dependencies("wrong_name"))
  })

  testthat::test_that("utils: state_dependencies - PASS - Correct function, one name", {
    suppressMessages({testthat::expect_message(SCpubr::state_dependencies("do_DimPlot"))})
  })

  testthat::test_that("utils: state_dependencies - PASS - Correct function, several names", {
    suppressMessages({testthat::expect_message(SCpubr::state_dependencies(c("do_DimPlot", "do_FeaturePlot")))})
  })

  testthat::test_that("utils: state_dependencies - PASS - Correct function, no parameters provided", {
    suppressMessages({testthat::expect_message(SCpubr::state_dependencies())})
  })

  # CHECK SEURAT

  # CHECK SUGGESTS
  testthat::test_that("utils: check_Seurat - FAIL - Not Seurat object", {

    testthat::expect_error(SCpubr:::check_Seurat("not a Seurat object"))
  })

  testthat::test_that("utils: check_suggests - PASS - Seurat object", {



    testthat::expect_silent(SCpubr:::check_Seurat(sample))
  })


  # CHECK COLORS
  testthat::test_that("utils: check_colors - FAIL - wrong color", {

    testthat::expect_error(SCpubr:::check_colors("not_a_color"))

  })

  testthat::test_that("utils: check_colors - FAIL - wrong color in a vector of colors", {

    testthat::expect_error(SCpubr:::check_colors(c("not_a_color", "red", "blue")))
  })

  testthat::test_that("utils: check_colors - PASS - One color", {

    testthat::expect_silent(SCpubr:::check_colors("red"))
  })

  testthat::test_that("utils: check_colors - PASS - Several colors", {

    testthat::expect_silent(SCpubr:::check_colors(c("red", "blue")))
  })




  # CHECK CONSISTENCY COLORS AND NAMES

  testthat::test_that("utils: check_consistency_colors_and_names - FAIL - more colors provided", {



    testthat::expect_error(SCpubr:::check_consistency_colors_and_names(sample = sample,
                                                                       colors = c("a" = "red", "b" = "blue"),
                                                                       grouping_variable = "orig.ident"))
  })

  testthat::test_that("utils: check_consistency_colors_and_names - FAIL - names of colors not matching", {



    testthat::expect_error(SCpubr:::check_consistency_colors_and_names(sample = sample,
                                                                       colors = c("a" = "red"),
                                                                       grouping_variable = "orig.ident"))
  })

  testthat::test_that("utils: check_consistency_colors_and_names - FAIL - less colors provided", {



    testthat::expect_error(SCpubr:::check_consistency_colors_and_names(sample = sample,
                                                                       colors = c("1" = "red"),
                                                                       grouping_variable = "seurat_clusters"))
  })

  testthat::test_that("utils: check_consistency_colors_and_names - PASS - Colors matching", {



    testthat::expect_silent(SCpubr:::check_consistency_colors_and_names(sample = sample,
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



    testthat::expect_silent(SCpubr:::check_consistency_colors_and_names(sample = sample,
                                                                        colors = c("Cell" = "red"),
                                                                        grouping_variable = "orig.ident"))
  })


  # GENERATE COLOR SCALE
  testthat::test_that("utils: generate_color_scale - PASS - equal length of output", {

    names_use <- c("a", "b", "c")
    colors <- colorspace::qualitative_hcl(length(names_use), palette = "Dark 3")
    testthat::expect_length(colors, length(names_use))
  })

  # COMPUTE SCALES LIMITS

  testthat::test_that("utils: compute_scale_limits - PASS - using a gene", {



    output <- SCpubr:::compute_scale_limits(sample = sample,
                                            feature = "EPC1")
    testthat::expect_length(output, 2)
  })

  testthat::test_that("utils: compute_scale_limits - PASS - using a metadata variable", {



    output <- SCpubr:::compute_scale_limits(sample = sample,
                                            feature = "orig.ident")
    testthat::expect_length(output, 2)
  })

  testthat::test_that("utils: compute_scale_limits - PASS - using dimensional reduction variable", {



    output <- SCpubr:::compute_scale_limits(sample = sample,
                                            feature = "PC_1")
    testthat::expect_length(output, 2)
  })

  # CHECK FEATURE

  testthat::test_that("utils: check_feature - FAIL - using the wrong gene", {



    testthat::expect_error(check_feature(sample = sample,
                                         features = "NOTEPC1"))
  })

  testthat::test_that("utils: check_feature - FAIL - using the wrong metadata", {



    testthat::expect_error(SCpubr:::check_feature(sample = sample,
                                                  features = "oris.ident"))
  })

  testthat::test_that("utils: check_feature - FAIL - using the wrong dimensional reduction variable", {



    testthat::expect_error(SCpubr:::check_feature(sample = sample,
                                                  features = "UMAP_38"))
  })

  testthat::test_that("utils: check_feature - FAIL - all features failing while in permissive mode", {



    testthat::expect_error(SCpubr:::check_feature(sample = sample,
                                                  features = c("NOTEPC1", "UMAP_38"),
                                                  permissive = TRUE))
  })

  testthat::test_that("utils: check_feature - WARNING - using one wrong gene and one good", {



    testthat::expect_warning(SCpubr:::check_feature(sample = sample,
                                                    features = c("NOTEPC1", "EPC1"),
                                                    permissive = TRUE))
  })

  testthat::test_that("utils: check_feature - WARNING - using one wrong metadata variable and one good", {



    testthat::expect_warning(SCpubr:::check_feature(sample = sample,
                                                    features = c("oris.ident", "orig.ident"),
                                                    permissive = TRUE))
  })

  testthat::test_that("utils: check_feature - WARNING - using one wrong dimensional reduction variable and one good", {



    testthat::expect_warning(SCpubr:::check_feature(sample = sample,
                                                    features = c("UMAP_38", "PC_1"),
                                                    permissive = TRUE))
  })

  testthat::test_that("utils: check_feature - PASS - dump reduction names", {



    dim_names <- SCpubr:::check_feature(sample = sample,
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
      features <- SCpubr:::check_feature(sample = sample,
                                         features = c("PC_1", "PC_99"),
                                         permissive = TRUE)
      testthat::expect_length(features, 1)
    })
  })

  testthat::test_that("utils: check_feature - PASS - permissive check length of output when both permissive and dump_reduction_names are present.", {



    output <- SCpubr:::check_feature(sample = sample,
                                     features = c("PC_1"),
                                     dump_reduction_names = TRUE,
                                     permissive = TRUE)
    testthat::expect_length(output, 2)
  })

  testthat::test_that("utils: check_feature - ERROR - using the wrong enforcer", {



    testthat::expect_error(SCpubr:::check_feature(sample = sample,
                                                  features = c("EPC1"),
                                                  enforce_check = "Gene",
                                                  enforce_parameter = "group.by"))
  })

  testthat::test_that("utils: check_feature - ERROR - using the wrong feature for the selected enforcer", {



    testthat::expect_error(SCpubr:::check_feature(sample = sample,
                                                  features = c("EPC1"),
                                                  enforce_check = "reductions",
                                                  enforce_parameter = "group.by"))
  })


  # REMOVE NOT FOUND FEATURES
  testthat::test_that("utils: remove_not_found_features - PASS - 0 features removed - character", {

    features <- c("a", "b")
    not_found_features <- ""
    output <- SCpubr:::remove_not_found_features(features = features, not_found_features = not_found_features)
    testthat::expect_length(output, 2)
    testthat::expect_type(output, "character")
  })

  testthat::test_that("utils: remove_not_found_features - PASS - 1 features removed - character", {

    features <- c("a", "b")
    not_found_features <- "a"
    output <- SCpubr:::remove_not_found_features(features = features, not_found_features = not_found_features)
    testthat::expect_length(output, 1)
    testthat::expect_type(output, "character")
  })

  testthat::test_that("utils: remove_not_found_features - PASS - 2 features removed - character", {

    features <- c("a", "b")
    not_found_features <- c("a", "b")
    output <- SCpubr:::remove_not_found_features(features = features, not_found_features = not_found_features)
    testthat::expect_length(output, 0)
    testthat::expect_type(output, "character")
  })

  testthat::test_that("utils: remove_not_found_features - PASS - 0 features removed - list", {

    features <- list("A" = c("a"),
                     "B" = c("b"))
    not_found_features <- ""
    output <- SCpubr:::remove_not_found_features(features = features, not_found_features = not_found_features)
    testthat::expect_length(output$A, 1)
    testthat::expect_length(output$B, 1)
    testthat::expect_type(output, "list")
  })

  testthat::test_that("utils: remove_not_found_features - PASS - 1 features removed - list", {

    features <- list("A" = c("a"),
                     "B" = c("b"))
    not_found_features <- "a"
    output <- SCpubr:::remove_not_found_features(features = features, not_found_features = not_found_features)
    testthat::expect_length(output$A, 0)
    testthat::expect_length(output$B, 1)
    testthat::expect_type(output, "list")
  })

  testthat::test_that("utils: remove_not_found_features - PASS - 2 features removed - list", {

    features <- list("A" = c("a"),
                     "B" = c("b"))
    not_found_features <- c("a", "b")
    output <- SCpubr:::remove_not_found_features(features = features, not_found_features = not_found_features)
    testthat::expect_length(output$A, 0)
    testthat::expect_length(output$B, 0)
    testthat::expect_type(output, "list")
  })


  # REMOVE DUPLICATED FEATURES

  testthat::test_that("utils: remove_duplicated_features - WARNING - having duplicated features - character", {

    features <- c("a", "a")
    testthat::expect_warning(SCpubr:::remove_duplicated_features(features))
    output <- suppressWarnings({SCpubr:::remove_duplicated_features(features)})
    testthat::expect_type(output, "character")
  })

  testthat::test_that("utils: remove_duplicated_features - WARNING - having duplicated features across lists - list", {

    features <- list("A" = c("a"),
                     "B" = c("a"))
    testthat::expect_warning(SCpubr:::remove_duplicated_features(features))
    output <- suppressWarnings({SCpubr:::remove_duplicated_features(features)})
    testthat::expect_type(output, "list")
  })

  testthat::test_that("utils: remove_duplicated_features - WARNING - having duplicated features within lists - list", {

    features <- list("A" = c("a", "a"),
                     "B" = c("b"))
    testthat::expect_warning(SCpubr:::remove_duplicated_features(features))
    output <- suppressWarnings({SCpubr:::remove_duplicated_features(features)})
    testthat::expect_type(output, "list")
  })

  testthat::test_that("utils: remove_duplicated_features - WARNING - having duplicated features across and between lists - list", {

    features <- list("A" = c("a", "a"),
                     "B" = c("a"))
    suppressWarnings({testthat::expect_warning(SCpubr:::remove_duplicated_features(features))})
    output <- suppressWarnings({SCpubr:::remove_duplicated_features(features)})
    testthat::expect_type(output, "list")

  })


  # CHECK IDENTITY

  testthat::test_that("utils: check_identity - FAIL - wrong identity", {



    testthat::expect_error(SCpubr:::check_identity(sample, "wrong_identity"))
  })

  testthat::test_that("utils: check_identity - PASS - right identity", {



    testthat::expect_silent(SCpubr:::check_identity(sample, "0"))
  })


  # CHECK AND SET REDUCTION

  testthat::test_that("utils: check_and_set_reduction - FAIL - no reductions", {



    test <- sample
    test@reductions[["pca"]] <- NULL
    test@reductions[["umap"]] <- NULL
    testthat::expect_error(SCpubr:::check_and_set_reduction(sample = test, reduction = "umap"))
  })

  testthat::test_that("utils: check_and_set_reduction - FAIL - wrong reductions", {



    testthat::expect_error(SCpubr:::check_and_set_reduction(sample = sample, reduction = "wrong_reduction"))
  })

  testthat::test_that("utils: check_and_set_reduction - PASS - null reduction, check that the output is the last computed reduction", {



    output <- SCpubr:::check_and_set_reduction(sample = sample)
    last_reduction <- names(sample@reductions)[length(names(sample@reductions))]
    testthat::expect_identical(output, last_reduction)
  })

  testthat::test_that("utils: check_and_set_reduction - PASS - provide a reduction", {



    output <- SCpubr:::check_and_set_reduction(sample = sample, reduction = "umap")
    reduction_check <- "umap"
    testthat::expect_identical(output, reduction_check)
  })

  testthat::test_that("utils: check_and_set_reduction - PASS - umap not in reductions", {



    sample@reductions$umap <- NULL
    output <- SCpubr:::check_and_set_reduction(sample = sample)
    reduction_check <- "pca"
    testthat::expect_identical(output, reduction_check)
  })


  # CHECK AND SET DIMENSIONS

  testthat::test_that("utils: check_and_set_dimensions - FAIL - dims not being a pair of values", {



    testthat::expect_error(SCpubr:::check_and_set_dimensions(sample = sample, reduction = "umap", dims = "wrong_input"))
  })

  testthat::test_that("utils: check_and_set_dimensions - FAIL - dims not being a pair of integers", {



    testthat::expect_error(SCpubr:::check_and_set_dimensions(sample = sample, reduction = "umap", dims = c(1, "wrong_input")))
  })

  testthat::test_that("utils: check_and_set_dimensions - FAIL - dims not being in the available list of dims", {



    testthat::expect_error(SCpubr:::check_and_set_dimensions(sample = sample, reduction = "umap", dims = c(1, 20)))
  })

  testthat::test_that("utils: check_and_set_dimensions - FAIL - reduction only having 1 dim", {



    test <- sample
    obj <- Seurat::CreateDimReducObject(test@reductions$umap[[]][, "UMAP_1", drop = FALSE], key = "UMAP_", assay = "SCT")
    test@reductions$umap <- obj
    testthat::expect_error(SCpubr:::check_and_set_dimensions(sample = test, reduction = "umap", dims = c(1, 2)))
  })

  testthat::test_that("utils: check_and_set_dimensions - PASS - NULL parameters", {



    output <- SCpubr:::check_and_set_dimensions(sample = sample)
    testthat::expect_identical(output, c(1, 2))
  })

  testthat::test_that("utils: check_and_set_dimensions - PASS - NULL dimension but provided dims", {



    output <- SCpubr:::check_and_set_dimensions(sample = sample, dims = c(2, 1))
    testthat::expect_identical(output, c(2, 1))
  })

  testthat::test_that("utils: check_and_set_dimensions - PASS - provided dimension and dims", {



    output <- SCpubr:::check_and_set_dimensions(sample = sample, reduction = "pca", dims = c(20, 11))
    testthat::expect_identical(output, c(20, 11))
  })


  # CHECK AND SET ASSAY
  testthat::test_that("utils: check_and_set_assay - FAIL - wrong assay type", {



    testthat::expect_error(SCpubr:::check_and_set_assay(sample = sample, assay = FALSE))
  })

  testthat::test_that("utils: check_and_set_assay - FAIL - no assays in object", {



    test <- sample
    test@assays$RNA <- NULL
    test@assays$SCT <- NULL
    testthat::expect_error(SCpubr:::check_and_set_assay(sample = test))
  })

  testthat::test_that("utils: check_and_set_assay - FAIL - assay not present", {



    testthat::expect_error(SCpubr:::check_and_set_assay(sample = sample, assay = "ATAC"))
  })

  testthat::test_that("utils: check_and_set_assay - PASS - null parameters", {



    output <- SCpubr:::check_and_set_assay(sample = sample)
    testthat::expect_identical(output$assay, Seurat::DefaultAssay(sample))
  })

  testthat::test_that("utils: check_and_set_assay - PASS - providing assay", {



    output <- SCpubr:::check_and_set_assay(sample = sample, assay = "SCT")
    testthat::expect_identical(output$assay, "SCT")
  })

  testthat::test_that("utils: check_and_set_assay - PASS - providing non defaultassay", {



    sample@assays$RNA <- sample@assays$SCT
    output <- SCpubr:::check_and_set_assay(sample = sample, assay = "RNA")
    testthat::expect_identical(output$assay, "RNA")
  })


  # CHECK TYPE

  testthat::test_that("utils: check_type - FAIL - wrong type", {

    parameters <- c("first" = 1,
                    "second" = 2,
                    "third" = "a")
    testthat::expect_error(SCpubr:::check_type(parameters = parameters, required_type = "numeric", test_function = is.numeric))
  })

  testthat::test_that("utils: check_type - PASS - numeric", {

    parameters <- c("first" = 1,
                    "second" = 2)
    testthat::expect_silent(SCpubr:::check_type(parameters = parameters, required_type = "numeric", test_function = is.numeric))
  })

  testthat::test_that("utils: check_type - PASS - numeric with NULL", {

    parameters <- c("first" = 1,
                    "second" = 2,
                    "third" = NULL)
    testthat::expect_silent(SCpubr:::check_type(parameters = parameters, required_type = "numeric", test_function = is.numeric))
  })

  testthat::test_that("utils: check_type - PASS - character", {

    parameters <- c("first" = "a",
                    "second" = "b")
    testthat::expect_silent(SCpubr:::check_type(parameters = parameters, required_type = "character", test_function = is.character))
  })

  testthat::test_that("utils: check_type - PASS - character with NULL", {

    parameters <- c("first" = "a",
                    "second" = "b",
                    "third" = NULL)
    testthat::expect_silent(SCpubr:::check_type(parameters = parameters, required_type = "character", test_function = is.character))
  })

  testthat::test_that("utils: check_type - PASS - logical", {

    parameters <- c("first" = TRUE,
                    "second" = FALSE)
    testthat::expect_silent(SCpubr:::check_type(parameters = parameters, required_type = "logical", test_function = is.logical))
  })

  testthat::test_that("utils: check_type - PASS - logical with NULL", {

    parameters <- c("first" = TRUE,
                    "second" = FALSE,
                    "third" = NULL)
    testthat::expect_silent(SCpubr:::check_type(parameters = parameters, required_type = "logical", test_function = is.logical))
  })

  testthat::test_that("utils: check_type - PASS - list", {

    parameters <- c("first" = list(),
                    "second" = list())
    testthat::expect_silent(SCpubr:::check_type(parameters = parameters, required_type = "list", test_function = is.list))
  })

  testthat::test_that("utils: check_type - PASS - list with NULL", {

    parameters <- c("first" = list(),
                    "second" = list(),
                    "third" = NULL)
    testthat::expect_silent(SCpubr:::check_type(parameters = parameters, required_type = "list", test_function = is.list))
  })


  # CHECK AND SET THE SLOT

  testthat::test_that("utils: check_and_set_slot - FAIL - wrong slot", {

    testthat::expect_error(SCpubr:::check_and_set_slot("wrong_slot"))
  })

  testthat::test_that("utils: check_and_set_slot - PASS - counts", {

    output <- SCpubr:::check_and_set_slot("counts")
    testthat::expect_identical(output, "counts")
  })

  testthat::test_that("utils: check_and_set_slot - PASS - data", {

    output <- SCpubr:::check_and_set_slot("data")
    testthat::expect_identical(output, "data")
  })

  testthat::test_that("utils: check_and_set_slot - PASS - scale.data", {

    output <- SCpubr:::check_and_set_slot("scale.data")
    testthat::expect_identical(output, "scale.data")
  })


  # CHECK LIMITS
  testthat::test_that("utils: check_and_set_slot - FAIL - wrong limit", {



    testthat::expect_error(SCpubr:::check_limits(sample = sample, feature = "EPC1", value_name = "scale.end", value = 30))
  })

  testthat::test_that("utils: check_and_set_slot - PASS - good limit", {



    testthat::expect_silent(SCpubr:::check_limits(sample = sample, feature = "EPC1", value_name = "scale.end", value = 2))
  })


  # COMPUTE FACTOR LEVELS

  testthat::test_that("utils: compute_factor_levels - FAIL - wrong position", {



    testthat::expect_error(SCpubr:::compute_factor_levels(sample = sample, feature = "seurat_clusters", position = "upper"))
  })

  testthat::test_that("utils: compute_factor_levels - PASS - order.by and group.by", {

    testthat::expect_type(SCpubr:::compute_factor_levels(sample = sample,
                                                         feature = "seurat_clusters",
                                                         position = "fill",
                                                         group.by = NULL), "character")

    testthat::expect_type(SCpubr:::compute_factor_levels(sample = sample,
                                                         feature = "seurat_clusters",
                                                         position = "fill",
                                                         group.by = "orig.ident",
                                                         order = TRUE), "character")

    testthat::expect_type(SCpubr:::compute_factor_levels(sample = sample,
                                                         feature = "seurat_clusters",
                                                         position = "fill",
                                                         group.by = "orig.ident"), "character")

    testthat::expect_type(SCpubr:::compute_factor_levels(sample = sample,
                                                         feature = "EPC1",
                                                         position = "fill",
                                                         group.by = "orig.ident",
                                                         order = TRUE), "character")

    testthat::expect_type(SCpubr:::compute_factor_levels(sample = sample,
                                                         feature = "EPC1",
                                                         position = "fill",
                                                         group.by = "orig.ident"), "character")

    testthat::expect_type(SCpubr:::compute_factor_levels(sample = sample,
                                                         feature = "seurat_clusters",
                                                         position = "stack",
                                                         group.by = "orig.ident",
                                                         order = TRUE), "character")

    testthat::expect_type(SCpubr:::compute_factor_levels(sample = sample,
                                                         feature = "seurat_clusters",
                                                         position = "stack",
                                                         group.by = "orig.ident"), "character")

    testthat::expect_type(SCpubr:::compute_factor_levels(sample = sample,
                                                         feature = "EPC1",
                                                         position = "stack",
                                                         group.by = "orig.ident",
                                                         order = TRUE), "character")

    testthat::expect_type(SCpubr:::compute_factor_levels(sample = sample,
                                                         feature = "EPC1",
                                                         position = "stack",
                                                         group.by = "orig.ident"), "character")
  })



  # CHECK VIRIDIS COLOR MAP

  testthat::test_that("utils: check_viridis_color_map - FAIL - wrong color map", {

    testthat::expect_error(SCpubr:::check_viridis_color_map("wrong_color_map"))
  })

  testthat::test_that("utils: check_viridis_color_map - PASS - using turbo with verbose = F", {

    testthat::expect_silent(SCpubr:::check_viridis_color_map("turbo"))
  })


  # CHECK LENGTH

  testthat::test_that("utils: check_length - FAIL - distinct length", {

    vector_parameters <- c(1, 2)
    vector_features <- c(1)
    parameters_name <- "A"
    features_name <- "B"
    testthat::expect_error(SCpubr:::check_length(vector_of_parameters = vector_parameters,
                                                 vector_of_features = vector_features,
                                                 parameters_name = parameters_name,
                                                 features_name = features_name))
  })

  testthat::test_that("utils: check_length - PASS - correct length", {

    vector_parameters <- c(1, 2)
    vector_features <- c(1, 2)
    parameters_name <- "A"
    features_name <- "B"
    testthat::expect_silent(SCpubr:::check_length(vector_of_parameters = vector_parameters,
                                                  vector_of_features = vector_features,
                                                  parameters_name = parameters_name,
                                                  features_name = features_name))
  })


  # USE DATASET
  testthat::test_that("utils: use_dataset - PASS - checks", {
    testthat::skip_on_cran()
    output <- SCpubr:::use_dataset()
    testthat::expect_length(colnames(output), 180)
    testthat::expect_length(rownames(output), 7761)
    testthat::expect_identical(Seurat::Assays(output), c("RNA", "SCT"))
    testthat::expect_identical(Seurat::Reductions(output), c("pca", "umap"))
  })


  # ADD SCALE
  testthat::test_that("utils: add_scale - PASS - checks", {



    p <- do_FeaturePlot(sample, features = "EPC1")
    output <- SCpubr:::add_scale(p = p, scale = "color", function_use = ggplot2::scale_color_viridis_b())
    testthat::expect_true("ggplot" %in% class(output))
  })


  # COMPUTE BARPLOT ANNOTATION

  testthat::test_that("utils: compute_barplot_annotation - PASS - checks", {



    out <- SCpubr:::compute_barplot_annotation(sample = sample, group.by = "seurat_clusters", annotation = "orig.ident")
    testthat::expect_true("tbl" %in% class(out))
  })


  # HEATMAP INNER

  testthat::test_that("utils: heatmap_inner - PASS - checks", {


    data <- data.frame("A" = c(0.0012, 0.0012, 0.0013, 0.0014),
                       "B" = c(0.0014, 0.0013, 0.0012, 0.0011))
    out <- SCpubr:::heatmap_inner(as.matrix(data), zeros_are_white = TRUE)
    testthat::expect_true("Legends" %in% class(out$legend))
    testthat::expect_true("Heatmap" %in% class(out$heatmap))

    out <- SCpubr:::heatmap_inner(as.matrix(data), zeros_are_white = FALSE)
    testthat::expect_true("Legends" %in% class(out$legend))
    testthat::expect_true("Heatmap" %in% class(out$heatmap))

    data <- data.frame("A" = c(0.0012, 0.0012, 0.0013, 0.0014),
                       "B" = c(0.0014, 0.0013, 0.0012, 0.0011))
    out <- SCpubr:::heatmap_inner(as.matrix(data), data_range = "only_pos", zeros_are_white = TRUE)
    testthat::expect_true("Legends" %in% class(out$legend))
    testthat::expect_true("Heatmap" %in% class(out$heatmap))

    out <- SCpubr:::heatmap_inner(as.matrix(data), data_range = "only_pos", zeros_are_white = FALSE)
    testthat::expect_true("Legends" %in% class(out$legend))
    testthat::expect_true("Heatmap" %in% class(out$heatmap))

    data <- data.frame("A" = c(-0.0012, -0.0012, -0.0013, -0.0014),
                       "B" = c(-0.0014, -0.0013, -0.0012, -0.0011))
    out <- SCpubr:::heatmap_inner(as.matrix(data), data_range = "only_neg", zeros_are_white = TRUE)
    testthat::expect_true("Legends" %in% class(out$legend))
    testthat::expect_true("Heatmap" %in% class(out$heatmap))

    out <- SCpubr:::heatmap_inner(as.matrix(data), data_range = "only_neg", zeros_are_white = FALSE)
    testthat::expect_true("Legends" %in% class(out$legend))
    testthat::expect_true("Heatmap" %in% class(out$heatmap))

    `%>%`<- purrr::`%>%`
    data <- as.matrix({
      sample@meta.data %>% dplyr::select(c(seurat_clusters)) %>% dplyr::group_by(seurat_clusters) %>% dplyr::summarise(n = dplyr::n()) %>% dplyr::mutate(test = seq(from = 1,
                                                                                                                                                                    to = 90,
                                                                                                                                                                    by = 10),
                                                                                                                                                         zeros = 0) %>% dplyr::select(-seurat_clusters)
    })
    testthat::expect_error(heatmap_inner(data, data_range = "both", range.data = 10))
    testthat::expect_error(heatmap_inner(matrix(c(1, 1))))
    out <- SCpubr:::heatmap_inner(data)
    testthat::expect_true("Legends" %in% class(out$legend))
    testthat::expect_true("Heatmap" %in% class(out$heatmap))

    out <- SCpubr:::heatmap_inner(data, symmetrical_scale = TRUE)
    testthat::expect_true("Legends" %in% class(out$legend))
    testthat::expect_true("Heatmap" %in% class(out$heatmap))

    out <- SCpubr:::heatmap_inner(data, range.data = c(-10, 10), data_range = "both", symmetrical_scale = TRUE)
    testthat::expect_true("Legends" %in% class(out$legend))
    testthat::expect_true("Heatmap" %in% class(out$heatmap))

    out <- SCpubr:::heatmap_inner(data,
                                  outlier.data = TRUE,
                                  outlier.down.label = "A",
                                  outlier.up.label = "B",
                                  range.data = c(0, 50),
                                  data_range = "both")
    testthat::expect_true("Legends" %in% class(out$legend))
    testthat::expect_true("Heatmap" %in% class(out$heatmap))

    out <- SCpubr:::heatmap_inner(data,
                                  outlier.data = TRUE,
                                  range.data = c(0, 50),
                                  data_range = "both")
    testthat::expect_true("Legends" %in% class(out$legend))
    testthat::expect_true("Heatmap" %in% class(out$heatmap))

    out <- SCpubr:::heatmap_inner(data,
                                  outlier.data = TRUE,
                                  outlier.down.label = "A",
                                  outlier.up.label = "B",
                                  range.data = 50,
                                  data_range = "only_pos")
    testthat::expect_true("Legends" %in% class(out$legend))
    testthat::expect_true("Heatmap" %in% class(out$heatmap))

    out <- SCpubr:::heatmap_inner(data,
                                  outlier.data = TRUE,
                                  range.data = 50,
                                  data_range = "only_pos")
    testthat::expect_true("Legends" %in% class(out$legend))
    testthat::expect_true("Heatmap" %in% class(out$heatmap))

    out <- SCpubr:::heatmap_inner(data,
                                  outlier.data = TRUE,
                                  range.data = 50,
                                  data_range = "only_pos",
                                  use_viridis = TRUE,
                                  zeros_are_white = FALSE,
                                  viridis_color_map = "G",
                                  viridis_direction = -1)
    testthat::expect_true("Legends" %in% class(out$legend))
    testthat::expect_true("Heatmap" %in% class(out$heatmap))

    out <- SCpubr:::heatmap_inner(data, colors.use = c("red", "yellow"))
    testthat::expect_true("Legends" %in% class(out$legend))
    testthat::expect_true("Heatmap" %in% class(out$heatmap))

    data.modified <- as.matrix({
      sample@meta.data %>% dplyr::select(c(seurat_clusters)) %>% dplyr::group_by(seurat_clusters) %>% dplyr::summarise(n = dplyr::n()) %>% dplyr::mutate(test = seq(from = -1,
                                                                                                                                                                    to = -90,
                                                                                                                                                                    by = -10),
                                                                                                                                                         zeros = 0,
                                                                                                                                                         n = -20) %>% dplyr::select(-seurat_clusters)
    })
    out <- SCpubr:::heatmap_inner(data.modified,
                                  outlier.data = TRUE,
                                  outlier.down.label = "A",
                                  outlier.up.label = "B",
                                  range.data = -50,
                                  data_range = "only_neg")
    testthat::expect_true("Legends" %in% class(out$legend))
    testthat::expect_true("Heatmap" %in% class(out$heatmap))

    out <- SCpubr:::heatmap_inner(data.modified,
                                  outlier.data = TRUE,
                                  range.data = -50,
                                  data_range = "only_neg")
    testthat::expect_true("Legends" %in% class(out$legend))
    testthat::expect_true("Heatmap" %in% class(out$heatmap))

    out <- SCpubr:::heatmap_inner(data.modified,
                                  outlier.data = TRUE,
                                  range.data = -50,
                                  data_range = "only_neg",
                                  zeros_are_white = TRUE)
    testthat::expect_true("Legends" %in% class(out$legend))
    testthat::expect_true("Heatmap" %in% class(out$heatmap))

    out <- SCpubr:::heatmap_inner(data.modified, data_range = "only_neg")
    testthat::expect_true("Legends" %in% class(out$legend))
    testthat::expect_true("Heatmap" %in% class(out$heatmap))

    out <- SCpubr:::heatmap_inner(data.modified, data_range = "only_neg", outlier.data = TRUE, range.data = -15)
    testthat::expect_true("Legends" %in% class(out$legend))
    testthat::expect_true("Heatmap" %in% class(out$heatmap))

    testthat::expect_error(heatmap_inner(data.modified, data_range = "only_pos"))
    testthat::expect_error(heatmap_inner(data, data_range = "only_neg"))

    out <- SCpubr:::heatmap_inner(data, data_range = "only_pos")
    testthat::expect_true("Legends" %in% class(out$legend))
    testthat::expect_true("Heatmap" %in% class(out$heatmap))

    out <- SCpubr:::heatmap_inner(data, data_range = "only_pos", outlier.data = TRUE, range.data = 15)
    testthat::expect_true("Legends" %in% class(out$legend))
    testthat::expect_true("Heatmap" %in% class(out$heatmap))

    out <- SCpubr:::heatmap_inner(data, outlier.data = TRUE, range.data = c(0, 15))
    testthat::expect_true("Legends" %in% class(out$legend))
    testthat::expect_true("Heatmap" %in% class(out$heatmap))

    out <- SCpubr:::heatmap_inner(data, column_title = "test")
    testthat::expect_true("Legends" %in% class(out$legend))
    testthat::expect_true("Heatmap" %in% class(out$heatmap))

    out <- SCpubr:::heatmap_inner(data, row_title = "test")
    testthat::expect_true("Legends" %in% class(out$legend))
    testthat::expect_true("Heatmap" %in% class(out$heatmap))

    out <- SCpubr:::heatmap_inner(data, row_names_side = "left")
    testthat::expect_true("Legends" %in% class(out$legend))
    testthat::expect_true("Heatmap" %in% class(out$heatmap))

    out <- SCpubr:::heatmap_inner(data, row_names_side = "right")
    testthat::expect_true("Legends" %in% class(out$legend))
    testthat::expect_true("Heatmap" %in% class(out$heatmap))

    out <- SCpubr:::heatmap_inner(data, column_names_side = "bottom")
    testthat::expect_true("Legends" %in% class(out$legend))
    testthat::expect_true("Heatmap" %in% class(out$heatmap))

    out <- SCpubr:::heatmap_inner(data, column_names_side = "top")
    testthat::expect_true("Legends" %in% class(out$legend))
    testthat::expect_true("Heatmap" %in% class(out$heatmap))

    out <- SCpubr:::heatmap_inner(data, cluster_columns = TRUE)
    testthat::expect_true("Legends" %in% class(out$legend))
    testthat::expect_true("Heatmap" %in% class(out$heatmap))

    out <- SCpubr:::heatmap_inner(data, cluster_columns = FALSE)
    testthat::expect_true("Legends" %in% class(out$legend))
    testthat::expect_true("Heatmap" %in% class(out$heatmap))

    out <- SCpubr:::heatmap_inner(data, cluster_rows = TRUE)
    testthat::expect_true("Legends" %in% class(out$legend))
    testthat::expect_true("Heatmap" %in% class(out$heatmap))

    out <- SCpubr:::heatmap_inner(data, cluster_rows = FALSE)
    testthat::expect_true("Legends" %in% class(out$legend))
    testthat::expect_true("Heatmap" %in% class(out$heatmap))

    out <- SCpubr:::heatmap_inner(data, border = TRUE)
    testthat::expect_true("Legends" %in% class(out$legend))
    testthat::expect_true("Heatmap" %in% class(out$heatmap))

    out <- SCpubr:::heatmap_inner(data, border = FALSE)
    testthat::expect_true("Legends" %in% class(out$legend))
    testthat::expect_true("Heatmap" %in% class(out$heatmap))

    out <- SCpubr:::heatmap_inner(data, row_dendogram = TRUE)
    testthat::expect_true("Legends" %in% class(out$legend))
    testthat::expect_true("Heatmap" %in% class(out$heatmap))

    out <- SCpubr:::heatmap_inner(data, row_dendogram = FALSE)
    testthat::expect_true("Legends" %in% class(out$legend))
    testthat::expect_true("Heatmap" %in% class(out$heatmap))

    out <- SCpubr:::heatmap_inner(data, column_dendogram = TRUE)
    testthat::expect_true("Legends" %in% class(out$legend))
    testthat::expect_true("Heatmap" %in% class(out$heatmap))

    out <- SCpubr:::heatmap_inner(data, column_dendogram = FALSE)
    testthat::expect_true("Legends" %in% class(out$legend))
    testthat::expect_true("Heatmap" %in% class(out$heatmap))

    out <- SCpubr:::heatmap_inner(data, column_title_side = "top")
    testthat::expect_true("Legends" %in% class(out$legend))
    testthat::expect_true("Heatmap" %in% class(out$heatmap))

    out <- SCpubr:::heatmap_inner(data, column_title_side = "bottom")
    testthat::expect_true("Legends" %in% class(out$legend))
    testthat::expect_true("Heatmap" %in% class(out$heatmap))

    out <- SCpubr:::heatmap_inner(data, row_title_side = "left")
    testthat::expect_true("Legends" %in% class(out$legend))
    testthat::expect_true("Heatmap" %in% class(out$heatmap))

    out <- SCpubr:::heatmap_inner(data, row_title_side = "right")
    testthat::expect_true("Legends" %in% class(out$legend))
    testthat::expect_true("Heatmap" %in% class(out$heatmap))

    out <- SCpubr:::heatmap_inner(data, outlier.data = TRUE, range.data = c(0, 15))
    testthat::expect_true("Legends" %in% class(out$legend))
    testthat::expect_true("Heatmap" %in% class(out$heatmap))

    sample <- SCpubr:::compute_enrichment_scores(sample, input_gene_list = "EPC1", nbin = 1, ctrl = 10)
    data <- as.matrix({
      sample@meta.data %>% dplyr::select(c(orig.ident, Input)) %>% dplyr::group_by(orig.ident) %>% dplyr::summarise(n = mean(Input)) %>% dplyr::pull(n)
    })
    testthat::expect_error({SCpubr:::heatmap_inner(data)})

    data <- as.matrix({
      sample@meta.data %>% dplyr::select(c(seurat_clusters)) %>% dplyr::group_by(seurat_clusters) %>% dplyr::summarise(n = dplyr::n()) %>% dplyr::pull(n)
    })

    data[2, ] <- 22
    obj <- ComplexHeatmap::HeatmapAnnotation(orig.ident = data, col = list(orig.ident = c("20" = "red", "22" = "green")), which = "row")
    out <- SCpubr:::heatmap_inner(data, row_annotation = obj, row_annotation_side = "right")
    testthat::expect_true("Legends" %in% class(out$legend))
    testthat::expect_true("Heatmap" %in% class(out$heatmap))

    out <- SCpubr:::heatmap_inner(data, row_annotation = obj, row_annotation_side = "left")
    testthat::expect_true("Legends" %in% class(out$legend))
    testthat::expect_true("Heatmap" %in% class(out$heatmap))

    obj <- ComplexHeatmap::HeatmapAnnotation(orig.ident = data, col = list(orig.ident = c("20" = "red")), which = "column")
    out <- SCpubr:::heatmap_inner(t(data), column_annotation = obj, column_annotation_side = "top")
    testthat::expect_true("Legends" %in% class(out$legend))
    testthat::expect_true("Heatmap" %in% class(out$heatmap))

    out <- SCpubr:::heatmap_inner(t(data), column_annotation = obj, column_annotation_side = "bottom")
    testthat::expect_true("Legends" %in% class(out$legend))
    testthat::expect_true("Heatmap" %in% class(out$heatmap))

  })


  # MODIFY STRING

  testthat::test_that("utils: modify_string - PASS - checks", {

    output <- SCpubr:::modify_string("This is a string to cut")
    testthat::expect_type(output, "character")
  })


  # COMPUTE ENRICHMENT SCORES

  testthat::test_that("utils: compute_enrichment_scores - PASS - checks", {



    output <- SCpubr:::compute_enrichment_scores(sample = sample, input_gene_list = list("test" = c("EPC1")), nbin = 1, ctrl = 10)
    testthat::expect_true("Seurat" %in% class(output))
    testthat::expect_true("test" %in% colnames(output@meta.data))

    output <- SCpubr:::compute_enrichment_scores(sample = sample, input_gene_list = list("test" = c("EPC1")), verbose = TRUE, nbin = 1, ctrl = 10)
    testthat::expect_true("Seurat" %in% class(output))
    testthat::expect_true("test" %in% colnames(output@meta.data))

    output <- SCpubr:::compute_enrichment_scores(sample = sample, input_gene_list = c("EPC1"), nbin = 1, ctrl = 10)
    testthat::expect_true("Seurat" %in% class(output))
    testthat::expect_true("Input" %in% colnames(output@meta.data))
  })

  # GET DATA COLUMN

  testthat::test_that("utils: get data column - PASS ", {



    data <- SCpubr:::get_data_column(sample = sample, feature = "EPC1", assay = "SCT", slot = "data")
    testthat::expect_true("data.frame" %in% class(data))
    testthat::expect_true("feature" %in% colnames(data))

    data <- SCpubr:::get_data_column(sample = sample, feature = "nCount_RNA", assay = "SCT", slot = "data")
    testthat::expect_true("data.frame" %in% class(data))
    testthat::expect_true("feature" %in% colnames(data))

    data <- SCpubr:::get_data_column(sample = sample, feature = "PC_1", assay = "SCT", slot = "data")
    testthat::expect_true("data.frame" %in% class(data))
    testthat::expect_true("feature" %in% colnames(data))
  })

  # CHECK PARAMETERS
  testthat::test_that("utils: check parameters - FAIL ", {

    testthat::expect_error({SCpubr:::check_parameters(parameter = -2, parameter_name = "viridis_direction")})
    testthat::expect_error({SCpubr:::check_parameters(parameter = "ERROR", parameter_name = "viridis_color_map")})
  })
}







