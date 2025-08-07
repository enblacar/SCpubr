if(base::isFALSE(dep_check[["do_LigandReceptorPlot"]])){
  testthat::test_that("do_LigandReceptorPlot: CRAN essentials", {
    p <- SCpubr::do_LigandReceptorPlot(liana_output = liana_output, verbose = FALSE)
    testthat::expect_true(ggplot2::is_ggplot(p))
  })

  testthat::test_that("do_LigandReceptorPlot: PASS - from output", {
    testthat::skip_on_cran()
    
    suppressMessages({testthat::expect_message({p <- SCpubr::do_LigandReceptorPlot(liana_output = liana_output, verbose = TRUE)})})
    testthat::expect_true(ggplot2::is_ggplot(p))
    
    p <- SCpubr::do_LigandReceptorPlot(liana_output = liana_output, verbose = FALSE, sort.by = "A", invert_specificity = TRUE, invert_magnitude = TRUE, sorting.type.specificity = "ascending", sorting.type.magnitude = "ascending")
    testthat::expect_true(ggplot2::is_ggplot(p))
    
    suppressMessages(p <- SCpubr::do_LigandReceptorPlot(liana_output = liana_output, sort.by = "A", invert_specificity = TRUE, invert_magnitude = TRUE, sorting.type.specificity = "ascending", sorting.type.magnitude = "ascending", verbose = TRUE))
    testthat::expect_true(ggplot2::is_ggplot(p))
    
    p <- SCpubr::do_LigandReceptorPlot(liana_output = liana_output, verbose = FALSE, sort.by = "A", invert_specificity = TRUE, invert_magnitude = TRUE, sorting.type.specificity = "descending", sorting.type.magnitude = "ascending")
    testthat::expect_true(ggplot2::is_ggplot(p))
    
    p <- SCpubr::do_LigandReceptorPlot(liana_output = liana_output, verbose = FALSE, sort.by = "A", invert_specificity = TRUE, invert_magnitude = TRUE, sorting.type.specificity = "ascending", sorting.type.magnitude = "descending")
    testthat::expect_true(ggplot2::is_ggplot(p))
    
    p <- SCpubr::do_LigandReceptorPlot(liana_output = liana_output, verbose = FALSE, sort.by = "A", invert_specificity = TRUE, invert_magnitude = TRUE, sorting.type.specificity = "descending", sorting.type.magnitude = "descending")
    testthat::expect_true(ggplot2::is_ggplot(p))
    
    
    p <- SCpubr::do_LigandReceptorPlot(liana_output = liana_output, verbose = FALSE, sort.by = "A", invert_specificity = TRUE, invert_magnitude = FALSE, sorting.type.specificity = "ascending", sorting.type.magnitude = "ascending")
    testthat::expect_true(ggplot2::is_ggplot(p))
    
    p <- SCpubr::do_LigandReceptorPlot(liana_output = liana_output, verbose = FALSE, sort.by = "A", invert_specificity = TRUE, invert_magnitude = FALSE, sorting.type.specificity = "descending", sorting.type.magnitude = "ascending")
    testthat::expect_true(ggplot2::is_ggplot(p))
    
    p <- SCpubr::do_LigandReceptorPlot(liana_output = liana_output, verbose = FALSE, sort.by = "A", invert_specificity = TRUE, invert_magnitude = FALSE, sorting.type.specificity = "ascending", sorting.type.magnitude = "descending")
    testthat::expect_true(ggplot2::is_ggplot(p))
    
    p <- SCpubr::do_LigandReceptorPlot(liana_output = liana_output, verbose = FALSE, sort.by = "A", invert_specificity = TRUE, invert_magnitude = FALSE, sorting.type.specificity = "descending", sorting.type.magnitude = "descending")
    testthat::expect_true(ggplot2::is_ggplot(p))
    
    
    p <- SCpubr::do_LigandReceptorPlot(liana_output = liana_output, verbose = FALSE, sort.by = "A", invert_specificity = FALSE, invert_magnitude = TRUE, sorting.type.specificity = "ascending", sorting.type.magnitude = "ascending")
    testthat::expect_true(ggplot2::is_ggplot(p))
    
    p <- SCpubr::do_LigandReceptorPlot(liana_output = liana_output, verbose = FALSE, sort.by = "A", invert_specificity = FALSE, invert_magnitude = TRUE, sorting.type.specificity = "descending", sorting.type.magnitude = "ascending")
    testthat::expect_true(ggplot2::is_ggplot(p))
    
    p <- SCpubr::do_LigandReceptorPlot(liana_output = liana_output, verbose = FALSE, sort.by = "A", invert_specificity = FALSE, invert_magnitude = TRUE, sorting.type.specificity = "ascending", sorting.type.magnitude = "descending")
    testthat::expect_true(ggplot2::is_ggplot(p))
    
    p <- SCpubr::do_LigandReceptorPlot(liana_output = liana_output, verbose = FALSE, sort.by = "A", invert_specificity = FALSE, invert_magnitude = TRUE, sorting.type.specificity = "descending", sorting.type.magnitude = "descending")
    testthat::expect_true(ggplot2::is_ggplot(p))
    
    
    p <- SCpubr::do_LigandReceptorPlot(liana_output = liana_output, verbose = FALSE, sort.by = "A", invert_specificity = FALSE, invert_magnitude = FALSE, sorting.type.specificity = "ascending", sorting.type.magnitude = "ascending")
    testthat::expect_true(ggplot2::is_ggplot(p))
    
    p <- SCpubr::do_LigandReceptorPlot(liana_output = liana_output, verbose = FALSE, sort.by = "A", invert_specificity = FALSE, invert_magnitude = FALSE, sorting.type.specificity = "descending", sorting.type.magnitude = "ascending")
    testthat::expect_true(ggplot2::is_ggplot(p))
    
    p <- SCpubr::do_LigandReceptorPlot(liana_output = liana_output, verbose = FALSE, sort.by = "A", invert_specificity = FALSE, invert_magnitude = FALSE, sorting.type.specificity = "ascending", sorting.type.magnitude = "descending")
    testthat::expect_true(ggplot2::is_ggplot(p))
    
    p <- SCpubr::do_LigandReceptorPlot(liana_output = liana_output, verbose = FALSE, sort.by = "A", invert_specificity = FALSE, invert_magnitude = FALSE, sorting.type.specificity = "descending", sorting.type.magnitude = "descending")
    testthat::expect_true(ggplot2::is_ggplot(p))
    
    
    
    
    
    
    
    p <- SCpubr::do_LigandReceptorPlot(liana_output = liana_output, verbose = FALSE, sort.by = "B", invert_specificity = TRUE, invert_magnitude = TRUE, sorting.type.specificity = "ascending", sorting.type.magnitude = "ascending")
    testthat::expect_true(ggplot2::is_ggplot(p))
    
    suppressMessages(p <- SCpubr::do_LigandReceptorPlot(liana_output = liana_output, sort.by = "B", invert_specificity = TRUE, invert_magnitude = TRUE, sorting.type.specificity = "ascending", sorting.type.magnitude = "ascending", verbose = TRUE))
    testthat::expect_true(ggplot2::is_ggplot(p))
    
    p <- SCpubr::do_LigandReceptorPlot(liana_output = liana_output, verbose = FALSE, sort.by = "B", invert_specificity = TRUE, invert_magnitude = TRUE, sorting.type.specificity = "descending", sorting.type.magnitude = "ascending")
    testthat::expect_true(ggplot2::is_ggplot(p))
    
    p <- SCpubr::do_LigandReceptorPlot(liana_output = liana_output, verbose = FALSE, sort.by = "B", invert_specificity = TRUE, invert_magnitude = TRUE, sorting.type.specificity = "ascending", sorting.type.magnitude = "descending")
    testthat::expect_true(ggplot2::is_ggplot(p))
    
    p <- SCpubr::do_LigandReceptorPlot(liana_output = liana_output, verbose = FALSE, sort.by = "B", invert_specificity = TRUE, invert_magnitude = TRUE, sorting.type.specificity = "descending", sorting.type.magnitude = "descending")
    testthat::expect_true(ggplot2::is_ggplot(p))
    
    
    p <- SCpubr::do_LigandReceptorPlot(liana_output = liana_output, verbose = FALSE, sort.by = "B", invert_specificity = TRUE, invert_magnitude = FALSE, sorting.type.specificity = "ascending", sorting.type.magnitude = "ascending")
    testthat::expect_true(ggplot2::is_ggplot(p))
    
    p <- SCpubr::do_LigandReceptorPlot(liana_output = liana_output, verbose = FALSE, sort.by = "B", invert_specificity = TRUE, invert_magnitude = FALSE, sorting.type.specificity = "descending", sorting.type.magnitude = "ascending")
    testthat::expect_true(ggplot2::is_ggplot(p))
    
    p <- SCpubr::do_LigandReceptorPlot(liana_output = liana_output, verbose = FALSE, sort.by = "B", invert_specificity = TRUE, invert_magnitude = FALSE, sorting.type.specificity = "ascending", sorting.type.magnitude = "descending")
    testthat::expect_true(ggplot2::is_ggplot(p))
    
    p <- SCpubr::do_LigandReceptorPlot(liana_output = liana_output, verbose = FALSE, sort.by = "B", invert_specificity = TRUE, invert_magnitude = FALSE, sorting.type.specificity = "descending", sorting.type.magnitude = "descending")
    testthat::expect_true(ggplot2::is_ggplot(p))
    
    
    p <- SCpubr::do_LigandReceptorPlot(liana_output = liana_output, verbose = FALSE, sort.by = "B", invert_specificity = FALSE, invert_magnitude = TRUE, sorting.type.specificity = "ascending", sorting.type.magnitude = "ascending")
    testthat::expect_true(ggplot2::is_ggplot(p))
    
    p <- SCpubr::do_LigandReceptorPlot(liana_output = liana_output, verbose = FALSE, sort.by = "B", invert_specificity = FALSE, invert_magnitude = TRUE, sorting.type.specificity = "descending", sorting.type.magnitude = "ascending")
    testthat::expect_true(ggplot2::is_ggplot(p))
    
    p <- SCpubr::do_LigandReceptorPlot(liana_output = liana_output, verbose = FALSE, sort.by = "B", invert_specificity = FALSE, invert_magnitude = TRUE, sorting.type.specificity = "ascending", sorting.type.magnitude = "descending")
    testthat::expect_true(ggplot2::is_ggplot(p))
    
    p <- SCpubr::do_LigandReceptorPlot(liana_output = liana_output, verbose = FALSE, sort.by = "B", invert_specificity = FALSE, invert_magnitude = TRUE, sorting.type.specificity = "descending", sorting.type.magnitude = "descending")
    testthat::expect_true(ggplot2::is_ggplot(p))
    
    
    p <- SCpubr::do_LigandReceptorPlot(liana_output = liana_output, verbose = FALSE, sort.by = "B", invert_specificity = FALSE, invert_magnitude = FALSE, sorting.type.specificity = "ascending", sorting.type.magnitude = "ascending")
    testthat::expect_true(ggplot2::is_ggplot(p))
    
    p <- SCpubr::do_LigandReceptorPlot(liana_output = liana_output, verbose = FALSE, sort.by = "B", invert_specificity = FALSE, invert_magnitude = FALSE, sorting.type.specificity = "descending", sorting.type.magnitude = "ascending")
    testthat::expect_true(ggplot2::is_ggplot(p))
    
    p <- SCpubr::do_LigandReceptorPlot(liana_output = liana_output, verbose = FALSE, sort.by = "B", invert_specificity = FALSE, invert_magnitude = FALSE, sorting.type.specificity = "ascending", sorting.type.magnitude = "descending")
    testthat::expect_true(ggplot2::is_ggplot(p))
    
    p <- SCpubr::do_LigandReceptorPlot(liana_output = liana_output, verbose = FALSE, sort.by = "B", invert_specificity = FALSE, invert_magnitude = FALSE, sorting.type.specificity = "descending", sorting.type.magnitude = "descending")
    testthat::expect_true(ggplot2::is_ggplot(p))
    
    
    
    p <- SCpubr::do_LigandReceptorPlot(liana_output = liana_output, verbose = FALSE, sort.by = "C", invert_specificity = TRUE, invert_magnitude = TRUE, sorting.type.specificity = "ascending", sorting.type.magnitude = "ascending")
    testthat::expect_true(ggplot2::is_ggplot(p))
    
    suppressMessages(p <- SCpubr::do_LigandReceptorPlot(liana_output = liana_output, sort.by = "C", invert_specificity = TRUE, invert_magnitude = TRUE, sorting.type.specificity = "ascending", sorting.type.magnitude = "ascending", verbose = TRUE))
    testthat::expect_true(ggplot2::is_ggplot(p))
    
    p <- SCpubr::do_LigandReceptorPlot(liana_output = liana_output, verbose = FALSE, sort.by = "C", invert_specificity = TRUE, invert_magnitude = TRUE, sorting.type.specificity = "descending", sorting.type.magnitude = "ascending")
    testthat::expect_true(ggplot2::is_ggplot(p))
    
    p <- SCpubr::do_LigandReceptorPlot(liana_output = liana_output, verbose = FALSE, sort.by = "C", invert_specificity = TRUE, invert_magnitude = TRUE, sorting.type.specificity = "ascending", sorting.type.magnitude = "descending")
    testthat::expect_true(ggplot2::is_ggplot(p))
    
    p <- SCpubr::do_LigandReceptorPlot(liana_output = liana_output, verbose = FALSE, sort.by = "C", invert_specificity = TRUE, invert_magnitude = TRUE, sorting.type.specificity = "descending", sorting.type.magnitude = "descending")
    testthat::expect_true(ggplot2::is_ggplot(p))
    
    
    p <- SCpubr::do_LigandReceptorPlot(liana_output = liana_output, verbose = FALSE, sort.by = "C", invert_specificity = TRUE, invert_magnitude = FALSE, sorting.type.specificity = "ascending", sorting.type.magnitude = "ascending")
    testthat::expect_true(ggplot2::is_ggplot(p))
    
    p <- SCpubr::do_LigandReceptorPlot(liana_output = liana_output, verbose = FALSE, sort.by = "C", invert_specificity = TRUE, invert_magnitude = FALSE, sorting.type.specificity = "descending", sorting.type.magnitude = "ascending")
    testthat::expect_true(ggplot2::is_ggplot(p))
    
    p <- SCpubr::do_LigandReceptorPlot(liana_output = liana_output, verbose = FALSE, sort.by = "C", invert_specificity = TRUE, invert_magnitude = FALSE, sorting.type.specificity = "ascending", sorting.type.magnitude = "descending")
    testthat::expect_true(ggplot2::is_ggplot(p))
    
    p <- SCpubr::do_LigandReceptorPlot(liana_output = liana_output, verbose = FALSE, sort.by = "C", invert_specificity = TRUE, invert_magnitude = FALSE, sorting.type.specificity = "descending", sorting.type.magnitude = "descending")
    testthat::expect_true(ggplot2::is_ggplot(p))
    
    
    p <- SCpubr::do_LigandReceptorPlot(liana_output = liana_output, verbose = FALSE, sort.by = "C", invert_specificity = FALSE, invert_magnitude = TRUE, sorting.type.specificity = "ascending", sorting.type.magnitude = "ascending")
    testthat::expect_true(ggplot2::is_ggplot(p))
    
    p <- SCpubr::do_LigandReceptorPlot(liana_output = liana_output, verbose = FALSE, sort.by = "C", invert_specificity = FALSE, invert_magnitude = TRUE, sorting.type.specificity = "descending", sorting.type.magnitude = "ascending")
    testthat::expect_true(ggplot2::is_ggplot(p))
    
    p <- SCpubr::do_LigandReceptorPlot(liana_output = liana_output, verbose = FALSE, sort.by = "C", invert_specificity = FALSE, invert_magnitude = TRUE, sorting.type.specificity = "ascending", sorting.type.magnitude = "descending")
    testthat::expect_true(ggplot2::is_ggplot(p))
    
    p <- SCpubr::do_LigandReceptorPlot(liana_output = liana_output, verbose = FALSE, sort.by = "C", invert_specificity = FALSE, invert_magnitude = TRUE, sorting.type.specificity = "descending", sorting.type.magnitude = "descending")
    testthat::expect_true(ggplot2::is_ggplot(p))
    
    
    p <- SCpubr::do_LigandReceptorPlot(liana_output = liana_output, verbose = FALSE, sort.by = "C", invert_specificity = FALSE, invert_magnitude = FALSE, sorting.type.specificity = "ascending", sorting.type.magnitude = "ascending")
    testthat::expect_true(ggplot2::is_ggplot(p))
    
    p <- SCpubr::do_LigandReceptorPlot(liana_output = liana_output, verbose = FALSE, sort.by = "C", invert_specificity = FALSE, invert_magnitude = FALSE, sorting.type.specificity = "descending", sorting.type.magnitude = "ascending")
    testthat::expect_true(ggplot2::is_ggplot(p))
    
    p <- SCpubr::do_LigandReceptorPlot(liana_output = liana_output, verbose = FALSE, sort.by = "C", invert_specificity = FALSE, invert_magnitude = FALSE, sorting.type.specificity = "ascending", sorting.type.magnitude = "descending")
    testthat::expect_true(ggplot2::is_ggplot(p))
    
    p <- SCpubr::do_LigandReceptorPlot(liana_output = liana_output, verbose = FALSE, sort.by = "C", invert_specificity = FALSE, invert_magnitude = FALSE, sorting.type.specificity = "descending", sorting.type.magnitude = "descending")
    testthat::expect_true(ggplot2::is_ggplot(p))
    
    
    
    
    p <- SCpubr::do_LigandReceptorPlot(liana_output = liana_output, verbose = FALSE, sort.by = "D", invert_specificity = TRUE, invert_magnitude = TRUE, sorting.type.specificity = "ascending", sorting.type.magnitude = "ascending")
    testthat::expect_true(ggplot2::is_ggplot(p))
    
    suppressMessages(p <- SCpubr::do_LigandReceptorPlot(liana_output = liana_output, sort.by = "D", invert_specificity = TRUE, invert_magnitude = TRUE, sorting.type.specificity = "ascending", sorting.type.magnitude = "ascending", verbose = TRUE))
    testthat::expect_true(ggplot2::is_ggplot(p))
    
    p <- SCpubr::do_LigandReceptorPlot(liana_output = liana_output, verbose = FALSE, sort.by = "D", invert_specificity = TRUE, invert_magnitude = TRUE, sorting.type.specificity = "descending", sorting.type.magnitude = "ascending")
    testthat::expect_true(ggplot2::is_ggplot(p))
    
    p <- SCpubr::do_LigandReceptorPlot(liana_output = liana_output, verbose = FALSE, sort.by = "D", invert_specificity = TRUE, invert_magnitude = TRUE, sorting.type.specificity = "ascending", sorting.type.magnitude = "descending")
    testthat::expect_true(ggplot2::is_ggplot(p))
    
    p <- SCpubr::do_LigandReceptorPlot(liana_output = liana_output, verbose = FALSE, sort.by = "D", invert_specificity = TRUE, invert_magnitude = TRUE, sorting.type.specificity = "descending", sorting.type.magnitude = "descending")
    testthat::expect_true(ggplot2::is_ggplot(p))
    
    
    p <- SCpubr::do_LigandReceptorPlot(liana_output = liana_output, verbose = FALSE, sort.by = "D", invert_specificity = TRUE, invert_magnitude = FALSE, sorting.type.specificity = "ascending", sorting.type.magnitude = "ascending")
    testthat::expect_true(ggplot2::is_ggplot(p))
    
    p <- SCpubr::do_LigandReceptorPlot(liana_output = liana_output, verbose = FALSE, sort.by = "D", invert_specificity = TRUE, invert_magnitude = FALSE, sorting.type.specificity = "descending", sorting.type.magnitude = "ascending")
    testthat::expect_true(ggplot2::is_ggplot(p))
    
    p <- SCpubr::do_LigandReceptorPlot(liana_output = liana_output, verbose = FALSE, sort.by = "D", invert_specificity = TRUE, invert_magnitude = FALSE, sorting.type.specificity = "ascending", sorting.type.magnitude = "descending")
    testthat::expect_true(ggplot2::is_ggplot(p))
    
    p <- SCpubr::do_LigandReceptorPlot(liana_output = liana_output, verbose = FALSE, sort.by = "D", invert_specificity = TRUE, invert_magnitude = FALSE, sorting.type.specificity = "descending", sorting.type.magnitude = "descending")
    testthat::expect_true(ggplot2::is_ggplot(p))
    
    
    p <- SCpubr::do_LigandReceptorPlot(liana_output = liana_output, verbose = FALSE, sort.by = "D", invert_specificity = FALSE, invert_magnitude = TRUE, sorting.type.specificity = "ascending", sorting.type.magnitude = "ascending")
    testthat::expect_true(ggplot2::is_ggplot(p))
    
    p <- SCpubr::do_LigandReceptorPlot(liana_output = liana_output, verbose = FALSE, sort.by = "D", invert_specificity = FALSE, invert_magnitude = TRUE, sorting.type.specificity = "descending", sorting.type.magnitude = "ascending")
    testthat::expect_true(ggplot2::is_ggplot(p))
    
    p <- SCpubr::do_LigandReceptorPlot(liana_output = liana_output, verbose = FALSE, sort.by = "D", invert_specificity = FALSE, invert_magnitude = TRUE, sorting.type.specificity = "ascending", sorting.type.magnitude = "descending")
    testthat::expect_true(ggplot2::is_ggplot(p))
    
    p <- SCpubr::do_LigandReceptorPlot(liana_output = liana_output, verbose = FALSE, sort.by = "D", invert_specificity = FALSE, invert_magnitude = TRUE, sorting.type.specificity = "descending", sorting.type.magnitude = "descending")
    testthat::expect_true(ggplot2::is_ggplot(p))
    
    
    p <- SCpubr::do_LigandReceptorPlot(liana_output = liana_output, verbose = FALSE, sort.by = "D", invert_specificity = FALSE, invert_magnitude = FALSE, sorting.type.specificity = "ascending", sorting.type.magnitude = "ascending")
    testthat::expect_true(ggplot2::is_ggplot(p))
    
    p <- SCpubr::do_LigandReceptorPlot(liana_output = liana_output, verbose = FALSE, sort.by = "D", invert_specificity = FALSE, invert_magnitude = FALSE, sorting.type.specificity = "descending", sorting.type.magnitude = "ascending")
    testthat::expect_true(ggplot2::is_ggplot(p))
    
    p <- SCpubr::do_LigandReceptorPlot(liana_output = liana_output, verbose = FALSE, sort.by = "D", invert_specificity = FALSE, invert_magnitude = FALSE, sorting.type.specificity = "ascending", sorting.type.magnitude = "descending")
    testthat::expect_true(ggplot2::is_ggplot(p))
    
    p <- SCpubr::do_LigandReceptorPlot(liana_output = liana_output, verbose = FALSE, sort.by = "D", invert_specificity = FALSE, invert_magnitude = FALSE, sorting.type.specificity = "descending", sorting.type.magnitude = "descending")
    testthat::expect_true(ggplot2::is_ggplot(p))
    
    
    
    
    p <- SCpubr::do_LigandReceptorPlot(liana_output = liana_output, verbose = FALSE, sort.by = "E", invert_specificity = TRUE, invert_magnitude = TRUE, sorting.type.specificity = "ascending", sorting.type.magnitude = "ascending")
    testthat::expect_true(ggplot2::is_ggplot(p))
    
    p <- SCpubr::do_LigandReceptorPlot(liana_output = liana_output, verbose = FALSE, sort.by = "E", invert_specificity = TRUE, invert_magnitude = TRUE, sorting.type.specificity = "descending", sorting.type.magnitude = "ascending")
    testthat::expect_true(ggplot2::is_ggplot(p))
    
    p <- SCpubr::do_LigandReceptorPlot(liana_output = liana_output, verbose = FALSE, sort.by = "E", invert_specificity = TRUE, invert_magnitude = TRUE, sorting.type.specificity = "ascending", sorting.type.magnitude = "descending")
    testthat::expect_true(ggplot2::is_ggplot(p))
    
    p <- SCpubr::do_LigandReceptorPlot(liana_output = liana_output, verbose = FALSE, sort.by = "E", invert_specificity = TRUE, invert_magnitude = TRUE, sorting.type.specificity = "descending", sorting.type.magnitude = "descending")
    testthat::expect_true(ggplot2::is_ggplot(p))
    
    
    p <- SCpubr::do_LigandReceptorPlot(liana_output = liana_output, verbose = FALSE, sort.by = "E", invert_specificity = TRUE, invert_magnitude = FALSE, sorting.type.specificity = "ascending", sorting.type.magnitude = "ascending")
    testthat::expect_true(ggplot2::is_ggplot(p))
    
    suppressMessages(p <- SCpubr::do_LigandReceptorPlot(liana_output = liana_output,  sort.by = "E", invert_specificity = TRUE, invert_magnitude = FALSE, sorting.type.specificity = "ascending", sorting.type.magnitude = "ascending", verbose = TRUE))
    testthat::expect_true(ggplot2::is_ggplot(p))
    
    p <- SCpubr::do_LigandReceptorPlot(liana_output = liana_output, verbose = FALSE, sort.by = "E", invert_specificity = TRUE, invert_magnitude = FALSE, sorting.type.specificity = "descending", sorting.type.magnitude = "ascending")
    testthat::expect_true(ggplot2::is_ggplot(p))
    
    p <- SCpubr::do_LigandReceptorPlot(liana_output = liana_output, verbose = FALSE, sort.by = "E", invert_specificity = TRUE, invert_magnitude = FALSE, sorting.type.specificity = "ascending", sorting.type.magnitude = "descending")
    testthat::expect_true(ggplot2::is_ggplot(p))
    
    p <- SCpubr::do_LigandReceptorPlot(liana_output = liana_output, verbose = FALSE, sort.by = "E", invert_specificity = TRUE, invert_magnitude = FALSE, sorting.type.specificity = "descending", sorting.type.magnitude = "descending")
    testthat::expect_true(ggplot2::is_ggplot(p))
    
    
    p <- SCpubr::do_LigandReceptorPlot(liana_output = liana_output, verbose = FALSE, sort.by = "E", invert_specificity = FALSE, invert_magnitude = TRUE, sorting.type.specificity = "ascending", sorting.type.magnitude = "ascending")
    testthat::expect_true(ggplot2::is_ggplot(p))
    
    p <- SCpubr::do_LigandReceptorPlot(liana_output = liana_output, verbose = FALSE, sort.by = "E", invert_specificity = FALSE, invert_magnitude = TRUE, sorting.type.specificity = "descending", sorting.type.magnitude = "ascending")
    testthat::expect_true(ggplot2::is_ggplot(p))
    
    p <- SCpubr::do_LigandReceptorPlot(liana_output = liana_output, verbose = FALSE, sort.by = "E", invert_specificity = FALSE, invert_magnitude = TRUE, sorting.type.specificity = "ascending", sorting.type.magnitude = "descending")
    testthat::expect_true(ggplot2::is_ggplot(p))
    
    p <- SCpubr::do_LigandReceptorPlot(liana_output = liana_output, verbose = FALSE, sort.by = "E", invert_specificity = FALSE, invert_magnitude = TRUE, sorting.type.specificity = "descending", sorting.type.magnitude = "descending")
    testthat::expect_true(ggplot2::is_ggplot(p))
    
    
    p <- SCpubr::do_LigandReceptorPlot(liana_output = liana_output, verbose = FALSE, sort.by = "E", invert_specificity = FALSE, invert_magnitude = FALSE, sorting.type.specificity = "ascending", sorting.type.magnitude = "ascending")
    testthat::expect_true(ggplot2::is_ggplot(p))
    
    p <- SCpubr::do_LigandReceptorPlot(liana_output = liana_output, verbose = FALSE, sort.by = "E", invert_specificity = FALSE, invert_magnitude = FALSE, sorting.type.specificity = "descending", sorting.type.magnitude = "ascending")
    testthat::expect_true(ggplot2::is_ggplot(p))
    
    p <- SCpubr::do_LigandReceptorPlot(liana_output = liana_output, verbose = FALSE, sort.by = "E", invert_specificity = FALSE, invert_magnitude = FALSE, sorting.type.specificity = "ascending", sorting.type.magnitude = "descending")
    testthat::expect_true(ggplot2::is_ggplot(p))
    
    p <- SCpubr::do_LigandReceptorPlot(liana_output = liana_output, verbose = FALSE, sort.by = "E", invert_specificity = FALSE, invert_magnitude = FALSE, sorting.type.specificity = "descending", sorting.type.magnitude = "descending")
    testthat::expect_true(ggplot2::is_ggplot(p))
    
    p <- SCpubr::do_LigandReceptorPlot(liana_output = liana_output, verbose = FALSE, sort.by = "E", invert_specificity = FALSE, invert_magnitude = FALSE, sorting.type.specificity = "descending", sorting.type.magnitude = "descending")
    testthat::expect_true(ggplot2::is_ggplot(p))
    
    out <- SCpubr::do_LigandReceptorPlot(liana_output = liana_output, verbose = FALSE, sort.by = "E", invert_specificity = FALSE, invert_magnitude = FALSE, sorting.type.specificity = "descending", sorting.type.magnitude = "descending", return_interactions = TRUE)
    testthat::expect_type(out, "list")
    
    
    
    p <- SCpubr::do_LigandReceptorPlot(liana_output = liana_output, dot_border = TRUE, use_viridis = TRUE, viridis.direction = 1, verbose = FALSE)
    testthat::expect_true(ggplot2::is_ggplot(p))
    
    p <- SCpubr::do_LigandReceptorPlot(liana_output = liana_output, dot_border = TRUE, use_viridis = TRUE, viridis.direction = -1, verbose = FALSE)
    testthat::expect_true(ggplot2::is_ggplot(p))
    
    p <- SCpubr::do_LigandReceptorPlot(liana_output = liana_output, dot_border = TRUE, use_viridis = FALSE, sequential.direction = -1, verbose = FALSE)
    testthat::expect_true(ggplot2::is_ggplot(p))
    
    p <- SCpubr::do_LigandReceptorPlot(liana_output = liana_output, dot_border = TRUE, use_viridis = FALSE, sequential.direction = 1, verbose = FALSE)
    testthat::expect_true(ggplot2::is_ggplot(p))
    
    p <- SCpubr::do_LigandReceptorPlot(liana_output = liana_output, dot_border = FALSE, use_viridis = TRUE, viridis.direction = 1, verbose = FALSE)
    testthat::expect_true(ggplot2::is_ggplot(p))
    
    p <- SCpubr::do_LigandReceptorPlot(liana_output = liana_output, dot_border = FALSE, use_viridis = TRUE, viridis.direction = -1, verbose = FALSE)
    testthat::expect_true(ggplot2::is_ggplot(p))
    
    p <- SCpubr::do_LigandReceptorPlot(liana_output = liana_output, dot_border = FALSE, use_viridis = FALSE, sequential.direction = -1, verbose = FALSE)
    testthat::expect_true(ggplot2::is_ggplot(p))
    
    p <- SCpubr::do_LigandReceptorPlot(liana_output = liana_output, dot_border = FALSE, use_viridis = FALSE, sequential.direction = 1, verbose = FALSE)
    testthat::expect_true(ggplot2::is_ggplot(p))
    
    p <- SCpubr::do_LigandReceptorPlot(liana_output = liana_output, plot.grid = TRUE, verbose = FALSE)
    testthat::expect_true(ggplot2::is_ggplot(p))
    
    p <- SCpubr::do_LigandReceptorPlot(liana_output = liana_output, plot.grid = FALSE, verbose = FALSE)
    testthat::expect_true(ggplot2::is_ggplot(p))
    
    p <- SCpubr::do_LigandReceptorPlot(liana_output = liana_output, plot.grid = TRUE, dot_border = FALSE, verbose = FALSE)
    testthat::expect_true(ggplot2::is_ggplot(p))

    p <- SCpubr::do_LigandReceptorPlot(liana_output = liana_output, verbose = FALSE)
    testthat::expect_true(ggplot2::is_ggplot(p))

    p <- SCpubr::do_LigandReceptorPlot(liana_output = liana_output, 
                                       keep_source = c("NK", "B"),
                                       keep_target = "CD8 T", verbose = FALSE)
    testthat::expect_true(ggplot2::is_ggplot(p))
  })

  testthat::test_that("do_LigandReceptorPlot: PASS - from output different n", {
    testthat::skip_on_cran()
    p <- SCpubr::do_LigandReceptorPlot(liana_output = liana_output,
                                       top_interactions = 50, verbose = FALSE)
    testthat::expect_true(ggplot2::is_ggplot(p))
  })


  testthat::test_that("do_LigandReceptorPlot: PASS - split.by", {
    testthat::skip_on_cran()
    p <- SCpubr::do_LigandReceptorPlot(liana_output = liana_output,
                                       split.by = "ligand.complex", verbose = FALSE)
    testthat::expect_true(ggplot2::is_ggplot(p))

    p <- SCpubr::do_LigandReceptorPlot(liana_output = liana_output,
                                       split.by = "receptor.complex", verbose = FALSE)
    testthat::expect_true(ggplot2::is_ggplot(p))
  })

  testthat::test_that("do_LigandReceptorPlot: PASS - from output, angle ", {
    testthat::skip_on_cran()
    p <- SCpubr::do_LigandReceptorPlot(liana_output = liana_output,
                                       axis.text.x.angle = 0, verbose = FALSE)
    testthat::expect_true(ggplot2::is_ggplot(p))

    p <- SCpubr::do_LigandReceptorPlot(liana_output = liana_output,
                                       axis.text.x.angle = 45, verbose = FALSE)
    testthat::expect_true(ggplot2::is_ggplot(p))

    p <- SCpubr::do_LigandReceptorPlot(liana_output = liana_output,
                                       axis.text.x.angle = 90, verbose = FALSE)
    testthat::expect_true(ggplot2::is_ggplot(p))
  })


  testthat::test_that("do_LigandReceptorPlot: PASS - from output legend.position", {
    testthat::skip_on_cran()
    p <- SCpubr::do_LigandReceptorPlot(liana_output = liana_output,
                                       top_interactions = 50,
                                       legend.position = "bottom", verbose = FALSE)
    testthat::expect_true(ggplot2::is_ggplot(p))

    p <- SCpubr::do_LigandReceptorPlot(liana_output = liana_output,
                                       top_interactions = 50,
                                       legend.position = "right", verbose = FALSE)
    testthat::expect_true(ggplot2::is_ggplot(p))
  })


  

  testthat::test_that("do_LigandReceptorPlot: PASS - sort interactions", {
    testthat::skip_on_cran()
    p <- SCpubr::do_LigandReceptorPlot(liana_output = liana_output,
                                       top_interactions = 50,
                                       sort_interactions_alphabetically =  TRUE, verbose = FALSE)
    testthat::expect_true(ggplot2::is_ggplot(p))

    p <- SCpubr::do_LigandReceptorPlot(liana_output = liana_output,
                                       top_interactions = 50,
                                       sort_interactions_alphabetically =  FALSE, verbose = FALSE)
    testthat::expect_true(ggplot2::is_ggplot(p))
  })

  testthat::test_that("do_LigandReceptorPlot: FAIL - wrong parameters", {
    testthat::skip_on_cran()
    testthat::expect_error({SCpubr::do_LigandReceptorPlot(liana_output = liana_output,
                                                          font.type = "wrong", verbose = FALSE)})

    testthat::expect_error({SCpubr::do_LigandReceptorPlot(liana_output = liana_output,
                                                          legend.type = "wrong", verbose = FALSE)})

    testthat::expect_error({SCpubr::do_LigandReceptorPlot(liana_output = liana_output,
                                                          axis.text.x.angle = 10, verbose = FALSE)})

    testthat::expect_error({SCpubr::do_LigandReceptorPlot(liana_output = liana_output,
                                                          font.type = "wrong", verbose = FALSE)})

    testthat::expect_error({SCpubr::do_LigandReceptorPlot(liana_output = liana_output,
                                                          legend.position = "wrong", verbose = FALSE)})

    testthat::expect_error({SCpubr::do_LigandReceptorPlot(liana_output = liana_output,
                                                          grid.type = "wrong", verbose = FALSE)})

    testthat::expect_error({SCpubr::do_LigandReceptorPlot(liana_output = liana_output,
                                                          split.by = "wrong", verbose = FALSE)})

  })
}

