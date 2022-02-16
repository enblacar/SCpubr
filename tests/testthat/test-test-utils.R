sample <- testing_data


# CHECK COLORS
testthat::test_that("utils: check_colors - FAIL - wrong color", {
  testthat::expect_error(check_colors("not_a_color"))

})

testthat::test_that("utils: check_colors - FAIL - wrong color in a vector of colors", {
  testthat::expect_error(check_colors(c("not_a_color", "red", "blue")))
})



# CHECK CONSISTENCY COLORS AND NAMES

testthat::test_that("utils: check_colors - FAIL - unequal number of colors", {
  testthat::expect_error(check_consistency_colors_and_names(sample = sample,
                                                            colors = c("a" = "red", "b" = "blue"),
                                                            grouping_variable = "orig.ident"))
})

testthat::test_that("utils: check_colors - FAIL - names of colors not matching", {
  testthat::expect_error(check_consistency_colors_and_names(sample = sample,
                                                            colors = c("a" = "red"),
                                                            grouping_variable = "orig.ident"))
})


# GENERATE COLOR SCALE
testthat::test_that("utils: generate_color_scale - PASS - equal length of output", {
  names_use <- c("a", "b", "c")
  colors <- colortools::setColors("#457b9d", length(names_use))
  testthat::expect_length(colors, length(names_use))
})



