# Declare .data as global variable to avoid R CMD checks.
utils::globalVariables(c(".data",
                       "..x..",
                       "..quantile..",
                       "..ecdf..",
                       "summarise"))
