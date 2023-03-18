# Declare .data as global variable to avoid R CMD checks.
utils::globalVariables(c(".data",
                        ".env",
                        ".",
                       "x",
                       "quantile",
                       "ecdf",
                       "summarise",
                       "stratum"))
