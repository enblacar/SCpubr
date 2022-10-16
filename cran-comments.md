This update is meant to fix urgent change requests from the CRAN team. 
In it, I have made sure that all packages on Suggests are used conditionally. 
This applies to the R code, to the examples and to the tests. Whenever the
suggested packages are not available, no code will be run for the R scripts,
a message will be displayed for the examples and no tests will be run.

## R CMD check results
There were no ERRORs, WARNINGs or NOTEs. 

## R CMD check results `_R_CHECK_DEPENDS_ONLY_` = TRUE
There were no ERRORs, WARNINGs or NOTEs.

## Downstream dependencies
There are currently no downstream dependencies for this package.
