This is the first release of SCpubr.
## Resubmission 4
This is a fourth resubmission. In this version I have:

* Modified examples so that `\dontrun` is removed, even for internal functions, as it was not necessary. The later two cases was an oversight from my side. This should be fixed now.

R CMD CHECK runs without ERRORs, WARNINGs and NOTEs.

## Resubmission 3
This is a third resubmission. In this version I have:

* Modified examples so that `\dontrun` is removed. In the cases in which the runtime is higher than 5s (measured with `system.time()` and reported also by `R CMD check` and `devtools::check()`), I used `\donttest` instead.
* Removed the use of `T` and `F` in examples, tests and scripts. Resorted to use `TRUE` and `FALSE` instead. 
* Moved `sysdata.rda` to `inst/extdata` so that it can be used for both examples and tests. 

R CMD CHECK runs without ERRORs, WARNINGs and NOTEs.

## Resubmission 2
This is a second resubmission. In this version I have:

* Modified Description field in DESCRIPTION so that it does not contain "This package".
* Modified the problematic link in `README.md` that was moved to another address. If this problem still persists, it will be helpful to get further feedback on how to correctly solve it.
* Drastically reduced the package size to be under 5 MB.

R CMD CHECK runs without ERRORs and WARNINGs.

After applying the changes, there is still one NOTE:
N  checking installed package size ... 
     installed size is  6.5Mb
     sub-directories of 1Mb or more:
       R   4.1Mb

Which is not present when `devtools::check()` is run.

## Resubmission
This is a resubmission. In this version I have:

* Converted the DESCRIPTION title to title case.

* Removed the package name from Description field in DESCRIPTION.
# Reformatted the Author@R entry so that my Name and family names are used without parameters and assigned to myself the role of "aut" to the existing "cre".
* Reformatted Title and Description so that it states that the focus is on Single Cell **Transcriptomics** data.
* Removed `library(monocle3)` call in setup.R in `tests/testthat/` as it was a leftover from the tests for a function that I removed as it still relied on a Github package stated previously in `Remotes:`. 
* Adapted the documentation in `SCpubr::do_CellularStatesPlot()` so that the DOI are used in the context of the `\doi` R Markdown function.

## R CMD check results
There were no ERRORs or WARNINGs. 

There was 1 NOTE (using the `Check` button in RStudio with --as-cran argument):

N  checking installed package size ...
     installed size is 15.7Mb
     sub-directories of 1Mb or more:
       R  13.2Mb

However, when run with `devtools::check(), this note turns into (this uses R CMD Check with arguments --no-manual --as-cran):
N  checking installed package size ...
     installed size is  9.1Mb
     sub-directories of 1Mb or more:
       R   8.8Mb

Which is more accurate to the file size I get from the command line.
  
Nevertheless, this note arises from the mock datasets needed to run the tests for this package. As this package needs single-cell datasets to work, producing a minimal example that runs without compromising the quality of the tests exceeds 1 Mb. This size is not expected to grow significantly larger in future releases. 

## Downstream dependencies
There are currently no downstream dependencies for this package.
