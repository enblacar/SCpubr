This is the first release of SCpubr.

## R CMD check results
There were no ERRORs or WARNINGs. 

There was 1 NOTE (using the `Check` button in RStudio with --as-cran argument):

N  checking installed package size ...
     installed size is 15.7Mb
     sub-directories of 1Mb or more:
       R  13.2Mb

However, when run with `devtools::check(), this note turns into:
N  checking installed package size ...
     installed size is  9.1Mb
     sub-directories of 1Mb or more:
       R   8.8Mb

Which is more accurate to the file size I get from the command line.
  
Nevertheless, this note arises from the mock datasets needed to run the tests for this package. As this package needs single-cell datasets to work, producing a minimal example that runs without compromising the quality of the tests exceeds 1 Mb. This size is not expected to grow significantly larger in future releases. 

## Downstream dependencies
There are currently no downstream dependencies for this package.
