This is the first release of SCpubr.

## R CMD check results
There were no ERRORs or WARNINGs. 

There was 1 NOTE:

* checking installed package size ... NOTE
    installed size is 17.7Mb
    sub-directories of 1Mb or more:
      help   2.9Mb
      R     13.3Mb
  
  This note arises from the mock datasets needed to run the tests for this package. As this package needs single-cell datasets to work, producing a minimal example that runs 
  without compromising the quality of the tests exceeds 1 Mb. This size is not expected to grow significantly larger in future releases. 

## Downstream dependencies
There are currently no downstream dependencies for this package.
