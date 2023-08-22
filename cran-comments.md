# Submission version 2.0.1
Fixed Startup messages and loading time at startup.

The package no longer calls for any mirrors to check package versions. I have removed every instance of this.

While `devtools::check` returns no warnings, `R CMD CHECK` does return:
─  checking package dependencies ...Warning: unable to access index for repository https://bioconductor.org/packages/3.15/bioc/src/contrib: (3s)
     cannot open URL 'https://bioconductor.org/packages/3.15/bioc/src/contrib/PACKAGES'
   Warning: unable to access index for repository https://bioconductor.org/packages/3.15/data/annotation/src/contrib:
     cannot open URL 'https://bioconductor.org/packages/3.15/data/annotation/src/contrib/PACKAGES'
   Warning: unable to access index for repository https://bioconductor.org/packages/3.15/data/experiment/src/contrib:
     cannot open URL 'https://bioconductor.org/packages/3.15/data/experiment/src/contrib/PACKAGES'
     
For which I am unable to locate the exact reason why it only appears there. Would it be possible to receive guidance on this if this
turns out to be a real problem?


## `devtools` R CMD check results
There were no ERRORs, WARNINGs or NOTEs. 

## `devtools` R CMD check results `_R_CHECK_DEPENDS_ONLY_` = TRUE
There were no ERRORs, WARNINGs or NOTEs.

## Downstream dependencies
There are currently no downstream dependencies for this package.
