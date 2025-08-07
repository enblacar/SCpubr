# Submission version 3.0.0
Major update, including standarization of function names on the context of future submission to scientific journal.

## `devtools` R CMD check results
There were no ERRORs, WARNINGs or NOTEs. 

## `devtools` R CMD check results `_R_CHECK_DEPENDS_ONLY_` = TRUE
There were no ERRORs, WARNINGs or NOTEs.

## Using `devtools::test()`
When using first in the session, two tests return a deprecation warning:

The `slot` argument of `FetchData()` is deprecated as of SeuratObject 5.0.0.
i Please use the `layer` argument instead.
i The deprecated feature was likely used in the Seurat package.
  Please report the issue at <https://github.com/satijalab/seurat/issues>.
Backtrace:
    ▆
 1. └─SCpubr::do_CellularStatesPlot(...) at test-do_CellularStatesPlot.R:21:5
 2.   └─SCpubr::do_FeaturePlot(...) at SCpubr/R/do_CellularStatesPlot.R:594:9
 3.     └─Seurat::FeaturePlot(...) at SCpubr/R/do_FeaturePlot.R:413:7
 4.       ├─SeuratObject::FetchData(...)
 5.       └─SeuratObject:::FetchData.Seurat(...)
 6.         └─SeuratObject::.Deprecate(...)

This is a known issue the solution of which is waiting for a fix in Seurat's future update.
Seurat's team is informed about this.

## Downstream dependencies
There are currently no downstream dependencies for this package.
