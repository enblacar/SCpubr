# SCpubr v1.0.0
- Removed `SCpubr::save_Plot()` function to align with CRAN policies that the package should not write to the file system. The code is still available in the [v1.0.0-dev-stable release](https://github.com/enblacar/SCpubr/releases/tag/v1.0.0-dev-stable) .
- Removed `SCpubr::do_LigandReceptorPlot()`, `SCpubr::do_SankeyPlot()` and `SCpubr::do_PseudotimePlot()` to align with CRAN policies and make it possible to publish the package. These functions can still be accessed in the [v1.0.0-dev-stable release](https://github.com/enblacar/SCpubr/releases/tag/v1.0.0-dev-stable).
- Removed `SCpubr::do_PseudotimePlot()` for the reason above and because the dependency `Matrix.utils` was removed from CRAN on *09-10-2022*. 

# SCpubr 0.0.0.9000

- Added a `NEWS.md` file to track changes to the package.
- Prepare package for submission to CRAN.
