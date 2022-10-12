# SCpubr v1.0.0
- Removed `SCpubr::save_Plot()` function to align with CRAN policies that the package should not write to the file system. The code can be requested in scpubr@gmail.com, if somebody still wants it or is still available in the [v.0.0.0.9000 release](https://github.com/enblacar/SCpubr/releases/tag/v.0.0.0.9000) in the [GitHub](https://github.com/enblacar/SCpubr).
- Removed `SCpubr::do_LigandReceptorPlot()`, `SCpubr::do_SankeyPlot()` to align with CRAN policies and make it possible to publish the package. These functions can still be accessed in the [v.0.0.0.9000 release](https://github.com/enblacar/SCpubr/releases/tag/v.0.0.0.9000) in the [GitHub](https://github.com/enblacar/SCpubr).
- Removed `SCpubr::do_PseudotimePlot()` for the reason above and because the dependency `Matrix.utils` was removed from CRAN on *09-10-2022*. 

# SCpubr 0.0.0.9000

- Added a `NEWS.md` file to track changes to the package.
- Prepare package for submission to CRAN.
