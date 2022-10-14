# SCpubr v1.0.0-dev-stable
- Same as v1.0.0, but with all the functions that do not pass CRAN checks. These functions are: `SCpubr::save_Plot()` `SCpubr::do_LigandReceptorPlot()` and `SCpubr::do_SankeyPlot()`.
 
This version can be obtained in the [v1.0.0-dev-stable release](https://github.com/enblacar/SCpubr/releases/tag/v1.0.0-dev-stable) release. 

# SCpubr v1.0.0
- Modified internal checks so that the functions that do not use `Seurat` do not require this to be installed. This is just for the very side case in which somebody downloads the package just for the `SCpubr::do_ColorPalette()` function.
- Removed the option to use `individual.titles`, 'individual.subtitles` and `individual.captions` in `SCpubr::do_NebulosaPlot()` as the benefit of such parameters did not surpass the problems the code was causing. The feature might come back in the future, once fully optimized.
- Removed `SCpubr::save_Plot()` function to align with CRAN policies that the package should not write to the file system. The code is still available in the v0.0.0.0.9000 release in Github.
- Removed `SCpubr::do_LigandReceptorPlot()`, `SCpubr::do_SankeyPlot()` and `SCpubr::do_PseudotimePlot()` to align with CRAN policies and make it possible to publish the package. These functions can still be accessed in the v0.0.0.0.9000 release in Github.
- Removed `SCpubr::do_PseudotimePlot()` for the reason above and because the dependency `Matrix.utils` was removed from CRAN on *09-10-2022*. 

This version can be obtained in the [v1.0.0 release](https://github.com/enblacar/SCpubr/releases/tag/v1.0.0) release. 

# SCpubr 0.0.0.9000

- Added a `NEWS.md` file to track changes to the package.
- Prepare package for submission to CRAN.
