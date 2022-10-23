# SCpubr v1.0.1.9000 (current in-development version)
Take this a teaser of things to come in the next version.

- Modified the code for `SCpubr::do_CopyNumberVariantPlot` to also report results for the whole chromosome as well as for each chromosome arm.
- Fixed potential bugs in `SCpubr::do_FeaturePlot` when setting `enforce_symmetry = TRUE`.
- Change color palette when using `enforce_symmetry = TRUE` to have the middle color as `grey95` instead of the previous one, which made middle values seem closer to the positive end of the scale.
- Include the `verbose` argument to `SCpubr::do_CopyNumberVariantPlot` to silence the messages when there are not enough genes in the chromosome to perform the analysis.
- Fixed bug in `SCpubr::do_GroupwiseDEPlot` in which the heatmap could not be computed.
- Added extra checks to ensure proper input in `SCpubr::do_GroupwiseDEPlot`. ([#3](https://github.com/enblacar/SCpubr/issues/3))
- Changed default value of `order` in `SCpubr::do_FeaturePlot()` from `TRUE` to `FALSE`.
- Added `min.cutoff` and `max.cutoff` parameters to `SCpubr::do_FeaturePlot()`. This allows to effectively subset the color scale to the values provided. Cells outside the range will be converted to the min or max values provided, respectively. ([#2](https://github.com/enblacar/SCpubr/issues/2))
- Fixed a typo that prevented labels to be bold in `SCpubr::do_DimPlot()` when cell borders are displayed.
- Fixed a typo that made the lines in `panel.grid.minor` to be displayed in `SCpubr::do_Ridgeplot()`.
- Added `split.by` functionality to `SCpubr::do_ViolinPlot()`. ([#4](https://github.com/enblacar/SCpubr/issues/4))
- Added `group.by` and `split.by` functionality to `SCpubr::do_DimPlot()`. ([#4](https://github.com/enblacar/SCpubr/issues/4)) 

# SCpubr v1.0.1

-   Rework on unit tests and examples so that it can pass CRAN R CMD Check without packages in Suggests. This is, to make sure all Suggested packages are used conditionally.

This version can be obtained in the [v1.0.1 release](https://github.com/enblacar/SCpubr/releases/tag/v1.0.1) release or **downloading it from CRAN** using:

``` r
install.packages("SCpubr")
```

## SCpubr v1.0.1-dev-stable

-   Same as v1.0.1, but with all the functions that do not pass CRAN checks. These functions are: `SCpubr::save_Plot()` `SCpubr::do_LigandReceptorPlot()` and `SCpubr::do_SankeyPlot()`.

This version can be obtained in the [v1.0.1-dev-stable release](https://github.com/enblacar/SCpubr/releases/tag/v1.0.1-dev-stable) release.

# SCpubr v1.0.0

-   Modified internal checks so that the functions that do not use `Seurat` do not require this to be installed. This is just for the very side case in which somebody downloads the package just for the `SCpubr::do_ColorPalette()` function.
-   Removed the option to use `individual.titles`, 'individual.subtitles`and`individual.captions`in`SCpubr::do_NebulosaPlot()\` as the benefit of such parameters did not surpass the problems the code was causing. The feature might come back in the future, once fully optimized.
-   Removed `SCpubr::save_Plot()` function to align with CRAN policies that the package should not write to the file system. The code is still available in the v0.0.0.0.9000 release in Github.
-   Removed `SCpubr::do_LigandReceptorPlot()`, `SCpubr::do_SankeyPlot()` and `SCpubr::do_PseudotimePlot()` to align with CRAN policies and make it possible to publish the package. These functions can still be accessed in the v0.0.0.0.9000 release in Github.
-   Removed `SCpubr::do_PseudotimePlot()` for the reason above and because the dependency `Matrix.utils` was removed from CRAN on *09-10-2022*.

This version can be obtained in the [v1.0.0 release](https://github.com/enblacar/SCpubr/releases/tag/v1.0.0) release.

## SCpubr v1.0.0-dev-stable

-   Same as v1.0.0, but with all the functions that do not pass CRAN checks. These functions are: `SCpubr::save_Plot()` `SCpubr::do_LigandReceptorPlot()` and `SCpubr::do_SankeyPlot()`.

This version can be obtained in the [v1.0.0-dev-stable release](https://github.com/enblacar/SCpubr/releases/tag/v1.0.0-dev-stable) release.

# SCpubr 0.0.0.9000

-   Added a `NEWS.md` file to track changes to the package.
-   Prepare package for submission to CRAN.
