# SCpubr v1.0.2.9000
## General
-   Added `min.cutoff` and `max.cutoff` parameter to effectively subset the color scale and remove the effect of extreme outliers in all ComplexHeatmap-based functions.
-   Added `min.cutoff` and `max.cutoff` parameter to effectively subset the color scale and remove the effect of extreme outliers in all ggplot2-based functions susceptible to be biased by outliers.


## `SCpubr::do_DimPlot()`
-   Implemented a change in which when using `split.by` and `group.by` in combination, the cells colored on top of the UMAP also have a border.
-   Implemented a bug-fix in which when using `split.by` and `group.by` in combination, the extra new layers would not raster if `raster = TRUE`.
-   Implemented a bug-fix in which when using `split.by` and `group.by` in combination, no plots will appear if `ncol` is set.
-   Implemented a new feature to add density line contours using `plot_density_contour`.
-   Implemented the conditional use of `raster.dpi` to Seurat versions higher or equal to 4.1.0.

## `SCpubr::do_EnrichmentHeatmap()`
-   Implemented a bug fix for internal checks in the function.
-   Added `plot_FeaturePlots` and `plot_GeyserPlots` to also report the enrichment scores in a gene set-based manner.
-   Added `flavor` parameter, that accepts `Seurat` and `UCell` to allow for different enrichment scoring methods. It requires `R 4.2.0` to run.
-   Renamed `symmetrical_scale` to `enforce_symmetry` to have a greater coherence accross functions.

## `SCpubr::do_FeaturePlot()`
-   Implemented a new feature to add density line contours using `plot_density_contour`.
-   Implemented the conditional use of `raster.dpi` to Seurat versions higher or equal to 4.1.0.

## `SCpubr::do_GeyserPlot()`
-   Fixed bug in which internal parameter names made it to the X axis title.
-   Removed `color.by` implementation due to it being very buggy. This will be re-implemented in a future patch.

## `SCpubr::do_RidgePlot()`
-   Implemented a bug-fix in which using `assay = "RNA"` or, in fact, any other assay rather than `SCT` will result in an error.


## `SCpubr::do_ViolinPlot()`
-   Corrected a bug in which legend title when using `split.by` was an actual line of code.
-   Added `legend.title` parameter to control the title of the legend.

# SCpubr v1.0.2

## General changes

-   Change color palette when using `enforce_symmetry = TRUE` to have the middle color as `grey95` instead of the previous one, which made middle values seem closer to the positive end of the scale.
-   Modified internal structure of all functions to abide with [tidyselect v1.2.0 lifecycle modifications](https://tidyselect.r-lib.org/news/index.html#lifecycle-changes-1-2-0).
-   Modified `rotate_x_axis_labels` parameter in all functions that made use of it. Now, instead of accepting a `logical`, accepts a `numeric`: either `0`, `45` or `90`, corresponding to the degrees in which the X axis labels should be rotated. ([#5](https://github.com/enblacar/SCpubr/issues/5#issuecomment-1289203453))

## `SCpubr::do_CopyNumberVariantPlot`

-   Modified the code for `SCpubr::do_CopyNumberVariantPlot` to also report results for the whole chromosome as well as for each chromosome arm.
-   Include the `verbose` argument to `SCpubr::do_CopyNumberVariantPlot` to silence the messages when there are not enough genes in the chromosome to perform the analysis.

## `SCpubr::do_DimPlot()`
-   Fixed a typo that prevented labels to be bold in `SCpubr::do_DimPlot()` when cell borders are displayed.
-   Added `group.by` and `split.by` functionality to `SCpubr::do_DimPlot()`. ([#4](https://github.com/enblacar/SCpubr/issues/4))

## `SCpubr::do_DotPlot()`
-   Added ticks to axes.
-   Modified default colors to convey a better aesthetic. 

## `SCpubr::do_FeaturePlot()`
-   Fixed potential bugs in `SCpubr::do_FeaturePlot` when setting `enforce_symmetry = TRUE`.
-   Changed default value of `order` in `SCpubr::do_FeaturePlot()` from `TRUE` to `FALSE`.
-   Added `min.cutoff` and `max.cutoff` parameters to `SCpubr::do_FeaturePlot()`. This allows to effectively subset the color scale to the values provided. Cells outside the range will be converted to the min or max values provided, respectively. ([#2](https://github.com/enblacar/SCpubr/issues/2))

## `SCpubr::do_GeyserPlot()`
-   Added `flip` parameter.

## `SCpubr::do_GroupwiseDEPlot()`
-   Fixed bug in `SCpubr::do_GroupwiseDEPlot` in which the heatmap could not be computed. ([#3](https://github.com/enblacar/SCpubr/issues/3))
-   Added extra checks to ensure proper input in `SCpubr::do_GroupwiseDEPlot`. ([#3](https://github.com/enblacar/SCpubr/issues/3))

## `SCpubr::do_LigandReceptorPlot()` (development release)
- Changed parameter `x_labels_angle` to `rotate_x_axis_labels` to keep a consistent terminology.

## `SCpubr::do_RidgePlot()`
-   Fixed a typo that made the lines in `panel.grid.minor` to be displayed in `SCpubr::do_Ridgeplot()`.
-   Added `flip` parameter.

## `SCpubr::do_ViolinPlot()`
-   Added `split.by` functionality to `SCpubr::do_ViolinPlot()`. ([#4](https://github.com/enblacar/SCpubr/issues/4), [#5](https://github.com/enblacar/SCpubr/issues/5))
-   Added `flip` parameter.
-   Now multiple features can be queried ad the same time. ([#5](https://github.com/enblacar/SCpubr/issues/5#issuecomment-1289203453))
-   Changed `feature` parameter to `features`, to better reflect the multiple feature behavior.
-   Recreated `Seurat`'s `share.y.lims` behavior and set it to `share.y.lims` parameter. ([#5](https://github.com/enblacar/SCpubr/issues/5#issuecomment-1289203453))

## SCpubr v1.0.2-dev-stable

-   Same as v1.0.2, but with all the functions that do not pass CRAN checks. These functions are: `SCpubr::save_Plot()` `SCpubr::do_LigandReceptorPlot()` and `SCpubr::do_SankeyPlot()`.

This version can be obtained in the [v1.0.2-dev-stable release](https://github.com/enblacar/SCpubr/releases/tag/v1.0.2-dev-stable) release.

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
