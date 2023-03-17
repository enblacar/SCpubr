# SCpubr v2.0.0 (to be released)

**THIS IS A SNEAK PEEK OF THE UPDATE TO COME TO HAVE PEOPLE WARNED ABOUT THE CHANGES THAT WILL LIKELY BREAK THE EXISTING CODE**

**The update will still take some time (date of writing: 15.03.2023)**

This major update focus on a complete re-implementation of all heatmap-based functions into `ggplot2` instead of `ComplexHeatmap`. This will lead to many of the existing code to break. Hence, the major update.

## Parameter renaming

-   Changed `viridis_color_map` to `viridis.palette`.
-   Changed `viridis_direction` to `viridis.direction`.
-   Changed `sequential_direction` to `sequential.direction`.

## New functions

-   `SCpubr::do_MetadataPlot()` to generate metadata heatmaps with ease both from Seurat object or from a data frame.
-   `SCpubr::do_SCExpressionHeatmap()` to generate heatmaps of expression of genes across all cells in the dataset.
-   `SCpubr::do_SCEnrichementHeatmap()` to generate heatmaps of enrichment of genes across all cells in the dataset.
-   `SCpubr::do_AffinityAnalysisPlot()` to assess the affinity of gene sets to subset of cells in the Seurat objects using the weighted means algorithms from [DecoupleR](https://saezlab.github.io/decoupleR/).

## General

-   Now when using `min.cutoff` or `max.cutoff`, the legend will show that the min/max value is higher/lower than the one provided, if such value appeared originally in the legend breaks. This potentially interacts with `enforce_symmetry`.
-   Added `number.breaks` parameter to control the number of breaks in the legend of ggplot2-based plots. It will not always work, as the function will try to fit the breaks accordingly. But still, will give some range of freedom to the user.
-   Removed `colorsteps` from `legend.type` parameters as it was prone to generate unintended bugs.
-   Changed default values from `min.cutoff` and `max.cutoff` from `NULL` to `NA`.
-   Implemented `diverging.palette` parameter in all plots that have a symmetrical color scale to help selecting other possible color scales for the plot.
-   Implemented `sequential.palette` parameter in all plots that have a continuous, non-symmetrical color scale to help selecting other possible color scales for the plot, in the case the user does not want to use viridis color scales.
-   Re-implemented all heatmap-based plots in ggplot2. This means that this is a **breaking change**. All heatmap-specific parameters have been replaced with the overarching parameters that are used across functions. This decision was taking after a lot of thought, but ultimately, having all plots rely on ggplot2 makes it way more compatible to work with them together, to debug, and to further implement new ideas.

## `SCpubr::do_BeeSwarmPlot`

-   Added `order` parameter to reorder the groups based on the median rank.

## `SCpubr::do_BoxPlot`

-   Change the reordering of boxplots based on the median rather than the mean.
-   Added `na.rm` to `geom_boxplot` to avoid unnecessary warnings when introducing NAs as part of the data.
-   Fixe bug in which `order` would not work if `NAs` are in the data.

## `SCpubr::do_CorrelationPlot`

-   Added parameter to fix a bug in which viridis scales did not apply due to the lack of the parameter.
-   Added `min.cutoff` and `max.cutoff` parameter to add cutoffs to the scales.
-   Added `mode = "jaccard"` to compute a correlation matrix of a list of gene sets based on jaccard similarity.
-   Added `use_viridis`, `sequential.palette` and `sequential_direction` and `diverging.palette` to control color palettes.

## `SCpubr::do_CopyNumberVariantPlot()`
-  Removed the option to compute Feature and Geyser plots.
-  Instead, a new paramerter `return_object` has been added to return the Seurat object with a new assay containing the CNV scores per cell on the `data` slot of the `CNV_scores` assay. 
-  The main output visualization is now a heatmap with the averaged scores by chromosome and groups and also by chromosome arms.

## `SCpubr::do_DimPlot`

-   Modified underlying code to correctly display borders around cells when `cells.highlight` or `idents.hightlight` are used. Also removed the "Not selected" item from the legend when doing so, as it was redundant.

## `SCpubr::do_DotPlot`

-   Removed `split.by` parameter in favor or the higher consistency and proper functionality accross parameters. Will probably come in the future, implemented outside of the umbrella of Seurat.

## `SCpubr::do_EnrichmentHeatmap`

-   Removed options to plot FeaturePlots, GeyserPlots, ViolinPlots, etc. - together with its related parameters. For the sake of simplicity in the function and its use, the user can get the Seurat object back with `return_object = TRUE` and plot the enrichment scores separately.

## `SCpubr::do_FeaturePlot`

-   Modified underlying code to show a border around selected cells when using `split.by`, `cells.hightlight` and `idents.highlight`.
-   Added parameter `border.density` to reduce the amount of extra cells drawn on the background to generate the borders. This will be a number between 0 and 1 corresponding to the quantile of the distribution of density of the points in the scatterplot drawn in the background. The lower the value, the harder it will be to keep a border around all cells, while it will significantly reduce the overall weight of the plot object.
-   Added parameter `group.by`, that allows to plot a big dot in the center of each group designated by `group.by` and thus allowing to locate easily where each identity is in the FeaturePlot. Also, plots a legend matching the color of the dots. This can be tweaked with additional parameters such as:
    -   `group.by.show.dots` to controlw hether these dots are plotted or not (to allow only plotting colored borders around cells - see below).
    -   `group.by.dot.size` to control the size of the introduced dots.
    -   `group.by.cell_border` to plot another contour-like border which also displays the color coding of the clusters designated by `group.by`, to signal the reach of each cluster. However, this basically signals the cluster the cells in the periphery of the cell blobs belong to. Take that into account.
    -   `group.by.cell_borders.alpha` controls the alpha of the new cell borders.
    -   `group.by.legend` controls the legend title of the new legend.

## `SCpubr::do_FunctionalAnnotationPlot`

-   Removed the tree plots as they proved to behave inconsistently accross datasets and the quality of visualizations were compromised.
-   Removed the option to plot the bar plots and dot plots in the sake of a more simplified, streamlined plot generation.
-   The option to return the result matrix using `return_matrix` is added, so that the user can use it to compute further analysis or visualizations.

## `SCpubr::do_GroupwiseDEPlot`

-   Added `set_min_expression_to_zero` parameter to set the viridis color scale to 0 in the expression heatmap.

## `SCpubr::do_PathwayActivityPlot`

-   Removed the option to plot geyser and feature plots to simplify the use (and computational time) of the function.
-   Introduced `return_object` parameter that returns the Seurat object with the new assay to use for other plotting purposes (such as Geyser and Feature plots).
-   Removed options to plot FeaturePlots, GeyserPlots - together with its related parameters. For the sake of simplicity in the function and its use, the user can get the Seurat object back with `return_object = TRUE` and plot the scores separately.

## `SCpubr::do_TFActivityPlot()`

-   Removed the option to plot geyser and feature plots to simplify the use (and computational time) of the function.
-   Introduced `return_object` parameter that returns the Seurat object with the new assay to use for other plotting purposes (such as Geyser and Feature plots).
-   Removed options to plot FeaturePlots, GeyserPlots - together with its related parameters. For the sake of simplicity in the function and its use, the user can get the Seurat object back with `return_object = TRUE` and plot the scores separately.

# SCpubr v1.1.2

More hotfixes in unit tests to comply with CRAN checks.

# SCpubr v1.1.1

Hotfixes in unit tests to comply with CRAN checks.

# SCpubr v1.1.0

## General

-   Increased the cell size of all heatmap-based functions from 5 to 8.
-   Decreased the thickness of frame and ticks of all ggplot2-based continuous legends to retrieve a similar behavior as in previous versions of ggplot2, as with the new update, the overall thickness of the frame and ticks increased, probably due to the changes related to `element_line`,
-   Added five new functions: `do_AlluvialPlot()`, `do_AzimuthAnalysisPlot()`, `do_ExpressionHeatmap()`, `do_GroupedGOTermPlot()` and `do_FunctionalAnnotationPlot()`.
-   Added `legend.ncol`, `legend.nrow`, `legend.title` and `legend.byrow` to as many functions as possible to further customize legend appearance.

## `SCpubr::do_BeeSwarmPlot()`

-   Added `min.cutoff` and `max.cutoff` parameter.
-   Added ticks to the plot, that were missing.
-   Added missing axes titles.
-   Added `viridis_direction` parameter to control how the continuous color scale is formed.
-   Added `return_object` parameter to return the Seurat object with the enrichment scores computed.
-   Added BoxPlots, BeeSwarmPlots and ViolinPlots to the possible outputs the user can choose from.
-   Make `legend.position` conditional of whether `continuous_feature` is set to TRUE. If it is false, legend is not displayed unless the user specifies otherwise.

## `SCpubr::do_BarPlot()`

-   Fixed a bug in which axes titles were not displaying correctly under certain combinations of `flip` and `split.by`.
-   Fixed a bug in which `x_lab` and `y_lab` would not rotate accordingly when using `flip = TRUE`.

## `SCpubr::do_BeeSwarmPlot()`

-   Adapted the code to the new 0.7.1 version of the package, thus deprecating the `groupOnX` parameter of `geom_quarirandom`. This will likely affect users with a lower version.
-   A warning has been placed for the users in lower versions of the need to upgrade to 0.7.1.
-   This changes are subject to the new behaviors/deprecations of ggplot2 and ggplot2.

## `SCpubr::do_BoxPlot()`

-   Set `assay` to NULL and will default to the default assay in the seurat object.

## `SCpubr::do_CellularStatesPlot()`

-   Fixed a bug that prevented FeaturePlots to have symmetrical axes with respect to the main plot.

## `SCpubr::do_CorrelationPlot()`

-   Added `viridis_direction` parameter.

## `SCpubr::do_DimPlot()`

-   Fixed a bug in which the legend title will not show up in regular basic plots even though the parameter `legend.title` was used.
-   Completely reformatted the way `split.by` works, so that now only one legend is displayed for the whole group and cells have border.
-   Added `label.size` and `label.box` parameters for further customize the appearance of the plot when using `label = TRUE`.
-   Changed `repel` to `FALSE` by default.

## `SCpubr::do_EnrichmentHeatmap()`

-   Fixed a bug in the code that prevented the feature plots and the geyser plots to be computed if the input gene list was not a named list of genes.
-   Added `flavor = "AUCell"`, that lets the user compute AUCell scoring of the gene sets.
-   Added the option to query multiple `group.by` parameters at the same time.
-   Fixed a bug in the code that prevented multiple outputs with different values of `group.by` to be returned properly, leading to the last value of `group.by` replacing all the rest.

## `SCpubr::do_FeaturePlot()`

-   Added `label`, `label.size` and `label.color` parameter to reproduce the same behavior as in `Seurat::FeaturePlot()`.

## `SCpubr::do_GroupwiseDEPlot()`

-   Set `assay` to NULL and will default to the default assay in the seurat object.

## `SCpubr::do_LigandReceptorPlot()`

-   Added `arrange_interactions_by` to control how output interactions are arranged (either by aggregate_rank, specificity, magnitude or a combination of magnitude and specificity).
-   Added `sort_interactions_alphabetically` to control whether the output dotplot has the interactions ordered alphabetically or as they come in the original matrix (meaning, they follow the arrangement specified in `arrange_interactions_by`). (([liana's issue #72](https://github.com/saezlab/liana/issues/72)))

## `do_PathwayActivityPlot()`

-   Added a fix in which when `enforce_symmetry` is set to `FALSE`, then the color scale turns into a viridis-based one instead of a two-color gradient scale.

## `do_TFActivityPlot()`

-   Added a fix in which when `enforce_symmetry` is set to `FALSE`, then the color scale turns into a viridis-based one instead of a two-color gradient scale.

## `SCpubr::do_ViolinPlot()`

-   Fixed a bug in the code in which no different colors could be passed to `colors.use`.
-   Reduced default line width from 1 to 0.5.

# SCpubr v1.0.4

-   Hotfix for v1.0.3 in which `do_GeyserPlot` with categorical variables had a bug that mapped the legend to the continuous axis.

# SCpubr v1.0.3

## General changes

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
-   Renamed `symmetrical_scale` to `enforce_symmetry` to have a greater coherence across functions.

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

### SCpubr v.1.0.3-dev-stable

-   Same as v1.0.3, but with all the functions that do not pass CRAN checks. These functions are: `SCpubr::save_Plot()` `SCpubr::do_LigandReceptorPlot()` and `SCpubr::do_SankeyPlot()`.

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

-   Changed parameter `x_labels_angle` to `rotate_x_axis_labels` to keep a consistent terminology.

## `SCpubr::do_RidgePlot()`

-   Fixed a typo that made the lines in `panel.grid.minor` to be displayed in `SCpubr::do_Ridgeplot()`.
-   Added `flip` parameter.

## `SCpubr::do_ViolinPlot()`

-   Added `split.by` functionality to `SCpubr::do_ViolinPlot()`. ([#4](https://github.com/enblacar/SCpubr/issues/4), [#5](https://github.com/enblacar/SCpubr/issues/5))
-   Added `flip` parameter.
-   Now multiple features can be queried ad the same time. ([#5](https://github.com/enblacar/SCpubr/issues/5#issuecomment-1289203453))
-   Changed `feature` parameter to `features`, to better reflect the multiple feature behavior.
-   Recreated `Seurat`'s `share.y.lims` behavior and set it to `share.y.lims` parameter. ([#5](https://github.com/enblacar/SCpubr/issues/5#issuecomment-1289203453))

### SCpubr v1.0.2-dev-stable

-   Same as v1.0.2, but with all the functions that do not pass CRAN checks. These functions are: `SCpubr::save_Plot()` `SCpubr::do_LigandReceptorPlot()` and `SCpubr::do_SankeyPlot()`.

# SCpubr v1.0.1

-   Rework on unit tests and examples so that it can pass CRAN R CMD Check without packages in Suggests. This is, to make sure all Suggested packages are used conditionally.

## SCpubr v1.0.1-dev-stable

-   Same as v1.0.1, but with all the functions that do not pass CRAN checks. These functions are: `SCpubr::save_Plot()` `SCpubr::do_LigandReceptorPlot()` and `SCpubr::do_SankeyPlot()`.

# SCpubr v1.0.0

-   Modified internal checks so that the functions that do not use `Seurat` do not require this to be installed. This is just for the very side case in which somebody downloads the package just for the `SCpubr::do_ColorPalette()` function.
-   Removed the option to use `individual.titles`, 'individual.subtitles`and`individual.captions`in`SCpubr::do_NebulosaPlot()\` as the benefit of such parameters did not surpass the problems the code was causing. The feature might come back in the future, once fully optimized.
-   Removed `SCpubr::save_Plot()` function to align with CRAN policies that the package should not write to the file system. The code is still available in the v0.0.0.0.9000 release in Github.
-   Removed `SCpubr::do_LigandReceptorPlot()`, `SCpubr::do_SankeyPlot()` and `SCpubr::do_PseudotimePlot()` to align with CRAN policies and make it possible to publish the package. These functions can still be accessed in the v0.0.0.0.9000 release in Github.
-   Removed `SCpubr::do_PseudotimePlot()` for the reason above and because the dependency `Matrix.utils` was removed from CRAN on *09-10-2022*.

## SCpubr v1.0.0-dev-stable

-   Same as v1.0.0, but with all the functions that do not pass CRAN checks. These functions are: `SCpubr::save_Plot()` `SCpubr::do_LigandReceptorPlot()` and `SCpubr::do_SankeyPlot()`.

# SCpubr 0.0.0.9000

-   Added a `NEWS.md` file to track changes to the package.
-   Prepare package for submission to CRAN.
