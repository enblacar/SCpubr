# SCpubr v2.0.0 

This major update focus on a complete re-implementation of all heatmap-based functions into `ggplot2` instead of `ComplexHeatmap`. This will lead to many of the existing code to break. The trade-off between the difficulty of debug, expand and maintain the existing heatmap-based functions with regards to the capabilities ComplexHeatmap offers with regards to ggplot2 was not worthy.

All heatmap-specific parameters have been replaced with the overarching parameters that are used across functions. This decision was taking after a lot of thought, but ultimately, having all plots rely on ggplot2 makes it way more compatible to work with them together, to debug, and to further implement new ideas.

Many (except a few selected cases) of the functions that returned list of different plots have been modified to return a single (and most important/relevant) plot and the option to return the Seurat object with the data generated added to it has been implemented so that the user can still generate plots with it. This goes in line with the fact that having so many interconnected functions made it very difficult to expand on them, if needed, as the downstream effects will cascade to other functions as well.

## Parameter renaming

-   Changed `viridis_color_map` to `viridis.palette`.
-   Changed `viridis_direction` to `viridis.direction`.
-   Changed `sequential_direction` to `sequential.direction`.
-   Changed `rotate_x_axis_labels` to `axis.text.x.angle`.
-   Changed `rotate_strip_text` to `strip.text.angle`.

## New functions (available on the development build for extended texting)

-   `SCpubr::do_MetadataPlot()` to generate metadata heatmaps with ease both from Seurat object or from a data frame. Will be first released as part of the `development version` and then released in CRAN as part of future updates. The idea is to gather feedback from users before officially releasing it.
-   `SCpubr::do_SCExpressionHeatmap()` to generate heatmaps of expression of genes across all cells in the dataset. Will be first released as part of the `development version` and then released in CRAN as part of future updates. The idea is to gather feedback from users before officially releasing it.
-   `SCpubr::do_SCEnrichementHeatmap()` to generate heatmaps of enrichment of genes across all cells in the dataset. Will be first released as part of the `development version` and then released in CRAN as part of future updates. The idea is to gather feedback from users before officially releasing it.
-   `SCpubr::do_AffinityAnalysisPlot()` to assess the affinity of gene sets to subset of cells in the Seurat objects using the weighted means algorithms from DecoupleR. Will be first released as part of the `development version` and then released in CRAN as part of future updates. The idea is to gather feedback from users before officially releasing it.
-   `SCpubr::do_LoadingsPlot()` to generate a summary heatmap of the PCA loadings (top and bottom scored genes for each PC) together with a expression heatmap of the same genes. Will be first released as part of the `development version` and then released in CRAN as part of future updates. The idea is to gather feedback from users before officially releasing it.
-   `SCpubr::do_DiffusionMapPlot()` to analyze the output of a diffusion map analysis on the context of enrichment in gene sets used for the generation of the diffusion map. Will be first released as part of the `development version` and then released in CRAN as part of future updates. The idea is to gather feedback from users before officially releasing it.
-   `SCpubr::check_dependencies()` to generate a per-function summary of the needed packages to run the function. The report has been enhanced with `cli` package and now clearly illustrates what is missing to run the function.

## Removed functions

-   `SCpubr::do_SankeyPlot()` has been removed and replaced by `SCpubr::do_AlluvialPlot()`, which is present in the official CRAN versions.
-   `SCpubr::do_PseudotimePlot()` has been removed indefinitely until a better, revamped, state-of-the-art version is generated.
-   `SCpubr::do_AzimuthAnalysisPlot()` has been removed as the output can be accomplished by a combination of the current functions in `SCpubr`. A vignette will be added to reproduce the same analysis.

## General

-   Now when using `min.cutoff` or `max.cutoff`, the legend will show that the min/max value is higher/lower than the one provided, if such value appeared originally in the legend breaks. This potentially interacts with `enforce_symmetry`.
-   Added `number.breaks` parameter to control the number of breaks in the legend of ggplot2-based plots. It will not always work, as the function will try to fit the breaks accordingly. But still, will give some range of freedom to the user.
-   Removed `colorsteps` from `legend.type` parameters as it was prone to generate unintended bugs.
-   Changed default values from `min.cutoff` and `max.cutoff` from `NULL` to `NA`.
-   Implemented `diverging.palette` parameter in all plots that have a symmetrical color scale to help selecting other possible color scales for the plot.
-   Implemented `sequential.palette` parameter in all plots that have a continuous, non-symmetrical color scale to help selecting other possible color scales for the plot, in the case the user does not want to use viridis color scales.
-   Renamed `SCpubr::state_dependencies()` to `SCpubr::check_dependencies()`.
-   Renewed printed messages at startup and while running functions using `cli` package.
-   Added the complete control of the font style of plot titles, subtitles, captions, axis titles, axis text, legend titles and legend text. For this, the following parameters have been added to all ggplot2-based functions:
    -   `plot.title.face`: To control the style of the **title**.
    -   `plot.subtitle.face`: To control the style of the **subtitle**.
    -   `plot.caption.face`: To control the style of the **caption**.
    -   `axis.title.face`: To control the style of the **axis title**.
    -   `axis.text.face`: To control the style of the **axis text**.
    -   `legend.title.face`: To control the style of the **legend title**.
    -   `legend.text.face`: To control the style of the **legend text**.
-   Changed default font style for legend text from `bold` to `plain`. 
-   Changed default font style for axis text from `bold` to `plain`.
-   When using `plot.axes = TRUE` parameter in `SCpubr::do_DimPlot()`, `SCpubr::do_FeaturePlot()` and `SCpubr::do_NebulosaPlot()`, now the entirety of the X and Y axis is removed, titles included.
-   Remove plot margin padding in `SCpubr::do_DimPlot()`, `SCpubr::do_FeaturePlot()` and `SCpubr::do_NebulosaPlot()`. 

## `SCpubr::do_AlluvialPlot`
-   Added `sequential.palette` and `sequential.direction` parameters.

## `SCpubr::do_BarPlot`

-   Added `facet.by` parameter to extra group the bars by a third metadata variable.
-   Added `facet.by.direction` parameter to decide in which direction the facets are drawn.
-   Added `order.by` to reorder the bars when using `position = fill` based on a value in `group.by`.
-   Limited the possible interactions from `group.by`, `split.by` and `order.by` to those that make sense to plot. For instance, a bar plot using `group.by` and `position = fill` but not using `split.by ` resulted in bars of equal lenght with only one value per group of proportion `1`.
-   Set default value of `plot.grid` to `FALSE`.
-   Added parameter `add.n` to display the total count on top when `position = fill`.
-   Added parameter `add.n.face` to control the appearance of the text displayed.
-   Added parameter `add.n.expand` to control the range of values in the Y axis. This has to be minimum 0 and maximum at least 1. This is set in order to tweak the limits so that the labels fit when `flip = TRUE`.

## `SCpubr::do_BeeSwarmPlot`

-   Added `order` parameter to reorder the groups based on the median rank.

## `SCpubr::do_BoxPlot`

-   Changed the reordering of boxplots based on the median rather than the mean.
-   Added `na.rm` to `geom_boxplot` to avoid unnecessary warnings when introducing NAs as part of the data.
-   Fixed a bug in which `order` would not work if `NAs` are in the data.
-   Changed default value of `boxplot.linewidth` from `1` to `0.5`.
-   Fixed a bug in which when using a combination of `group.by` and `split.by`, the package would check that the colors provided to `colors.use` need to match the values in `group.by` and not `split.by`.

## `SCpubr::do_CorrelationPlot`

-   Added parameter to fix a bug in which viridis scales did not apply due to the lack of the parameter.
-   Added `min.cutoff` and `max.cutoff` parameter to add cutoffs to the scales.
-   Added `mode = "jaccard"` to compute a correlation matrix of a list of gene sets based on jaccard similarity.
-   Added `use_viridis`, `sequential.palette` and `sequential_direction` and `diverging.palette` to control color palettes.
-   Added `cluster` parameter to toggle on/off the clustering of the rows and columns in the heatmap.
-   Added `remove.diagonal` parameter to toggle on/off the conversion of the diagonal in the correlation matrix to `NA`.
-   Fixed several issues with setting cutoffs for the color scale using `min.cutoff` and `max.cutoff`. 
-   Fixed an issue where `number.breaks` will not work in `mode = "jaccard"`.

## `SCpubr::do_CopyNumberVariantPlot()`

-   Removed the option to compute Feature and Geyser plots.
-   Instead, a new paramerter `return_object` has been added to return the Seurat object with a new assay containing the CNV scores per cell on the `data` slot of the `CNV_scores` assay.
-   The main output visualization is now a heatmap with the averaged scores by chromosome and groups and also by chromosome arms.

## `SCpubr::do_DimPlot`

-   Modified underlying code to correctly display borders around cells when `cells.highlight` or `idents.hightlight` or `idents.keep` are used. Also removed the "Not selected" item from the legend when doing so, as it was redundant.
-   Fixed a bug in which multiple legend would appear when using a combination of `group.by` and `split.by`, given that the individual UMAPs would not have the same number of entities to plot and color.

## `SCpubr::do_DotPlot`

-   Added `scale` parameter to allow for the data to be scaled or not scaled.
-   Removed `split.by` parameter in favor or the higher consistency and proper functionality accross parameters. Will probably come in the future, implemented outside of the umbrella of Seurat.
-   Renamed parameter `cluster.idents` to `cluster`.
-   Removed the limitation of `flip` when `features` was a list of genes. Now any combination of `flip` and `features` is possible.

## `SCpubr::do_EnrichmentHeatmap`

-   Removed options to plot FeaturePlots, GeyserPlots, ViolinPlots, etc. - together with its related parameters. For the sake of simplicity in the function and its use, the user can get the Seurat object back with `return_object = TRUE` and plot the enrichment scores separately, that are stored as a new Assay.
-   Removed `return_matrix` parameter as the scores can now be retrieved from the Seurat object as an assay.
-   Enforcing the use of `named lists` as input for the function.
-   Added `cluster` parameter to allow for clustering of rows and columns.
-   Added `groups.order` to allow for specifically arrange the groups defined by `group.by` with a given order.
-   Added `features.order` to allow for specifically arrange the gene sets defined by `input_gene_list`.

## `SCpubr::do_ExpressionHeatmap`

-   Added `cluster` parameter to allow for clustering of rows and columns.
-   Added `groups.order` to allow for specifically arrange the groups defined by `group.by` with a given order.
-   Added `features.order` to allow for specifically arrange the features defined by `features`.

## `SCpubr::do_FeaturePlot`

-   Modified underlying code to show a border around selected cells when using `split.by`, `cells.hightlight` and `idents.highlight`.
-   Added parameter `border.density` to reduce the amount of extra cells drawn on the background to generate the borders. This will be a number between 0 and 1 corresponding to the quantile of the distribution of density of the points in the scatterplot drawn in the background. The lower the value, the harder it will be to keep a border around all cells, while it will significantly reduce the overall weight of the plot object.
-   Added parameter `group.by`, that allows to plot a big dot in the center of each group designated by `group.by` and thus allowing to locate easily where each identity is in the FeaturePlot. Also, plots a legend matching the color of the dots. This can be tweaked with additional parameters such as:
-   `group.by.show.dots` to controlw hether these dots are plotted or not (to allow only plotting colored borders around cells - see below).
-   `group.by.dot.size` to control the size of the introduced dots.
-   `group.by.cell_border` to plot another contour-like border which also displays the color coding of the clusters designated by `group.by`, to signal the reach of each cluster. However, this basically signals the cluster the cells in the periphery of the cell blobs belong to. Take that into account.
-   `group.by.cell_borders.alpha` controls the alpha of the new cell borders.
-   `group.by.legend` controls the legend title of the new legend.
-   Renamed `split.by.idents` to `idents.keep` to better synergize with the parameter in `SCpubr::do_DimPlot`. Only works when `split.by` is used.

## `SCpubr::do_FunctionalAnnotationPlot`

-   Removed the tree plots as they proved to behave inconsistently across datasets and the quality of visualizations were compromised.
-   Removed the option to plot the bar plots and dot plots in the sake of a more simplified, streamlined plot generation.
-   The option to return the result matrix using `return_matrix` is added, so that the user can use it to compute further analysis or visualizations.

## `SCpubr::do_FunctionalAnnotationPlot`
-   Renamed `order_by_mean` to `order`.
-   Ordering using `order = TRUE` now is done based on the median instead of the mean.

## `SCpubr::do_LigandReceptorPlot()`
-   Modified the accepted input so that only the result of `liana::liana_aggregate()` is taken into account.
-   Removed `arrange_interactions_by` as now the function only accepts the output of `liana::liana_aggregate()`.
-   Added a `sort.by` parameter instead to select how the output of `liana::liana_aggregate()` should be ordered prior the subset by `top_interactions`. Five modes are available:
    - `A`: Orders the output by `specificity`.
    - `B`: Orders the output by `magnitude`.
    - `C`: Orders the output by `specificity` then `magnitude`. This prioritizes the `specificity` column.
    - `D`: Orders the output by `magnitude` then `specificity`. This prioritizes the `magnitude` column.
    - `E`: Orders the output by `specificity` and `magnitude` providing equal weights to both columns.
-   Removed `flip` parameter as the output was prone to errors.
-   Removed parameter `compute_ChordDiagrams` and added `return_interactions`. This parameter returns two tibbles that can be used alongside `SCpubr::do_ChordDiagramPlot` to plot the diagrams. 
-   Now the filtering applied by using `keep_source` and `keep_target` takes place before subsetting for the top N interactions defined by `top_interactions`. This ensures that, if the user wants to focus on a given interaction subset, we retrieve the most important interactions for the subset.
-   Added `magnitude` and `specificity` columns to allow the user to choose which variables to use for plotting.
-   Added `sorting.type.magnitude` and `sorting.type.specificity` to allow the user to choose how the columns are sorted prior plotting.
-   Added `invert_magnitude` and `invert_specificity` to allow the user to choose how the data is displayed for columns that tend to 0. Inverting performs a `-log10` transformation on the column.
-   Added a `verbose` parameter and set it to `TRUE` by default to inform the user of the arrangements taking place in the output of `liana::liana_aggregate()` prior plotting.

## `SCpubr::do_PathwayActivityPlot()`

-   Removed the option to plot geyser and feature plots to simplify the use (and computational time) of the function.
-   Introduced `return_object` parameter that returns the Seurat object with the new assay to use for other plotting purposes (such as Geyser and Feature plots).
-   Removed options to plot FeaturePlots, GeyserPlots - together with its related parameters. For the sake of simplicity in the function and its use, the user can get the Seurat object back with `return_object = TRUE` and plot the scores separately.
-   Added `slot` parameter to decide whether to plot scale data or not.
-   Fixed bug in which after setting `enforce_symmetry = FALSE` the color palette used was `diverging.palette` instead.

## `SCpubr::do_TFActivityPlot()`

-   Removed the option to plot geyser and feature plots to simplify the use (and computational time) of the function.
-   Introduced `return_object` parameter that returns the Seurat object with the new assay to use for other plotting purposes (such as Geyser and Feature plots).
-   Removed options to plot FeaturePlots, GeyserPlots - together with its related parameters. For the sake of simplicity in the function and its use, the user can get the Seurat object back with `return_object = TRUE` and plot the scores separately.
-   Added `slot` parameter to decide whether to plot scale data or not.
-   Fixed bug in which after setting `enforce_symmetry = FALSE` the color palette used was `diverging.palette` instead.

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
