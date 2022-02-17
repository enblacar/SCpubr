# Bee Swarm plots

This one is a very interesting plot. It stems from the idea that we can order (rank) the cells in a given variable. This variable has to be a continuous variable, for a better representation. The order goes from lowest to maximum value. Then, the cells are grouped into any other variable of interest and displayed in a scatter plot fashion. This is achieved by using the [ggbeeswarm package](https://github.com/eclarke/ggbeeswarm). The idea of using the `ggbeeswarm::geom_quasirandom()` geometry provided by this package and implement it for single-cell analyses came from [this tutorial from the Broad Institute](https://broadinstitute.github.io/2019_scWorkshop/functional-pseudotime-analysis.html#diffusion-map-pseudotime).

## Using categorical variables
Let's say we want to focus on how much each cluster is driven by the PC_1 and PC_2. The first thought is to just use `SCpubr::do_Dimplot()` to plot the PCA embedding instead of the UMAP. We also query PC_3 and PC_4 to have a not-so-clear example.



```r
p1 <- SCpubr::do_DimPlot(sample, reduction = "pca", label = T, legend = F, dims = c(1, 2)) 
p2 <- SCpubr::do_DimPlot(sample, reduction = "pca", label = T, legend = F, dims = c(3, 4)) 

p1 | p2
```

With this, we get right away a decent overview. Clusters 0, 5, 7 and 8 separate on PC_1 from the rest. However, in many cases this will not be clear, such as the image on the right. This is where Bee Swarm plots come in handy. This is implemented in `SCpubr::do_BeeSwarmPlot()`. This function needs the user to provide:
- The variable to rank to `feature_to_rank`.
- The groups to divide the plot into to `group.by`.
- Whether the output should be colored with a categorical or continuous scale, with `continuous_feature`.


```r
p1 <- SCpubr::do_DimPlot(sample, reduction = "pca", label = T, legend = F, dims = c(1, 2))
p2 <- SCpubr::do_BeeSwarmPlot(sample = sample, feature_to_rank = "PC_1", group.by = "seurat_clusters", continuous_feature = F)
p3 <- SCpubr::do_DimPlot(sample, reduction = "pca", label = T, legend = F, dims = c(3, 4)) 
p4 <- SCpubr::do_BeeSwarmPlot(sample = sample, feature_to_rank = "PC_4", group.by = "seurat_clusters", continuous_feature = F)

(p1 | p2) / (p3 | p4)
```

Here, we have selected PC_1 and PC_4. We can observe how the X axis of the Bee Swarm plot displays the ordering (rank) of all of the cells across the selected feature. Focusing on PC_1, we can see that cluster 0 is completely shifted to the right on PC_1, with is nicely displayed in the Bee Swarm plot by having all of the cells also ranked high (the higher the rank, the bigger the "value" of the feature to rank, in this case, the PC_1 value). In the case of PC_4, the Bee Swarm plot nicely shows which clusters lay on the upper, lower or middle part of the PC_4.

A very important thing to note in these kind of plots is that no cells will have the same rank. This is, imagine a scenario like PC_4, but we artificially remove clusters 0, 3, 5, 7, 8, 9, leaving only those forming a "straight line" in PC_4. The nature of this plot will also separate the remaining clusters:


```r
# Clusters to exclude.
clusters_exclude <- c("0", "3", "5", "7", "8", "9")

# Keep the original coloring.
cols.use <- colortools::setColors("steelblue", length(levels(sample)))
names(cols.use) <- levels(sample)

p1 <- SCpubr::do_DimPlot(sample[, !(sample$seurat_clusters %in% clusters_exclude)], reduction = "pca", label = T, legend = F, dims = c(3, 4), colors.use = cols.use) 
p2 <- SCpubr::do_BeeSwarmPlot(sample = sample[, !(sample$seurat_clusters %in% clusters_exclude)], feature_to_rank = "PC_4", group.by = "seurat_clusters", continuous_feature = F, colors.use = cols.use)

p1 | p2
```


See, we still clearly see two groups, formed by clusters 1 and 2, and clusters 4 and 6. We could even remove clusters 1 and 2 and still see a similar effect.


```r
# Clusters to exclude.
clusters_exclude <- c("0", "1", "2", "3", "5", "7", "8", "9")

# Keep the original coloring.
cols.use <- colortools::setColors("steelblue", length(levels(sample)))
names(cols.use) <- levels(sample)

p1 <- SCpubr::do_DimPlot(sample[, !(sample$seurat_clusters %in% clusters_exclude)], reduction = "pca", label = T, legend = F, dims = c(3, 4), colors.use = cols.use) 
p2 <- SCpubr::do_BeeSwarmPlot(sample = sample[, !(sample$seurat_clusters %in% clusters_exclude)], feature_to_rank = "PC_4", group.by = "seurat_clusters", continuous_feature = F, colors.use = cols.use)

p1 | p2
```


As can be seen here, both clusters now span all X axis. The cells have still ranked, therefore showing a cloud of dots. With this, we would just want that, as with any data visualization technique, each plot comes with a set of benefits and caveats. This visualization suffers from trying to plot highly similar values. Therefore, it is key to **understand the nature of the variable you want to rank** beforehand. 

## Using continuous variables

There are also scenarios in which we want to rank the cells to a continuous variable, but instead of showing colors for each group (which is anyway depicted in the Y axis), we want to introduce a continuous color scale. This is specially interesting to assess enrichment of clusters towards a given set of features. 



```r
# Set up list of a genes to compute enrichment. Let's use a monocyte signature.
genes.use <- c("CD14", "LYZ")

# Compute enrichment and rename the output.
sample <- Seurat::AddModuleScore(sample, features = genes.use, name = "Monocyte_signature")
sample$Monocyte_signature <- sample$Monocyte_signature1
sample$Monocyte_signature1 <- NULL

p1 <- SCpubr::do_DimPlot(sample, label = T, legend = F)
p2 <- SCpubr::do_FeaturePlot(sample, features = "Monocyte_signature") 
p3 <- SCpubr::do_BeeSwarmPlot(sample, feature_to_rank = "Monocyte_signature", group.by = "seurat_clusters", continuous_feature = T)
p1 | p2 | p3
```

By using this combination of figures, we can also assess that the monocyte signature seems to be predominantly enriched in clusters 0 and 7.

## Modify color maps for continuous variables
Same as in `SCpubr::do_FeaturePlot()`, it is also change the color map of the plot to one of the eight possible ones defined in [viridis](https://cran.r-project.org/web/packages/viridis/vignettes/intro-to-viridis.html). This is achieved by using `viridis_color_map` parameter and the color map name or code (capital letter). Options are:

- A - magma color map.
- B - inferno color map.
- C - plasma color map.
- D - viridis color map.
- E - cividis color map.
- F - rocket color map.
- G - mako color map.
- H - turbo  color map.


```r
p1 <- SCpubr::do_BeeSwarmPlot(sample = sample, feature_to_rank = "Monocyte_signature", group.by = "seurat_clusters", continuous_feature = TRUE, fontsize = 10, viridis_color_map = "A", plot.title = "Magma")
p2 <- SCpubr::do_BeeSwarmPlot(sample = sample, feature_to_rank = "Monocyte_signature", group.by = "seurat_clusters", continuous_feature = TRUE, fontsize = 10, viridis_color_map = "B", plot.title = "Inferno")
p3 <- SCpubr::do_BeeSwarmPlot(sample = sample, feature_to_rank = "Monocyte_signature", group.by = "seurat_clusters", continuous_feature = TRUE, fontsize = 10, viridis_color_map = "C", plot.title = "Plasma")
p4 <- SCpubr::do_BeeSwarmPlot(sample = sample, feature_to_rank = "Monocyte_signature", group.by = "seurat_clusters", continuous_feature = TRUE, fontsize = 10, viridis_color_map = "D", plot.title = "Viridis")
p5 <- SCpubr::do_BeeSwarmPlot(sample = sample, feature_to_rank = "Monocyte_signature", group.by = "seurat_clusters", continuous_feature = TRUE, fontsize = 10, viridis_color_map = "E", plot.title = "Cividis")
p6 <- SCpubr::do_BeeSwarmPlot(sample = sample, feature_to_rank = "Monocyte_signature", group.by = "seurat_clusters", continuous_feature = TRUE, fontsize = 10, viridis_color_map = "F", plot.title = "Rocket")
p7 <- SCpubr::do_BeeSwarmPlot(sample = sample, feature_to_rank = "Monocyte_signature", group.by = "seurat_clusters", continuous_feature = TRUE, fontsize = 10, viridis_color_map = "G", plot.title = "Mako")
p8 <- SCpubr::do_BeeSwarmPlot(sample = sample, feature_to_rank = "Monocyte_signature", group.by = "seurat_clusters", continuous_feature = TRUE, fontsize = 10, viridis_color_map = "H", plot.title = "Turbo")

p <- patchwork::wrap_plots(list(p1, p2, p3, p4, p5, p6, p7, p8), ncol = 2, byrow = TRUE)
p
```

