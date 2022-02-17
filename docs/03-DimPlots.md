# Dim plots

DimPlots are, probably, one of the most iconic visualizations from `Seurat`. It allows the user to visualize the cells in a dimensional reduction embedding such as `PCA` or `UMAP`. The cells can be, then, colored by any desired groups. In short, this visualization allows the user to plot **any kind of categorical data** onto the cells in the dimesional reduction embedding. This is achieved by using `Seurat::DimPlot()` funtion:


```r
Seurat::DimPlot(sample)
```

Overall, this is a pretty neat visualization, but there are quite some changes that one would like to implement. For instance, shuffling the cells so that there is no overlap of cells just due to the cluster names.


```r
Seurat::DimPlot(sample, shuffle = T)
```

Furthermore, one would think about the *need* of the axes. If, by consensus, UMAPs are shown plotting the first UMAP component on the X axis and the second on the Y axis, then showing them becomes redundant, specially when one can not truly rely on the numbers shown by the scales. 


```r
Seurat::DimPlot(sample, shuffle = T) + Seurat::NoAxes()
```

Right now, we can observe a couple of things. First, is that the dot size is rather small. Let's set it to 0.5.


```r
Seurat::DimPlot(sample, shuffle = T, pt.size = 0.5) + Seurat::NoAxes()
```

Still, the legend seems rather small. Let's increase it's font size and set it to bold so that it can better read.


```r
Seurat::DimPlot(sample, shuffle = T, pt.size = 0.5) + 
  Seurat::NoAxes() +
  ggplot2::theme(legend.text = ggplot2::element_text(size = 16, face = "bold"),
                 legend.title = ggplot2::element_text(size = 16, face = "bold"))
```

We would also like to add a title to our plot, to best describe it.


```r
Seurat::DimPlot(sample, shuffle = T, pt.size = 0.5) + 
  Seurat::NoAxes() +
  ggplot2::ggtitle("My awesome SC dataset") +
  ggplot2::theme(legend.text = ggplot2::element_text(size = 16, face = "bold"),
                 legend.title = ggplot2::element_text(size = 16, face = "bold"))
```

And, naturally, we would like to increase the font size of the title and put it in bold and centered.


```r
Seurat::DimPlot(sample, shuffle = T, pt.size = 0.5) + 
  Seurat::NoAxes() +
  ggplot2::ggtitle("My awesome SC dataset") +
  ggplot2::theme(plot.title = ggplot2::element_text(size = 20, face = "bold", hjust = 0.5),
                 legend.text = ggplot2::element_text(size = 16, face = "bold"),
                 legend.title = ggplot2::element_text(size = 16, face = "bold"))
```

Now, we would like to modify the color palette. This palette has too bright colors, and we want something more toned down.


```r
num_clusters <- length(unique(sample$seurat_clusters))
color_scale <- colortools::setColors("#457b9d", num_clusters)
names(color_scale) <- sort(unique(sample$seurat_clusters))

Seurat::DimPlot(sample, shuffle = T, pt.size = 0.5, cols = color_scale) + 
  Seurat::NoAxes() +
  ggplot2::ggtitle("My awesome SC dataset") +
  ggplot2::theme(plot.title = ggplot2::element_text(size = 20, face = "bold", hjust = 0.5),
                 legend.text = ggplot2::element_text(size = 16, face = "bold"),
                 legend.title = ggplot2::element_text(size = 16, face = "bold"))
```

The legend on the right side seems off, what if we were to have long cluster names? It would inevitable take a lot of space from the actual plot. Let's better put it on the bottom.


```r
num_clusters <- length(unique(sample$seurat_clusters))
color_scale <- colortools::setColors("#457b9d", num_clusters)
names(color_scale) <- sort(unique(sample$seurat_clusters))

Seurat::DimPlot(sample, shuffle = T, pt.size = 0.5, cols = color_scale) + 
  ggpubr::theme_pubr(legend = "bottom") +
  Seurat::NoAxes() +
  ggplot2::ggtitle("My awesome SC dataset") +
  ggplot2::theme(plot.title = ggplot2::element_text(size = 20, face = "bold", hjust = 0.5),
                 legend.text = ggplot2::element_text(size = 16, face = "bold"),
                 legend.title = ggplot2::element_text(size = 16, face = "bold"))
```

Still, there are too many columns in the legend. Let's rearrange it into four columns.


```r
num_clusters <- length(unique(sample$seurat_clusters))
color_scale <- colortools::setColors("#457b9d", num_clusters)
names(color_scale) <- sort(unique(sample$seurat_clusters))

Seurat::DimPlot(sample, shuffle = T, pt.size = 0.5, cols = color_scale) + 
  ggpubr::theme_pubr(legend = "bottom") +
  Seurat::NoAxes() +
  ggplot2::ggtitle("My awesome SC dataset") +
  ggplot2::theme(plot.title = ggplot2::element_text(size = 20, face = "bold", hjust = 0.5),
                 legend.text = ggplot2::element_text(size = 16, face = "bold"),
                 legend.title = ggplot2::element_text(size = 16, face = "bold")) + 
  ggplot2::guides(color = ggplot2::guide_legend(ncol = 4,
                                                byrow = F))
```

But now the legend icon sizes are too small! We have to fix this.


```r
num_clusters <- length(unique(sample$seurat_clusters))
color_scale <- colortools::setColors("#457b9d", num_clusters)
names(color_scale) <- sort(unique(sample$seurat_clusters))

Seurat::DimPlot(sample, shuffle = T, pt.size = 0.5, cols = color_scale) + 
  ggpubr::theme_pubr(legend = "bottom") +
  Seurat::NoAxes() +
  ggplot2::ggtitle("My awesome SC dataset") +
  ggplot2::theme(plot.title = ggplot2::element_text(size = 20, face = "bold", hjust = 0.5),
                 legend.text = ggplot2::element_text(size = 16, face = "bold"),
                 legend.title = ggplot2::element_text(size = 16, face = "bold")) + 
  ggplot2::guides(color = ggplot2::guide_legend(ncol = 4,
                                                byrow = F,
                                                override.aes = list(size = 4)))
```

As of now, this plot looks much better and polished than the default counterpart. This, is the setting ground for `SCpubr::do_DimPlot()`. 

## Regular DimPlots

This is the default output from `SCpubr::do_DimPlot()`. 


```r
SCpubr::do_DimPlot(sample)
```

We can add a title with the `plot.title` parameter.


```r
SCpubr::do_DimPlot(sample, plot.title = "My awesome SC data set")
```

We can change the legend location and number of columns with `legend.position` and `legend.ncol`.


```r
SCpubr::do_DimPlot(sample, 
                   plot.title = "My awesome SC data set", 
                   legend.position = "left", 
                   legend.ncol = 4)
```


## Highlighting cells

One of the nice features of `Seurat::DimPlot()` is the possibility of highlighting a certain group of cells in the DimPlot. This is achieved by using the `cells.highligh` parameter. This is how the default plot looks like.


```r
# Select 1000 random cells out of clusters 1, 5 and 7.
cells.use <- sample(colnames(sample[, sample$seurat_clusters %in% c("1", "5", "7")]), 1000)
Seurat::DimPlot(sample, cells.highlight = cells.use)
```

This is how SCpubr returns this plot. For this, the same parameter has to be set up.


```r
# Select 1000 random cells out of clusters 1, 5 and 7.
cells.use <- sample(colnames(sample[, sample$seurat_clusters %in% c("1", "5", "7")]), 1000)
SCpubr::do_DimPlot(sample, cells.highlight = cells.use)
```

By default, the size of all cells in `SCpubr::do_DimPlot()` is the same. However, the size of the highlighted dots can be modified with the parameter `sizes.highlight` from Seurat.


```r
# Select 1000 random cells out of clusters 1, 5 and 7.
cells.use <- sample(colnames(sample[, sample$seurat_clusters %in% c("1", "5", "7")]), 1000)
SCpubr::do_DimPlot(sample, cells.highlight = cells.use, sizes.highlight = 2)
```

## Splitting by a category

Another useful paramter of `Seurat::DimPlot` is `split.by`, which allows you to split your DimPlot into multiple panels, each one containing a different unique value of the metadata variable you have provided to the argument. One can understand this as using the `group.by` parameter and then splitting the resulting DimPlot into different panels. In this example, we are going to use the different clusters as an example This is how it looks by default:


```r
# Using ncol = 5 to maintain some of the proportions. 
Seurat::DimPlot(sample, split.by = "seurat_clusters", ncol = 5)
```
As can be observed, this plots accomplish the task of separating the cells into each panel, but the approach followed actually makes interpretation difficult. Clusters such as Cluster 9, with fewer cells, tell pretty much nothing. Not knowing how the original UMAP looked like is a major downside of this approach. This is where `SCpubr` focus. Instead of using `Seurat`'s `split.by` parameter, it generates as many plots as unique values in the metadata to split the plot by, but uses `cells.highlight` instead, which keeps the rest of cells greyed out. This is how it looks:


```r
# Using ncol = 5 to maintain some of the proportions.
# Using legend = F to remove unwanted repeated legends.
SCpubr::do_DimPlot(sample, split.by = "seurat_clusters", ncol = 5, legend = F)
```

This way, we can see that clusters such as Cluster 7 are way more disperse than the rest, accounting not only for standalone groups of cells but also blending in other bigger clusters. Actually, the user might want to change the color of the highlighted cells in this split DimPlot. This is achieved by using `colors.split` parameter and providing either a color name recognized by `ggplot2` or (recommended) a HEX code.



```r
# Using ncol = 5 to maintain some of the proportions.
# Using legend = F to remove unwanted repeated legends.
SCpubr::do_DimPlot(sample, split.by = "seurat_clusters", ncol = 5, legend = F, colors.split = "black")
```



Furthermore, one might also want to color each cluster by the original color. This can be achieved by using the argument `colorss.split`, either providing a named vector of each cluster (or metadata variable unique value) as names and color hex codes as values or `TRUE`, thus resorting to the default `SCpubr` categorical coloring. 


```r
# Using ncol = 5 to maintain some of the proportions.
# Using legend = F to remove unwanted repeated legends.
SCpubr::do_DimPlot(sample, split.by = "seurat_clusters", ncol = 5, legend = F, colors.split = TRUE)
```






