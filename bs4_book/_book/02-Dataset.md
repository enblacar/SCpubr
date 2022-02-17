# Dataset {-}

Through this manual we are going to use a [publicly available dataset](https://www.10xgenomics.com/resources/datasets/10k-human-pbmcs-3-ht-v3-1-chromium-x-3-1-high) containing 10K raw cells. The following code is used to generate a Seurat object ready for plotting.


```r
counts_path <- "path_to_count_matrix"

# Path count matrix.
counts <- Seurat::Read10X(counts_path)
# Create Seurat object.
sample <- Seurat::CreateSeuratObject(counts = counts, project = "10K_pbmc")
# Compute percentage of mithochondrial RNA.
sample <- Seurat::PercentageFeatureSet(sample, pattern = "^MT-", col.name = "percent.mt")
# Compute QC.
mask1 <- sample$nCount_RNA >= 1000
mask2 <- sample$nFeature_RNA >= 500
mask3 <- sample$percent.mt <= 20
mask <- mask1 & mask2 & mask3
sample <- sample[, mask]
# Normalize.
sample <- Seurat::SCTransform(sample)

# Dimensional reduction.
sample <- Seurat::RunPCA(sample)
sample <- Seurat::RunUMAP(sample, dims = 1:30)
# Find clusters.
sample <- Seurat::FindNeighbors(sample, dims = 1:30)
sample <- Seurat::FindClusters(sample, resolution = 0.2)
```



