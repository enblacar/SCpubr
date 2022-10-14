<!-- badges: start -->
[![R-CMD-check](https://github.com/enblacar/SCpubr/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/enblacar/SCpubr/actions/workflows/R-CMD-check.yaml)
[![Code Coverage](https://codecov.io/gh/enblacar/SCpubr/branch/main/graph/badge.svg?token=HK7JB08VFD)](https://app.codecov.io/gh/enblacar/SCpubr/)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)


<!-- badges: end -->
  
# SCpubr
<span class="border-0"><img src="man/figures/logo.png" class="mx-auto d-block" alt="" style="box-shadow: none; width: 100%"/></span>

This package aims to provide a streamlined way of generating publication ready plots for known **S**ingle-**C**ell visualizations in a "**pub**lication **r**eady" format (**SCpubr**). This is, the aim is to automatically generate plots with the highest quality possible, that can be used right away or with minimal modifications for a research article. 


For further information read the [publication](https://www.biorxiv.org/content/10.1101/2022.02.28.482303v1).

For installation and tutorials consult the [reference manual](https://enblacar.github.io/SCpubr-book/).

**SCpubr is still under development. A near future CRAN release of v1.0.0 is planned.**

# Installation

**SCpubr** can be installed via:

```r
# From CRAN:
# Future sumission to CRAN.

# From GitHub.
if(!requireNamespace("devtools", quietly = TRUE)){
  install.packages("devtools") # If not installed.
}

## Latest stable development version.
devtools::install_github("enblacar/SCpubr", ref = "v1.0.0-dev-stable")

## Latest release
devtools::install_github("enblacar/SCpubr", ref = "v1.0.0")
```

By default, dependencies should not be installed. In order to access all functions in the package, the following
packages should also be installed:

```r
# Install CRAN packages.
cran_packages <- c("circlize",
                   "colorspace",
                   "dplyr",
                   "forcats",
                   "ggbeeswarm",
                   "ggdist",
                   "ggExtra",
                   "ggplot2",
                   "ggplotify",
                   "ggrastr",
                   "ggrepel",
                   "ggridges",
                   "ggsignif",
                   "grDevices",
                   "grid",
                   "magrittr",
                   "patchwork",
                   "pbapply",
                   "plyr",
                   "rlang",
                   "scales",
                   "scattermore",
                   "Seurat",
                   "stats",
                   "stringr",
                   "svglite",
                   "tibble",
                   "tidyr",
                   "viridis")

install.packages(cran_packages)

# Install bioconductor packages.
bioconductor_packages <- c("ComplexHeatmap",
                           "infercnv",
                           "Nebulosa")
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(bioconductor_packages)

# Install github packages.
github_packages <- c("ggsankey",
                     "liana",
                     "monocle3")

if (!requireNamespace("remotes", quietly = TRUE))
    install.packages("remotes")

remotes::install_github(github_packages)
```

# Updates
`SCpubr` is an active package currently aiming to improve and add new functionalities.

Keep track of our new updates in the [NEWS page](https://github.com/enblacar/SCpubr/blob/master/NEWS.md)!

# Citation
To cite `SCpubr` in your publications, please use: 

```
Blanco-Carmona, E. Generating publication ready visualizations 
for Single Cell transcriptomics using SCpubr. bioRxiv (2022) 
doi:10.1101/2022.02.28.482303.
```

# Contact
scpubr@gmail.com

