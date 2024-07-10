# SCpubr <a href="https://github.com/enblacar/SCpubr"><img src="man/figures/SCpubr_Logo512.png" align="right" height="200" width="200"/></a>
<!-- badges: start -->
[![CRAN status](https://www.r-pkg.org/badges/version/SCpubr)](https://CRAN.R-project.org/package=SCpubr)
[![R-CMD-check](https://github.com/enblacar/SCpubr/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/enblacar/SCpubr/actions/workflows/R-CMD-check.yaml)
[![Code Coverage](https://codecov.io/gh/enblacar/SCpubr/branch/main/graph/badge.svg?token=HK7JB08VFD)](https://app.codecov.io/gh/enblacar/SCpubr/)
[![CodeFactor](https://www.codefactor.io/repository/github/enblacar/scpubr/badge/main)](https://www.codefactor.io/repository/github/enblacar/scpubr/overview/main)
[![metacran downloads](https://cranlogs.r-pkg.org/badges/grand-total/SCpubr)](https://cran.r-project.org/package=SCpubr)
<!-- badges: end -->


This package aims to provide a streamlined way of generating publication ready plots for known **S**ingle-**C**ell visualizations in a "**pub**lication **r**eady" format (**SCpubr**). This is, the aim is to automatically generate plots with the highest quality possible, that can be used right away or with minimal modifications for a research article. 

## Installation

**SCpubr** can be installed via:

```r
# From CRAN - Official release:
install.packages("SCpubr")

# From GitHub - Latest stable development version:
if(!requireNamespace("devtools", quietly = TRUE)){
  install.packages("devtools") # If not installed.
}

devtools::install_github("enblacar/SCpubr", ref = "v2.0.0-dev-stable")
```

For tutorials, consult the [reference manual](https://enblacar.github.io/SCpubr-book/).


## Updates
`SCpubr` is an active package currently aiming to improve and add new functionalities.

Keep track of our new updates in the [NEWS page](https://github.com/enblacar/SCpubr/blob/master/NEWS.md)!

## Citation
To cite `SCpubr` in your publications, please use: 

```
Blanco-Carmona, E. Generating publication ready visualizations 
for Single Cell transcriptomics using SCpubr. bioRxiv (2022) 
doi:10.1101/2022.02.28.482303.
```

## Contact
`scpubr@gmail.com`
