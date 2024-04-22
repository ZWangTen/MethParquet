
<!-- README.md is generated from README.Rmd. Please edit that file -->

# MethParquet

Manipulate and analyze methylation data using Apache Parquet

MethParquet leverages the columnar format Parquet for efficient DNA
methylation data analysis. Epigenome-wide association analysis,
methylation risk score calculation and testing, as well as CpG subset
extraction are streamlined via the Methlist object containing connection
to methylation Parquet directory along with CpG and subject annotation
data.

<img
src="https://github.com/ZWangTen/MethParquet/blob/main/man/Figure/MethParquet.png"
width="700" height="580" />

## Installation

MethParquet requires fucntions from Bioconductor packages
[limma](https://bioconductor.org/packages/release/bioc/html/limma.html)
and
[GENESIS](https://bioconductor.org/packages/release/bioc/html/GENESIS.html).
Please make sure they are installed before installing MethParquet.

``` r
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("limma")
BiocManager::install("GENESIS")

devtools::install_github("ZWangTen/MethParquet",build_vignettes = TRUE)

# To view Vignette for association analysis
browseVignettes('MethParquet')
```

## Citation
