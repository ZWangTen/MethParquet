
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
width="700" height="580" /> After Methlist construction,
`dev_meth_score` function can be used to calculate the MRS as a weighted
sum of an individual’s DNAm values and pre-calculated weights associated
with the phenotype.

EWAS using MethParquet can be conducted on the full set or a subset of
DNAm data, through linear regression (`lm_ewas_outcome`), robust linear
regression (`rlm_ewas_outcome`), generalized or mixed linear regression
(`ewas_meth_exposure`, after fitting a null model using `NullModel`).

`cpg_extract` allows users to extract subset of DNAm data by gene or CpG
name, chromosome number, and row indices without loading all the data.

To read specific examples please refer to Vignette.

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

# To view Vignette for Methparquet
browseVignettes('MethParquet')
```

# Citation

Wang, Z., Cassidy, M., Wallace, D. A., & Sofer, T. (2024). MethParquet:
An R package for Rapid and efficient DNA methylation Association
Analysis adopting Apache Parquet. *Bioinformatics*.
<https://doi.org/10.1093/bioinformatics/btae410>
