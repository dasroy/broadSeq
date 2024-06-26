
<!-- README.md is generated from README.Rmd. Please edit that file -->

# broadSeq

<!-- badges: start -->
<!-- badges: end -->

The goal of broadSeq is to do easily RNA-seq data analysis with multiple
methods (usually which needs many different input formats). Here the
user will provide the expression data as a
*[SummarizedExperiment](https://bioconductor.org/packages/3.17/SummarizedExperiment)*
object and will get results from different methods. This function
oriented package will give user freedom to develop customized and
reproducible workflow. Additionally, it will also help to quickly
evaluate different methods easily.

![](man/figures/README-unnamed-chunk-2-1.png)

##### This is under active development. Please report if you find any bug.

## Installation

``` r
if (!require("BiocManager"))
    install.packages("BiocManager")
BiocManager::install("broadSeq")
```

You can install the development version of broadSeq from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("dasroy/broadSeq")
```

For package documentation

``` r
devtools::install_github("dasroy/broadSeq", build_vignettes = TRUE)
browseVignettes(package = "broadSeq")
```
