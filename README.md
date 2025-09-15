
# OmniNorm

<!-- badges: start -->

[![GitHub
issues](https://img.shields.io/github/issues/fgualdr/OmniNorm)](https://github.com/fgualdr/OmniNorm/issues)
[![GitHub
pulls](https://img.shields.io/github/issues-pr/fgualdr/OmniNorm)](https://github.com/fgualdr/OmniNorm/pulls)
[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
[![R-CMD-check-bioc](https://github.com/fgualdr/OmniNorm/actions/workflows/R-CMD-check-bioc.yaml/badge.svg)](https://github.com/fgualdr/OmniNorm/actions/workflows/R-CMD-check-bioc.yaml)
[![R-CMD-check](https://github.com/fgualdr/OmniNorm/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/fgualdr/OmniNorm/actions/workflows/R-CMD-check.yaml)
<!-- badges: end --> 

OmniNorm is an R library that utilizes a mixed Skewed distribution to normalize allegedly any numerical (matrix-arranged) dataset. It has been developed to address the normalization of heavily unbalanced distributions where classical approaches such as median, quantile, etc., normalization fail.

The library has been tested on multi-omics datasets including RNA-seq, ChIP-seq, ATAC-seq, and proteomic-based datasets, as well as heavily degraded assays such as CETSA-MS.

The fundamental idea behind OmniNorm is that it is not always accurate to assume that given a perturbation, the number of changing observations is the minority, and that the number of up or down changes is equally distributed. This observation implies that instead of symmetrically distributed changes, the distributions of changes will display various degrees of skewness, with the unperturbed population being the only one with a quasi-normal distribution.

A more detailed description will soon be added to the repository, as a short article is under preparation.

## Installation instructions

Get the latest stable `R` release from
[CRAN](http://cran.r-project.org/). Then install `OmniNorm`
using devtool:

``` r
require(devtools)
install_github("fgualdr/OmniNorm")
```
