# OmniNorm

<!-- badges: start -->
[![GitHub issues](https://img.shields.io/github/issues/fgualdr/OmniNorm)](https://github.com/fgualdr/OmniNorm/issues)
[![GitHub pulls](https://img.shields.io/github/issues-pr/fgualdr/OmniNorm)](https://github.com/fgualdr/OmniNorm/pulls)
[![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
[![R-CMD-check](https://github.com/fgualdr/OmniNorm/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/fgualdr/OmniNorm/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

**OmniNorm** is an R package for robust normalization of numerical matrices using **mixtures of skewed distributions**. It is designed for complex and unbalanced datasets in which standard normalization assumptions may not hold, including bulk and single-cell omics and noisy or degraded assays.

## Motivation

Many normalization approaches rely, explicitly or implicitly, on assumptions such as:

- most measured features being unchanged across conditions;
- approximately balanced increases and decreases across the measured feature space.

These assumptions can be violated when biological or technical perturbations affect a large fraction of measured features or produce strongly asymmetric changes.

**OmniNorm** addresses this problem by modeling pairwise log-ratio distributions using **skewed mixture models** and estimating scaling factors from the inferred invariant component of the data.

The approach is intended to provide robust normalization in datasets characterized by asymmetric biological changes, heterogeneous distributions, high noise, or sparsity.

## Applications

OmniNorm can be applied to numerical matrices arising from different omics technologies, including:

- **Bulk RNA-seq** — transcriptional profiling
- **Single-cell omics** — sparse and heterogeneous molecular measurements
- **ChIP-seq** — chromatin-associated signal
- **ATAC-seq** — chromatin accessibility
- **Proteomics** — protein abundance measurements
- **CETSA-MS** — thermal stability proteomics

## Installation

Install the latest development version directly from GitHub:

```r
# install.packages("devtools")
devtools::install_github("fgualdr/OmniNorm")
```

## Citation

If OmniNorm materially contributes to analyses reported in a scientific publication, please cite the software.

Until a designated OmniNorm publication is available, please cite:

> Gualdrini F. *OmniNorm: Normalization of Skewed Numerical Datasets Using Skewed Mixture Distributions*. R software package, version 1.0.0. https://github.com/fgualdr/OmniNorm

Citation metadata are also provided in [`CITATION.cff`](CITATION.cff).

## License

OmniNorm is provided for **non-commercial academic, scientific, educational, and research use**.

Scientific results generated using OmniNorm may be published, subject to the attribution and citation requirements specified in the license.

Commercial use, including incorporation into commercial products or services or use primarily for commercial advantage, requires separate written permission from the copyright holder.

Modified versions must be clearly identified as modified and must not be represented as the original OmniNorm software or methodology.

See [`LICENSE`](LICENSE) for the complete terms.

Copyright © 2025–2026 Francesco Gualdrini.
