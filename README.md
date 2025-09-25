# OmniNorm

**OmniNorm** is an R package for robust normalization of numerical matrices using **mixtures of skewed distributions**. It is designed to handle complex, unbalanced datasets where standard normalization methods (e.g., median, quantile) fail — including **single-cell omics**, **bulk multi-omics**, and **noisy or degraded assays**.

---

## 🔍 Motivation

In many biological experiments, it is often assumed that:

- Most features are unchanged across conditions,
- Up- and down-regulated features are balanced.

However, **real-world data rarely follows these assumptions**. Perturbations may affect a large proportion of features, or introduce **directionally biased (skewed)** changes. In such cases, classical normalization can introduce significant artifacts.

**OmniNorm** addresses this by modeling pairwise log-ratio distributions using **skewed mixture models**, providing robust scaling across diverse datasets — including high-noise, high-sparsity data like single-cell experiments.

---

## 🧪 Applications

OmniNorm has been tested and validated on a wide range of omics datasets, including:

- 🧬 **Bulk RNA-seq** – transcriptional profiling
- 🧫 **Single-cell techniques** – sparse and noisy
- 🧬 **ChIP-seq** – chromatin-bound protein occupancy
- 🧬 **ATAC-seq** – chromatin accessibility profiling
- 🧊 **Proteomics** – protein abundance quantification
- 🧊 **CETSA-MS** – thermal shift-based proteomics under degradation

---

## ⚙️ Installation

You can install the latest development version from GitHub using:

```r
# install.packages("devtools")
devtools::install_github("fgualdr/OmniNorm")
