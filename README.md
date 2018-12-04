
<!-- README.md is generated from README.Rmd. Please edit that file -->
RSoptSC
=======

<img src="man/figures/logo.svg" align="center" alt="" width="720" />

We present SoptSC, a similarity matrix optimization method for single-cell data analysis, which performs clustering, pseudotemporal ordering, lineage and marker gene identification from a cell similarity matrix. SoptSC also presents a new function: cell-cell signaling network inference, enabling the reconstruction of complex lineage relationships and associated feedback/feedforward interactions. As we show by application to several datasets, SoptSC can predict the number of clusters and the initial state unsupervised, and outperforms current methods for inference of clusters and pseudotime from single-cell data. Learn more at `vignette("Joost_et_al")`.

Installation
------------

``` r
# Download the zip from GitHub and unzip
# 
devtools::install("path to unzipped file")
```

Features
--------

-   Regularized network inference on single-cell gene expression.

-   Unsupervised clustering

-   Semi-supervised or unsupervised pseudotime by root cell/cluster inference and graph embedding
