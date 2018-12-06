
<!-- README.md is generated from README.Rmd. Please edit that file -->
RSoptSC
=======

<img src="man/figures/logo.svg" width="60%" />

RSoptSC is an R package for SoptSC: performing analysis and prediction on single-cell RNA sequencing (scRNA-seq) data. SoptSC can infer cell-cell communication between single cells, and completes multiple single-cell analysis tasks from a single coherent framework. This enables the reconstruction of complex lineage relationships alongside prediction of specific feedback or feedforward interactions.

-   Demo SoptSC at `vignette("Joost_et_al")`

-   Read the SoptSC preprint on [bioRxiv](https://www.biorxiv.org/content/early/2018/05/12/168922)

Installation
------------

``` r
install.packages("devtools")
library(devtools)
install_github("mkarikom/RSoptSC")
```

Features
--------

-   Regularized network inference on single-cell gene expression.

-   Unsupervised clustering

-   Semi-supervised or unsupervised pseudotime by root cell/cluster inference and graph embedding
