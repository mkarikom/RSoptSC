
<!-- README.md is generated from README.Rmd. Please edit that file -->
RSoptSC
=======

<img src="man/figures/logo.svg" width="60%" />

RSoptSC is an R package for SoptSC: performing analysis and prediction on single-cell RNA sequencing (scRNA-seq) data. SoptSC can infer cell-cell communication between single cells, and completes multiple single-cell analysis tasks from a single coherent framework. This enables the reconstruction of complex lineage relationships alongside prediction of specific feedback or feedforward interactions.

SoptSC is also available as a MATLAB package [here](https://github.com/wangshuxiong/soptsc).

-   Demo SoptSC at [`vignette("Joost_et_al")`](https://mkarikom.github.io/RSoptSC/articles/Joost_et_al.html)

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

-   Inference of cell-cell communication between single cells

-   Integration of multiple analyses in a consistent mathematical framework: clustering, marker genes, pseudotime, and lineage inference

-   Cell-cell similarity matrix construction to improve clustering

-   NMF-based marker gene identification

-   Prediction of the number of clusters present (via eigengap properties of the similarity matrix)

-   Prediction of the initial cell in pseudotime

Documentation
-------------

Full details of RSoptSC and examples are available [here](https://mkarikom.github.io/RSoptSC).
