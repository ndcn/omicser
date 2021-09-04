
<!-- README.md is generated from README.Rmd. Please edit that file -->

# omicser

<!-- badges: start -->

[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
<!-- badges: end -->

NDCN -Omics Browser (omicser TBC) CZI, NDCN -Omics Browser Data
exploration playground for general -omics sharing. Developed as part of
the NDCN Open Science initiative Data Science Pilot March-Sept 202.

# Introduction

The goal is to make a shiny app which contains a simple omics data
browser. A variety of data formats can be configured to browse. Check
the `Quickstart` guide.

# Installation

This shiny app is build inside an R package and can be installed with
`devtools::install_github("ergonyc/omicser")`. It is crucial to set up a
python environment for reticulate which includes scanpy.

# Usage

## Locally

You can use it locally from RStudio by:

    library(omicser)
    launchApp()

## Shiny server

To run it from a shiny server you have to make sure the package is
installed and then you only have to create a file called `app.R` in the
folder where you want to run it from. The file `app.R` should only
contain:

    omicser::launchApp()

# Example data

I added 2 files (positive and negative mode) as example data. After
installing the package you can find them in the folder `extdata`. In the
repository you can find them in `inst/extdata`.

## Installation

You can install the current version of omicser from giuthub with:

`devtools::install_github("omicser")`

but first run this::

``` r
install.packages("devtools")
install.packages("BiocManager")
BiocManager::install("SingleCellExperiment")

install.packages("reticulate")

require(reticulate)
reticulate::install_miniconda()
reticulate::conda_create("omxr",python_version = 3.9)
reticulate::conda_install(envname = "omxr", packages = "scanpy")
reticulate::conda_install(envname="omxr",
                          channel = "conda-forge",
                          packages = c("leidenalg") )
```
