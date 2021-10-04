
<!-- README.md is generated from README.Rmd. Please edit that file -->

# NDCN -Omics Browser (name TBC)

<!-- badges: start -->

[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
<!-- badges: end -->

## Introduction

The goal of this package is to allow exploration of -omics data using an
interactive “playground” (e.g., Shiny application), accessible through
an interactive web browser.

> *Note*: although the NDCN browser is not yet officially named, there
> are references to the browser/package as `omicser` or `omxr`.

There are two types of users for this package:

-   **Curators**: data analysts who configure the app and prepare
    datasets for viewing. These individuals should be fairly proficient
    with R coding, including familiarity with building Shiny apps and
    knowledgeable about -omics data structures.
-   **Viewers**: scientists who are interested in exploring -omics data
    using the application. They do not require extensive knowledge about
    R, but may need to have it installed if they are running the app on
    their own computer.

The diagram below summarizes the steps required to use this package:

![Overview of data visualization using the
omicser](man/figures/README-omicser-overview.png)

Support for the following *post-QC* -omics data types has been tested:

-   *transcript*-omics
-   *prote*-omics
-   *metabol*-omics
-   *lipid*-omics

The package is able to provide the following outputs:

-   images visualizing data
-   subsets of data
-   web application packaged with data that can be shared with
    collaborators, so they can run the browser on their own computer

## Installation and usage for Viewers

Once a curator has prepared the data and configured the application, the
package can be installed with R:

``` r
devtools::install_github("ergonyc/omicser")
library(omicser)
```

And then then the browser can be launched:

``` r
omicser::run_app()
```

For more information on launching the web application with your data,
please consult with the Curator who created it.

## Installation and usage for Curators

This package was developed to create a
[`{shiny}`](https://shiny.rstudio.com/) application using  
[`{golem}`](https://github.com/ThinkR-open/golem), “an opinionated
framework for building production-grade shiny applications.” It requires
a number of additional R and Python packages. The package creates the
Shiny browser app as well as several functions for curating and
packaging the data being curated.

Documentation to walk through the process of installing necessary
software, creating your application, and browsing your data is available
in the [Quickstart
guide](https://ergonyc.github.io/omicser/articles/00_quickstart.html):

-   [**1-Set up
    environment**](https://ergonyc.github.io/omicser/articles/01_environment_setup.html):
    ensure underlying tools/packages from R and python are available
-   [**2-Install
    package**](https://ergonyc.github.io/omicser/articles/02_install.html):
    includes data curation functions and software that runs the browser
    application
-   [**3-Curate
    data**](https://ergonyc.github.io/omicser/articles/03_data_curation.html):
    the *-omics* data curated into a *database*
-   [**4-Configure
    application**](https://ergonyc.github.io/omicser/articles/04_configuration.html):
    connecting the *-omics database* to the browser app
-   [**5-Browse
    data**](https://ergonyc.github.io/omicser/articles/05_browsing.html):
    explore the data, test hypotheses, create visualizations
-   [**6-Share
    application**](https://ergonyc.github.io/omicser/articles/06_sharing.html):
    perform final curation steps that allow you to share the browser
    with collaborators (Viewers)

Please view the [vignettes directory of the GitHub
repository](https://github.com/ergonyc/omicser/tree/main/vignettes) for
this project to view additional files used in the Quickstart.

## Contributing to this project

This project began as part of the CZI NDCN Open Science Initiative Data
Science Pilot (March-Sept 2021).

Please see our [Contributing guidelines](CONTRIBUTING.md) for more
information on any of the following:

-   reporting bugs and problems using the software
-   requesting additional features and functionality
-   submitting your own code or documentation to become a part of this
    project
-   sharing information about how you’ve used this package in your own
    work

Anyone interacting with our project is expected to follow our [Code of
Conduct](CODE_OF_CONDUCT.md).

### Roadmap

Please view our [issues](https://github.com/ergonyc/omicser/issues),
especially those tagged as [high
priority](https://github.com/ergonyc/omicser/issues?q=is%3Aopen+is%3Aissue+label%3A%22high+priority%22),
for more information about immediate plans for development.
