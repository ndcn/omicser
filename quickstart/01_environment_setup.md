---
title: "Environment Setup"
author: "andy henrie"
date: "9/14/2021"
output:
  html_document: 
    toc: yes
    keep_md: yes
  md_document:
    variant: markdown_github
---




# Environment Setup

This package has been developed with a recent (4.1.0) version of R, shiny 1.6.0, and RStudio v1.4.1717.  

```bash
R --version
```

## Overview

In addition to a recently fresh R we will also need a few more things:

1. python 3.9, which we will leverage via the `reticulate` R package. 
2. R `devtools`, so we can install the `omicser` and additional dependencies, such as `reticulate`


## conda python 

The current version of the app looks for a conda environment so we need to have miniconda or anaconda installed.   For example for MacOS we can install miniconda like this:


```bash
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-MacOSX-x86_64.sh -O ~/miniconda.sh
bash ~/miniconda.sh -b -p $HOME/miniconda
```

Now we can use conda to install a python 3.9 environment to use with `reticulate`.   I like to use a command line to create the environmenbt, but if you prefer you can also use R. (I'll illustrate this below) 
For this example we've called the environment `omxr`, and installed `scanpy` and `leidenalg` which will make sure we have everything we need on the python side.

```bash
conda create --name omxr python=3.9 scanpy                                                                             
conda install -c conda-forge leidenalg                                                                                  
conda activate omxr
```


## R packages
For the sake of simplicity and clarity lets just install `devtools`, and `reticulate`.  When we install the `omicser` browser package we can leverage the dependencies to get everything else we might need.



```r
install.packages("devtools")
install.packages("reticulate")
```

Now if you want to use R to create your environment we just need to use the following `reticulate` calls:


```r
YOUR_CONDA_ENV <- "omxr"
require(reticulate)
reticulate::install_miniconda() 
reticulate::conda_create(YOUR_CONDA_ENV,python_version = 3.9)
reticulate::conda_install(envname = YOUR_CONDA_ENV, packages = "scanpy")
reticulate::conda_install(envname= YOUR_CONDA_ENV,
                          channel = "conda-forge",
                          packages = c("leidenalg") )
```


## Next -> Installation

link <check 02_install.Rmd>
