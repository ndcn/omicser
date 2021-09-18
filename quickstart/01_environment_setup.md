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





#  Environment Setup

This package has been developed with a recent (4.1.0) version of R, shiny 1.6.0, and RStudio v1.4.1717.  

```bash
R --version
```


In addition to a recently fresh R we will also need a few more things:

1. python 3.9, which we will leverage via the `reticulate` R package. 
2. R `devtools`, so we can install the `omicser` and additional dependencies, such as `reticulate`


### Conda Python 

The browser app looks for a conda environment so first we need to have miniconda or anaconda installed. (Please file a bug or a pull request to begin the process of supporting other python environments.)  If you don't already have anaconda or mini-conda, its easy to  [install](https://conda.io/projects/conda/en/latest/user-guide/install/index.html).

Or, if you prefer to do everything in R, it can be installed via `reticulate`. 

The browser uses python to manage the data on the back end an we require v3.9. 

### R packages
For the sake of simplicity and clarity we can just install `devtools`, and `reticulate`.  When we install the `omicser` browser package we can leverage the dependencies to get everything else we might need.



```r

install.packages("devtools")
install.packages("reticulate")
```

Now, to use R to create your environment we just need to use the following `reticulate` calls:


```r
YOUR_CONDA_ENV <- "omxr"
PYTHON_VERSION <- 3.9

require(reticulate)
reticulate::install_miniconda() 
reticulate::conda_create(YOUR_CONDA_ENV, python_version = PYTHON_VERSION)
```

Now that we have conda and our environment available we need to install two packages to do everything contained in the example scripts.

(Warning: if you have bioconductor configured as a channel it may cause problems.  Forcing install from the "conda-forge" channel is preferred)


```r
YOUR_CONDA_ENV <- "omxr"
reticulate::conda_install(envname = YOUR_CONDA_ENV, packages = "scanpy")
reticulate::conda_install(envname= YOUR_CONDA_ENV,
                          channel = "conda-forge",
                          packages = c("leidenalg") )
```

And we are now ready for [Installation](02_install.md). 


### ETC 
#### Alternate Command Line Instructions

For example for MacOS we can install miniconda like this:


```bash
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-MacOSX-x86_64.sh -O ~/miniconda.sh
bash ~/miniconda.sh -b -p $HOME/miniconda
```

Now we can use conda to install a python 3.9 environment to use with `reticulate`. 
For this example we've called the environment `omxr`, and installed `scanpy` and `leidenalg` which will make sure we have everything we need on the python side.

```bash
conda create --name omxr python=3.9 scanpy                                                                             
conda install -c conda-forge leidenalg                                                                                  
conda activate omxr
```

