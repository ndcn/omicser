
<!-- README.md is generated from README.Rmd. Please edit that file -->

# Quickstart: NDCN -Omics Browser

# TODO:

-   add links between teh headings and to the files referenced.

## Overview

These *quickstart* materials are intended to quickly help a user /
curator get the NDCN Omics browser up and running. Contained in this
directory are a sequence of vignettes which detail each step of:

1.  [Environment Setup](01_environment_setup.md)
2.  [Installation](02_install.md)
3.  [Data Curation](03_data_curation.md)
4.  [Configuration](04_configuration.md)
5.  [Browsing](05_browsing.md)

There is also a subdirectory of [`examples/`](examples/) which have some
R-scripts used to curate and execute example datasets, and *stubs* to
some directories holding the curated *DATABASE* examples and
configuration files (e.g. [`omicser_options.yml`](omicser_options.yml))

Below the entire cycle of vignettes is copied into a single tutorial.

## Quickstart Tutorial

### Environment Setup

This package has been developed with a recent (4.1.0) version of R,
shiny 1.6.0, and RStudio v1.4.1717.

``` bash
R --version
```

In addition to a recently fresh R we will also need a few more things:

1.  python 3.9, which we will leverage via the `reticulate` R package.
2.  R `devtools`, so we can install the `omicser` and additional
    dependencies, such as `reticulate`

#### conda python

The current version of the app looks for a conda environment so we need
to have miniconda or anaconda installed. For example for MacOS we can
install miniconda like this:

``` bash
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-MacOSX-x86_64.sh -O ~/miniconda.sh
bash ~/miniconda.sh -b -p $HOME/miniconda
```

Now we can use conda to install a python 3.9 environment to use with
`reticulate`. I like to use a command line to create the environmenbt,
but if you prefer you can also use R. (I’ll illustrate this below) For
this example we’ve called the environment `omxr`, and installed `scanpy`
and `leidenalg` which will make sure we have everything we need on the
python side.

``` bash
conda create --name omxr python=3.9 scanpy                                                                             
conda install -c conda-forge leidenalg                                                                                  
conda activate omxr
```

#### R packages

For the sake of simplicity and clarity lets just install `devtools`, and
`reticulate`. When we install the `omicser` browser package we can
leverage the dependencies to get everything else we might need.

``` r
install.packages("devtools")
install.packages("reticulate")
```

Now if you want to use R to create your environment we just need to use
the following `reticulate` calls:

``` r
YOUR_CONDA_ENV <- "omxr"
require(reticulate)
reticulate::install_miniconda() 
reticulate::conda_create(YOUR_CONDA_ENV,python_version = 3.9)
reticulate::conda_install(envname = YOUR_CONDA_ENV, packages = "scanpy")
reticulate::conda_install(envname= YOUR_CONDA_ENV,
                          channel = "conda-forge",
                          packages = c("leidenalg") )
```

And we are now reatdy for Installation.

### NDCN Browser Installation

Now that we’ve set up our enviroment, installing the NDCN browser with
`devtools::install_github` is easy.

``` r
devtools::install_github("ergonyc/omicser")
```

Thats it! Now comes the crucial (and fun!) part: curating your data.

### Data Curation and Preparation

The most crucial step is *curating* your data into a database that can
be loaded by the browser. As the *curator* you will have the
responsibility to make some choices about what and how the data can be
seen.  
This process results in specifying the location of the data, and
creating files which are formatted for the browser.

As an example there is the `curate_domenico_stem_cell.R` script in the
/examples subdirectory, which illustrates the steps and some of the
choices required to curate your data into a browsable *database*.

The process as illustrated in the example can be broken into the
following steps:

1.  provenance & meta data setup - define the meta-data and context for
    the dataset
2.  raw data ingest - translate the outputs of your QC to our database
    format 3a. define helper function - any helpers you need to read and
    process the outputs of your QC 3b. load the raw data
3.  pack into the browser data format - pack into the anndata structure
    (scanpy/python)
4.  post-processing - compute relavent marginal quantities, define
    additional annotation and grouping variables, etc. 5a. dimension
    reduction - compute and cluster if needed
5.  differential expression tables - compute and/or formate existing
    tables
6.  write database - write the files to the database location

#### database path

The first step will be to make a folder to contain your data. We will
call this folder and contents the *DATABASE* . Each *DATABASE* folder
that one might want to load into the browser should be in the same path.
e.g. the with the repositories `quickstart/` path there is a `test_db`
folder which contains several sub-folder *DATABASES*

#### `anndata` Schema

The anndata scheme requires us to define three pieces of data:

1.  DATA: a matrix - e.g. transcriptomics - count matrix of cells X
    genes.
2.  FEATURE METADATA - a table of `omic` annotation - e.g. gene/protein
    names, families, “highly variable genes”, “is marker gene”
3.  SAMPLE METADATA: a table of sample annotations - e.g. cell types,
    sample \#, batch ID sex, experiemntal condition, etc.

More info here
<https://cran.r-project.org/web/packages/anndata/readme/README.html>

![Anndata scheme](../inst/anndata_for_r.png)

#### Loading the data

In addition to these three pieces differential expression tables need to
be pre-computed.

Here are some of the key helper functions and he section they fall into.

Here is an example of loading three data files and then pakaging them
with a “helper function” (defined in section \#2 of the example curation
script) into the `data_list` which will be used by the next stage.

#### `setup_database()`

The `omicser::setup_database()` function packages the separate tables –
DATA matrix, omic FEATURE METADATA annotations, and SAMPLE METADATA –
into the anndata object. This function can also take the name of a
seurat object file.

``` r
ad <- omicser::setup_database(database_name=DB_NAME,
                              db_path=DS_ROOT_PATH,
                              data_in=data_list,
                              db_meta=NULL ,
                              re_pack=TRUE)
```

Once we have packed the data into the `anndata` object we can leverage
`scanpy` and the `reticulate` python backend to do dimension reduction
and clustering. Although this stage is not nescessary, it demonstrates
the bridge to `cellxgene`.

``` r
#==== 5-a. dimension reduction - PCA / umap ========================================================
sc <- import("scanpy")
# scanpy pre-processing - sc$pp
sc$pp$pca(ad)
sc$pp$neighbors(ad)
sc$tl$leiden(ad)
sc$tl$umap(ad)
```

#### Differential Expression Tables

This is probably the trickiest part of the curation Fortunately we have
some helper functions to help us compute them. Often – especially for
*prote-* omic databases – differential expressions are computed as the
output of the instrumentataion by a commercial software. These
algorithms leverage bespoke statistics, so it will be best to reformat
those tables.

##### DE Table Schema

Most proteomic, metabelomic and lipidomic data will have differential
calculations at the output of the instrumentation (which leverages know
statastical assumptions of the quantifications) we can also use scanpys
tools to compute differential expression. The diff\_exp tables will be
needed for volcano plots either way.

The differential expression table has these fields:

-   `group` - the comparison {names}V{reference}
-   `names` - what are we comparing?
-   `obs_name` - name of the meta data variable
-   `test_type` - what statistic are we using
-   `reference` - the denomenator. or the condition we are comparing
    expressions values to
-   `comp_type` - grpVref or grpVrest. rest is all other conditions
-   `logfoldchanges` - log2(name/reference)
-   `scores` - statistic score
-   `pvals` - pvalues from the stats test. e.g. t-test
-   `pvals_adj` - adjusted pvalue (Q)
-   `versus` - label which we will choose in the browser

##### omicser::compute\_de\_table()

In [`R/pre_process_helpers.R`](../R/pre_process_helpers.R) theres a
function which leverages `scanpy` and the `anndata` format we have
packed to do differential expression. We just need to pass a few
quantites and it returns a properly formatted differential expression
table.

parameters:

-   `ad` - the anndata object
-   `comp_types` - what kind of comparisons? there are two types
    -   “allVrest” which takes each of our experimental conditions in
        turn and compares against the “rest” of the data.
    -   “{a}V{b}” or “firstgroupVsecondgroup” which compares the
        experimental condition “firstgroup” against “secondgroup”
-   `test_types` - statistical tests. See the examples or `scanpy`
    documentation for which test types are available.
-   `obs_names` -name of the ad$obs column defining the comparision
    groups
-   `sc` - the scanpy data object we imported with `reticulate`

Here’s an example which computes a differential expression with a
`wilcoxon` test of significance for each *disease* with respect to the
rest of the distease states, AND for each *cell\_type* with-respect-to
the “rest of” the *cell\_type* s.

``` r
sc <- reticulate::import("scanpy")
test_types <- c('wilcoxon')
comp_types <- c('allVrest')
obs_names <- c('disease','cell_type')
diff_exp <- omicser::compute_de_table(ad,comp_types, test_types, obs_names,sc)
```

#### Save the data to the *DATABASE*

Finally we write this the anndata data object to our
*database\_DATABASE* folder. In the examples contained in the
`quickstart/` folder we defined `DS_ROOT_PATH <- "test_db"`, and
`DB_NAME <- "domenico_stem_cell"` .

``` r
ad$write_h5ad(filename=file.path(DS_ROOT_PATH,DB_NAME,"db_data.h5ad"))
```

In the end each DATABASE folder should now these three files: 1.
`db_data.h5ad` - the anndata object  
2. `db_de_table.rds` - differntial expression table 3, `db_meta.yml` -
list of database ‘meta’ information

You might also want to save some *intermediate* files such as in the
[examples](examples/curate_pbmc3k.R) which also generate:
`normalized_data.h5ad`,`core_data.h5ad`,and `norm_data_plus_dr.h5ad`.

We are almost there. In the next section we will cover the configuration
which will add a `db_config.yml` files to the *DATABASE* directory, and
create an [`omicser_options.yml`](omicser_options.yml) in the directory
where the app will be executed.

### Configuration

#### WIP

### Browsing

#### WIP

### Sharing

#### WIP
