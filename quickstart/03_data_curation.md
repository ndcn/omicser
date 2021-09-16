---
title: "Data Curation"
author: "andy henrie"
date: "9/14/2021"
output:
  html_document: 
    keep_md: yes
    toc: true
  md_document:
    variant: markdown_github
---


# Data curation and preparation

We set up our environment and installed the `omicser` package, and now we are getting into the most cruical step: DATA CURATION

## Overview

The setup of the browser is 5 steps:

1. Environment Setup: underlying tools/packages from R and pyhon <quickstart/01_environment_setup.Rmd>
2. Install: creates the browser and curation functions as an R package from github <quickstart/02_install.Rmd>
3. Data Curation: the _-omics_ data curated into a _database_  <quickstart/03_data_curation.Rmd>
4. Configuration: connecting the _-omics database_ to the browser app <quickstart/04_configurationn.Rmd>
5. Browse: explore the data!  <quickstart/05_browse.Rmd>


## Data curation and preparation

The most crucial step is _curating_ your data into a database that can be loaded by the browser.  As the _curator_ you will have the responsibility to make some choices about what and how the data can be seen.  
This process results in specifying the location of the data, and creating files which are formatted for the browser.

As an example there is the `curate_domenico_stem_cell.R` script in the /examples subdirectory, which illustrates the steps and some of the choices required to curate your data into a browsable _database_. 

The process as illustrated in the example can be broken into the following steps:

1. provenance & meta data setup - define the meta-data and context for the dataset
3. raw data ingest - translate the outputs of your QC to our database format
  3a. define helper function - any helpers you need to read and process the outputs of your QC
  3b. load the raw data
4. pack into the browser data format - pack into the anndata structure (scanpy/python)
5. post-processing - compute relavent marginal quantities, define additional annotation and grouping variables, etc.
  5a. dimension reduction - compute and cluster if needed
6. differential expression tables - compute and/or formate existing tables
7. write database - write the files to the database location

## database path 

The first step will be to make a folder to contain your data.  We will call this folder and contents the _DATABASE_ .   Each _DATABASE_ folder that one might want to load into the browser should be in the same path.  e.g. the with the repositories `quickstart/` path there is a `test_db` folder which contains several sub-folder _DATABASES_

              
## anndata schema


The anndata scheme requires us to define three pieces of data:


1. DATA: a matrix - e.g. transcriptomics - count matrix of cells X genes. (WARNING: make sure its a proper matrix, e.g.: `as.matrix(data)`)
2. FEATURE METADATA - a table of `omic` annotation - e.g. gene/protein names, families, "highly variable genes", "is marker gene"
3. SAMPLE METADATA: a table of sample annotations - e.g.  cell types, sample #, batch ID sex, experiemntal condition, etc.

More info here https://cran.r-project.org/web/packages/anndata/readme/README.html 

![Anndata scheme](/Users/ahenrie/Projects/NDCN_dev/omicser/inst/anndata_for_r.png)


## loading the data

In addition to these three pieces differential expression tables need to be pre-computed.

Here are some of the key helper functions and he section they fall into.

Here is an example of loading three data files and then pakaging them with a "helper function" (defined in section #2 of the example curation script) into the `data_list` which will be used by the next stage.

### `setup_database()`
The `omicser::setup_database()` function packages the separate tables -- DATA matrix, omic FEATURE METADATA annotations, and SAMPLE METADATA -- into the anndata object.  This function can also take the name of a seurat object file.  



```r
#==== 4. pack into anndata =========================================================================
ad <- omicser::setup_database(database_name=DB_NAME,
                              db_path=DS_ROOT_PATH,
                              data_in=data_list,
                              db_meta=NULL ,
                              re_pack=TRUE)
```


Once we have packed the data into the `anndata` object we can leverage `scanpy` and the `reticulate` python backend to do dimension reduction and clustering. 
Although this stage is not nescessary, it demonstrates the bridge to `cellxgene`.


```r
#==== 5-a. dimension reduction - PCA / umap ========================================================
sc <- import("scanpy")
# scanpy pre-processing - sc$pp
sc$pp$pca(ad)
sc$pp$neighbors(ad)
sc$tl$leiden(ad)
sc$tl$umap(ad)
```

## Differential Expression Tables 
This is probably the trickiest part of the curation  Fortunately we have some helper functions to help us compute them.  Often -- especially for _prote-_ omic databases -- differential expressions are computed as the output of the instrumentataion by a commercial software.  These algorithms leverage bespoke statistics, so it will be best to reformat those tables. 



### DE Table Schema

Most proteomic, metabelomic and lipidomic data will have differential calculations at the output of the instrumentation (which leverages know statastical assumptions of the quantifications) we can also use scanpys tools to compute differential expression.  The diff_exp tables will be needed for volcano plots either way.

The differential expression table has these fields:

   - group - the comparison   {names}V{reference}
   - names - what are we comparing?
   - obs_name  - name of the meta data variable
   - test_type - what statistic are we using
   - reference - the denomenator. or the condition we are comparing expressions values to
   - comp_type - grpVref or grpVrest. rest is all other conditions
   - logfoldchanges - log2(name/reference)
   - scores - statistic score
   - pvals - pvalues from the stats test. e.g. t-test
   - pvals_adj - adjusted pvalue (Q)
   - versus - label which we will choose in the browser
   
   
   
### omicser::compute_de_table()

In `R/pre_process_helpers.R` theres a function which leverages `scanpy` and the `anndata` format we have packed to do differential expression.  We just need to pass a few quantites and it returns a properly formatted differential expression table.

parameters:

  - `ad` - the anndata object
  - `comp_types` - what kind of comparisons? there are two types
    - "allVrest" which takes each of our experimental conditions in turn and compares against the "rest" of the data.
    - "{a}V{b}" or "firstgroupVsecondgroup" which compares the experimental condition "firstgroup" against "secondgroup"
  - `test_types` - statistical tests.  See the examples or `scanpy` documentation for which test types are available.
  - `obs_names` -name of the ad$obs column defining the comparision groups
  - `sc` - the scanpy data object we imported with `reticulate`

Here's an example which computes a differential expression with a `wilcoxon` test of significance for each _disease_ with respect to the rest of the distease states, AND for each _cell_type_ with-respect-to the "rest of" the _cell_type_ s. 


```r
sc <- reticulate::import("scanpy")
test_types <- c('wilcoxon')
comp_types <- c('allVrest')
obs_names <- c('disease','cell_type')
diff_exp <- omicser::compute_de_table(ad,comp_types, test_types, obs_names,sc)
```



## Save the data to the _DATABASE_

Finally we write this the anndata data object to our _database_DATABASE_ folder.   In the examples contained in the `quickstart/` folder we defined `DS_ROOT_PATH <- "test_db"`, and `DB_NAME <- "domenico_stem_cell"` .

```r
#==== 8. write data file to load  =========================================================================
ad$write_h5ad(filename=file.path(DS_ROOT_PATH,DB_NAME,"db_data.h5ad"))
```

In the end each DATABASE folder should now these three files:
1. `db_data.h5ad` - the anndata object          
2. `db_de_table.rds` - differntial expression table
3, `db_meta.yml` - list of database 'meta' information

You might also want to save some _intermediate_ files such as in the examples which also generate: `normalized_data.h5ad`,`core_data.h5ad`,and `norm_data_plus_dr.h5ad`.


## Next -> "configure"
We are almost there.  In the next section we will cover the configuration which will add a `db_config.yml` files to the _DATABASE_ directory, and create an `omicser_options.yml` in the directory where the app will be executed.

<LINK TO 04_configuration.Rmd>
