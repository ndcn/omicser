---
title: "Quickstart"
author: "andy henrie"
date: "9/14/2021"
output:
  html_document: 
    keep_md: yes
    toc: true
  md_document:
    variant: markdown_github
---



## Overview

In order to use the NDCN ‘omics’ browser we need to set a few things up.  This should be four simple steps:

1. Setup the environment: The underlying tools/packages from R and python
2. Create the browser itself: which is an R package available on github
3. Curate a _database_ of ‘omics’ data for browsing
4. Configure the _database_ for the browser, and finally
5. Run the browser!

The following materials should provide a simple set of steps to accomplish this and start browsing ‘omics’!
This document assumes software installation and environment configuration documented in the README.  
    
## Data curation and preparation

The most crucial step is _curating_ your data into a database that can be loaded by the browser.  As the _curator_ you will have the responsibility to make some choices about what and how the data can be seen.  
This process results in specifying the location of the data, and creating files which are formatted for the browser.

As an example there is the `curate_domenico_stem_cell.R` script in the /examples subdirectory, which illustrates the steps and some of the choices required to curate your data into a browsable _database_. 

The process as illustrated in the example can be broken into the following steps:

1. provenance & meta data setup - define the meta-data and context for the dataset
2. helper function definition - any helpers you need
3. raw data ingest - load the raw data from whatever format they live in
4. pack into the browser data format - pack into the anndata structure (scanpy/python)
5. post-processing - compute relavent marginal quantities, define additional annotation and grouping variables, etc.
  5a. dimension reduction - compute and cluster if needed
6. differential expression tables - compute and/or formate existing tables
7. write database - write the files to the database location


The first step will be to make the .yml file that will let the browser know what/where the data will be.  Your options are to edit the "omxr_options.yml" directly or make a new one. e.g. for the stem cell proteomics example:


```r
dataset_names <- list(
  "Domenico DIA" = "domenico_stem_cell"
)
# where do our data files live. # WARNING do not use the ~ alias for Home
# MUST BE FULL or RELATIVE PATH... ~ will cause loading to
ds_root_path = 'test_db'

# python environment
conda_environment = 'omxr'

omicser_options <- list(dataset_names=dataset_names,
                        ds_root_path=ds_root_path,
                        conda_environment=conda_environment)


require(configr)
configr::write.config(config.dat = omicser_options, file.path = "omxr_options.yml",
                      write.type = "yaml", indent = 4)
```

<<<<<<< HEAD
The anndata scheme requires us to define three main tables:
- the omic data matrix - e.g. transcriptomics - count matrix of cells X genes.
- omic annotation data - e.g. gene families, "highly variable genes", "is marker gene"
- sample meta data - e.g.  sex, experiemntal condition, etc.
=======
The anndata scheme requires us to define three pieces of data:
>>>>>>> 352d709 (clean up ingest / config logic, step 1)

1. DATA: a matrix - e.g. transcriptomics - count matrix of cells X genes.
2. FEATURE METADATA - a table of `omic` annotation - e.g. gene/protein names, families, "highly variable genes", "is marker gene"
3. SAMPLE METADATA: a table of sample annotations - e.g.  cell types, sample #, batch ID sex, experiemntal condition, etc.

More info here https://cran.r-project.org/web/packages/anndata/readme/README.html 

![Anndata scheme](/Users/ahenrie/Projects/NDCN_dev/omicser/inst/anndata_for_r.png)

In addition to these three pieces differential expression tables need to be pre-computed.

Here are some of the key helper functions and he section they fall into.

Here is an example of loading three data files and then pakaging them with a "helper function" (defined in section #2 of the example curation script) into the `data_list` which will be used by the next stage.


```r
#==== 3. load data -========================================================================================
matrix_data_file <- "20210524_093609_170805_aging_against_SC_merged_all_lib_2_Report.xls"
# candidate table without filter
annot_de_file <- "170805_aging_against_SC_merged_all_lib_2_candidates.xls"
# condition setup
conditions_table_file <- "170805_aging_against_SC_merged_all_lib_2_ConditionSetup.xls"
data_list <- prep_DIA_files(matrix_data_file,annot_de_file,conditions_table_file,RAW_DIR)
# save diff expression data for later...
diff_exp <- data_list$de
# saveRDS(diff_exp, file = file.path(DB_DIR, "diff_expr_table.rds"))
saveRDS(diff_exp, file.path(DS_ROOT_PATH,DB_NAME, "diff_expr_table.rds"))
```

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


Most proteomic, metabelomic and lipidomic data will have differential calculations at the output of the instrumentation (which leverages know statastical assumptions of the quantifications) we can also use scanpys tools to compute differentalial expression.  The diff_exp tables will be needed for volcano plots either way.


```r
#==== 6. differential expression =====================================================================
test_types <- c('wilcoxon','t-test_overestim_var')
comp_types <- c("grpVrest")
obs_names <- c('disease','cell_type')
diff_exp <- omicser::compute_de_table(ad,comp_types, test_types, obs_names)
```

Finally we write this the anndata file to our database folder as specified in the  omxr_options.yml.


```r
#==== 8. write data file to load  =========================================================================
ad$write_h5ad(filename=file.path(DS_ROOT_PATH,DB_NAME,"omxr_data.h5ad"))
```

Please refer to the full curation script for more context. 

TODO:  more information on the diff_exp data schema expected by the browser.


## Configuration
This involves executing a few “ingestor” helper functions, and a few choices by the ‘curator’ to tell the browser where to find the data.  

Edit the .yml or better yet include it in the curation script.  e.g. part 7

Finally we need to define the configuration.  Most of these fields *could* be inferred from the anndata file, but this is where curation is important.  Lets choose the most reasonable quantities _only_.  


```r
#==== 7. create configs =========================================================================
# differentials  #if we care we need to explicitly state. defaults will be the order...
conf_list <- list(
  x_obs = c("Is.Reference","Condition","Replicate", "Label"),
  y_obs =  c("expr_var", "expr_mean", "expr_frac", "sample_ID", "leiden"), #MEASURES
  obs_groupby = c("Is.Reference","Condition","Replicate", "Label"),
  obs_subset = c("Is.Reference","Condition","Replicate", "Label"),

  x_var = character(0),
  y_var = c("expr_geomean", "expr_mean", "expr_var", "expr_frac" ),

  var_groupby = character(0),
  var_subset = character(0),

  diffs = list(diff_exp_comps = levels(factor(diff_exp$versus)),
               diff_exp_comp_type =  levels(factor(diff_exp$comp_type)),
               diff_exp_obs_name =  levels(factor(diff_exp$obs_name)),
               diff_exp_tests =  levels(factor(diff_exp$test_type))
  ),

  layers = c("X","raw","X_is_scaled_na_to_0","scaled","zro_na"),

  # Dimred
  dimreds = list(obsm = ad$obsm_keys(),
                 varm = ad$varm_keys()),

  # what ad$obs do we want to make default values for...
  # # should just pack according to UI?
  default_factors = c("Condition","Color","Replicate")

)
configr::write.config(config.dat = conf_list, file.path = file.path(DS_ROOT_PATH,DB_NAME,"config.yml" ),
                      write.type = "yaml", indent = 4)
```

## Browse!!

Assuming you have already loaded the omicser package, once the .yml files have been generated and teh data placed in the right directories you are good to browse!


```r
run_app(options = list(launch.browser = TRUE))
```

Note, that if you are using a development setup, you might want to run the run_dev script which will handle unloading / re-loading for you.
