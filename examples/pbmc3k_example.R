# Overview --------------
#### Create an app to browse PBMC3k data from 10X Genomics
require(reticulate) # runs python when curating data

#  Step 1: Set paths--------------
#  path where the omicser repo is cloned to (see install_scrip.R)
REPO_PATH <- getwd() # /path/to/cloned/omicser
OMICSER_RUN_DIR <- REPO_PATH

DB_NAME <- list("10x PBMC3k" = "pbmc3k") # name of database(s)
RAW_DATA_DIR <- file.path(OMICSER_RUN_DIR,"quickstart/raw_data",DB_NAME)

DB_ROOT_PATH <- file.path(OMICSER_RUN_DIR,"quickstart/test_db")

OMICSER_PYTHON <-  "pyenv_omxr"
# installation type (see install_script.R)
CLONED_OMICSER <- TRUE

# Step 2:  Load the omicser package
if (CLONED_OMICSER <- TRUE){
  golem::document_and_reload(pkg = REPO_PATH)
} else {
  require(omicser)
  #see install_script.R if not installed
}

# Step 3: Assert python back-end ----------------
#  for the curation we need to have scanpy
if (FALSE){  #you should already have installed miniconda and created the env
    reticulate::install_miniconda() #in case it is not already installed

    # full conda install
    # packages1 <- c("seaborn", "scikit-learn", "statsmodels", "numba", "pytables")
    # packages2 <- c("python-igraph", "leidenalg")
    #
    # reticulate::conda_create(OMICSER_PYTHON, python_version = 3.8,packages = packages1)
    # reticulate::conda_install(envname=OMICSER_PYTHON,
    #                         channel = "conda-forge",
    #                         packages = packages2 )
    # reticulate::conda_install(envname=OMICSER_PYTHON,
    #                           channel = "conda-forge",
    #                           pip = TRUE,
    #                           packages = "scanpy" )

    # simpler pip pypi install
    packages <- c("scanpy[leiden]")
    reticulate::conda_create(OMICSER_PYTHON, python_version = 3.8)
    reticulate::conda_install(envname=OMICSER_PYTHON,
                              channel = "conda-forge",
                              pip = TRUE,
                              packages =  packages )

}
reticulate::use_condaenv(condaenv = OMICSER_PYTHON,
                                conda = reticulate::conda_binary(),
                                required = TRUE)
# check that we have our python on deck
reticulate::py_discover_config()
if (!reticulate::py_module_available(module = "scanpy") ) { #    reticulate::conda_create(OMICSER_PYTHON, python_version = 3.8,packages = packages1)

  reticulate::conda_install(envname=OMICSER_PYTHON,
                                          packages = "scanpy[leiden]",
                                          pip = TRUE,
                                          conda=reticulate::conda_binary())
}


# Step 4:  get the data ---------------
# create directory structure for data and databases
DB_DIR = file.path(DB_ROOT_PATH,DB_NAME)
if (!dir.exists(DB_DIR)) {
  dir.create(DB_DIR)
}

if (!dir.exists(RAW_DATA_DIR)) {
  dir.create(RAW_DATA_DIR)
}
# change paths to make data manipulations easier
setwd(RAW_DATA_DIR)
# download data
data_file <- "http://cf.10xgenomics.com/samples/cell-exp/1.1.0/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz"
download.file(data_file, "pbmc3k_filtered_gene_bc_matrices.tar.gz")
# extract all downloaded files
untar("pbmc3k_filtered_gene_bc_matrices.tar.gz")
# compress individual files (required for scanpy)
tar("filtered_gene_bc_matrices/hg19/matrix.mtx.gz", files = "filtered_gene_bc_matrices/hg19/matrix.mtx", compression = "gzip")
tar("filtered_gene_bc_matrices/hg19/barcodes.tsv.gz", files = "filtered_gene_bc_matrices/hg19/barcodes.tsv", compression = "gzip")
tar("filtered_gene_bc_matrices/hg19/genes.tsv.gz", files = "filtered_gene_bc_matrices/hg19/genes.tsv", compression = "gzip")

# Step 5:  define for source helper functions ------------

# N/A



# Step 6:  pack data into AnnData format --------------
#### Data curation 2. Format and ingest raw data
# make scanpy functions available
sc <- reticulate::import("scanpy")

# load the PBMC dataset
adata <- sc$read_10x_mtx(
  # locate directory containing mtx file
  "filtered_gene_bc_matrices/hg19/",
  # use gene symbols for the variable names (variables-axis index)
  var_names='gene_symbols',
  # write a cache file for faster subsequent reading
  cache=TRUE)
# save as file
#
#

# reset working directory
setwd(OMICSER_RUN_DIR)

adata$write_h5ad(file.path(DB_ROOT_PATH, DB_NAME,"core_data.h5ad"))

# identify location of raw data
data_list <- list(object=file.path(DB_ROOT_PATH,DB_NAME,"core_data.h5ad"))

# create database formatted as AnnData
adata <- omicser::setup_database(database_name = DB_NAME,
                                 db_path = DB_ROOT_PATH,
                                 data_in = data_list,
                                 db_meta = NULL,
                                 re_pack = TRUE)




# Steps 7-9: CURATION
# Step 7: additional data processing ----
adata$var_names_make_unique()
# unnecessary if using `var_names='gene_ids'` in `sc.read_10x_mtx`

# filter data
sc$pp$filter_cells(adata, min_genes=200)
sc$pp$filter_genes(adata, min_cells=10)

# annotate the group of mitochondrial genes as 'mt'
adata$var['mt'] <- startsWith(adata$var_names,'MT-')
sc$pp$calculate_qc_metrics(adata, qc_vars=list('mt'), percent_top=NULL, log1p=FALSE, inplace=TRUE)

# filter data
adata <- adata[adata$obs$n_genes_by_counts < 2500, ]
adata <- adata[adata$obs$pct_counts_mt < 5, ]
sc$pp$normalize_total(adata, target_sum=1e4)
sc$pp$log1p(adata)
sc$pp$highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)

# transform data
#adata = adata[, adata$var$highly_variable]
sc$pp$regress_out(adata, list('total_counts', 'pct_counts_mt'))
sc$pp$scale(adata, max_value=10)

# choose top 40 genes by variance across dataset as "targets"
adata$var$var_rank <- order(adata$var$dispersions_norm)
target_omics <- adata$var_names[which(adata$var$var_rank <= 40)]

# calculate deciles
adata$var$decile <- dplyr::ntile(adata$var$dispersions_norm, 10)
#raw <- ad$raw$to_adata()

# save intermediate database file
if (FALSE){
  adata$write_h5ad(filename=file.path(DB_ROOT_PATH,DB_NAME,"normalized_data.h5ad"))
}

#7-b. dimension reduction - PCA / umap
 #pca
sc$pp$pca(adata)
# compute neighbor graph
sc$pp$neighbors(adata)
## infer clusters
sc$tl$leiden(adata)
# compute umap
sc$tl$umap(adata)

# save intermediate database file
if (FALSE){
  adata$write_h5ad(filename=file.path(DB_ROOT_PATH,DB_NAME,"norm_data_plus_dr.h5ad"))
}

# Step 8: pre-compute differential expression
# identify stats
# see scanpy documentation for possible stat test choices
test_types <- c('wilcoxon')
comp_types <- c("grpVrest")
obs_names <- c('leiden')
# calculate DE
diff_exp <- omicser::compute_de_table(adata,comp_types, test_types, obs_names,sc)

### WARNING:  there's an overflow bug in the logfoldchange values for this dataset
### Might need to rescale?

# save intermediate database file
if (FALSE){
  adata$write_h5ad(filename=file.path(DB_ROOT_PATH,DB_NAME,"norm_data_with_de.h5ad"))
}

# save DE tables
saveRDS(diff_exp, file = file.path(DB_ROOT_PATH, DB_NAME, "db_de_table.rds"))


# Step 9: Write data files to database directory -------
# write final database
adata$write_h5ad(filename = file.path(DB_ROOT_PATH, DB_NAME, "db_data.h5ad"))

# Step 10:  configure browser ----
# load intermediate files if available
if (FALSE) {
  adata <- anndata::read_h5ad(filename=file.path(DB_ROOT_PATH,DB_NAME,"db_data.h5ad"))
  diff_exp <- readRDS( file = file.path(DB_ROOT_PATH,DB_NAME, "db_de_table.rds"))
}

# Step 9: Write data files to database directory -------
DB_DIR = file.path(DB_ROOT_PATH,DB_NAME)
if (!dir.exists(DB_DIR)) {
  dir.create(DB_DIR)
}

# save diff expression data
diff_exp <- data_list$de

saveRDS(diff_exp, file.path(DB_ROOT_PATH,DB_NAME, "db_de_table.rds"))

# write the anndata object
adata$write_h5ad(filename=file.path(DB_ROOT_PATH,DB_NAME,"db_data.h5ad"))

#reload
if (FALSE) adata <- anndata::read_h5ad(filename=file.path(DB_ROOT_PATH,DB_NAME,"db_data.h5ad"))



# Step 10:  configure browser ----

omic_type <- "transcript" #c("transcript","prote","metabol","lipid","other")
aggregate_by_default <- (if (omic_type=="transcript") TRUE else FALSE ) #e.g.  single cell
# choose top 40 proteins by variance across dataset as our "targets"
target_features <- adata$var_names[which(adata$var$var_rank <= 40)]
#if we care we need to explicitly state. defaults will be the order...
config_list <- list(
  # meta-tablel grouping "factors"
  group_obs = c("leiden"),
  group_var = c("decile","highly_variable"),

  # LAYERS
  # each layer needs a label/explanation
  layer_values = c("X","raw"),
  layer_names = c("norm-count","counts" ),

  # ANNOTATIONS / TARGETS
  # what adata$obs do we want to make default values for...
  # # should just pack according to UI?
  default_obs =  c("Condition","leiden"), #subset & ordering

  obs_annots = c( "leiden",
                  "n_genes","n_genes_by_counts","total_counts","total_counts_mt","pct_counts_mt"),

  default_var = character(0), #just use them in order as defined
  var_annots = c(
    "n_cells",
    "mt",
    "n_cells_by_counts",
    "mean_counts",
    "pct_dropout_by_counts",
    "total_counts",
    "highly_variable",
    "dispersions_norm",
    "decile"),


  target_features = target_features,
  feature_deets = c( "feature_name",
                     "gene_ids",
                     "n_cells",
                     "mt",
                     "n_cells_by_counts",
                     "mean_counts",
                     "pct_dropout_by_counts",
                     "total_counts",
                     "highly_variable",
                     "means",
                     "dispersions",
                     "dispersions_norm",
                     "mean",
                     "std",
                     "var_rank",
                     "decile" ),

  filter_feature = c("dispersions_norm"), #if null defaults to "fano_factor"

  xr_groupby = c("decile","highly_variable"),
  var_subset = c("decile","highly_variable"),

  # differential expression
  diffs = list( diff_exp_comps = levels(factor(diff_exp$versus)),
                diff_exp_obs_name =  levels(factor(diff_exp$obs_name)),
                diff_exp_tests =  levels(factor(diff_exp$test_type))
  ),

  # Dimension reduction (depricated)
  dimreds = list(obsm = adata$obsm_keys(),
                 varm = adata$varm_keys()),


  #meta info
  annotation_database =  NA,
  publication = "TBD",
  method = "bulk", # c("single-cell","bulk","other")
  omic_type = omic_type, #c("transcript","prote","metabol","lipid","other")
  aggregate_by_default = aggregate_by_default, #e.g.  single cell

  organism = "human",
  lab = "",
  source = "peripheral blood mononuclear cells (PBMCs)",
  annotation_database =  "",
  title = "pbmc3k",
  omic_type = omic_type,
  measurment = "normalized counts",
  pub = "10X Genomics",
  url = "https://support.10xgenomics.com/single-cell-gene-expression/datasets/1.1.0/pbmc3k",
  date = format(Sys.time(), "%a %b %d %X %Y")
)

omicser::write_db_conf(config_list,DB_NAME, db_root = DB_ROOT_PATH)


#### Configuration 2. Browser configuration

# Now move to the directory where you want to execute the omicser
# and make the omicser_options.yml

omicser_options <- list(database_names=DB_NAME,
                        db_root_path=DB_ROOT_PATH,
                        conda_environment=CONDA_ENV)

omicser::write_config(omicser_options,in_path = OMICSER_RUN_DIR )

#### Launch browser

omicser::run_app(options = list(launch.browser = TRUE))



