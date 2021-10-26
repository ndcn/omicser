#### Create an app to browse PBMC3k data from 10X Genomics ####

#### Environment setup ####

library(devtools) # to install omicser package
library(reticulate) # runs python when curating data
library(fs) # file management

# name of conda environment
CONDA_ENV <- "omxr"

# install and configure environment
reticulate::install_miniconda()
reticulate::conda_create(CONDA_ENV, python_version = 3.9)
reticulate::conda_install(envname=CONDA_ENV,
                          channel = "conda-forge",
                          packages = c("scanpy","leidenalg") )
# load conda environment
reticulate::use_condaenv(condaenv = CONDA_ENV, required = TRUE)

# install omicser package
devtools::install_github("ndcn/omicser")

# location where data, database, and browser files will be created
dir_create("omicser_test")
setwd("omicser_test/")
OMICSER_RUN_DIR <- getwd() # run directory; must be absolute path
DB_ROOT_PATH <- path(OMICSER_RUN_DIR,"test_db") # location of database
DB_NAME <- list("10x PBMC3k" = "pbmc3k") # name of database(s)
RAW_DIR <- path(OMICSER_RUN_DIR,"raw_data", DB_NAME) # location of data

#### Download and organize data ###

# create directory structure for data and databases
dir_create(c(OMICSER_RUN_DIR, path(DB_ROOT_PATH, DB_NAME), RAW_DIR))
# change paths to make data manipulations easier
setwd(RAW_DIR)
# download data
data_file <- "http://cf.10xgenomics.com/samples/cell-exp/1.1.0/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz"
download.file(data_file, "pbmc3k_filtered_gene_bc_matrices.tar.gz")
# extract all downloaded files
untar("pbmc3k_filtered_gene_bc_matrices.tar.gz")
# compress individual files (required for scanpy)
tar("filtered_gene_bc_matrices/hg19/matrix.mtx.gz", files = "filtered_gene_bc_matrices/hg19/matrix.mtx", compression = "gzip")
tar("filtered_gene_bc_matrices/hg19/barcodes.tsv.gz", files = "filtered_gene_bc_matrices/hg19/barcodes.tsv", compression = "gzip")
tar("filtered_gene_bc_matrices/hg19/genes.tsv.gz", files = "filtered_gene_bc_matrices/hg19/genes.tsv", compression = "gzip")
setwd(OMICSER_RUN_DIR)

#### Data curation: 1. Metadata documentation ####

# load browser app package (should also load dependencies)
library(omicser)

# aggregate metadata for dataset
db_meta <- list(
  organism = "human",
  lab = "",
  source = "peripheral blood mononuclear cells (PBMCs)",
  annotation_database =  "",
  title = "pbmc3k",
  omic_type = "Transcriptomics",
  measurment = "normalized counts",
  pub = "10X Genomics",
  date = format(Sys.time(), "%a %b %d %X %Y"),
  url = "https://support.10xgenomics.com/single-cell-gene-expression/datasets/1.1.0/pbmc3k"
)

# write metadata configuration to db_meta.yml
write_db_meta(db_meta, DB_NAME, db_root = DB_ROOT_PATH)

#### Data curation 2. Format and ingest raw data ####

# make scanpy functions available
sc <- reticulate::import("scanpy")

# load the PBMC dataset
adata <- sc$read_10x_mtx(
  # locate directory containing mtx file
  "raw_data/pbmc3k/filtered_gene_bc_matrices/hg19/",
  # use gene symbols for the variable names (variables-axis index)
  var_names='gene_symbols',
  # write a cache file for faster subsequent reading
  cache=TRUE)
# save as file
adata$write_h5ad(path(DB_ROOT_PATH, DB_NAME,"core_data.h5ad"))

# identify location of raw data
data_list <- list(object=path(DB_ROOT_PATH,DB_NAME,"core_data.h5ad"))

# create database formatted as AnnData
adata <- omicser::setup_database(database_name = DB_NAME,
                                 db_path = DB_ROOT_PATH,
                                 data_in = data_list,
                                 db_meta = NULL,
                                 re_pack = TRUE)

#### Data curation 3. Post-processing ####

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

## Dimension reduction - PCA / umap
# run PCA
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

#### Data curation 4. Differential expression ####

# identify stats
test_types <- c('wilcoxon','t-test_overestim_var')
comp_types <- c("grpVrest")
obs_names <- c('leiden')
# calculate DE
diff_exp <- omicser::compute_de_table(adata,comp_types, test_types, obs_names,sc)

# save intermediate database file
if (FALSE){
  adata$write_h5ad(filename=file.path(DB_ROOT_PATH,DB_NAME,"norm_data_with_de.h5ad"))
}

# save DE tables
saveRDS(diff_exp, file = file.path(DB_ROOT_PATH, DB_NAME, "db_de_table.rds"))

#### Data curation 5. Write database ####

# write final database
adata$write_h5ad(filename = file.path(DB_ROOT_PATH, DB_NAME, "db_data.h5ad"))

#### Configuration 1. Database configuration ####

# load intermediate files if available
if (FALSE) {
  adata <- anndata::read_h5ad(filename=file.path(DB_ROOT_PATH,DB_NAME,"norm_data_with_de.h5ad"))
  diff_exp <- readRDS( file = file.path(DB_ROOT_PATH,DB_NAME, "db_de_table.rds"))
}

# set default values for ad$obs
default_factors <- c("leiden")

# differentials
# if we care we need to explicitly state. defaults will be the order...
config_list <- list(
  x_obs = c("leiden"),
  y_obs =  c('n_genes', 'n_genes_by_counts', 'total_counts', 'total_counts_mt', 'pct_counts_mt'), #MEASURES
  obs_groupby = c("leiden"),
  obs_subset = c("leiden"),

  x_var = c("decile","highly_variable"),
  y_var = c( 'n_cells', 'mt', 'n_cells_by_counts', 'mean_counts',
            'pct_dropout_by_counts', 'total_counts',
            'means', 'dispersions', 'dispersions_norm', 'mean', 'std'),

  var_groupby = c("decile","highly_variable"),
  var_subset = c("decile","highly_variable"),

  diffs = list(diff_exp_comps = levels(factor(diff_exp$versus)),
               diff_exp_comp_type =  levels(factor(diff_exp$comp_type)), #i don't think we need this
               diff_exp_obs_name =  levels(factor(diff_exp$obs_name)),
               diff_exp_tests =  levels(factor(diff_exp$test_type))
  ),

  layers = c("X","raw"),

  # dimension reduction
  dimreds = list(obsm = adata$obsm_keys(),
                 varm = adata$varm_keys()),

  # what ad$obs do we want to make default values for...
  # # should just pack according to UI?
  default_factors = default_factors,
  target_omics = target_omics,
  omic_details = c('gene_ids', 'n_cells', 'mt', 'n_cells_by_counts',
                   'mean_counts', 'pct_dropout_by_counts', 'total_counts',
                   'highly_variable', 'means', 'dispersions',
                   'dispersions_norm', 'mean', 'std', 'var_rank', 'decile')
)


# write db_config.yml
omicser::write_db_conf(config_list,DB_NAME, db_root = DB_ROOT_PATH)

#### Configuration 2. Browser configuration ####
# Now move to the directory where you want to execute the omicser
# and make the omicser_options.yml

omicser_options <- list(database_names=DB_NAME,
                        db_root_path=DB_ROOT_PATH,
                        conda_environment=CONDA_ENV)

omicser::write_config(omicser_options,in_path = OMICSER_RUN_DIR )

#### Launch browser ####

omicser::run_app(options = list(launch.browser = TRUE))

#==== create meta-data narrative for rendering in INGEST  =========================================================================

# ALSO create an `additional_info.Rmd` as in `inst/app/www/additional_info.Rmd`

# ALTERNATE SEURAT STYLE =================
#'
#' # FROM SEURAT::
#' #' Read 10X hdf5 file
#' #'
#' #' Read count matrix from 10X CellRanger hdf5 file.
#' #' This can be used to read both scATAC-seq and scRNA-seq matrices.
#' #'
#' #' @param filename Path to h5 file
#' #' @param use.names Label row names with feature names rather than ID numbers.
#' #' @param unique.features Make feature names unique (default TRUE)
#' #'
#' #' @return Returns a sparse matrix with rows and columns labeled. If multiple
#' #' genomes are present, returns a list of sparse matrices (one per genome).
#' #'
#' #' @export
#' #' @concept preprocessing
#' #'
#' Read10X_h5 <- function(filename, use.names = TRUE, unique.features = TRUE) {
#'   if (!requireNamespace('hdf5r', quietly = TRUE)) {
#'     stop("Please install hdf5r to read HDF5 files")
#'   }
#'   if (!file.exists(filename)) {
#'     stop("File not found")
#'   }
#'   infile <- hdf5r::H5File$new(filename = filename, mode = 'r')
#'   genomes <- names(x = infile)
#'   output <- list()
#'   if (hdf5r::existsGroup(infile, 'matrix')) {
#'     # cellranger version 3
#'     if (use.names) {
#'       feature_slot <- 'features/name'
#'     } else {
#'       feature_slot <- 'features/id'
#'     }
#'   } else {
#'     if (use.names) {
#'       feature_slot <- 'gene_names'
#'     } else {
#'       feature_slot <- 'genes'
#'     }
#'   }
#'   for (genome in genomes) {
#'     counts <- infile[[paste0(genome, '/data')]]
#'     indices <- infile[[paste0(genome, '/indices')]]
#'     indptr <- infile[[paste0(genome, '/indptr')]]
#'     shp <- infile[[paste0(genome, '/shape')]]
#'     features <- infile[[paste0(genome, '/', feature_slot)]][]
#'     barcodes <- infile[[paste0(genome, '/barcodes')]]
#'     sparse.mat <- sparseMatrix(
#'       i = indices[] + 1,
#'       p = indptr[],
#'       x = as.numeric(x = counts[]),
#'       dims = shp[],
#'       giveCsparse = FALSE
#'     )
#'     if (unique.features) {
#'       features <- make.unique(names = features)
#'     }
#'     rownames(x = sparse.mat) <- features
#'     colnames(x = sparse.mat) <- barcodes[]
#'     sparse.mat <- as(object = sparse.mat, Class = 'dgCMatrix')
#'     # Split v3 multimodal
#'     if (infile$exists(name = paste0(genome, '/features'))) {
#'       types <- infile[[paste0(genome, '/features/feature_type')]][]
#'       types.unique <- unique(x = types)
#'       if (length(x = types.unique) > 1) {
#'         message("Genome ", genome, " has multiple modalities, returning a list of matrices for this genome")
#'         sparse.mat <- sapply(
#'           X = types.unique,
#'           FUN = function(x) {
#'             return(sparse.mat[which(x = types == x), ])
#'           },
#'           simplify = FALSE,
#'           USE.NAMES = TRUE
#'         )
#'       }
#'     }
#'     output[[genome]] <- sparse.mat
#'   }
#'   infile$close_all()
#'   if (length(x = output) == 1) {
#'     return(output[[genome]])
#'   } else{
#'     return(output)
#'   }
#' }
# pbmc_data <- Read10X(data.dir = "raw_data/pbmc3k/filtered_gene_bc_matrices/hg19/")
