# curate pbmnc3k
#
#
## The data consist of *3k PBMCs from a Healthy Donor* and are freely available from 10x Genomics
# ([here](http://cf.10xgenomics.com/samples/cell-exp/1.1.0/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz)
# from this [webpage](https://support.10xgenomics.com/single-cell-gene-expression/datasets/1.1.0/pbmc3k)).
#
#  last line creates a directory for writing processed data.
#  We will work with the processed PBMC3K dataset. This is a dataset of around
#3000 peripheral blood mononuclear cells that was produced by 10x Genomics and
# has been processed as described in the scanpy PBMC3K tutorial.

#


DB_NAME <-  list("10x PBMC3k" = "pbmc3k")

require(omicser)
OMICSER_RUN_DIR <- file.path(getwd(),"examples")
setwd(OMICSER_RUN_DIR) # incase we have relative dir
# BOOTSTRAP the options we have already set up...
# NOTE: we are looking in the "quickstart" folder.  the default is to look for the config in with default getwd()



# if you already have an omicser_options.yml for other databases run the following lines
if (FALSE) {
  omicser_options <- omicser::get_config(in_path = OMICSER_RUN_DIR)
  CONDA_ENV <- omicser_options$conda_environment
  DB_ROOT_PATH <- omicser_options$db_root_path
  if (! (DB_NAME %in% omicser_options$database_names)){
    omicser_options$database_names <- c(omicser_options$database_names,DB_NAME)
    omicser::write_config(omicser_options,in_path = OMICSER_RUN_DIR )
  }


} else {
  CONDA_ENV <- "omxr"
  DB_ROOT_PATH <- file.path(OMICSER_RUN_DIR,"test_db") #or as `relative path` = "test_db"

  omicser_options <- list(
      conda_environment = CONDA_ENV,
      db_root_path = DB_ROOT_PATH,
      database_names = DB_NAME)

}



DB_DIR = file.path(DB_ROOT_PATH,DB_NAME)
if (!dir.exists(DB_DIR)) {
  dir.create(DB_DIR)
}


#==== 1. documentation / provenance ==============================================================
# TODO:  markdown file or html with some copy about the database
#  - lab, paper link/name
#  summarize results / data origin whatever

db_meta <- list(
  organism = 'human',
  lab = "",
  source = "blook monocular cells",
  annotation_database =  "",
  title = "pbmc3k",
  omic_type = "Transcriptomics",
  measurment = "normalized counts",
  pub = "TBD",
  date = format(Sys.time(), "%a %b %d %X %Y"),
  url = "https://support.10xgenomics.com/single-cell-gene-expression/datasets/1.1.0/pbmc3k"
)

write_db_meta(db_meta,DB_NAME, db_root = DB_ROOT_PATH)

# ALSO create an `additional_info.Rmd` as in `inst/app/www/additional_info.Rmd`
#==== 2. helper functions =================================================================================


#==== 3. load data -========================================================================================


RAW_DIR <- file.path(OMICSER_RUN_DIR,"raw_data", "pbmc3k")
if (!dir.exists(RAW_DIR))  dir.create(RAW_DIR)

data_file <- 'http://cf.10xgenomics.com/samples/cell-exp/1.1.0/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz'
download.file(data_file, file.path(RAW_DIR,'pbmc3k_filtered_gene_bc_matrices.tar.gz'))
system( "tar -xzf raw_data/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz" )


#
# The data consist of *3k PBMCs from a Healthy Donor* and are freely available from 10x Genomics
# ([here](http://cf.10xgenomics.com/samples/cell-exp/1.1.0/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz)
# from this [webpage](https://support.10xgenomics.com/single-cell-gene-expression/datasets/1.1.0/pbmc3k)).
# On a unix system, you can uncomment and run the following to download and unpack the data. The
#  last line creates a directory for writing processed data.
#f we will work with the processed PBMC3K dataset. This is a dataset of around
#3000 peripheral blood mononuclear cells that was produced by 10x Genomics and
# has been processed as described in the scanpy PBMC3K tutorial.
#
#https://github.com/theislab/scanpy-tutorials/blob/master/pbmc3k.ipynb
#
#  # !mkdir data
# !wget http://cf.10xgenomics.com/samples/cell-exp/1.1.0/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz -O raw_data/pbmc3k_filtered_gene_bc_matrices.tar.gz
# !cd data; tar -xzf pbmc3k_filtered_gene_bc_matrices.tar.gz
# !mkdir write
#
sc <- reticulate::import("scanpy")

# Load the PBMC dataset
adata = sc$read_10x_mtx(
  'raw_data/pbmc3k/filtered_gene_bc_matrices/hg19/',  # the directory with the `.mtx` file
  var_names='gene_symbols',                # use gene symbols for the variable names (variables-axis index)
  cache=TRUE)                              # write a cache file for faster subsequent reading

adata$write_h5ad(filename=file.path(DB_ROOT_PATH,DB_NAME,"core_data.h5ad"))

########### this is filtered but not many fields to "browse"
# data_file <- 'https://raw.githubusercontent.com/chanzuckerberg/cellxgene/main/example-dataset/pbmc3k.h5ad'
# download.file(data_file, file.path(RAW_DIR,'pbmc3k.h5ad'))


#==== 4. pack into anndata =========================================================================

data_list <- list(object=file.path(DB_ROOT_PATH,DB_NAME,"core_data.h5ad"))

adata <- omicser::setup_database(database_name=DB_NAME,
                              db_path=DB_ROOT_PATH,
                              data_in=data_list,
                              db_meta=NULL ,
                              re_pack=TRUE)

#adata <- anndata::read_h5ad(file.path(RAW_DIR,'pbmc3k.h5ad'))



# sc <- reticulate::import("scanpy")
# ad = sc$datasets$pbmc3k_processed()

#==== 5. post processing =========================================================================               --

### WARNING:  but if you re-read the object from the newly written file vs. keep using...
#

adata$var_names_make_unique()  # this is unnecessary if using `var_names='gene_ids'` in `sc.read_10x_mtx`

sc$pp$filter_cells(adata, min_genes=200)
sc$pp$filter_genes(adata, min_cells=10)

adata$var['mt'] = startsWith(adata$var_names,'MT-')  # annotate the group of mitochondrial genes as 'mt'
sc$pp$calculate_qc_metrics(adata, qc_vars=list('mt'), percent_top=NULL, log1p=FALSE, inplace=TRUE)


adata = adata[adata$obs$n_genes_by_counts < 2500, ]
adata = adata[adata$obs$pct_counts_mt < 5, ]
sc$pp$normalize_total(adata, target_sum=1e4)
sc$pp$log1p(adata)
sc$pp$highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)

adata$raw = adata$copy()


#adata = adata[, adata$var$highly_variable]
sc$pp$regress_out(adata, list('total_counts', 'pct_counts_mt'))
sc$pp$scale(adata, max_value=10)

adata$var$var_rank <- order(adata$var$dispersions_norm)
# choose top 40 genes by variance across dataset as our "targets"
target_omics <- adata$var_names[which(adata$var$var_rank <= 40)]


adata$var$decile <- dplyr::ntile(adata$var$dispersions_norm, 10)
#raw <- ad$raw$to_adata()


if (FALSE){ # save an intermediate file (incase we want to revert...)
  adata$write_h5ad(filename=file.path(DB_ROOT_PATH,DB_NAME,"normalized_data.h5ad"))
}
#==== 5-a. dimension reduction - PCA / umap  =========================================================================
# ## Step 2: Normalize with a very vanilla recipe
# normalized_data = sc$pp$recipe_seurat(raw, copy=TRUE)

## Step 3: Do some basic preprocessing to run PCA and compute the neighbor graph
sc$pp$pca(adata)
sc$pp$neighbors(adata)

## Step 4: Infer clusters with the Louvain algorithm
#sc$tl$louvain(ad_tmp)
sc$tl$leiden(adata)
## Step 5: Compute tsne and umap embeddings
sc$tl$umap(adata)
#ad_tmp$layers <- list(non_regressed=ad$X) #list('count'=layers)

if (FALSE){ # save an intermediate file (incase we want to revert...)
  adata$write_h5ad(filename=file.path(DB_ROOT_PATH,DB_NAME,"norm_data_plus_dr.h5ad"))
}

#==== 6. differential expression =========================================================================


test_types <- c('wilcoxon','t-test_overestim_var')


comp_types <- c("grpVrest")
obs_names <- c('leiden')
diff_exp <- omicser::compute_de_table(adata,comp_types, test_types, obs_names,sc)

if (FALSE){ # save an intermediate file (incase we want to revert...)
  adata$write_h5ad(filename=file.path(DB_ROOT_PATH,DB_NAME,"norm_data_with_de.h5ad"))
}

saveRDS(diff_exp, file = file.path(DB_ROOT_PATH,DB_NAME, "db_de_table.rds"))



#==== 7. create configs =========================================================================
if (FALSE) { # load intermediate files if available
  adata <- anndata::read_h5ad(filename=file.path(DB_ROOT_PATH,DB_NAME,"norm_data_with_de.h5ad"))
  diff_exp <- readRDS( file = file.path(DB_ROOT_PATH,DB_NAME, "db_de_table.rds"))
}

# what ad$obs do we want to make default values for...
default_factors <- c("leiden")

# differentials  #if we care we need to explicitly state. defaults will be the order...
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

  # Dimred
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



omicser::write_db_conf(config_list,DB_NAME, db_root = DB_ROOT_PATH)

#==== 8. write data file to load  =========================================================================
adata$write_h5ad(filename=file.path(DB_ROOT_PATH,DB_NAME,"db_data.h5ad"))




#==== 9. create meta-data narrative for rendering in INGEST  =========================================================================


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
