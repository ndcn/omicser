# INFO ========================================================
##  Microglia scRNAseq Transcriptomics
##  - from Vilas:
#
##
#  formely `Vila_A` which was derived from a online text repositories
#
# library(configr)
# configr::read.config('config.yml')
#

#==== 0. preamble/setup =========================================================================
library(magrittr)
require(reticulate)
reticulate::use_condaenv(required = TRUE, condaenv = 'omxr')
require(Matrix)

require(anndata)
require(data.table)
require(configr)



#==== 1. documentation / provenance =========================================================================

db_meta <- list(
  organism = '',
  lab = "Menon",
  annotation_database = NA,
  title = "single-cell microglia data",
  omic_type = "Transcriptomics",
  measurment = "normalized counts",
  pub = "TBD",
  date = format(Sys.time(), "%a %b %d %X %Y")
)

RAW_DIR <- "ingest/Vilas_A"
DB_NAME = "vilas_microglia"

out_fn <-  file.path("data-raw",DB_NAME,"db_meta.yml")
write.config(config.dat = db_meta, file.path = out_fn,
             write.type = "yaml", indent = 4)


#==== 2. helper functions =========================================================================
# "helper" functions to prep your "raw" data go here or should be sourced here...

#==== 3. load data =========================================================================
RAW_DIR <- "ingest/Vilas_A"

# 3a. matrix data --------------------
#  #check for raw_data_matrix.rds and load if exist
if (file.exists(file.path(RAW_DIR, "raw_data_matrix.rds")) ){
  data_mat <- readRDS(data_mat, file = file.path(RAW_DIR, "raw_data_matrix.rds"))

  omics <- row.names(data_mat)
  sample_ID = col.names(data_mat)
  var_annot = data.frame(omics_name=omics,row.names = omics)

} else {
  # else read the csv
  # TODO:  use readr or vroom instread of fread sinnce we are not using data.tables at this stage
  matrix_data_file <- "41467_2020_19737_MOESM15_ESM.csv"
  #count_table <- read.csv(file=file.path(RAW_DIR, matrix_data_file), header=TRUE, sep=",", row.names=1)
  count_table <- fread(file=file.path(RAW_DIR, matrix_data_file))
  omics <- count_table$V1
  sample_ID = colnames(count_table)[-1]

  var_annot = data.frame(omics_name=omics,row.names = omics)

  data_mat <- as.matrix(count_table,rownames=1)
  data_mat <- as(data_mat, "dgTMatrix")
  #data_mat <- t(data_mat)  #transpose so samples are rows (as in ANNDATA format)
  # data1 <- as.matrix(data) #2.2GB
  # data2 <- as(data1, "dgeMatrix") #4.4GB
  # data3 <- as(data1, "dgTMatrix") #0.29GB
  saveRDS(data_mat, file = file.path(RAW_DIR, "raw_data_matrix.rds"))
}


# 3b. load annotation data --------------------
file_path <- file.path(RAW_DIR,"41467_2020_19737_MOESM17_ESM.csv")
#annots <- read_csv(file_path)
obs_meta_file <- "41467_2020_19737_MOESM17_ESM.csv"
obs_meta <- fread(file=file.path(RAW_DIR, obs_meta_file))
# add rowname for convenience (even though its not data.table kosher)
row.names(obs_meta) <- obs_meta$sample_id
obs_meta$cluster_name = paste0("Cluster_",obs_meta$cluster_label)



#==== 4. pack into anndata =========================================================================

DB_NAME = "vilas_microglia"
helper_function<-('data-raw/ingest_helpers.R')
source(helper_function)

data_list <- list(data_mat = data_mat,
                  obs_meta = obs_meta,
                  var_annot =var_annot,
                  omics = omics,
                  sample_ID = sample_ID,
                  etc = NULL)

ad <- setup_database(database_name=DB_NAME, data_in=data_list, db_meta=db_meta , re_pack=TRUE)
ad$write_h5ad(filename=file.path("data-raw",DB_NAME,"core_data.h5ad"))




#==== 5. post processing =========================================================================
require(reticulate)
reticulate::use_condaenv(required = TRUE, condaenv = 'omxr')
require(anndata)

DB_NAME = "vilas_microglia"
ad <- read_h5ad(file.path("data-raw",DB_NAME,"core_data.h5ad"))
ad

sc <- import("scanpy")

sc$pp$filter_cells(ad, min_genes=200)
sc$pp$filter_genes(ad, min_cells=3)

raw <- ad$copy()

# to conform with cellXgene scheme gene filtering is tricky.
#
# we need to retain the list / shape of the original genes so we must

sc$pp$normalize_per_cell(ad, counts_per_cell_after=1e4)
sc$pp$log1p(ad)  # already logg-ed?
sc$pp$highly_variable_genes(ad, min_mean=0.0125, max_mean=3, min_disp=0.5  )
#sc$pp$regress_out(ad_tmp, 'n_counts' )


ad$raw <- raw
ad$write_h5ad(filename=file.path("data-raw",DB_NAME,"normalized_data.h5ad"))

#==== 5-a. dimension reduction - PCA / umap  =========================================================================
# ## Step 2: Normalize with a very vanilla recipe
# normalized_data = sc$pp$recipe_seurat(raw, copy=TRUE)
## Step 3: Do some basic preprocessing to run PCA and compute the neighbor graph
sc$pp$pca(ad)
sc$pp$neighbors(ad)

## Step 4: Infer clusters with the Louvain algorithm
#sc$tl$louvain(ad_tmp)
sc$tl$leiden(ad)
## Step 5: Compute tsne and umap embeddings
sc$tl$umap(ad)
#ad_tmp$layers <- list(non_regressed=ad$X) #list('count'=layers)

ad$write_h5ad(filename=file.path("data-raw",DB_NAME,"norm_data_plus_dr.h5ad"))


#==== 6. differential expression =========================================================================

sc <- import("scanpy")
helper_function<-('data-raw/compute_de_table.R')
source(helper_function)

test_types <- c('wilcoxon','t-test_overestim_var')


comp_types <- c("grpVrest")

obs_names <- c('cluster_name','leiden')
diff_exp <- compute_de_table(ad,comp_types, test_types, obs_names)


ad$write_h5ad(filename=file.path("data-raw",DB_NAME,"norm_data_with_de.h5ad"))
saveRDS(diff_exp, file = file.path("data-raw",DB_NAME, "diff_expr_table.rds"))



#==== 7. create configs =========================================================================
DB_NAME = "vilas_microglia"

#ad <- read_h5ad(filename=file.path("data-raw",DB_NAME,"normalized_data.h5ad"))
ad <- read_h5ad(filename=file.path("data-raw",DB_NAME,"norm_data_with_de.h5ad"))
diff_exp <- readRDS( file = file.path("data-raw",DB_NAME, "diff_expr_table.rds"))


# differentials  #if we care we need to explicitly state. defaults will be the order...
conf_list <- list(
  x_obs = c('cluster_name', 'leiden', 'batch', 'color' ),
  y_obs = c("n_genes","n_counts"), #MEASURES
  obs_groupby = c('cluster_name', 'leiden', 'batch', 'color' ),
  obs_subset = c('cluster_name', 'leiden', 'batch', 'color' ),
  x_var = c("highly_variable"),
  y_var = c("n_cells",  "means" , "dispersions", "dispersions_norm" ),
  var_groupby = c("highly_variable"),
  var_subset = c("highly_variable"),  # NOTE:  <omic selector> is NOT in the data object so its not actually going to load

  layers = c("X","raw"),

  diffs = list(diff_exp_comps = levels(factor(diff_exp$versus)),
               diff_exp_comp_type =  levels(factor(diff_exp$comp_type)), #i don't think we need this
               diff_exp_obs_name =  levels(factor(diff_exp$obs_name)),
               diff_exp_tests =  levels(factor(diff_exp$test_type))
  ),

  # Dimred
  dimreds = list(obsm = c('X_pca', 'X_tsne'),
                 varm = c('PCs')),

  # what ad$obs do we want to make default values for...
  # # should just pack according to UI?
  default_factors = c('cluster_name', 'leiden', 'batch')
)

configr::write.config(config.dat = conf_list, file.path = file.path("data-raw",DB_NAME,"config.yml" ),
                      write.type = "yaml", indent = 4)


#==== 8. write data file to load  =========================================================================
ad$write_h5ad(filename=file.path("data-raw",DB_NAME,"omxr_data.h5ad"))




#==== DONE- =========================================================================

ad <- read_h5ad(filename=file.path("data-raw",DB_NAME,"omxr_data.h5ad"))
conf_list_out <- configr::read.config( file.path("data-raw",DB_NAME,"config.yml" ) )
