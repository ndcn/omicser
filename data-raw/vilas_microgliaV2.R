#######################################################################
#######################################################################
##
##  Microglia scRNAseq Transcriptomics
##  - from Vilas:
#
##
#######################################################################
#######################################################################
#  formely `Vilas` which was derived from a online text repositories
#
# library(configr)
# configr::read.config('config.yml')
#

# ------------------------------------------
# 0. preamble/setup -------------------------
# ------------------------------------------
library(magrittr)
require(reticulate)
reticulate::use_condaenv(required = TRUE, condaenv = 'omxr')
require(Matrix)

require(anndata)
# create the folder to contain the raw data
DB_NAME = "vilas_microglia"
DB_DIR = file.path("data-raw",DB_NAME)
if (!dir.exists(DB_DIR)) {
  dir.create(DB_DIR)
}

# ------------------------------------------
# 1. documentation / provenance ------------
# ------------------------------------------

db_meta <- list(
  organizm = '',
  lab = "Menon",
  annotation_database = NA
)


# ------------------------------------------
# 2. helper functions ----------------------
# ------------------------------------------


# ------------------------------------------
# 3. load data -----------------------------
# ------------------------------------------
RAW_DIR <- "ingest/Vilas_A"

file_path <- file.path(RAW_DIR,"41467_2020_19737_MOESM15_ESM.csv")
count_table <- read.csv(file=file_path, header=TRUE, sep=",", row.names=1)
#features <- data.frame("genes" = row.names(count_table))
features <- row.names(count_table)
countm <- as.matrix(count_table)
countm <- as(countm, "dgTMatrix")
#counts <- t(countm)  #transpose so samples are rows (as in ANNDATA format)

# data <- read.csv(file=file_path, header=TRUE, sep=",", row.names=1)
# data1 <- as.matrix(data) #2.2GB
# data2 <- as(data1, "dgeMatrix") #4.4GB
# data3 <- as(data1, "dgTMatrix") #0.29GB


# 2. load annotation data --------------------
file_path <- file.path(RAW_DIR,"41467_2020_19737_MOESM17_ESM.csv")
#annots <- read_csv(file_path)
obs_annots <- read.csv(file=file_path, header=TRUE, sep=",", row.names=NULL)
#obs_annots$obs_names <- obs$sample_id
obs_annots$cluster_name = paste0("Cluster_",obs_annots$cluster_label)
obs_annots$cluster_id = obs_annots$cluster_label
row.names(obs_annots) <- obs_annots$sample_id


# pack into X, obs, var, meta?
#
# scemaa


# ------------------------------------------
# 4. pack into anndata                    --
# ------------------------------------------


# scemaa
#c("object", "data_mat","obs_meta","var_annot","omics","sample_ID","etc")
data_list <- list(data_mat = (countm),
                  obs_meta = obs_annots,
                  var_annot = data.frame(features,row.names = features),
                  omics = features,
                  sample_ID = colnames(count_table),
                  etc = NULL)


DB_NAME = "vilas_microglia"

helper_function<-('R/fct_ingestor.R')
source(helper_function)

ad <- setup_database(database_name=DB_NAME, data_in=data_list, db_meta=NULL , re_pack=TRUE)
ad$write_h5ad(filename=file.path("data-raw",DB_NAME,"core_data.h5ad"))



# ------------------------------------------
# 5. post processing                      --
# ------------------------------------------
require(anndata)
require(reticulate)
DB_NAME = "vilas_microglia"


ad <- read_h5ad(file.path("data-raw",DB_NAME,"core_data.h5ad"))
ad
ad_tmp <- ad
raw <- ad$copy()

sc <- import("scanpy")

#outvals <- sc$pp$filter_cells(ad, min_genes=200, inplace=FALSE)

sc$pp$filter_cells(ad, min_genes=200)
sc$pp$filter_genes(ad, min_cells=3)
ad_tmp <- ad$copy()


# to conform with cellXgene scheme gene filtering is tricky.
#
# we need to retain the list / shape of the original genes so we must
# replace zeros/nulls with NA

# outvals <- sc$pp$filter_genes(ad, min_cells=2, inplace=FALSE)

sc$pp$normalize_per_cell(ad, counts_per_cell_after=1e4)
sc$pp$log1p(ad)
sc$pp$highly_variable_genes(ad, min_mean=0.0125, max_mean=3, min_disp=0.5  )
#sc$pp$regress_out(ad_tmp, 'n_counts' )

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

ad$raw <- ad_tmp
#ad_tmp$layers <- list(non_regressed=ad$X) #list('count'=layers)


ad$write_h5ad(filename=file.path("data-raw",DB_NAME,"processed_data.h5ad"))




#
# sc$pp$normalize_total(adata, target_sum=1e4)
#
# normpc_adata <- sc$pp$normalize_per_cell(ad,copy=TRUE)
# norm_adata <- sc$pp$normalize_total(ad,copy=TRUE)
#
# sc$pp$filter_cells(adata, min_genes=200)
# sc$pp$filter_genes(adata, min_cells=2)
#
# #z_adata <- sc$pp$scale(ad,copy=TRUE)
# # log_adata <- sc$pp$log1p(ad,copy=TRUE)
# #
# # sc$pp$normalize_total(ad)
# # sc$pp$log1p(ad)
# # sc$pp$recipe_seurat(ad, log=TRUE)
# # recipe_zheng17
#
# sc$pp$filter_genes(ad, min_counts=1)         # only consider genes with more than 1 count
# sc$pp$normalize_per_cell(                       # normalize with total UMI count per cell
#   ad, key_n_counts='n_counts_all'
# )
# # filter_result = sc$pp$filter_genes_dispersion(  # select highly-variable genes
# #   ad$X, flavor='cell_ranger', n_top_genes=3000, log=FALSE
# # )
# #adata = adata[:, filter_result.gene_subset]     # subset the genes
# sc$pp$normalize_per_cell(ad)                 # renormalize after filtering
# sc$pp$log1p(ad)                      # log transform: adata.X = log(adata.X + 1)
# sc$pp$scale(ad)                              # scale to unit variance and shift to zero mean
#
# #ad$write_csvs(file.path("data-raw",DB_NAME,"whoa.csv"))


test_types <- c('wilcoxon','t-test_overestim_var')
comp_types <- c("allVrest")



helper_function<-('data-raw/compute_de_table.R')
source(helper_function)

diff_exp <- compute_de_table(ad,comp_types, test_types, obs_name = c('cluster_name'))

#ad$raw <- raw

# put the logvals in layers of ad
# copy the whole thing and replace X to copy the uns to ad
ad$write_h5ad(filename=file.path("data-raw",DB_NAME,"core_data_plus_de.h5ad"))

db_prefix = "de"
saveRDS(diff_exp, file = paste0("data-raw/",DB_NAME, "/", db_prefix, "_table.rds"))


# ------------------------------------------
# 6. dimension reduction - PCA / umap    --
# ------------------------------------------
require(anndata)
require(reticulate)

DB_NAME = "vilas_microglia"
DB_DIR = file.path("data-raw",DB_NAME)
#RAW_DIR <- "ingest/Vilas_B"

ad <- read_h5ad(file.path(DB_DIR,"core_data_plus_de.h5ad"))





# ------------------------------------------
# 7 . create config and default files                   --
# ------------------------------------------
require(anndata)
require(reticulate)

DB_NAME = "vilas_microglia"
DB_DIR = file.path("data-raw",DB_NAME)
#RAW_DIR <- "ingest/Vilas_B"

ad <- read_h5ad(file.path(DB_DIR,"core_data_plus_de.h5ad"))
db_prefix = "de"
diff_exp = readRDS( file = paste0(DB_DIR, "/", db_prefix, "_table.rds"))





###

# measures
#  This ordering is the "default"
measures <- list(obs = c("n_genes"),
                 var = ad$var_keys())

# differentials  #if we care we need to explicitly state. defaults will be the order...
diffs <- list(diff_exp_groups =  levels(factor(diff_exp$group)),
              diff_exp_comp_type =  levels(factor(diff_exp$comp_type)),
              diff_exp_obs_name =  levels(factor(diff_exp$obs_name)),
              diff_exp_tests =  levels(factor(diff_exp$test_type)))


# what ad$obs do we want to make default values for...
# # should just pack according to UI?
default_factors <- c("cluster_name","batch","color")



# Dimred
dimreds <- list(varm = NA,
                obsm = NA)


helper_function<-('data-raw/create_config_table.R')

source(helper_function)

db_prefix = "omxr"
conf_and_def <- create_config_table(ad,
                                    measures,
                                    diffs,
                                    dimreds,
                                    default_factors,
                                    db_prefix= db_prefix,
                                    db_dir = DB_DIR)


vilas_microglia_conf = readRDS( paste0(DB_DIR,"/",db_prefix,"_conf.rds") )
vilas_microglia_def = readRDS( paste0(DB_DIR,"/",db_prefix,"_def.rds") )
vilas_microglia_omics = readRDS( paste0(DB_DIR,"/",db_prefix,"_omics.rds") )
vilas_microglia_meta = readRDS( paste0(DB_DIR,"/",db_prefix,"_meta.rds") )


# usethis::use_data(vilas_microglia_conf, vilas_microglia_def, vilas_microglia_omics, vilas_microglia_meta, overwrite = TRUE)

