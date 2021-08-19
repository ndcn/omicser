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
#
# Golem options (Inactive) ---------------------------
# Set options here
# options(golem.app.prod = FALSE) # TRUE = production mode, FALSE = development mode
# # Detach all loaded packages and clean your environment
# golem::detach_all_attached()
# # rm(list=ls(all.names = TRUE))
# # Document and reload your package
# golem::document_and_reload()

# ------------------------------------------
# 0. preamble/setup -------------------------
# ------------------------------------------

require(reticulate)
reticulate::use_condaenv(required = TRUE, condaenv = 'sc39')
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

# TODO:  markdown file or html with some copy about the database
#  - lab, paper link/name
#  summarize results / data origin whatever

organism <- ""
lab <- "Menon"
annotation_database <- "NA"


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

# ------------------------------------------
# 4. pack into anndata                    --
# ------------------------------------------
X <- Matrix::t(countm)  #transpose so rows and colums correspond to anndata expectation
obs <- obs_annots
var_ <- data.frame(features,row.names = features)


db_prefix = "core_data"
saveRDS(X, file = paste0(DB_DIR, "/", db_prefix, "_X.rds"))
saveRDS(obs, file = paste0(DB_DIR, "/", db_prefix, "_obs.rds"))
saveRDS(var_, file = paste0(DB_DIR, "/", db_prefix, "_var.rds"))


ad <- AnnData(
  X = X,
  obs = obs,
  var = as.data.frame(var_),
)
ad


ad$write_h5ad(filename=file.path(DB_DIR,"core_data.h5ad"))



# ------------------------------------------
# 5. post processing                      --
# ------------------------------------------
require(anndata)
require(reticulate)
DB_NAME = "vilas_microglia"
DB_DIR = file.path("data-raw",DB_NAME)
RAW_DIR <- "ingest/Vilas_A"

db_prefix = "core_data"
X = readRDS( file = paste0(DB_DIR, "/", db_prefix, "_X.rds"))
obs = readRDS(  file = paste0(DB_DIR, "/", db_prefix, "_obs.rds"))
var_ = readRDS( file = paste0(DB_DIR, "/", db_prefix, "_var.rds"))


ad <- read_h5ad(file.path(DB_DIR,"core_data.h5ad"))
ad


sc <- import("scanpy")
#
# #normpc_adata <- sc$pp$normalize_per_cell(ad,copy=TRUE)
# #norm_adata <- sc$pp$normalize_total(ad,copy=TRUE)
#z_adata <- sc$pp$scale(ad,copy=TRUE)
# log_adata <- sc$pp$log1p(ad,copy=TRUE)

#sc$pp$normalize_total(ad)
sc$pp$recipe_seurat(ad)
#sc$pp$log1p(ad)




test_types <- c('wilcoxon','t-test_overestim_var')
comp_types <- c("allVrest")


helper_function<-('data-raw/compute_de_table.R')
source(helper_function)

diff_exp <- compute_de_table(ad,comp_types, test_types, obs_name = c('cluster_name'))


# put the logvals in layers of ad
# copy the whole thing and replace X to copy the uns to ad
ad$write_h5ad(filename=file.path(DB_DIR,"core_data_plus_de.h5ad"))

db_prefix = "de"
saveRDS(diff_exp, file = paste0(DB_DIR, "/", db_prefix, "_table.rds"))


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


usethis::use_data(vilas_microglia_conf, vilas_microglia_def, vilas_microglia_omics, vilas_microglia_meta, overwrite = TRUE)

