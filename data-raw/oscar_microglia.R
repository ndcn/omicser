#######################################################################
#######################################################################
##
##  Microglia scRNAseq Transcriptomics
##  - from Oscar:
#
##
#######################################################################
#######################################################################
#  formely `oscar_A` which was derived from a Seurat file
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
DB_NAME = "oscar_microglia"
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
lab <- "Harari"
annotation_database <- "NA"


# ------------------------------------------
# 2. helper functions ----------------------
# ------------------------------------------



# ------------------------------------------
# 3. load data -----------------------------
# ------------------------------------------
RAW_DIR <- "ingest/Oscar_microglia1"


file_path <- file.path(RAW_DIR,"microglia_matrix.mtx")
counts <- Matrix::readMM(file_path) #DGTMATRIX .34GB

# META-DATA
file_path <- file.path(RAW_DIR,"microglia_meta.csv")
meta_data <- read.csv(file_path, sep = ",", header = TRUE, stringsAsFactors = FALSE)

# FEATURES
file_path <- file.path(RAW_DIR,"microglia_features.tsv")
features <- read.csv(file_path, sep = "\t", header = FALSE, stringsAsFactors = FALSE)
# features2 <- vroom::vroom(file_path,delim = "\t")
colnames(features) <- "genes"

# BARCODES
file_path <- file.path(RAW_DIR,"microglia_barcodes.tsv")
barcodes <- read.csv(file_path, sep = "\t", header = FALSE, stringsAsFactors = FALSE)
# barcodes2 <- vroom::vroom(file_path,delim = "\t")

# assert that meta_data$barcode == barcodes
all(meta_data$barcode == barcodes)
meta_data$sample_id <- paste0("sample_",1:length(meta_data$barcode))


dimnames(counts) <- list(
  features$genes,  # row names
  meta_data$sample_id # column names / or barcodes$V1
)

# fix the "status" entries..
meta_data$Status <- tolower(meta_data$Status)



meta_data$cluster_label <- paste0("Cluster_", meta_data$Cluster)
all_clusters <- unique(meta_data$cluster_label[order(meta_data$Cluster)])


# counts2 <- as(counts, "dgCMatrix") #0.257FGB
# counts3 <- as(counts, "dgeMatrix") #2.3GB

#save(data_matrix,file=out_file_path)


# ------------------------------------------
# 4. pack into anndata                    --
# ------------------------------------------
X <- t(counts)  #transpose so rows and colums correspond to anndata expectation
obs <- meta_data
var_ <- features


db_prefix = "core_data"
saveRDS(X, file = paste0(DB_DIR, "/", db_prefix, "_X.rds"))
saveRDS(obs, file = paste0(DB_DIR, "/", db_prefix, "_obs.rds"))
saveRDS(var_, file = paste0(DB_DIR, "/", db_prefix, "_var.rds"))


ad <- AnnData(
  X = X,
  obs = obs,
  var = var_,
)

ad




ad$write_h5ad(filename=file.path(DB_DIR,"core_data.h5ad"))





# ------------------------------------------
# 5. post processing                      --
# ------------------------------------------
require(anndata)
require(reticulate)
DB_NAME = "oscar_microglia"
DB_DIR = file.path("data-raw",DB_NAME)
RAW_DIR <- "ingest/Oscar_microglia1"

db_prefix = "core_data"
X = readRDS( file = paste0(DB_DIR, "/", db_prefix, "_X.rds"))
obs = readRDS(  file = paste0(DB_DIR, "/", db_prefix, "_obs.rds"))
var_ = readRDS( file = paste0(DB_DIR, "/", db_prefix, "_var.rds"))


ad <- read_h5ad(file.path(DB_DIR,"core_data.h5ad"))
ad


sc <- import("scanpy")


# normalize??
sc$pp$normalize_total(ad)
sc$pp$log1p(ad)

#
# #normpc_adata <- sc$pp$normalize_per_cell(ad,copy=TRUE)
# z_adata <- sc$pp$scale(ad,copy=TRUE)
# log_adata <- sc$pp$log1p(ad,copy=TRUE)



test_types <- c('wilcoxon','t-test_overestim_var')

diff_exp <- data.frame()
#log_adata <- sc$pp$log1p(ad,copy=TRUE)
condition_key = 'cluster_label'
comp_type <- "allVrest"
test_type <- test_types[1]
key <- paste0(test_type,"_", comp_type)
sc$tl$rank_genes_groups(ad, condition_key, method=test_type, key_added = key)
de_table <- sc$get$rank_genes_groups_df(ad, NULL, key=key)
de_table$condition <- comp_type
de_table$test_type <- test_type
diff_exp <- dplyr::bind_rows(diff_exp, de_table)# for this dataset its straightforward to do all comparisons...


test_type <- test_types[2]
key <- paste0(test_type,"_", comp_type)
sc$tl$rank_genes_groups(ad, condition_key, method=test_type, key_added = key)
de_table <- sc$get$rank_genes_groups_df(ad, NULL, key=key)
de_table$condition <- comp_type
de_table$test_type <- test_type
diff_exp <- dplyr::bind_rows(diff_exp, de_table)# for this dataset its straightforward to do all comparisons...




# put the logvals in layers of ad
# copy the whole thing and replace X to copy the uns to ad



ad$write_h5ad(filename=file.path(DB_DIR,"core_data_plus_de.h5ad"))



# also need to pack the diff_exp1 and diff_exp2 into easy to deal wiht tables for volcanos...


db_prefix = "de"
saveRDS(diff_exp, file = paste0(DB_DIR, "/", db_prefix, "_table.rds"))


# ------------------------------------------
# 6. create config and default files                   --
# ------------------------------------------

#################
##################
##################


#usethis::use_data(oscar_A_obsm,oscar_A_varm,overwrite = TRUE)


# OBSERVABLES
#
marginal_measures <- c("fracexp", "meanexp", "sdexp", "mulogexp", "sdlogexp")
observables <- list(obs = marginal_measures,
                    var = marginal_measures,
                    layers = NA,
                    raw = c("X"),
                    obsm = NA) # this might not even be possible


# COMPARABLES
comparables <- list(varm = names(varm),
                    obsm = NA)
# Dimred
dimreds <- list(varm = NA,
                obsm = NA)


#######################################################################
#######################################################################
##
##  Create data tables / h5
##  1. config (XXX_conf.rds)  : structure of the database
##  2. defaults (XXX_def.rds) : the startup parameters
##  3. meta_data (XXX_meta.rds)  : data about the observations (obs)
##  4. omics  (XXX_omics.rds):  data about the  (var)
##  3. gene_expression (xxx_gexpr.h5):  count matrix of observ X var (X)
##
#######################################################################
#######################################################################

# X <- omicser::oscar_A_X
# obs <- omicser::oscar_A_obs
# var_ <- omicser::oscar_A_var
# obsm <- omicser::oscar_A_obsm
# varm <- omicser::oscar_A_varm
# uns <- omicser::oscar_A_uns


#TODO:  pack into a list or anndata structure for simplicity...

require("data.table")
## TODO:  check
##      - do we need default_1, 2 and multi?
##

#helper_functions<-('data-raw/data_prep_functions.R')
helper_function<-('data-raw/make_ingest_file_primitives.R')

source(helper_function)

db_dir = "data-raw"
db_prefix <- "Oscar_A_"

make_ingest_file_primitives(X,obs,var_,obsm,varm,uns, layers,
                            observables, comparables, dimreds,
                            default_omic = NA, default_dimred = NA, meta_to_include = NA,
                            gex.assay = NA, gex.slot = c("data", "scale.data", "counts"),
                            gene_mapping = FALSE, db_prefix = db_prefix, db_dir = db_dir,
                            chunk_size = 500,  legend_cols = 4,
                            max_levels_ui = 50)
#
# make_ingest_file_primitives(X,obs,var_,obsm=obsm, varm=varm,
#                              uns=uns, gex.assay = NA, gex.slot = c("data", "scale.data", "counts"),
#                              gene_mapping = FALSE, db_prefix = db_prefix, db_dir = "data-raw",
#                                         default_omics1 = NA, default_omics2 = NA, default_multi = NA,
#                                         default_dimred = NA, chunk_size = 500, meta_to_include = NA, legend_cols = 4,
#                                         max_levels_ui = 50)

#oscar_A_conf = readRDS(file.path(db_dir,"test1conf.rds"))
oscar_A_conf = readRDS( paste0(db_dir,"/",db_prefix,"conf.rds") )

# defaults:  list of meta1, meta2, omics1, omics2, omics (list of 10). dimred, grp1, grp2
oscar_A_def = readRDS( paste0(db_dir,"/",db_prefix,"def.rds") )

# list of vars )e/g/ 3000 genes with counts?
oscar_A_omics = readRDS( paste0(db_dir,"/",db_prefix,"omics.rds") )
# use this sorted one to resort everything before packing into anndata

oscar_A_meta = readRDS( paste0(db_dir,"/",db_prefix,"meta.rds") )
oscar_A_X = readRDS( paste0(db_dir,"/",db_prefix,"X.rds") )
oscar_A_obs = readRDS( paste0(db_dir,"/",db_prefix,"obs.rds") )
oscar_A_obsm = readRDS( paste0(db_dir,"/",db_prefix,"obsm.rds") )
oscar_A_var = readRDS( paste0(db_dir,"/",db_prefix,"var.rds") )
oscar_A_varm = readRDS( paste0(db_dir,"/",db_prefix,"varm.rds") )
oscar_A_uns = readRDS( paste0(db_dir,"/",db_prefix,"uns.rds") )
oscar_A_layers = readRDS( paste0(db_dir,"/",db_prefix,"layers.rds") )



# usethis::use_data(oscar_A_X,oscar_A_var,oscar_A_obs, overwrite = TRUE)
# usethis::use_data(oscar_A_obsm,oscar_A_varm,oscar_A_layers,oscar_A_uns, overwrite = TRUE)
usethis::use_data(oscar_A_conf, oscar_A_def, oscar_A_omics, oscar_A_meta, overwrite = TRUE)

#######################################################################
#######################################################################
##
##  Create some "differential tables" of type:
##  1. aggregated over "observations"
##  2. aggregated over "features"
##  3. aggretated over features and observations.
##
#######################################################################
#######################################################################

#######################################################################
#######################################################################
##
##  ANNDATA example
##
#######################################################################
#######################################################################

library(anndata)

#
# X <- omicser::oscar_A_X
# obs <- omicser::oscar_A_obs
# var_ <- omicser::oscar_A_var
# obsm <- omicser::oscar_A_obsm
# varm <- omicser::oscar_A_varm
# uns <- omicser::oscar_A_varm
# layers <- omicser::oscar_A_layers
#
#
#
# a_X <- omicser::oscar_A_X
# a_obs <- omicser::oscar_A_obs
# a_var_ <- omicser::oscar_A_var
# a_obsm <- omicser::oscar_A_obsm
# a_varm <- omicser::oscar_A_varm
# a_uns <- omicser::oscar_A_varm
# a_layers <- omicser::oscar_A_layers
#
#
# X <- omicser::Domenico_A_X
# obs <- omicser::Domenico_A_obs
# var_ <- omicser::Domenico_A_var
# obsm <- omicser::Domenico_A_obsm
# varm <- omicser::Domenico_A_varm
# uns <- omicser::Domenico_A_uns
#
#

X <- oscar_A_X
obs <- oscar_A_obs
var_ <- oscar_A_var
obsm <- oscar_A_obsm
varm <- oscar_A_varm
uns <- oscar_A_uns
layers <- oscar_A_layers

ad <- AnnData(
  X = X,
  obs = obs,
  var = var_,
  layers = layers,
  obsm = obsm,
  varm = varm,
  uns = uns
)

ad

#write_h5ad(anndata = ad, filename = file.path(db_dir,"data-raw/Vilas_A.h5ad"))
# anndata R wrapper is broken.. .invoke python
#
ad$write_h5ad(filename="data-raw/oscar_A.h5ad")

#######################################################################
#######################################################################
##
##  Create some "differential tables" of type:
##  1. aggregated over "observations"
##  2. aggregated over "features"
##  3. aggretated over features and observations.
##
#######################################################################
#######################################################################

