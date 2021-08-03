## code to prepare `vilas_B` dataset goes here

# #
#         Here’s a link to our single-cell microglia data set as a
#         Seurat object, with counts tables, normalized counts tables,
#         UMAP coordinates, and cell metadata in the table. The UMAP
#         coordinates and clusters were generated with a previous v
#         ersion of Seurat, with obsolete normalization and clustering
#         routines. However, for visualization and to kick the tires,
#         I hope this is a good starting data set. This is also what
#         Chris Sifuentes has been playing with in cellxgene, so down
#         the road it could also be a good data set if we want to explore
#         cross-connectivity between visualization and cellxgene options.

library("Seurat")


# ingest script -----------------------
DATA_DIR <- "ingest"
DB_NAME <- "Vilas_B"

data_file <- "microglia_data_with_counts_RNA_SCT.rda"

file_name <- file.path(DB_NAME, data_file)
file_path <- file.path(DATA_DIR, file_name)

load(file_path) # microglia_data


# UPDATE OBJECT
new_microglia_data <- UpdateSeuratObject(object = microglia_data)
# # OUTPUT:
#   Validating object structure
#   Updating object slots
#   Ensuring keys are in the proper strucutre
#   Ensuring feature names don't have underscores or pipes
#   Object representation is consistent with the most current Seurat version

saveRDS(new_microglia_data, file = file.path(DATA_DIR, DB_NAME, "new_microglia_data.rds"))

# convert to ANNDATA
#
# if (!requireNamespace("remotes", quietly = TRUE)) {
#   install.packages("remotes")
# }
# remotes::install_github("mojaveazure/seurat-disk")
# devtools::install_github('satijalab/seurat-data')



library("Seurat")

DATA_DIR <- "ingest"
DB_NAME <- "Vilas_B"
library(SeuratData)
library(SeuratDisk)

out_file_name <- "new_microglia_data.h5Seurat"
out_file_path <- file.path(DATA_DIR, DB_NAME, out_file_name)

new_microglia_data <- readRDS(file.path(DATA_DIR, DB_NAME, "new_microglia_data.rds"))


SaveH5Seurat(new_microglia_data, filename = out_file_path, overwrite = TRUE)
# Creating h5Seurat file for version 3.1.5.9900
# Adding counts for RNA
# Adding data for RNA
# No variable features found for RNA
# No feature-level metadata found for RNA
# Adding counts for SCT
# Adding data for SCT
# Adding scale.data for SCT
# Adding variable features for SCT
# Adding feature-level metadata for SCT
# Adding counts for counts
# Adding data for counts
# No variable features found for counts
# No feature-level metadata found for counts
# Adding cell embeddings for pca
# Adding loadings for pca
# No projected loadings for pca
# Adding standard deviations for pca
# No JackStraw data for pca
# Adding cell embeddings for tsne
# No loadings for tsne
# No projected loadings for tsne
# No standard deviations for tsne
# No JackStraw data for tsne
#
#

Convert(out_file_path, dest = "h5ad",overwrite = TRUE)
######     from https://mojaveazure.github.io/seurat-disk/reference/Convert.html
######---- Assay data
######
###### X will be filled with scale.data if scale.data is present; otherwise,
######        it will be filled with data
###### var will be filled with meta.features only for the features present
######        in X; for example, if X is filled with scale.data, then var
######        will contain only features that have been scaled
###### raw.X will be filled with data if X is filled with scale.data;
######        otherwise, it will be filled with counts. If counts is not
######        present, then raw will not be filled
###### raw.var will be filled with meta.features with the features present
######        in raw.X; if raw.X is not filled, then raw.var will not be filled
######
######



# Validating h5Seurat file
# Adding scale.data from SCT as X  <-----
# Transferring meta.features to var
# Adding data from SCT as raw
# Transfering meta.features to raw/var
# Transfering meta.data to obs
# Adding dimensional reduction information for pca
# Adding feature loadings for pca
# Adding dimensional reduction information for tsne
#
##  Are there arguments to add to make sure that the counts, and "non-scaled" values are written to the file?

require(anndata)
new_file_path <- gsub(".h5Seurat", ".h5ad", out_file_path)
ad <- read_h5ad(new_file_path)

# NOTE:  raw counts are NOT contained in the anndata file... they could be added back with layers from the seurat file *if* needed

# now lets add uns entries for the varm and obsm variables....
#
uns <- list(frac_expres=colnames(varm$frac_expres),
            mean_z = colnames(varm$mean_z),
            mean_log = colnames(varm$mean_log),
            mean_z_log = colnames(varm$mean_z_log))



# OBSERVABLES
#
observables <- list(obs = c("nCount_SCT","nFeature_SCT"),
                    var = names(ad$var)
                    layers = NULL,
                    raw = c("X"))


# COMPARABLES
comparables <- list(varm = names$varm,
                    obsm = NULL)
# Dimred
dimreds <- list(varm = c('X_pca', 'X_tsne'),
                    obsm = c('PCs'))


####################################
####################################
##
##   read from .h5ad
##
####################################
####################################

helper_function<-('data-raw/make_ingest_file_primitives.R')

source(helper_function)


DATA_DIR <- "ingest"
DB_NAME <- "Vilas_B"


file_name <- file.path(DB_NAME, "new_microglia_data.h5ad")
file_path <- file.path(DATA_DIR, file_name)

# ui_conf = create_config(file_path,
#   meta_to_include = NA, legend_cols = 4,
#   max_levels = 50
# )
# inp_obj = file_path
#
# db_prefix = "Vilas_B_"
# db_dir = "data-raw"
# gene.mapping = FALSE
# default_omics1 = NA
# default_omics2 = NA
# default_multi = NA
# default_dimred = NA
# chunk_size = 500
# #gex.assay = "X" # X.raw" # X is filled with scale.data / data instead of counts... counts are in raw.X
# #gex.slot = "scale.data"  # only for Seurat.
#
# vilasb_conf = make_ingest_files_h5ad(inp_obj, ui_conf,
#   db_prefix = "Vilas_B_", db_dir = "data-raw",
#   gene.mapping = FALSE,
#   default_omics1 = NA, default_omics2 = NA, default_multi = NA,
#   default_dimred = NA, chunk_size = 500
# )

# extract_primitives_from_ad <- function(ad){
#   X <- ad$X
#   obs <- ad$obs
#   obsm <- ad$obsm
#   var_names <- ad$var_names
#
#   var_ <- ad$var
#   varm <- ad$varm
#   var_keys <- ad$var_keys()
#   uns <- ad$uns
#   layers <- ad$obs
#   rawX <- ad$raw$X
#   rawvar <- ad$raw$var
# }
#
# X <- ad$X
# obs <- ad$obs
# varm <- ad$obs

## TODO:  check if
## varm & obsm are backwards...
##

X <- ad$X
obs <- ad$obs
obsm <- ad$obsm
var_names <- ad$var_names

var_ <- ad$var
varm <- ad$varm
var_keys <- ad$var_keys()
uns <- ad$uns
layers <- ad$obs
rawX <- ad$raw$X
rawvar <- ad$raw$var


db_dir = "data-raw"
db_prefix <- "Vilas_B_"
make_ingest_file_primitives(X,obs,var_,obsm=obsm, varm=varm,
                            uns=uns, gex.assay = NA, gex.slot = c("data", "scale.data", "counts"),
                            gene_mapping = FALSE, db_prefix = db_prefix, db_dir = "data-raw",
                            default_omics1 = NA, default_omics2 = NA, default_multi = NA,
                            default_dimred = NA, chunk_size = 500, meta_to_include = NA, legend_cols = 4,
                            max_levels_ui = 50)

#vilas_B_conf = readRDS(file.path(db_dir,"test1conf.rds"))
vilas_B_conf = readRDS( paste0(db_dir,"/",db_prefix,"conf.rds") )

# defaults:  list of meta1, meta2, omics1, omics2, omics (list of 10). dimred, grp1, grp2
vilas_B_def = readRDS( paste0(db_dir,"/",db_prefix,"def.rds") )

# list of vars )e/g/ 3000 genes with counts?
vilas_B_omics = readRDS( paste0(db_dir,"/",db_prefix,"omics.rds") )
# use this sorted one to resort everything before packing into anndata

vilas_B_meta = readRDS( paste0(db_dir,"/",db_prefix,"meta.rds") )
vilas_B_X = readRDS( paste0(db_dir,"/",db_prefix,"X.rds") )
vilas_B_obs = readRDS( paste0(db_dir,"/",db_prefix,"obs.rds") )
vilas_B_obsm = readRDS( paste0(db_dir,"/",db_prefix,"obsm.rds") )
vilas_B_var = readRDS( paste0(db_dir,"/",db_prefix,"var.rds") )
vilas_B_varm = readRDS( paste0(db_dir,"/",db_prefix,"varm.rds") )
vilas_B_uns = readRDS( paste0(db_dir,"/",db_prefix,"uns.rds") )
vilas_B_layers = readRDS( paste0(db_dir,"/",db_prefix,"layers.rds") )



usethis::use_data(vilas_B_X,vilas_B_var,vilas_B_obs,
                  vilas_B_obsm,vilas_B_varm,vilas_B_layers,vilas_B_uns,
                  vilas_B_conf, vilas_B_def, vilas_B_omics, vilas_B_meta, overwrite = TRUE)

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


# X <- omicser::vilas_B_X
# obs <- omicser::vilas_B_obs
# var_ <- omicser::vilas_B_var
# obsm <- omicser::vilas_B_obsm
# varm <- omicser::vilas_B_varm
# uns <- omicser::vilas_B_varm
# layers <- omicser::vilas_B_layers
X <- vilas_B_X
obs <- vilas_B_obs
var_ <- vilas_B_var
obsm <- vilas_B_obsm
varm <- vilas_B_varm
uns <- vilas_B_uns
layers <- vilas_B_layers

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
ad$write_h5ad(filename="data-raw/vilas_B.h5ad")






new_file_path <- gsub(".h5Seurat", ".h5ad", out_file_path)
ad_og <- read_h5ad(new_file_path)

ad_og

ad_og$write_h5ad(filename="data-raw/vilas_Bog.h5ad")


# source("data-raw/domenico_A.R",echo = FALSE)

#
#
# ####################################
# ####################################
# ##
# ## read from .h5ad directly... (depricated)
# ##
# ####################################
# ####################################
# library("Seurat")
#
#
# DATA_DIR <- "ingest"
# DB_NAME <- "Vilas_B"
# library(SeuratData)
# library(SeuratDisk)
#
# out_file_name <- "new_microglia_data.h5Seurat"
# out_file_path <- file.path(DATA_DIR, DB_NAME, out_file_name)
#
# new_microglia_data <- readRDS(file.path(DATA_DIR, DB_NAME, "new_microglia_data.rds"))
#
#
#
#
#
# ui_conf2 = create_config(new_microglia_data,
#                         meta_to_include = NA, legend_cols = 4,
#                         max_levels = 50
# )
#
# db_prefix = "Vilas_B2_"
# db_dir = "data-raw"
# gene.mapping = FALSE
# default_omics1 = NA
# default_omics2 = NA
# default_multi = NA
# default_dimred = NA
# chunk_size = 500
# gex.assay = "SCT"  # this makes the data match the anndata version ,
# gex.slot = "scale.data" #  scale.data only includes 3000 genes instead of 18913 in counts/data
#
# vilasb_conf = make_ingest_files_etc(new_microglia_data, ui_conf2,
#                                      db_prefix = "Vilas_B2_", db_dir = "data-raw", gex.assay = gex.assay, gex.slot = gex.slot,
#                                      gene.mapping = FALSE,
#                                      default_omics1 = NA, default_omics2 = NA, default_multi = NA,
#                                      default_dimred = NA, chunk_size = 500
#                                     )
#
#
# vilas_B2_conf = readRDS(file.path(db_dir,"vilas_B2_conf.rds"))
# vilas_B2_def = readRDS(file.path(db_dir,"vilas_B2_def.rds"))
# vilas_B2_omics = readRDS(file.path(db_dir,"vilas_B2_omics.rds"))
# vilas_B2_meta = readRDS(file.path(db_dir,"vilas_B2_meta.rds"))
#
# file.copy(from=file.path(db_dir,"Vilas_B2_gexpr.h5"), to=file.path("data","Vilas_B2_gexpr.h5"),
#           overwrite = TRUE, recursive = FALSE,
#           copy.mode = TRUE)
#
#
# usethis::use_data(vilas_B2_conf, vilas_B2_def, vilas_B2_omics, vilas_B2_meta, overwrite = TRUE)
#
