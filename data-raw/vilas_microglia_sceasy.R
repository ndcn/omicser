#######################################################################
#######################################################################
##
##  Microglia scRNAseq Transcriptomics
##  - from Vilas:
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
##
#######################################################################
#######################################################################
#  formely `vilas_B` which was derived from a Seurat file
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
require("Seurat")

# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install(c("LoomExperiment", "SingleCellExperiment"))
#devtools::install_github("cellgeni/sceasy")
require("sceasy")

require(reticulate)
reticulate::use_condaenv(required = TRUE, condaenv = 'sc39')
require(anndata)

# create the folder to contain the raw data
DB_NAME = "vilas_microglia_sceasy"
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
RAW_DIR <- "ingest/Vilas_B"
# Seurat Data object
data_file <- "microglia_data_with_counts_RNA_SCT.rda"
load(file.path(RAW_DIR,data_file) ) # microglia_data

# 3a. update seurat object -----------------------------
microglia_data_updated <- UpdateSeuratObject(object = microglia_data)



#> OUTPUT:
#   Validating object structure
#   Updating object slots
#   Ensuring keys are in the proper strucutre
#   Ensuring feature names don't have underscores or pipes
#   Object representation is consistent with the most current Seurat version

saveRDS(microglia_data_updated,
        file = file.path(RAW_DIR, "microglia_data_seu_new.rds"))



loadRDS(file = file.path(RAW_DIR, "microglia_data_seu_new.rds"))




ad <- sceasy::convertFormat(microglia_data_updated, from="seurat", to="anndata",
                      outFile = NULL,
                      assay = 'SCT',
                      main_layer = 'data',
                      transfer_layers = c('data', 'counts', 'scale.data')
                      )

raw <- sceasy::convertFormat(microglia_data_updated, from="seurat", to="anndata",
                            outFile = NULL,
                            assay = 'RNA',
                            main_layer = 'counts',
                            transfer_layers = NULL)


ad$raw <- raw

ad$write_h5ad(filename=file.path(DB_DIR,"core_data_plus_de.h5ad"))


outFile=file.path(DB_DIR,"core_data.h5ad")

sceasy::convertFormat(microglia_data_updated, from="seurat", to="anndata",
                      outFile=file.path(DB_DIR,"core_data.h5ad"))



sceasy::convertFormat(microglia_data_updated, from="seurat", to="anndata",
                      outFile=file.path(DB_DIR,"core_data.h5ad"))

ad <- read_h5ad(file.path(DB_DIR,"core_data.h5ad"))



# ------------------------------------------
# 4. pack into anndata   (redundant?)      --
# ------------------------------------------0


var_ <- ad$var
obs <- ad$obs
X <- ad$X
ad_raw <- ad$raw
obsm <- ad$obsm
varm <- ad$varm

db_prefix = "core_data"
saveRDS(X, file = paste0(DB_DIR, "/", db_prefix, "_X.rds"))
saveRDS(obs, file = paste0(DB_DIR, "/", db_prefix, "_obs.rds"))
saveRDS(var_, file = paste0(DB_DIR, "/", db_prefix, "_var.rds"))

layers <- ad$layers$get('counts')

ad_ <- AnnData(
  X = X,
  obs = obs,
  var = var_,
  raw = ad_raw,
  obsm = obsm,
  varm=varm,
  layers =  list(counts=ad$layers$get('counts')) #list('count'=layers)
)

ad_

ad_$write_h5ad(filename=file.path(DB_DIR,"recast.h5ad"))



#ad$write_h5ad(filename=file.path(DB_DIR,"core_data.h5ad"))

# ------------------------------------------
# 5. post processing                      --
# ------------------------------------------
require(anndata)
require(reticulate)
DB_NAME = "vilas_microglia_seu"
DB_DIR = file.path("data-raw",DB_NAME)
RAW_DIR <- "ingest/Vilas_B"

db_prefix = "core_data"
X = readRDS( file = paste0(DB_DIR, "/", db_prefix, "_X.rds"))
obs = readRDS(  file = paste0(DB_DIR, "/", db_prefix, "_obs.rds"))
var_ = readRDS( file = paste0(DB_DIR, "/", db_prefix, "_var.rds"))


ad <- read_h5ad(file.path(DB_DIR,"core_data.h5ad"))
ad


sc <- import("scanpy")


test_types <- c('wilcoxon','t-test_overestim_var')
comp_types <- c("allVrest")


helper_function<-('data-raw/compute_de_table.R')
source(helper_function)

diff_exp <- compute_de_table(ad,comp_types, test_types)


# put the logvals in layers of ad
# copy the whole thing and replace X to copy the uns to ad

ad$write_h5ad(filename=file.path(DB_DIR,"core_data_plus_de.h5ad"))

# also need to pack the diff_exp1 and diff_exp2 into easy to deal wiht tables for volcanos...


db_prefix = "de"
saveRDS(diff_exp, file = paste0(DB_DIR, "/", db_prefix, "_table.rds"))


# ------------------------------------------
# 6. dimension reduction - PCA / umap    --
# ------------------------------------------
require(anndata)
require(reticulate)

DB_NAME = "vilas_microglia_seu"
DB_DIR = file.path("data-raw",DB_NAME)
#RAW_DIR <- "ingest/Vilas_B"

ad <- read_h5ad(file.path(DB_DIR,"core_data_plus_de.h5ad"))




# ------------------------------------------
# 7 . create config and default files                   --
# ------------------------------------------
require(anndata)
require(reticulate)

DB_NAME = "vilas_microglia_seu"
DB_DIR = file.path("data-raw",DB_NAME)
#RAW_DIR <- "ingest/Vilas_B"

ad <- read_h5ad(file.path(DB_DIR,"core_data_plus_de.h5ad"))
db_prefix = "de"
diff_exp = readRDS( file = paste0(DB_DIR, "/", db_prefix, "_table.rds"))




###

# measures
#  This ordering is the "default"
measures <- list(obs = c("nCount_SCT","nFeature_SCT","nCount_RNA","nFeature_RNA"),
                 var = ad$var_keys())
# [1] "sct.detection_rate"    "sct.gmean"             "sct.variance"
# [4] "sct.residual_mean"     "sct.residual_variance" "sct.variable"

# differentials  #if we care we need to explicitly state. defaults will be the order...
diffs <- list(diff_exp_groups =  levels(factor(diff_exp$group)),
              diff_exp_comp_type =  levels(factor(diff_exp$comp_type)),
              diff_exp_obs_name =  levels(factor(diff_exp$obs_name)),
              diff_exp_tests =  levels(factor(diff_exp$test_type)))

# Dimred
dimreds <- list(obsm = c('X_pca', 'X_tsne'),
                    varm = c('PCs'))

# what ad$obs do we want to make default values for...
# # should just pack according to UI?
default_factors <- c("tissue","disease","cell_type")




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


vilas_microglia_seu_conf = readRDS( paste0(DB_DIR,"/",db_prefix,"_conf.rds") )
vilas_microglia_seu_def = readRDS( paste0(DB_DIR,"/",db_prefix,"_def.rds") )
vilas_microglia_seu_omics = readRDS( paste0(DB_DIR,"/",db_prefix,"_omics.rds") )
vilas_microglia_seu_meta = readRDS( paste0(DB_DIR,"/",db_prefix,"_meta.rds") )


usethis::use_data(vilas_microglia_seu_conf, vilas_microglia_seu_def, vilas_microglia_seu_omics, vilas_microglia_seu_meta, overwrite = TRUE)


