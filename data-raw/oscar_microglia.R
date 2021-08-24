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
reticulate::use_condaenv(required = TRUE, condaenv = 'omxr')
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

countm <- as(counts, "dgTMatrix")

# ------------------------------------------
# 4. pack into anndata                    --
# ------------------------------------------
X <- Matrix::t(counts)  #transpose so rows and colums correspond to anndata expectation
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
# for some reason we lost the var_names and col_names
ad$var_names <- colnames(X)
ad$obs_names <- rownames(X)



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


# #normpc_adata <- sc$pp$normalize_per_cell(ad,copy=TRUE)
# #norm_adata <- sc$pp$normalize_total(ad,copy=TRUE)
#z_adata <- sc$pp$scale(ad,copy=TRUE)
# log_adata <- sc$pp$log1p(ad,copy=TRUE)

sc$pp$recipe_seurat(ad)
#sc$pp$log1p(ad)

sc$pp$normalize_total(ad)


test_types <- c('wilcoxon','t-test_overestim_var')
comp_types <- c("grpVrest")


helper_function<-('data-raw/compute_de_table.R')
source(helper_function)

diff_exp <- compute_de_table(ad,comp_types, test_types, obs_name = c('Status','cluster_label'))


#
# diff_exp <- data.frame()
# #log_adata <- sc$pp$log1p(ad,copy=TRUE)
# condition_key = 'cluster_label'
# comp_type <- "allVrest"
# test_type <- test_types[1]
# key <- paste0(test_type,"_", comp_type)
# sc$tl$rank_genes_groups(ad, condition_key, method=test_type, key_added = key)
# de_table <- sc$get$rank_genes_groups_df(ad, NULL, key=key)
# de_table$condition <- comp_type
# de_table$test_type <- test_type
# diff_exp <- dplyr::bind_rows(diff_exp, de_table)# for this dataset its straightforward to do all comparisons...
#
#
# test_type <- test_types[2]
# key <- paste0(test_type,"_", comp_type)
# sc$tl$rank_genes_groups(ad, condition_key, method=test_type, key_added = key)
# de_table <- sc$get$rank_genes_groups_df(ad, NULL, key=key)
# de_table$condition <- comp_type
# de_table$test_type <- test_type
# diff_exp <- dplyr::bind_rows(diff_exp, de_table)# for this dataset its straightforward to do all comparisons...
#
#


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

DB_NAME = "oscar_microglia"
DB_DIR = file.path("data-raw",DB_NAME)
#RAW_DIR <- "ingest/Vilas_B"

ad <- read_h5ad(file.path(DB_DIR,"core_data_plus_de.h5ad"))





# ------------------------------------------
# 7 . create config and default files                   --
# ------------------------------------------
require(anndata)
require(reticulate)

DB_NAME = "oscar_microglia"
DB_DIR = file.path("data-raw",DB_NAME)
#RAW_DIR <- "ingest/Vilas_B"

ad <- read_h5ad(file.path(DB_DIR,"core_data_plus_de.h5ad"))
db_prefix = "de"
diff_exp = readRDS( file = paste0(DB_DIR, "/", db_prefix, "_table.rds"))



# OBSERVABLES
#
measures <- list(obs = c("n_genes","nCount_RNA", "nFeature_RNA"," percent.mito","n_genes"),
                 var = ad$var_keys()[2:4])

# differentials  #if we care we need to explicitly state. defaults will be the order...
diffs <- list(diff_exp_groups =  levels(factor(diff_exp$group)),
              diff_exp_comp_type =  levels(factor(diff_exp$comp_type)),
              diff_exp_obs_name =  levels(factor(diff_exp$obs_name)),
              diff_exp_tests =  levels(factor(diff_exp$test_type)))

default_factors <- c("Status","cluster_label","Gender")


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


oscar_microglia_conf = readRDS( paste0(DB_DIR,"/",db_prefix,"_conf.rds") )
oscar_microglia_def = readRDS( paste0(DB_DIR,"/",db_prefix,"_def.rds") )
oscar_microglia_omics = readRDS( paste0(DB_DIR,"/",db_prefix,"_omics.rds") )
oscar_microglia_meta = readRDS( paste0(DB_DIR,"/",db_prefix,"_meta.rds") )


#usethis::use_data(oscar_microglia_conf, oscar_microglia_def, oscar_microglia_omics, oscar_microglia_meta, overwrite = TRUE)

