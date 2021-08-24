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
DB_NAME = "olah_et_al"
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


omxr_pack_anndata <- function (data_in){

  #tools::file_path_sans_ext(data_in)
  if ( class(data_in)[1] == "list" ) {
    # multple files in c("object", "data_mat","obs_meta","var_annot","omics","sample_ID",etc")
    # data_mat - data matrix.  need to assert that this matrix object has omics and sample+ID names

    if (dim(data_in$data_mat)[2] != dim(data_in$var_annot)[1]) {
      X <-Matrix::t(data_in$data_mat)
    } else {
      X <- data_in$data_mat
    }
    if (is.null(dimnames(X)[1])) { rownames(X) <- data_in$sample_ID }
    if (is.null(dimnames(X)[2])) { colnames(X) <- data_in$omics  }

    if(is.null(rownames(obs)) ){
      rownames(obs) <- data_in$sample_ID
    }
    # obs_meta - ensure that we are factored and sample_ID is first column
    obs <- data_in$obs_meta
    id_col <- colnames(obs)[1]
    if (all(obs[[id_col]] == data_in$sample_ID)) {
      obs <- obs %>% dplyr::rename(sample_ID=all_of(id_col))
    } else {
      obs <- obs %>% dplyr::mutate(sample_ID=data_in$sample_ID)
    }

    if(is.null(rownames(obs)) ){
      rownames(obs) <- data_in$sample_ID
    }

    # var_annot - ensure that "omics" is first column
    var_ <- data_in$var_annot
    omics_col <- colnames(var_)[1]

    if ( all(var_[[omics_col]] == data_in$omics) ) {
      var_ <- var_ %>% dplyr::rename(omics_name=all_of(omics_col))
    } else {
      var_ <- var_ %>% dplyr::mutate(omics_name=data_in$omics)
    }
    if(is.null(rownames(var_)) ){
      rownames(var_) <- data_in$omics
    }

    # etc goes into an uns entry
    #
    if ( class(data_in$uns)=="list" ){
      uns <- data_in$uns
    } else {
      uns <- list(etc=data_in$uns)
    }

    ad <- anndata::AnnData(
      X = X,
      obs = obs,
      var = var_,
      uns = uns
    )

  } else if (tolower(tools::file_ext(data_in)) == "rds") {

    data_in = readRDS( file = data_in )


    if(class(data_in)[1] == "Seurat"){
      # how stereotyped is this pattern?  check for Oscar...
      ad <- sceasy::convertFormat(data_in, from="seurat", to="anndata",
                                  outFile = NULL,
                                  assay = 'SCT',
                                  main_layer = 'data',
                                  transfer_layers = c('data', 'counts', 'scale.data')
      )

      raw <- sceasy::convertFormat(data_in, from="seurat", to="anndata",
                                   outFile = NULL,
                                   assay = 'RNA',
                                   main_layer = 'counts',
                                   transfer_layers = NULL)


      ad$raw <- raw


    } else if (class(data_in)[1] == "SingleCellExperiment") {
      print("SingleCellExperiment not enabled")
      ad <- NULL

    } else if ("data.frame" %in% class(data_in)) {
      # could _everything be in a dataframe???
      # yes... lipidomic... strip off first two columns?
      print("enable this for lipidomics? not enabled")
      ad <- NULL

    }

  } else if (tolower(tools::file_ext(data_in)) == "h5ad") {
    ad <- anndata::read_h5ad(data_in)


  } else if (tolower(tools::file_ext(data_in)) == "loom"){
    print("loom loading not enabled")
    ad <- NULL
  }


  return(ad)
}


setup_database <- function(database_name, data_in, db_meta , re_pack=TRUE){
  #LOAD & PACK into ANNDATA
  ##if data_in contains filenames they must be the full path (i.e. RAW_DIR inlcuded)

  DB_NAME <- database_name


  require(reticulate)
  reticulate::use_condaenv(required = TRUE, condaenv = 'omxr')
  require(anndata)

  DB_DIR = file.path("data-raw",DB_NAME)
  if (!dir.exists(DB_DIR)) {
    dir.create(DB_DIR)
  }

  # # db_meta e.g.
  #   organism <- ""
  #   lab <- "Menon"
  #   annotation_database <- "NA"
  # check to see what level of data we were given in our data_list
  #
  if (length(data_in) & names(data_in)=="object") {
    data_in <- unlist(data_in)
  }

  # check if we have it or are forcing
  if ( !file.exists(paste0(DB_DIR,"/core_data.h5ad"))
       | re_pack ) {
    #create it
    # sub-functions to deal with what kind of data we have...
    ad <- omxr_pack_anndata(data_in)

  } else {
    #load it
    ad <- anndata::read_h5ad(file.path(DB_DIR,"core_data.h5ad"))
  }

  return(ad)
}




# ------------------------------------------
# 3. load data -----------------------------
# ------------------------------------------
RAW_DIR <- "ingest/Vilas_C"

SupplementaryData14.csv
SupplementaryData15.csv

cell_barcodes.csv

file_path <- file.path(RAW_DIR,"SupplementaryData14.csv")
counts <- read.csv(file=file_path, header=TRUE, sep=",", row.names=1)
#features <- data.frame("genes" = row.names(count_table))
features <- row.names(counts)
countm <- as.matrix(counts)
countm <- as(countm, "dgTMatrix")
#counts <- t(countm)  #transpose so samples are rows (as in ANNDATA format)

# data <- read.csv(file=file_path, header=TRUE, sep=",", row.names=1)
# data1 <- as.matrix(data) #2.2GB
# data2 <- as(data1, "dgeMatrix") #4.4GB
# data3 <- as(data1, "dgTMatrix") #0.29GB


# 2. load annotation data --------------------
file_path <- file.path(RAW_DIR,"SupplementaryData15.csv")
#annots <- read_csv(file_path)
obs_annots <- read.csv(file=file_path, header=TRUE, sep=",", row.names=NULL)
#obs_annots$obs_names <- obs$sample_id
obs_annots$cluster_name = paste0("Cluster_",obs_annots$cluster_label)
obs_annots$cluster_id = obs_annots$cluster_label
row.names(obs_annots) <- obs_annots$sample_id

file_path <- file.path(RAW_DIR,"cell_barcodes.csv")
barcodes <- read.csv(file=file_path, header=TRUE, sep=",", row.names=NULL)

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

raw <- ad$copy()

sc <- import("scanpy")
#
normpc_adata <- sc$pp$normalize_per_cell(ad,copy=TRUE)
norm_adata <- sc$pp$normalize_total(ad,copy=TRUE)

sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)

#z_adata <- sc$pp$scale(ad,copy=TRUE)
# log_adata <- sc$pp$log1p(ad,copy=TRUE)
#
# sc$pp$normalize_total(ad)
# sc$pp$log1p(ad)
# sc$pp$recipe_seurat(ad, log=TRUE)
# recipe_zheng17

sc$pp$filter_genes(ad, min_counts=1)         # only consider genes with more than 1 count
sc$pp$normalize_per_cell(                       # normalize with total UMI count per cell
  ad, key_n_counts='n_counts_all'
)
# filter_result = sc$pp$filter_genes_dispersion(  # select highly-variable genes
#   ad$X, flavor='cell_ranger', n_top_genes=3000, log=FALSE
# )
#adata = adata[:, filter_result.gene_subset]     # subset the genes
sc$pp$normalize_per_cell(ad)                 # renormalize after filtering
sc$pp$log1p(ad)                      # log transform: adata.X = log(adata.X + 1)
sc$pp$scale(ad)                              # scale to unit variance and shift to zero mean

#ad$write_csvs(file.path("data-raw",DB_NAME,"whoa.csv"))


test_types <- c('wilcoxon','t-test_overestim_var')
comp_types <- c("allVrest")


helper_function<-('data-raw/compute_de_table.R')
source(helper_function)

diff_exp <- compute_de_table(ad,comp_types, test_types, obs_name = c('cluster_name'))

ad$raw <- raw

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

