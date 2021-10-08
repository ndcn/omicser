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


DB_NAME <-  list("Vilas Microglia (sceasy)" = "vilas_microglia_sceasy")



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


write_db_meta(db_meta,DB_NAME, db_root = DB_ROOT_PATH)



#==== 2. helper functions =========================================================================
# "helper" functions to prep your "raw" data go here or should be sourced here...

prep_vilas_mg_files <- function(data_file,path_root){

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

  # WIP  -------------------------
  exp_group <- raw_table$Group
  exp_group[exp_group==""] <- "CONTROL"  #i'm not sure this is true

  meta_cols <- c("Name","Group")
  obs_meta <- raw_table[,c("Name","Group")]
  obs_meta$Group <- exp_group #recode
  rownames(obs_meta) <- raw_table$Name


  dat_mat <- raw_table[,!c("Name","Group")]
  dat_mat <- as.matrix( dat_mat)
  rownames(dat_mat) <- obs_meta$Name

  class(dat_mat) <- "numeric"


  # make var_annots
  lipids=colnames(dat_mat)
  var_annot <- as.data.frame(lipids)
  rownames(var_annot) <- var_annot$lipids

  raw <- dat_mat
  # zero out NA: is this tricky because its sparse?
  dat_mat[which(is.na(dat_mat), arr.ind = TRUE)] <- 0

  ####  marginals on un-scaled data with zeroed NA
  var_annot$mean <- colMeans(raw,na.rm = TRUE)
  var_annot$var <- matrixStats::colVars(raw,na.rm = TRUE)

  obs_meta$var <- matrixStats::rowVars(raw,na.rm = TRUE)
  obs_meta$mean <-rowMeans(raw,na.rm = TRUE)
  # experiment with the scaled (including)

  excess_zero_conc <- ( colSums(dat_mat==0) > 2/3*dim(dat_mat)[1] )
  var_annot$excess_zero_conc <- excess_zero_conc
  # var annotation
  #poor_conc <- which(colSums(conc_mat==0) > 2/3*dim(conc_mat)[1])


  data_list <- list(data_mat = dat_mat,
                    obs_meta = obs_meta,
                    var_annot = var_annot,
                    omics = rownames(var_annot),
                    sample_ID = rownames(obs_meta),
                    etc = NULL,
                    raw = raw)

  return(data_list)
}






#==== 3. load data =========================================================================
RAW_DIR <- "ingest/Vilas_A"


ad <- setup_database(database_name=DB_NAME, data_in=data_list, db_meta=db_meta , re_pack=TRUE)
ad$write_h5ad(filename=file.path("data-raw",DB_NAME,"core_data.h5ad"))


# TODO:::  fix this here:::

#==== 4. pack into anndata =========================================================================

conc_csv_name <- "example_species_conc.csv"
conc_dat_list <- prep_lipidizer_files(conc_csv_name,RAW_DIR)
# b. load composition data --------------------

conc <- omicser::setup_database(database_name=DB_NAME,
                                db_path=DB_ROOT_PATH,
                                data_in=conc_dat_list,
                                db_meta=NULL ,
                                re_pack=TRUE)



#==== 5. post processing =========================================================================
sc <- reticulate::import("scanpy")
# re-read the anndata file to instantiate the anndata in R wrappers

raw <- ad$raw$copy()

#outvals <- sc$pp$filter_cells(ad, min_genes=200, inplace=FALSE)

sc$pp$filter_cells(ad, min_genes=200)
sc$pp$filter_genes(ad, min_cells=20)


rawX <- raw$X

# FIX RAW
colnames(rawX) <- raw$var_names

rawX <- rawX[which(raw$obs_names %in% ad$obs_names),(raw$var_names %in% ad$var_names)]
raw$obs_names <- raw$obs_names[which(raw$obs_names %in% ad$obs_names)]



newraw <- anndata::AnnData(
  X=rawX,
  var = data.frame(names=colnames(rawX)),
  obs = data.frame(sampleID=rownames(ad$X))

)
ad$raw <- newraw

# to conform with cellXgene scheme gene filtering is tricky.
#
# we need to retain the list / shape of the original genes so we must
# replace zeros/nulls with NA

# outvals <- sc$pp$filter_genes(ad, min_cells=2, inplace=FALSE)

sc$pp$normalize_per_cell(ad, counts_per_cell_after=1e4)
sc$pp$log1p(ad)  # already logg-ed?
sc$pp$highly_variable_genes(ad, min_mean=0.0125, max_mean=3, min_disp=0.5  )
#sc$pp$regress_out(ad_tmp, 'n_counts' )

#sc$pp$highly_variable_genes(ad,n_top_genes=40)


#  don't know how to make this work....
#sc$pp$highly_variable_genes(ad,n_top_genes=40)
ad$var$var_rank <- order(ad$var$sct.residual_variance)
# choose top 40 genes by variance across dataset as our "targets"
target_omics <- ad$var_names[which(ad$var$var_rank <= 40)]


ad$var$decile <- dplyr::ntile(ad$var$sct.residual_variance, 10)



if (FALSE){ # save an intermediate file (incase we want to revert...)
  ad$write_h5ad(filename=file.path(DB_ROOT_PATH,DB_NAME,"normalized_data.h5ad"))
}
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

if (FALSE){ # save an intermediate file (incase we want to revert...)
  ad$write_h5ad(filename=file.path(DB_ROOT_PATH,DB_NAME,"norm_data_plus_dr.h5ad"))
}
#==== 6. differential expression =========================================================================

sc <- import("scanpy")

test_types <- c('wilcoxon','t-test_overestim_var')


comp_types <- c("grpVrest")
obs_names <- c('disease','cell_type')
diff_exp <- omicser::compute_de_table(ad,comp_types, test_types, obs_names,sc)


comp_types <- c("Mild_Cognitive_ImpairmentVAlzheimer_Disease",
                "Mild_Cognitive_ImpairmentVTemporal_Lobe_Epilepsy",
                "Alzheimer_DiseaseVTemporal_Lobe_Epilepsy")
obs_names <- c('disease')
diff_exp2 <- omicser::compute_de_table(ad,comp_types, test_types, obs_names,sc)

diff_exp <- dplyr::bind_rows(diff_exp, diff_exp2)
#diff_exp <- rbind(diff_exp,diff_exp2)
if (FALSE){ # save an intermediate file (incase we want to revert...)
 ad$write_h5ad(filename=file.path(DB_ROOT_PATH,DB_NAME,"norm_data_with_de.h5ad"))
}
saveRDS(diff_exp, file = file.path(DB_ROOT_PATH,DB_NAME, "db_de_table.rds"))

#==== 7. create configs =========================================================================
DB_NAME = "vilas_microglia_sceasy"
if (FALSE) { # load intermediate files if available
  ad <- read_h5ad(filename=file.path("data-raw",DB_NAME,"normalized_data.h5ad"))
  ad <- read_h5ad(filename=file.path("data-raw",DB_NAME,"norm_data_with_de.h5ad"))
  diff_exp <- readRDS( file = file.path("data-raw",DB_NAME, "db_de_table.rds"))
}

# differentials  #if we care we need to explicitly state. defaults will be the order...
config_list <- list(
  x_obs = c("tissue", "disease", "cell_type", "sex", "leiden"),
  y_obs = c("nCount_RNA","nFeature_RNA","nCount_SCT","nFeature_SCT", "n_genes","n_counts"), #MEASURES
  obs_groupby = c("tissue", "disease", "cell_type", "sex", "leiden"),
  obs_subset = c("tissue", "disease", "cell_type", "sex", "leiden"),

  x_var = c("highly_variable","decile"),
  y_var = c("sct.detection_rate", "sct.gmean", "sct.variance","sct.residual_mean","sct.residual_variance", "sct.variable",
            "n_cells",  "means" , "dispersions", "dispersions_norm" ),
  var_groupby = c("highly_variable","decile"),
  var_subset = c("highly_variable","decile"),  # NOTE:  <omic selector> is NOT in the data object so its not actually going to load

  layers = c("X","raw","counts"),

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
  default_factors = c("cell_type","disease","tissue"),
  target_omics = target_omics,
  omic_details = c("highly_variable",
                   "decile",
                   "sct.detection_rate",
                   "sct.gmean",
                   "sct.variance",
                   "sct.residual_variance",
                   "sct.variable",
                   "n_cells")

)

omicser::write_db_conf(config_list,DB_NAME, db_root = DB_ROOT_PATH)

#==== 8. write data file to load  =========================================================================

ad$write_h5ad(filename=file.path(DB_ROOT_PATH,DB_NAME,"db_data.h5ad"))



#==== 9. create meta-data narrative for rendering in INGEST  =========================================================================


# ALSO create an `additional_info.Rmd` as in `inst/app/www/additional_info.Rmd`



