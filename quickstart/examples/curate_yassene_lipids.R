#==== INFO ========================================================
## code to prepare `yassene_example` dataset goes here
##
##  "WHAT IS THE CONTEXT?" lipidomics composition data
##  - data from yassenne / ricoo
##
#  formely `yassene_A` (conc and compos)

#==== 0. preamble/setup ==================================================

DEV_OMICSER <- TRUE


if (DEV_OMICSER){
  # this should be a full path... e.g. ~/Projects/NDCN_dev/omicser
  # but for github, we will set relative to the repo BASE
  REPO_PATH <- "/Users/ahenrie/Projects/NDCN_dev/omicser"
  OMICSER_RUN_DIR <- file.path(REPO_PATH,"quickstart")
  golem::document_and_reload(pkg = REPO_PATH)
} else {

  require(omicser)
  OMICSER_RUN_DIR <- file.path(REPO_PATH,"quickstart")

}

# BOOTSTRAP the options we have already set up...
# NOTE: we are looking in the "quickstart" folder.  the default is to look for the config in with default getwd()
omicser_options <- omicser::get_config(in_path = OMICSER_RUN_DIR)


CONDA_ENV <- omicser_options$conda_environment
DB_ROOT_PATH <- omicser_options$db_root_path

DB_NAME <-  list("Yassene Lipid concentraions & compositions" ="yassene_lipid")

if (! (DB_NAME %in% omicser_options$database_names)){
  omicser_options$database_names <- c(omicser_options$database_names,DB_NAME)
  omicser::write_config(omicser_options,in_path = OMICSER_RUN_DIR )
}


#DB_NAME = "yassene_lipid"
DB_DIR = file.path(DB_ROOT_PATH,DB_NAME)
if (!dir.exists(DB_DIR)) {
  dir.create(DB_DIR)
}


#==== 1. documentation / provenance ==============================================================
# TODO:  markdown file or html with some copy about the database
#  - lab, paper link/name
#  summarize results / data origin whatever
db_meta <- list(
  organism = 'mmusculus',
  lab = "Giera",
  annotation_database =  "TBD",
  title = "lipidizer lipidomics",
  omic_type = "Lipidomics",
  measurment = "concentration and composition",
  pub = "TBD",
  date = format(Sys.time(), "%a %b %d %X %Y")
)

write_db_meta(db_meta,DB_NAME, db_root = DB_ROOT_PATH)


#==== 2. helper functions =================================================================================


prep_lipidizer_files <- function(data_file,path_root){

  raw_table <- data.table::fread(file=file.path(path_root,data_file), header=TRUE)

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



# TODO:  define "get marker genes" or some such


#==== 3. load data -========================================================================================

RAW_DIR <- file.path(OMICSER_RUN_DIR,"raw_data", "Yassene_A")


# a. load concentration data --------------------

conc_csv_name <- "example_species_conc.csv"
conc_dat_list <- prep_lipidizer_files(conc_csv_name,RAW_DIR)
# b. load composition data --------------------
#
comp_csv_name <- "example_composition.csv"
comp_dat_list <- prep_lipidizer_files(comp_csv_name,RAW_DIR)


#==== 4. pack into anndata =========================================================================

# helper_function<-('data-raw/ingest_helpers.R')
# source(helper_function)



DB_NAME = "yassene_lipid"

conc <- omicser::setup_database(database_name=DB_NAME,
                                db_path=DB_ROOT_PATH,
                                data_in=conc_dat_list,
                                db_meta=NULL ,
                                re_pack=TRUE)

comp <- omicser::setup_database(database_name=DB_NAME,
                                db_path=DB_ROOT_PATH,
                                data_in=comp_dat_list,
                                db_meta=NULL ,
                                re_pack=TRUE)

ad<- conc$copy()
raw <- ad$copy()
raw$X <- conc_dat_list$raw
ad$raw <- raw
ad$layers <- list(composition=comp$X,
                  raw_comp=comp_dat_list$raw)


ad$write_h5ad(filename=file.path(DB_ROOT_PATH,DB_NAME,"core_data.h5ad"))



#==== 5. post processing =========================================================================               --
# use scanpy to do some scaling and calculations...
sc <- reticulate::import("scanpy")


# 5a. lasso regression to choose sig ----------------------------
#ZSCALE the matrix
dat_mat <- scale(ad$X) #including zeros
test_vals <- dat_mat
regr_group <- ad$obs$Group
g <- unique(regr_group)
# g
#
regr_group <- as.numeric(factor(ad$obs$Group))
ind_rem_group <- which(regr_group == which(table(regr_group) < 3))

# remove groups with less than 3 observations..
if (length(ind_rem_group) > 0) {
  test_vals <- test_vals[-ind_rem_group, ]
  regr_group <- regr_group[-ind_rem_group]
}

# choose the "significant" columns via lasso regression (glmnet)
set.seed(100)
cvfit <- glmnet::cv.glmnet(test_vals, regr_group, nlambda = 100, alpha = .8, family = "multinomial", type.multinomial = "grouped")
coef <- coef(cvfit, s = "lambda.min")
tmp <- as.matrix(coef$"1")
tmp1 <- tmp[which(tmp != 0)]
coef_names <- rownames(tmp)[which(tmp != 0)][-1]
ind_coef <- which(colnames(test_vals) %in% coef_names)

ad$var$sig_lasso_coef <- (colnames(test_vals) %in% coef_names)
ad$uns <- list(comp_lasso_coef=coef)


#  don't know how to make this work....
#sc$pp$highly_variable_genes(ad,n_top_genes=40)
ad$var$var_rank <- order(ad$var$var)
# choose top 40 proteins by variance across dataset as our "targets"
target_omics <- ad$var_names[which(ad$var$var_rank <= 40)]

ad$var$decile <- dplyr::ntile(ad$var$var, 10)


# save an intermediate file (incase we want to revert...)
ad$write_h5ad(filename=file.path(DB_ROOT_PATH,DB_NAME,"normalized_data.h5ad"))

#==== 5-a. dimension reduction - PCA / umap =========================================================================

## Step 3: Do some basic preprocessing to run PCA and compute the neighbor graphs
##

comp <- ad$copy()
comp$X <- ad$layers$get('composition')
sc$pp$pca(comp)
sc$pp$neighbors(comp)
## Step 4: Infer clusters with the leiden algorithm
sc$tl$leiden(comp)
## Step 5: Compute tsne and umap embeddings
sc$tl$umap(comp)


sc$pp$pca(ad)
sc$pp$neighbors(ad)
sc$tl$leiden(ad)
sc$tl$umap(ad)


ad$obsm$comp_X_pca <- comp$obsm$X_pca
ad$obsm$comp_X_umap <- comp$obsm$X_umap
ad$varm$comp_PCs <- comp$varm$PCs
ad$obsp$comp_distances <- comp$obsp$distances
ad$obsp$comp_connectivities <- comp$obsp$connectivities
# save an intermediate file (incase we want to revert...)
ad$write_h5ad(filename=file.path(DB_ROOT_PATH,DB_NAME,"norm_data_plus_dr.h5ad"))


#==== 6. differential expression  ======================================================================

test_types <- c('wilcoxon','t-test_overestim_var')

comp_types <- c("grpVrest")
obs_names <- c('Group','leiden')

diff_exp <- omicser::compute_de_table(ad,comp_types, test_types, obs_names,sc)


ad$write_h5ad(filename=file.path(DB_ROOT_PATH,DB_NAME,"norm_data_with_de.h5ad"))
saveRDS(diff_exp, file = file.path(DB_ROOT_PATH,DB_NAME, "db_de_table.rds"))

#==== 7. create configs =========================================================================
# what ad$obs do we want to make default values for...

config_list <- list(
  x_obs = c('Group','leiden' ),
  y_obs = c('var', 'mean'), #MEASURES
  obs_groupby = c('Group','leiden'),
  obs_subset = c('Group','leiden' ),
  x_var = c('excess_zero_conc', 'sig_lasso_coef',"decile"),
  y_var = c('mean', 'var'),
  var_groupby = c("highly_variable","decile"),
  var_subset = c("highly_variable","decile"),  # NOTE:  <omic selector> is NOT in the data object so its not actually going to load

  layers = c("X","raw",'composition','raw_comp'),

  diffs = list(diff_exp_comps = levels(factor(diff_exp$versus)),
               diff_exp_comp_type =  levels(factor(diff_exp$comp_type)), #i don't think we need this
               diff_exp_obs_name =  levels(factor(diff_exp$obs_name)),
               diff_exp_tests =  levels(factor(diff_exp$test_type))
  ),


  # Dimred
  dimreds = list(obsm = c('X_pca', 'X_umap', 'comp_X_pca', 'comp_X_umap'),
                 varm = c('PCs', 'comp_PCs')),

  # what ad$obs do we want to make default values for...
  # # should just pack according to UI?
  default_factors = c('Group', 'leiden', 'Group'),
  target_omics = target_omics,
  omic_details = c('omics_name', 'mean', 'var', 'excess_zero_conc', 'sig_lasso_coef', 'var_rank', 'decile')

)



omicser::write_db_conf(config_list,DB_NAME, db_root = DB_ROOT_PATH)

#==== 8. write data file to load  =========================================================================
ad$write_h5ad(filename=file.path(DB_ROOT_PATH,DB_NAME,"db_data.h5ad"))

