## code to prepare `yassene_example` dataset goes here
#######################################################################
#######################################################################
##
##  "WHAT IS THE CONTEXT?" lipidomics composition data
##  - data from yassenne / ricoo
##
#######################################################################
#######################################################################
#  formely `yassene_A` (conc and compos)
#


# ------------------------------------------
# 0. preamble/setup -------------------------
# ------------------------------------------
require(tidyverse)
require(reticulate)
reticulate::use_condaenv(required = TRUE, condaenv = 'omxr')
require(anndata)

require(glmnet)
require(readxl)
require(matrixStats)


# create the folder to contain the raw data
DB_NAME = "yassene_A_compos"
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
lab <- ""
annotation_database <- ""



# ------------------------------------------
# 2. helper functions ----------------------
# ------------------------------------------



# ------------------------------------------
# 3. load data -----------------------------
# ------------------------------------------
RAW_DIR <- "ingest/Yassene_A"



# 1. load composition data --------------------
csv_name <- "example_composition.csv"

file_name <- file.path(RAW_DIR,csv_name)
lipid_composition <- read.csv(file=file_name, header=TRUE, as.is = TRUE)

meta_cols <- c("Name","Group")
obs_meta <- lipid_composition[meta_cols]
comp_group <- lipid_composition$Group
comp_group[comp_group==""] <- "CONTROL"  #i'm not sure this is true
obs_meta$Group <- comp_group

rownames(obs_meta) <- lipid_composition$Name

# # convert to matrix
# comp_mat <- as.matrix( lipid_composition[ , !(names(lipid_composition) %in% meta_cols)])

comp_data <- lipid_composition %>% dplyr::select(!matches(meta_cols))
rownames(comp_data) <- obs_meta$Name


comp_mat <- as.matrix(comp_data)
rownames(comp_mat) <- obs_meta$Name

class(comp_mat) <- "numeric"  # converts to numbers / NAs
#ZSCALE the matrix
zcomp_raw <- scale(comp_mat)


comp_mat[which(is.na(comp_mat), arr.ind = TRUE)] <- 0
# FIND columns wiht >2/3 zeros for the "REMOVE ZEROS" MATRIX
excess_zero_comp <- ( colSums(comp_mat==0) > 2/3*dim(comp_mat)[1] )
#poor_comp <- which(colSums(comp_mat==0) > 2/3*dim(comp_mat)[1])

# ZSCALE the zeroed matrix
zcomp <- scale(comp_mat)

mu_comp <- colMeans(comp_mat)
var_comp <- colVars(comp_mat)

# obs_meta$mean <- rowMeans(comp_mat)
# obs_meta$var <- apply(comp_mat,1,var)

obs_meta$var_comp <- rowVars(comp_mat)
obs_meta$mean_comp <-rowMeans(comp_mat)

# experiment with the scaled (including)
regr_group <- obs$Group
g <- unique(regr_group)
# g
#
test_comp <- zcomp
regr_group <- as.numeric(factor(obs_meta$Group))

ind_rem_group <- which(regr_group == which(table(regr_group) < 3))

# remove groups with less than 3 observations..
if (length(ind_rem_group) > 0) {
  test_comp <- test_comp[-ind_rem_group, ]
  regr_group <- regr_group[-ind_rem_group]
}

dim(test_comp)

lipids=colnames(test_comp)
var_ <- as.data.frame(lipids)


# TODO: make these generic so its the same for -compos and -conc databases
var_$mean_comp <- mu_comp
var_$var_comp <- var_comp
var_$excess_zero_comp <- excess_zero_comp

# choose the "significant" columns via lasso regression (glmnet)
set.seed(100)
cvfit <- cv.glmnet(test_comp, regr_group, nlambda = 100, alpha = .8, family = "multinomial", type.multinomial = "grouped")
coef <- coef(cvfit, s = "lambda.min")
tmp <- as.matrix(coef$"1")
tmp1 <- tmp[which(tmp != 0)]
coef_names <- rownames(tmp)[which(tmp != 0)][-1]
ind_coef <- which(colnames(test_comp) %in% coef_names)

uns <- list(comp_lasso_coef=coef)

var_$comp_sig_lasso_coef <- (colnames(test_comp) %in% coef_names)


rownames(var_) <- var_$lipids



# ------------------------------------------
# 4. pack into anndata                    --
# ------------------------------------------

#

X <- comp_mat
obs <- obs_meta

db_prefix = "core_data"
saveRDS(X, file = paste0(DB_DIR, "/", db_prefix, "_X.rds"))
saveRDS(obs, file = paste0(DB_DIR, "/", db_prefix, "_obs.rds"))
saveRDS(var_, file = paste0(DB_DIR, "/", db_prefix, "_var.rds"))
saveRDS(uns, file = paste0(DB_DIR, "/", db_prefix, "_uns.rds"))


ad <- AnnData(
  X = X,
  obs = obs,
  var = var_,
  uns = uns
)

ad




ad$write_h5ad(filename=file.path(DB_DIR,"core_data.h5ad"))





# ------------------------------------------
# 5. post processing                      --
# ------------------------------------------
# use scanpy to compute:
#     - differential measures (volcano)
#     - marginal stats
#     - dimension reductions (PCA,UMAP, tSNE)


#> AnnData object with n_obs × n_vars = 2 × 3
#>     obs: 'group'

require(anndata)
require(reticulate)
DB_NAME = "yassene_A_compos"
DB_DIR = file.path("data-raw",DB_NAME)
RAW_DIR <- "ingest/Yassene_A"

db_prefix = "core_data"
X = readRDS( file = paste0(DB_DIR, "/", db_prefix, "_X.rds"))
obs = readRDS(  file = paste0(DB_DIR, "/", db_prefix, "_obs.rds"))
var_ = readRDS( file = paste0(DB_DIR, "/", db_prefix, "_var.rds"))




ad <- read_h5ad(file.path(DB_DIR,"core_data.h5ad"))
ad



sc <- import("scanpy")


# norm_adata <- sc$pp$normalize_total(ad,copy=TRUE)
# #normpc_adata <- sc$pp$normalize_per_cell(ad,copy=TRUE)
# z_adata <- sc$pp$scale(ad,copy=TRUE) #this does the same colmean as above
# log1p_adata <- sc$pp$log1p(ad,copy=TRUE)





test_types <- c('wilcoxon','t-test_overestim_var')
comp_types <- c("grpVrest")


helper_function<-('data-raw/compute_de_table.R')
source(helper_function)

diff_exp <- compute_de_table(ad,comp_types, test_types, obs_names = 'Group')

#
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

# not sure how to add layers...

# TODO:  add de_table into uns
#new_adata$uns$de_table <- diff_exp
ad$write_h5ad(filename=file.path(DB_DIR,"core_data_plus_de.h5ad"))

db_prefix = "de"
saveRDS(diff_exp, file = paste0(DB_DIR, "/", db_prefix, "_table.rds"))



# ------------------------------------------
# 6. dimension reduction - PCA / umap    --
# ------------------------------------------
require(anndata)
require(reticulate)

DB_NAME = "yassene_A_compos"
DB_DIR = file.path("data-raw",DB_NAME)
#RAW_DIR <- "ingest/Vilas_B"

ad <- read_h5ad(file.path(DB_DIR,"core_data_plus_de.h5ad"))




# ------------------------------------------
# 7 . create config and default files                   --
# ------------------------------------------
require(anndata)
require(reticulate)

DB_NAME = "yassene_A_compos"
DB_DIR = file.path("data-raw",DB_NAME)
#RAW_DIR <- "ingest/Vilas_B"

ad <- read_h5ad(file.path(DB_DIR,"core_data_plus_de.h5ad"))
db_prefix = "de"
diff_exp = readRDS( file = paste0(DB_DIR, "/", db_prefix, "_table.rds"))




###

# measures
#  This ordering is the "default"
measures <- list(obs = ad$obs_keys()[3:4],
                 var = ad$var_keys()[2:3])

# differentials  #if we care we need to explicitly state. defaults will be the order...
diffs <- list(diff_exp_groups =  levels(factor(diff_exp$group)),
              diff_exp_comp_type =  levels(factor(diff_exp$comp_type)),
              diff_exp_obs_name =  levels(factor(diff_exp$obs_name)),
              diff_exp_tests =  levels(factor(diff_exp$test_type)))

# Dimred
dimreds <- list(obsm = NA,
                varm = NA)

# what ad$obs do we want to make default values for...
# # should just pack according to UI?
default_factors <- c("Group")




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


yassene_A_compos_conf = readRDS( paste0(DB_DIR,"/",db_prefix,"_conf.rds") )
yassene_A_compos_def = readRDS( paste0(DB_DIR,"/",db_prefix,"_def.rds") )
yassene_A_compos_omics = readRDS( paste0(DB_DIR,"/",db_prefix,"_omics.rds") )
yassene_A_compos_meta = readRDS( paste0(DB_DIR,"/",db_prefix,"_meta.rds") )


#usethis::use_data(yassene_A_compos_conf, yassene_A_compos_def, yassene_A_compos_omics, yassene_A_compos_meta, overwrite = TRUE)




