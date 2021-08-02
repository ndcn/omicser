## code to prepare `yassene_example` dataset goes here


################################################
###
### SCRIPT 1
###
#################################################

require(glmnet)

DATA_DIR <- "ingest"
DB_NAME = "Yassene_A"

# 1. load composition data --------------------
csv_name <- "example_composition.csv"

file_name <- file.path(DATA_DIR,DB_NAME,csv_name)
lipid_composition <- read.csv(file=file_name, header=TRUE, as.is = TRUE)

meta_cols <- c("Name","Group")
comp_group <- lipid_composition$Group
obs <- lipid_composition[meta_cols]
rownames(obs) <- lipid_composition$Name
comp_mat <- as.matrix( lipid_composition[ , !(names(lipid_composition) %in% meta_cols)])
rownames(comp_mat) <- obs$Name

class(comp_mat) <- "numeric"  # converts to numbers / NAs
#ZSCALE the matrix
zcomp_raw <- scale(comp_mat)

comp_mat[which(is.na(comp_mat), arr.ind = TRUE)] <- 0
# ZSCALE the zeroed matrix
zcomp <- scale(comp_mat)
# FIND columns wiht >2/3 zeros for the "REMOVE ZEROS" MATRIX
poor_comp <- which(colSums(comp_mat==0) > 2/3*dim(comp_mat)[1])


# 2. load concentration data --------------------
csv_name <- "example_species_conc.csv"
file_name <- file.path(DATA_DIR,DB_NAME,csv_name)
lipid_concentration <- read.csv(file=file_name, header=TRUE, as.is = TRUE)
conc_group <- lipid_concentration$Group

# make sure meta is aligned between concentration and composition
all(conc_group == comp_group)
all(lipid_concentration$Name == lipid_composition$Name)

xx <- paste(obs[, 2], obs[, 1], sep = "_")


conc_mat <- as.matrix( lipid_concentration[ , !(names(lipid_concentration) %in% meta_cols)])
rownames(conc_mat) <- obs$Name

class(conc_mat) <- "numeric"
#ZSCALE the matrix
zconc_raw <- scale(conc_mat)

# zero out NA: is this tricky because its sparse?
conc_mat[which(is.na(conc_mat), arr.ind = TRUE)] <- 0

poor_conc <- which(colSums(conc_mat==0) > 2/3*dim(conc_mat)[1])

# ZSCALE the matrix
zconc <- scale(conc_mat)

# # remove_zeros
# comp_scl_no0 <- comp[, -ind_rem]
# conc_scl_no0 <- conc[, -ind_rem]
all(poor_comp==poor_conc)



# experiment with the scaled (including)
regr_group <- obs$Group
g <- unique(regr_group)
# g
#
test_comp <- zcomp
test_conc <- zconc
regr_group <- as.numeric(factor(obs$Group))

ind_rem_group <- which(regr_group == which(table(regr_group) < 3))

# remove groups with less than 3 observations..
if (length(ind_rem_group) > 0) {
  test_comp <- test_comp[-ind_rem_group, ]
  test_conc <- test_conc[-ind_rem_group, ]
  regr_group <- regr_group[-ind_rem_group]
}

dim(test_comp)

lipids=colnames(test_comp)
var_ <- as.data.frame(lipids)

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


set.seed(100)
cvfit <- cv.glmnet(test_conc, regr_group, nlambda = 100, alpha = .8, family = "multinomial", type.multinomial = "grouped")
coef <- coef(cvfit, s = "lambda.min")
tmp <- as.matrix(coef$"1")
tmp1 <- tmp[which(tmp != 0)]
coef_names <- rownames(tmp)[which(tmp != 0)][-1]
ind_coef <- which(colnames(test_conc) %in% coef_names)

uns$conc_lasso_coef <- coef

var_$conc_sig_lasso_coef <- (colnames(test_conc) %in% coef_names)
rownames(var_) <- var_$lipids




raw <- conc_mat
X <- zconc

layers <- list(comp = conc_mat)
layers$zcomp <- zcomp

obsm <- NULL
varm <- NULL

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

# X <- omicser::yassene_A_X
# obs <- omicser::yassene_A_obs
# var_ <- omicser::yassene_A_var
# obsm <- omicser::yyassene_A_obsm
# varm <- omicser::yassene_A_varm
# uns <-  omicser::yyassene_A_uns
# layers <-  omicser::yassene_A_layers


#TODO:  pack into a list or anndata structure for simplicity...

require("data.table")
## TODO:  check
##      - do we need default_1, 2 and multi?
##

helper_functions<-('data-raw/data_prep_functions.R')

source(helper_functions)

db_dir = "data-raw"
# ui_config <- make_ingest_file_primitives(X,obs,var_,obsm=NA,varm=NA,uns=NA,
#                                           gex.assay = NA, gex.slot = c("data", "scale.data", "counts"),
#                                           gene_mapping = FALSE, db_prefix = "Vilas_A", db_dir = db_dir,
#                                           default_omics1 = NA, default_omics2 = NA, default_multi = NA,
#                                           default_dimred = NA, chunk_size = 500, meta_to_include = NA, legend_cols = 4,
#                                           max_levels_ui = 50)

db_prefix <- "yassene_A_"
make_ingest_files_primitive(X,obs,var_,obsm=obsm, varm=varm,
                            uns=uns, layers = layers, gex.assay = NA, gex.slot = c("data", "scale.data", "counts"),
                            gene_mapping = FALSE, db_prefix = db_prefix, db_dir = "data-raw",
                            default_omics1 = NA, default_omics2 = NA, default_multi = NA,
                            default_dimred = NA, chunk_size = 500, meta_to_include = NA, legend_cols = 4,
                            max_levels_ui = 50)

#yassene_A_conf = readRDS(file.path(db_dir,"test1conf.rds"))
yassene_A_conf = readRDS( paste0(db_dir,"/",db_prefix,"conf.rds") )

# defaults:  list of meta1, meta2, omics1, omics2, omics (list of 10). dimred, grp1, grp2
yassene_A_def = readRDS( paste0(db_dir,"/",db_prefix,"def.rds") )

# list of vars )e/g/ 3000 genes with counts?
yassene_A_omics = readRDS( paste0(db_dir,"/",db_prefix,"omics.rds") )
# use this sorted one to resort everything before packing into anndata

yassene_A_meta = readRDS( paste0(db_dir,"/",db_prefix,"meta.rds") )
yassene_A_X = readRDS( paste0(db_dir,"/",db_prefix,"X.rds") )
yassene_A_obs = readRDS( paste0(db_dir,"/",db_prefix,"obs.rds") )
yassene_A_obsm = readRDS( paste0(db_dir,"/",db_prefix,"obsm.rds") )
yassene_A_var = readRDS( paste0(db_dir,"/",db_prefix,"var.rds") )
yassene_A_varm = readRDS( paste0(db_dir,"/",db_prefix,"varm.rds") )
yassene_A_uns = readRDS( paste0(db_dir,"/",db_prefix,"uns.rds") )
yassene_A_layers = readRDS( paste0(db_dir,"/",db_prefix,"layers.rds") )



usethis::use_data(yassene_A_X,yassene_A_var,yassene_A_obs,
                  yassene_A_obsm,yassene_A_varm,yassene_A_layers,yassene_A_uns,
                  yassene_A_conf, yassene_A_def, yassene_A_omics, yassene_A_meta, overwrite = TRUE)

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


# X <- omicser::yassene_A_X
# obs <- omicser::yassene_A_obs
# var_ <- omicser::yassene_A_var
# obsm <- omicser::yassene_A_obsm
# varm <- omicser::yassene_A_varm
# uns <- omicser::yassene_A_varm
# layers <- omicser::yassene_A_layers

X <- yassene_A_X
obs <- yassene_A_obs
var_ <- yassene_A_var
obsm <- yassene_A_obsm
varm <- yassene_A_varm
uns <-  yassene_A_uns
layers <-  yassene_A_layers


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


adraw <- AnnData(
  X = raw,
  var = var_
)

ad$raw <- adraw


#write_h5ad(anndata = ad, filename = file.path(db_dir,"data-raw/Vilas_A.h5ad"))
# anndata R wrapper is broken.. .invoke python
#
ad$write_h5ad(filename="data-raw/yassene_A.h5ad")
#> AnnData object with n_obs × n_vars = 2 × 3

