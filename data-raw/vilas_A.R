#######################################################################
#######################################################################
##
##  Vilas "A" dataset.
##  - info??
##
#######################################################################
#######################################################################

require(Matrix)

# ingest script -----------------------
DATA_DIR <- "/Users/ahenrie/Projects/NDCN_dev/omicser/ingest"
DATA_DIR <- "ingest"
DB_NAME = "Vilas_A"

# 1. load count data --------------------
csv_name <- "41467_2020_19737_MOESM15_ESM.csv"

file_name <- file.path(DB_NAME,csv_name)
file_path <- file.path(DATA_DIR,file_name)

#require(readr)
#count_table <- read_csv(file_path)  #readr returns a tibble

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
annnot_csv <- "41467_2020_19737_MOESM17_ESM.csv"

#file_name <- file.path(DB_NAME,annnot_csv)
file_path <- file.path(DATA_DIR,DB_NAME,annnot_csv)

#annots <- read_csv(file_path)
obs_annots <- read.csv(file=file_path, header=TRUE, sep=",", row.names=NULL)
#obs_annots$obs_names <- obs$sample_id
obs_annots$cluster_name = paste0("Cluster_",obs_annots$cluster_label)
obs_annots$cluster_id = obs_annots$cluster_label
row.names(obs_annots) <- obs_annots$sample_id

# pack into X, obs, var, meta?
#
X <- t(countm)[obs_annots$sample_id,features]
Vilas_A_X <- X
Vilas_A_obs <- obs_annots
Vilas_A_var <- data.frame(features,row.names = features)

# var_names <- row.names(var)
# obs_names <- row.names(obs)


#usethis::use_data(Vilas_A_X,Vilas_A_obs,Vilas_A_var,overwrite = TRUE)
#usethis::use_data(Vilas_A_var,overwrite = TRUE)




# ### annot_mat
# ###save generated files - note that the expression matrix is re-saved, in case any changes are made to it###
# # basedir <- "/Users/ahenrie/Projects/NDCN_dev/Vilas/shinydotplot/"
# # load(file=paste0(basedir,"generated_files/anno_expr_mat.rda"))
#
# X <- omicser::Vilas_A_X  #should this be reactive?
# obs <- omicser::Vilas_A_obs  #should this be reactive?
# var <- omicser::Vilas_A_var  #should this be reactive?


# 3. calculate fractional expression and differential levels in groups/clusters --------------------
X <- Vilas_A_X
obs <- Vilas_A_obs
var_ <- Vilas_A_var

# calculate the fractional expression for each gene within each group.

all_clusters <- unique(obs_annots$cluster_name[order(obs_annots$cluster_label)])

frac_exp_mat <- matrix(0, nrow = ncol(X), ncol = length(all_clusters) )

# create
rownames(frac_exp_mat) <- colnames(X)
colnames(frac_exp_mat) <- all_clusters

# grand stats over all observations
#mu_mat <- apply(X, 2, mean)

mu_mat <- colMeans(X)
sd_mat <- apply(X, 2, sd)
sd_mat <- sd_mat + 1e-100

log_X <- log(X + 1.)

# copy the empty matrices
mean_z_mat <- frac_exp_mat
mean_log_mat <- frac_exp_mat
mean_z_log_mat <- frac_exp_mat


Z <- scale(X)
#z_submat <- (submat-mu_mat)/ %*% diag(1 / sd_mat)
#z_log_mat <- ( log_X -  rep(mu_log_mat, each=nrow(X)) )  / rep(sd_log_mat, each=nrow(X))
Z_log <- scale(log_X)
mu_log <- colMeans(log_X)
sd_log <- apply(log_X,2,sd)

# i should just be able to normalize these....
# expressed <- (X > 0)

for (clust_i in all_clusters) {
  subs <- obs_annots$sample_id[ obs_annots$cluster_name == clust_i ]

  # is colmeans more performant?
  frac_exp_mat[, clust_i] <- colMeans(X[subs,]>0)
  mean_z_mat[, clust_i] <- colMeans(Z[subs,])
  mean_log_mat[, clust_i] <- colMeans(log_X[subs,] )
  mean_z_log_mat[, clust_i] <- colMeans(Z_log[subs,] )
  # #
  print(c(clust_i, " done"))
} # rows=genes, columns = clusters

# add an overall value...

frac_exp_mat <- cbind(frac_exp_mat, all=colMeans(X>0))

mean_z_mat <- cbind(mean_z_mat,mu=mu_mat,sd=sd_mat)

mean_log_mat <- cbind(mean_log_mat,mu=mu_log,sd=sd_log)

mean_z_log_mat <- cbind(mean_z_log_mat,mu=mu_log,sd=sd_log)


# # now stack the summary frac_mat, mean_z_mat, mean_log_mat,mean_z_log_mat into vars
# colnames(frac_mat) <- paste0("frac_exp_",all_clusters)
# colnames(mean_z_mat) <- paste0("mean_z_",all_clusters)
# colnames(mean_log_mat) <- paste0("mean_log10mu_",all_clusters)
# colnames(mean_z_log_mat) <- paste0("mean_zlog10_", all_clusters)
#
# new_vars <- cbind(frac_mat,mean_z_mat,mean_log_mat,mean_z_log_mat)
#
#

######3   these go in var... probably don't need varm... but they would be
######
#######varM for "matrix"... put things like dimreduction and 'PCs'
varm = list(
  frac_expres = frac_exp_mat,
  mean_z = mean_z_mat,
  mean_log = mean_log_mat,
  mean_z_log = mean_z_log_mat
)

obsm = NULL

Vilas_A_obsm <- obsm
Vilas_A_varm <- varm

uns <- list(frac_expres=colnames(varm$frac_expres),
            mean_z = colnames(varm$mean_z),
            mean_log = colnames(varm$mean_log),
            mean_z_log = colnames(varm$mean_z_log))
Vilas_A_uns <- uns

#usethis::use_data(Vilas_A_obsm,Vilas_A_varm,Vilas_A_uns,overwrite = TRUE)






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

# X <- omicser::Vilas_A_X
# obs <- omicser::Vilas_A_obs
# var_ <- omicser::Vilas_A_var
# obsm <- omicser::Vilas_A_obsm
# varm <- omicser::Vilas_A_varm
# uns <- omicser::Vilas_A_uns
# layers <- NULL


#TODO:  pack into a list or anndata structure for simplicity...

require("data.table")
## TODO:  check
##      - do we need default_1, 2 and multi?
##

#helper_functions<-('data-raw/data_prep_functions.R')
helper_function<-('data-raw/make_ingest_file_primitives.R')

source(helper_function)

db_dir = "data-raw"
db_prefix <- "Vilas_A_"
make_ingest_file_primitives(X,obs,var_,obsm=obsm, varm=varm,
                             uns=uns, gex.assay = NA, gex.slot = c("data", "scale.data", "counts"),
                             gene_mapping = FALSE, db_prefix = db_prefix, db_dir = "data-raw",
                                        default_omics1 = NA, default_omics2 = NA, default_multi = NA,
                                        default_dimred = NA, chunk_size = 500, meta_to_include = NA, legend_cols = 4,
                                        max_levels_ui = 50)

#vilas_A_conf = readRDS(file.path(db_dir,"test1conf.rds"))
vilas_A_conf = readRDS( paste0(db_dir,"/",db_prefix,"conf.rds") )

# defaults:  list of meta1, meta2, omics1, omics2, omics (list of 10). dimred, grp1, grp2
vilas_A_def = readRDS( paste0(db_dir,"/",db_prefix,"def.rds") )

# list of vars )e/g/ 3000 genes with counts?
vilas_A_omics = readRDS( paste0(db_dir,"/",db_prefix,"omics.rds") )
# use this sorted one to resort everything before packing into anndata

vilas_A_meta = readRDS( paste0(db_dir,"/",db_prefix,"meta.rds") )
vilas_A_X = readRDS( paste0(db_dir,"/",db_prefix,"X.rds") )
vilas_A_obs = readRDS( paste0(db_dir,"/",db_prefix,"obs.rds") )
vilas_A_obsm = readRDS( paste0(db_dir,"/",db_prefix,"obsm.rds") )
vilas_A_var = readRDS( paste0(db_dir,"/",db_prefix,"var.rds") )
vilas_A_varm = readRDS( paste0(db_dir,"/",db_prefix,"varm.rds") )
vilas_A_uns = readRDS( paste0(db_dir,"/",db_prefix,"uns.rds") )
vilas_A_layers = readRDS( paste0(db_dir,"/",db_prefix,"layers.rds") )



usethis::use_data(vilas_A_X,vilas_A_var,vilas_A_obs,
                  vilas_A_obsm,vilas_A_varm,vilas_A_layers,vilas_A_uns,
                  vilas_A_conf, vilas_A_def, vilas_A_omics, vilas_A_meta, overwrite = TRUE)

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
# X <- omicser::vilas_A_X
# obs <- omicser::vilas_A_obs
# var_ <- omicser::vilas_A_var
# obsm <- omicser::vilas_A_obsm
# varm <- omicser::vilas_A_varm
# uns <- omicser::vilas_A_varm
# layers <- omicser::vilas_A_layers
#
#
#
# a_X <- omicser::vilas_A_X
# a_obs <- omicser::vilas_A_obs
# a_var_ <- omicser::vilas_A_var
# a_obsm <- omicser::vilas_A_obsm
# a_varm <- omicser::vilas_A_varm
# a_uns <- omicser::vilas_A_varm
# a_layers <- omicser::vilas_A_layers
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

X <- vilas_A_X
obs <- vilas_A_obs
var_ <- vilas_A_var
obsm <- vilas_A_obsm
varm <- vilas_A_varm
uns <- vilas_A_uns
layers <-vilas_A_layers

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
ad$write_h5ad(filename="data-raw/vilas_A.h5ad")
#> AnnData object with n_obs × n_vars = 2 × 3
#>     obs: 'group'
#>     var: 'type'
#>     uns: 'a', 'b', 'c'
#>     obsm: 'ones', 'rand', 'zeros'
#>     varm: 'ones', 'rand', 'zeros'
#>     layers: 'spliced', 'unspliced'
source("data-raw/vilas_B.R",echo = FALSE)
