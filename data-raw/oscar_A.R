## code to prepare `oscar_microglia` dataset goes here

# not actually a good idea... better to just keep the annotation separate. as in ANNDATA format
require(Matrix)


# better script -------------------------

DATA_DIR <- "/Users/ahenrie/Projects/NDCN_dev/dbrowse/ingest"
DATA_DIR <- "ingest"
#basedir <- "/Users/ahenrie/Projects/NDCN_dev/Oscar"
DB_NAME = "Oscar_microglia1"

root <- file.path(DB_NAME,"microglia_")
file_path <- file.path(DATA_DIR,paste0(root,"matrix.mtx"))

counts <- Matrix::readMM(file_path) #DGTMATRIX .34GB

# META-DATA
#filenm <- "microglia_meta.csv"
file_path <- file.path(DATA_DIR,paste0(root,"meta.csv"))

#basedir <- "/Users/ahenrie/Projects/NDCN_dev/dbrowse/ingest/Oscar_microglia1/microglia_matrix.mtx"
meta_data <- read.csv(file_path, sep = ",", header = TRUE, stringsAsFactors = FALSE)
# meta_data2 <- vroom::vroom(file_path)

# FEATURES
#filenm <- "microglia_features.tsv"
file_path <- file.path(DATA_DIR,paste0(root,"features.tsv"))
features <- read.csv(file_path, sep = "\t", header = FALSE, stringsAsFactors = FALSE)
# features2 <- vroom::vroom(file_path,delim = "\t")
colnames(features) <- "genes"

# BARCODES
#filenm <- "microglia_barcodes.tsv"
file_path <- file.path(DATA_DIR,paste0(root,"barcodes.tsv"))
barcodes <- read.csv(file_path, sep = "\t", header = FALSE, stringsAsFactors = FALSE)
# barcodes2 <- vroom::vroom(file_path,delim = "\t")

# assert that meta_data$barcode == barcodes
all(meta_data$barcode == barcodes)

meta_data$sample_id <- paste0("sample_",1:length(meta_data$barcode))

# # Reload meta-data as factors for sumplicity
# meta_data$PA_DB_UID <- factor(meta_data$PA_DB_UID)
# meta_data$Sample_ID <- factor(meta_data$Sample_ID)
# meta_data$Gender <- factor(meta_data$Gender)
# meta_data$Status <- factor(meta_data$Status)
# meta_data$nAPOE <- factor(meta_data$nAPOE)
# meta_data$AOD <- factor(meta_data$AOD)
# meta_data$FBW <- factor(meta_data$FBW)
# meta_data$Cluster <- factor(meta_data$Cluster)


dimnames(counts) <- list(
  features$genes,  # row names
  meta_data$sample_id # column names / or barcodes$V1
)

X <- t(counts)  #transpose so rows and colums correspond to anndata expectation
obs <- meta_data
var_ <- features


# counts2 <- as(counts, "dgCMatrix") #0.257FGB
# counts3 <- as(counts, "dgeMatrix") #2.3GB

#save(data_matrix,file=out_file_path)


oscar_A_X <- X
oscar_A_obs <- obs

# usethis::use_data(oscar_A_X,oscar_A_obs,oscar_A_var)



X <- oscar_A_X
obs <- oscar_A_obs


##################################################v
##################################################
##################################################
#
#  calculated some marginal statistics and other things...
#
#
#  TODO:  use scanpy to add the marginal stats instead...
#  TODO:  conform to cellXgene schema
#
##################################################
##################################################
##################################################
##################################################
# 3. calculate fractional expression and differential levels in groups/clusters --------------------

##################
##################
# calculate the fractional expression for each gene within each group.

meta_data$cluster_label <- paste0("Cluster_", meta_data$Cluster)
all_clusters <- unique(meta_data$cluster_label[order(meta_data$Cluster)])

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



mean_z_mat <- frac_exp_mat
mean_log_mat <- frac_exp_mat
mean_z_log_mat <- frac_exp_mat


Z <- scale(X)
#z_submat <- (submat-mu_mat)/ %*% diag(1 / sd_mat)
#z_log_mat <- ( log_X -  rep(mu_log_mat, each=nrow(X)) )  / rep(sd_log_mat, each=nrow(X))
Z_log <- scale(log_X)
mu_log <- colMeans(log_X)
sd_log <- apply(log_X,2,sd)

fracexp <- colMeans(X>0)

# TODO:  Pack the overal mean, sd and logmean, logvar into var_
#
# copy the empty matrices
#

var_ <- cbind(var_,fracexp=fracexp,
              meanexp=mu_mat,
              sdexp=sd_mat,
              mulogexp=mu_log,
              sdlogexp=sd_log)

oscar_A_var <- var_



mu_matr <- rowMeans(X)
sd_matr <- apply(X, 1, sd)

mu_logr <- rowMeans(log_X)
sd_logr <- apply(log_X,1,sd)
fracexpr <- rowMeans(X>0)
obs <- cbind(obs,fracexp=fracexpr,
              meanexp=mu_matr,
              sdexp=sd_matr,
              mulogexp=mu_logr,
              sdlogexp=sd_logr)

oscar_A_obs <- obs

# i should just be able to normalize these....
# expressed <- (X > 0)

for (clust_i in all_clusters) {
  subs <- meta_data$sample_id[ meta_data$cluster_label == clust_i ]
  #subs <- meta_data$barcode[ meta_data$Cluster == clust_i ]

  # is colmeans more performant?
  frac_exp_mat[, clust_i] <- colMeans(X[subs,]>0)
  mean_z_mat[, clust_i] <- colMeans(Z[subs,])
  mean_log_mat[, clust_i] <- colMeans(log_X[subs,] )
  mean_z_log_mat[, clust_i] <- colMeans(Z_log[subs,] )
  # #
  print(c(clust_i, " done"))
} # rows=genes, columns = clusters

# add an overall value...

frac_exp_mat <- cbind(frac_exp_mat, all=fracexp)

mean_z_mat <- cbind(mean_z_mat,mu=mu_mat,sd=sd_mat)

mean_log_mat <- cbind(mean_log_mat,mu=mu_log,sd=sd_log)

mean_z_log_mat <- cbind(mean_z_log_mat,mu=mu_log,sd=sd_log)

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
layers <- NULL

oscar_A_obsm <- obsm
oscar_A_varm <- varm

uns <- list(frac_expres=colnames(varm$frac_expres),
            mean_z = colnames(varm$mean_z),
            mean_log = colnames(varm$mean_log),
            mean_z_log = colnames(varm$mean_z_log))

oscar_A_uns <- uns
oscar_A_var <- var_
oscar_A_layers <- layers

##################
##################
##################
##################
##################
##################
##################
##################
##################
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
