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


usethis::use_data(Vilas_A_X,Vilas_A_obs,Vilas_A_var,overwrite = TRUE)
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


# calculate the fractional expression for each gene within each group.

all_clusters <- unique(obs_annots$cluster_name[order(obs_annots$cluster_label)])

frac_exp_mat <- matrix(0, nrow = ncol(X), ncol = length(all_clusters) )

# create
rownames(frac_exp_mat) <- colnames(X)
colnames(frac_exp_mat) <- all_clusters

# grand stats over all observations
mu_mat <- apply(X, 2, mean)
sd_mat <- apply(X, 2, sd)


log_X <- log10(X + 1.0)


# mu_logf <- function(inval){
#   outval <- mean(log10(inval + 1))
#   return(outval)
# }
#
# sd_logf <- function(inval){
#   outval <- sd(log10(inval + 1))
#   return(outval)
# }
mu_log_mat <- apply(log_X, 2, mean)
sd_log_mat <- apply(log_X, 2, sd)

# mu_log_mat <- apply(X, 2, mu_logf)
# sd_log_mat <- apply(X, 2, sd_logf)

# copy the empty matrices
mean_z_mat <- frac_exp_mat
mean_log_mat <- frac_exp_mat
mean_z_log_mat <- frac_exp_mat


# zscore against our previously computed grand mu & sd
z_mat <-  (X-mu_mat) / rep(sd_mat, each=nrow(X))
#z_submat <- (submat-mu_mat)/ %*% diag(1 / sd_mat)
z_log_mat <- ( log_X -  rep(mu_log_mat, each=nrow(X)) )  / rep(sd_log_mat, each=nrow(X))

# i should just be able to normalize these....
expressed <- (X > 0)

for (clust_i in all_clusters) {
  subs <- obs_annots$sample_id[ obs_annots$cluster_label == clust_i ]
  # submat <- X[subs,]
  # log_submat <- log_X[subs,]
  # z_submat <- z_mat[subs,]
  # x_log_submat <- z_log_mat[subs,]
  #
  frac_exp_mat[, clust_i] <- apply(expressed[subs,], 2, mean)

  # #
  # z_submat <-  (submat-mu_mat) / rep(sd_mat, each=nrow(submat))
  # #z_submat <- (submat-mu_mat)/ %*% diag(1 / sd_mat)
  # z_log_submat <- ( log10(submat+1.) -  rep(mu_log_mat, each=nrow(submat)) )  / rep(sd_log_mat, each=nrow(submat))

  mean_z_mat[, clust_i] <- apply(z_mat[subs,] , 2, mean)
  #mean_log_mat[, clust_i] <- apply(submat , 2, mu_logf)
  mean_log_mat[, clust_i] <- apply(log_X[subs,] , 2, mean)
  mean_z_log_mat[, clust_i] <- apply(z_log_mat[subs,] , 2, mean)

  print(c(clust_i, " done"))
} # rows=genes, columns = clusters


# # now stack the summary frac_mat, mean_z_mat, mean_log_mat,mean_z_log_mat into vars
# colnames(frac_mat) <- paste0("frac_exp_",all_clusters)
# colnames(mean_z_mat) <- paste0("mean_z_",all_clusters)
# colnames(mean_log_mat) <- paste0("mean_log10mu_",all_clusters)
# colnames(mean_z_log_mat) <- paste0("mean_zlog10_", all_clusters)
#
# new_vars <- cbind(frac_mat,mean_z_mat,mean_log_mat,mean_z_log_mat)
#
#

varm = list(
  frac_expres = frac_exp_mat,
  mean_z = mean_z_mat,
  mean_log = mean_log_mat,
  mean_z_log = mean_z_log_mat
)

obsm = NULL

Vilas_A_obsm <- obsm
Vilas_A_varm <- varm

usethis::use_data(Vilas_A_obsm,Vilas_A_varm,overwrite = TRUE)


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


X <- omicser::Vilas_A_X
obs <- omicser::Vilas_A_obs
var_ <- omicser::Vilas_A_var
#obs <- obs_annots
#var_ <- data.frame(features,row.names = features)
obsm <- omicser::Vilas_A_obsm
varm <- omicser::Vilas_varm
uns <- list()

X <- omicser::Domenico_A_X
obs <- omicser::Domenico_A_obs
var_ <- omicser::Domenico_A_var
obsm <- omicser::Domenico_A_obsm
varm <- omicser::Domenico_A_varm
uns <- omicser::Domenico_A_uns


ad <- AnnData(
  X = X,
  obs = obs,
  var = var_,
  layers = list(),
  obsm = obsm,
  varm = varm,
  uns = list()
)

ad
#> AnnData object with n_obs × n_vars = 2 × 3
#>     obs: 'group'
#>     var: 'type'
#>     uns: 'a', 'b', 'c'
#>     obsm: 'ones', 'rand', 'zeros'
#>     varm: 'ones', 'rand', 'zeros'
#>     layers: 'spliced', 'unspliced'
