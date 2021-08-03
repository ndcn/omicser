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
  meta_data$barcode # column names / or barcodes$V1
)

X <- t(counts)  #transpose so rows and colums correspond to anndata expectation
obs <- meta_data
var_ <- features


# counts2 <- as(counts, "dgCMatrix") #0.257FGB
# counts3 <- as(counts, "dgeMatrix") #2.3GB

#save(data_matrix,file=out_file_path)


oscar_microglia_X <- X
oscar_microglia_obs <- obs
oscar_microglia_vars <- var_

usethis::use_data(oscar_microglia_X,oscar_microglia_obs,oscar_microglia_vars)



# 3. calculate fractional expression and differential levels in groups/clusters --------------------
X <- oscar_microglia_X
obs <- oscar_microglia_obs


# calculate the fractional expression for each gene within each group.

all_clusters <- unique(meta_data$Cluster[order(meta_data$Cluster)])


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

for (i in 1:length(all_clusters)) {
  clust_i <- all_clusters[i]
  subs <- meta_data$barcode[ meta_data$Cluster == clust_i ]
  # submat <- X[subs,]
  # log_submat <- log_X[subs,]
  # z_submat <- z_mat[subs,]
  # x_log_submat <- z_log_mat[subs,]
  #
  frac_exp_mat[, i] <- colMeans(expressed[subs,])

  mean_z_mat[, i] <- colMeans(z_mat[subs,])
  mean_log_mat[, i] <- colMeans(log_X[subs,] )
  mean_z_log_mat[, i] <- colMeans(z_log_mat[subs,] )
  # #
  # z_submat <-  (submat-mu_mat) / rep(sd_mat, each=nrow(submat))
  # #z_submat <- (submat-mu_mat)/ %*% diag(1 / sd_mat)
  # z_log_submat <- ( log10(submat+1.) -  rep(mu_log_mat, each=nrow(submat)) )  / rep(sd_log_mat, each=nrow(submat))



  # mean_z_mat[, i] <- apply(z_mat[subs,] , 2, mean)
  # #mean_log_mat[, clust_i] <- apply(submat , 2, mu_logf)
  # mean_log_mat[, i] <- apply(log_X[subs,] , 2, mean)
  # mean_z_log_mat[, i] <- apply(z_log_mat[subs,] , 2, mean)

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

oscar_microglia_obsm <- obsm
oscar_microglia_varm <- varm

usethis::use_data(oscar_microglia_obsm,oscar_microglia_varm,overwrite = TRUE)


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


X <- omicser::oscar_microglia_X
obs <- omicser::oscar_microglia_obs
var_ <- omicser::oscar_microglia_var
#obs <- obs_annots
#var_ <- data.frame(features,row.names = features)
obsm <- omicser::oscar_microglia_obsm
varm <- omicser::oscar_microglia_varm
uns <- list()
layers <- list()

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
#> AnnData object with n_obs × n_vars = 2 × 3
#>     obs: 'group'
#>     var: 'type'
#>     uns: 'a', 'b', 'c'
#>     obsm: 'ones', 'rand', 'zeros'
#>     varm: 'ones', 'rand', 'zeros'
#>     layers: 'spliced', 'unspliced'
