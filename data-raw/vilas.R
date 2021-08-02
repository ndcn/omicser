## code to prepare `vilas` dataset goes here
#install.packages("filesstrings")
#library(filesstrings) # should probably be using fs
#install.packages("fs")
#path_ext_remove(path)

require(Matrix)

#setwd("~/Projects/NDCN_dev/Vilas/shinydotplot")
basedir <- "/Users/ahenrie/Projects/NDCN_dev/Vilas/shinydotplot/"


# ingest script -----------------------
DATA_DIR <- "/Users/ahenrie/Projects/NDCN_dev/dbrowse/ingest"
DATA_DIR <- "ingest"
#basedir <- "/Users/ahenrie/Projects/NDCN_dev/Oscar"
DB_NAME = "Vilas_A"

# 1. load count data --------------------
csv_name <- "41467_2020_19737_MOESM15_ESM.csv"

file_name <- file.path(DB_NAME,csv_name)
file_path <- file.path(DATA_DIR,file_name)
read.csv <- read.csv(file=file_path, header=FALSE, sep=",",row.names )
count_table <- read.csv(file=file_path, header=FALSE, sep=",", row.names=1)
features <- data.frame("genes" = row.names(count_table))

countm <- as.matrix(count_table)
#counts <- as(counts, "dgeMatrix")
countm <- as(countm, "dgTMatrix")
counts <- t(countm)  #transpose so samples are rows (as in ANNDATA format)

# norm_data <- read.csv(file=file_path, header=TRUE, sep=",", row.names=1)
# norm_data1 <- as.matrix(norm_data) #2.2GB
# norm_data2 <- as(norm_data1, "dgeMatrix") #4.4GB
# norm_data3 <- as(norm_data1, "dgTMatrix") #0.29GB

# 2. load annotation data --------------------
annnot_csv <- "41467_2020_19737_MOESM17_ESM.csv"

file_name <- file.path(DB_NAME,annnot_csv)
file_path <- file.path(DATA_DIR,file_name)

annots <- read.csv(file=file_path, header=TRUE, sep=",", row.names=NULL)

# pack into X, obs, vars, meta?
X <- counts
obs <- annots
vars <- features



gene_names <- unique(vars$genes)



### annot_mat
###save generated files - note that the expression matrix is re-saved, in case any changes are made to it###
# basedir <- "/Users/ahenrie/Projects/NDCN_dev/Vilas/shinydotplot/"
# load(file=paste0(basedir,"generated_files/anno_expr_mat.rda"))

X <- omicser::Vilas_A_X  #should this be reactive?
obs <- omicser::Vilas_A_obs  #should this be reactive?
vars <- omicser::Vilas_A_vars  #should this be reactive?


# 3. calculate fractional expression and differential levels in groups/clusters --------------------
dat_mat <- X
###generate annotation matrix###
annot_mat <- data.frame(sample_name = obs$sample_id,
                        cluster_id = obs$cluster_label,
                        cluster_label = paste0("Cluster_",obs$cluster_label),
                        stringsAsFactors = F)

# calculate the fractional expression for each gene within each group.

all_clusters <- unique(annot_mat$cluster_label[order(annot_mat$cluster_id)])
frac_mat <- matrix(0, nrow = ncol(dat_mat), ncol = length(all_clusters) )

# create
rownames(frac_mat) <- colnames(dat_mat)
colnames(frac_mat) <- all_clusters

mu_mat <- apply(dat_mat, 2, mean)
sd_mat <- apply(dat_mat, 2, sd)


mu_logf <- function(inval){
            outval <- mean(log10(inval+1))
            return(outval)
            }

sd_logf <- function(inval){
            outval <- sd(log10(inval+1))
            return(outval)
          }

mu_log_mat <- apply(dat_mat, 2, mu_logf)
sd_log_mat <- apply(dat_mat, 2, sd_logf)

mean_z_mat <- frac_mat
mean_log_mat <- frac_mat
mean_z_log_mat <- frac_mat


for (clust_i in all_clusters) {
  submat <- dat_mat[annot_mat$sample_name[annot_mat$cluster_label == clust_i],]

  frac_mat[, clust_i] <- apply(submat > 0, 2, mean)

  z_submat <-  (submat-mu_mat) / rep(sd_mat, each=nrow(submat))
  #z_submat <- (submat-mu_mat)/ %*% diag(1 / sd_mat)
  z_log_submat <- (log10(submat+1.) -  rep(mu_log_mat, each=nrow(submat)))  / rep(sd_log_mat, each=nrow(submat))

  mean_z_mat[, clust_i] <- apply(z_submat , 2, mean)
  mean_log_mat[, clust_i] <- apply(submat , 2, mu_logf)
  mean_z_log_mat[, clust_i] <- apply(z_log_submat , 2, mean)

  print(c(clust_i, " done"))
} # rows=genes, columns = clusters

# now stack the summary frac_mat, mean_z_mat, mean_log_mat,mean_z_log_mat into vars
colnames(frac_mat) <- paste0("frac_exp_",all_clusters)
colnames(mean_z_mat) <- paste0("mean_z_",all_clusters)
colnames(mean_log_mat) <- paste0("mean_log10mu_",all_clusters)
colnames(mean_z_log_mat) <- paste0("mean_zlog10_", all_clusters)

new_vars <- cbind(frac_mat,mean_z_mat,mean_log_mat,mean_z_log_mat)

variables<-features[2:dim(features)[1],]



Vilas_A_X <- X
Vilas_A_obs <- obs
Vilas_A_vars <- new_vars

usethis::use_data(Vilas_A_X,Vilas_A_obs,Vilas_A_vars,overwrite = TRUE)


#######################################################################
#######################################################################
########################################################
#######################################################################
#######################################################################
#######################################################################

obs <- omicser::Vilas_A_obs
vars <- omicser::Vilas_A_vars
X <- omicser::Vilas_A_X



# do the subsets ANNDATA style... rows (features/cell-ID/Groups) x columns (genes/proteins)

# this is where we could subset the genes to speed things up...
gene_names <- colnames(dat_mat)#[match(tolower(genes), tolower(rownames(norm_data)))]
cluster_labels <- annot_mat$cluster_label

submat <- matrix(NA, nrow = length(gene_names), ncol = length(all_clusters) + 1)
colnames(submat) <- c(all_clusters, "50% expression size")
rownames(submat) <- gene_names
fracmat <- submat


# norms
norm_submat <- submat
norm_fracmat <- submat

# do each cluster relative to the rest.
for (clust_i in all_clusters) {

  clusts <- annot_mat$cluster_label[match(clust_i, annot_mat$cluster_label)]

  # annotmat indices
  in_clust <- annot_mat$cluster_label %in% clusts
  keepfeats <- which(in_clust)
  normfeats <- which(!in_clust)

  keepinds <- annot_mat$sample_name[keepfeats]
  norminds <- annot_mat$sample_name[normfeats]
  n_cells <- length(keepinds)
  n_normcells <- length(norminds)


  # SAVE THIS FOR computing pairs of clusters
  # n_cells <- c()
  # n_normcells <- c()
  # keepinds <- list()
  # norminds <- list()
  # # just in case the cells aren't both sorted we'll look the indicies up by name and intersect them
  # for (jj in clusts) {
  #   keepinds[[jj]] <- annot_mat$sample_name[intersect(which(cluster_labels == jj),keepfeats)]
  #   norminds[[jj]] <- annot_mat$sample_name[intersect(which(cluster_labels != jj),normfeats)]
  #   n_cells <- c(n_cells, length(keepinds[[jj]]))
  #   n_normcells <- c(n_normcells, length(norminds[[jj]]))
  # }
  # names(n_cells) <- clusts

  # for now just z-scale... we can pre-compute all the others later
  sc_allvals <- c(scale(allvals))
  sc2_outvals <- c(log10(allvals + 1))
  sc3_outvals <- c(scale(log10(allvals + 1)))

  for (ii in 1:length(gene_names)) {
    allvals <- dat_mat[keepfeats, gene_names[ii]]
    names(allvals) <- rownames(dat_mat)[keepfeats]

    # for now just z-scale... we can pre-compute all the others later
    sc_allvals <- c(scale(allvals))
    sc2_outvals <- c(log10(allvals + 1))
    sc3_outvals <- c(scale(log10(allvals + 1)))

    names(sc_allvals) <- names(allvals)


    #allvals2 is the "rest of the genes outside of the cluster"
    allvals2 <- dat_mat[normfeats, gene_names[ii]]
    names(allvals2) <- rownames(dat_mat)[normfeats]

    # for now just z-scale... we can pre-compute all the others later
    sc_allvals2 <- c(scale(allvals2))
    names(sc_allvals2) <- names(allvals2)

    submat[genes[ii]] <- mean(allvals2[keepinds)
    frac_mat[genes[ii]] <- sum(allvals[keepinds > 0) / numcells

    submat[genes[ii]] <- mean(allvals2[keepinds])
    frac_mat[genes[ii]] <- sum(allvals[keepinds > 0) / numcells


    # for (jj in clusts) {
    #   submat[genes[ii], jj] <- mean(allvals2[keepinds[[jj]]])
    #   frac_mat[genes[ii], jj] <- sum(allvals[keepinds[[jj]]] > 0) / numcells[jj]
    #
    #   submat[genes[ii], jj] <- mean(allvals2[keepinds[[jj]]])
    #   frac_mat[genes[ii], jj] <- sum(allvals[keepinds[[jj]]] > 0) / numcells[jj]
    #
    # }
  }

  frac_mat[, ncol(frac_mat)] <- 0.5


}
# get the subset containing our clusst




colnames(submat) <- c(groups, "50% expression size")
rownames(submat) <- genes



# helper functions
#############

my_scale <- function(invals,scaletype="zscore"){
  if (scaletype == "zscore") {
    outvals <- c(scale(invals))
  } else if (scaletype == "logcpm") {
    outvals <- c(log10(invals + 1))
  } else if (scaletype == "zscorelogcpm") {
    outvals <- c(scale(log10(invals + 1)))
  } else { # should never get here
    outvals <- c(scale(invals))
  }
  return(outvals)
}

set_plot_dimensions <- function(width_choice, height_choice) {
  options(repr.plot.width = width_choice, repr.plot.height = height_choice)
}

# Function to set plot colormap
my_lut_plot <- function(scaletype,submat){
  colinds <- c(submat)
  if (scaletype == "logcpm") {
    colinds[colinds > 5] <- 5
    colinds <- round(colinds * 200 / 5) + 1
    colvec <- colorRampPalette(c("white", "yellow", "red"))
    colvec <- colvec(201)
    colplot <- colvec[colinds]
  } else {
    colinds[colinds > 2] <- 2
    colinds[colinds < -2] <- -2
    colinds <- round((colinds + 2) * 200 / (4)) + 1
    colvec <- colorRampPalette(c("blue", "white", "red")) #don't need to do this EVERY time
    colvec <- colvec(201)
    colplot <- colvec[colinds]
  }
  outlist <- list("lut"=colvec,"lut_plot"=colplot)
  return(outlist)
}

#colorbar
my_color_bar <- function(scaletype="zsocre",lut, nticks = 11) {
  #lut = colvec
  if (scaletype == "zscore") {
    min = -2
    max = 2
    ylabval = "Mean z-score (over selected cells)"
    #color.bar(lut = colvec, min = -2, max=2, ylabval = "Mean z-score (over selected cells)")
  } else if (scaletype == "logcpm") {
    min = 0
    min = 5
    ylabval = "Log10(CPM+1)"
  } else if (scaletype == "zscorelogcpm") {
    min = -2
    max = 2
    ylabval = "Mean z-score of Log10(CPM+1) (over selected cells)"
  } else { # defensive... should never get here spawn error
    min = -2
    max=2
    ylabval = "Mean z-score (over selected cells)"
  }

  ticks = seq(min, max, len = nticks)
  scale <- (length(lut) - 1) / (max - min)

  # dev.new(width=1.75, height=5)
  plot(c(0, 10), c(min, max), type = "n", bty = "n", xaxt = "n", xlab = "", yaxt = "n", ylab = ylabval)
  axis(2, ticks, las = 1)
  for (i in 1:(length(lut) - 1)) {
    y <- (i - 1) / scale + min
    rect(0, y, 10, y + 1 / scale, col = lut[i], border = NA)
  }
}

### generate plotting matrix####

derive_dotplot  <- function( genes, clusts, scaletype, annot_mat, norm_data){

      inc <- 1 / 6
      #genes <- isolate(getgenelist())
      #clusts <- isolate(clustsetInput())

      # enforce some cluster niceness
      if (is.numeric(clusts[1])) {
        clusts <- annot_mat$cluster_label[match(clusts, annot_mat$cluster_id)]
      } else if (!(tolower(clusts[1]) %in% c("all", "")) ) {
        clusts <- annot_mat$cluster_label[match(clusts, annot_mat$cluster_label)]
      } else {
        clusts <- unique(annot_mat$cluster_label[order(annot_mat$cluster_id)])
      }

      genes <- rownames(norm_data)[match(tolower(genes), tolower(rownames(norm_data)))]


      # get the subset containing our clusst
      keepcols <- which(annot_mat$cluster_label %in% clusts)
      groups <- clusts



      testvec <- annot_mat$cluster_label

      submat <- matrix(NA, nrow = length(genes), ncol = length(groups) + 1)
      #, nrow = ncol(dat_mat)


      colnames(submat) <- c(groups, "50% expression size")
      rownames(submat) <- genes
      frac_mat <- submat
      numcells <- c()
      keepinds <- list()
      for (jj in groups) {
        keepinds[[jj]] <- annot_mat$sample_name[intersect(which(testvec == jj), keepcols)]
        numcells <- c(numcells, length(keepinds[[jj]]))
      }
      names(numcells) <- groups
      for (ii in 1:length(genes)) {
        allvals <- norm_data[genes[ii], keepcols]
        names(allvals) <- colnames(norm_data)[keepcols]

        allvals2 <- my_scale(allvals, scaletype)


        names(allvals2) <- names(allvals)
        for (jj in groups) {
          submat[genes[ii], jj] <- mean(allvals2[keepinds[[jj]]])
          frac_mat[genes[ii], jj] <- sum(allvals[keepinds[[jj]]] > 0) / numcells[jj]
        }
      }

      frac_mat[, ncol(frac_mat)] <- 0.5
      xpos <- rep(1:ncol(submat), each = nrow(submat))
      ypos <- rep(1:nrow(submat), ncol(submat))

      tmpvals <- my_lut_plot(scaletype,submat)
      lut_plot <- tmpvals$lut_plot
      lut <- tmpvals$lut

      sizemat <- c(frac_mat)

      # set up the plotting dataframe : dfx
      dfx <- data.frame(x = c(xpos), y = max(ypos) - c(ypos) + 1, sizeval = sizemat, colsplot = lut_plot)
      dfx <- dfx[which(!is.na(dfx$sizeval)), ]
      dfx <- dfx[dfx$sizeval > 0, ]

      set_plot_dimensions(width_choice = ncol(frac_mat), height_choice = nrow(frac_mat))

      layout(matrix(c(rep(1, 14), 3, 2), 2, 8, byrow = FALSE))
      par("mar" = c(12, 7, 14, 0))
      plot(c(1, ncol(submat)), c(1, nrow(submat)), pch = "", xlab = "", ylab = "", xaxt = "n", yaxt = "n", bty = "n")
      # submat[is.na(submat)]=0

      # horiz lines
      for (ii in 1:nrow(submat)) {
        lines(c(range(which(!is.na(submat[ii, ])))), rep(max(ypos) - ii, 2), col = "grey")
      }

      # vert lines
      for (jj in 1:(ncol(submat) - 1)) {
        upval <- max(dfx$y[dfx$x == jj])
        lines(c(jj, jj), c(upval, 0), col = "grey")
      }

      with(dfx,
           symbols(x = x, y = y,
                   circles = sqrt(sizeval),
                   inches = inc, ann = F,
                   bg = as.character(colsplot), fg = "black",
                   xlab = colnames(submat), add = T,
                   xlim = c(1, ncol(submat)),
                   ylim = c(1, nrow(submat)),
                   xaxt = "n", yaxt = "n"))
      axis(3, at = 1:ncol(submat), labels = c(colnames(submat)[-ncol(submat)], ""), las = 2, cex.axis = 1.2)
      axis(1, at = 1:ncol(submat), labels = c(paste0("n = ", numcells), colnames(submat)[ncol(submat)]), las = 2, cex.axis = 1.2)
      axis(2, at = 1:nrow(submat), label = rev(genes), las = 2, cex.axis = 1)

      my_color_bar(scaletype,colvec)

}

genes <- c("SNAP25,ENO2,SLC17A7,SLC17A6,GAD1,GAD2,SLC32A1,LAMP5,SST,CHODL,PVALB,VIP,CUX2,RORB,RBP4,GJA1,FGFR3,GFAP,OLIG1,OPALIN,PDGFRA,AIF1,TYROBP,NOSTRIN")
genes <-strsplit(genes, ",| |;|, ")[[1]]

clusts <- unique(annot_mat$cluster_label[order(annot_mat$cluster_id)])
clusts_i <- unique(annot_mat$cluster_id)

scaletype = "logcpm"
derive_dotplot( genes, clusts, scaletype, annot_mat, norm_data)


vilasA_data_table
vilasA_data_table

usethis::use_data(vilas, overwrite = TRUE)

#                                                                                                                                                                                        1L, 2L, 2L)), .Names = c("kwd1", "kwd2", "similarity"), class = "data.frame", row.names = c(NA, -4L))
# # Query the quantification_genes table to get expression and join with the Demographics table to get the sex
# salmon_query <- str_glue("
#                          SELECT participant_id, sex, sample_id, Name, TPM
#                          FROM `{BQ_TIER1_RELEASE_DATASET_STD}.Demographics`
#                          JOIN `{BQ_RNASEQ_RELEASE_DATASET_STD}.quantification_genes`
#                          USING (participant_id)
#                          WHERE (Name like 'ENSG00000229807.%'
#                          OR Name like 'ENSG00000183878.%')
#                          ORDER BY participant_id, sample_id, Name")
# print(salmon_query)
#
# read.csv(file='test_dataframe.csv', header=TRUE, sep=","



usethis::use_data(vilas, overwrite = TRUE)



library(anndata)
#
# ad <- AnnData(
#   X = matrix(1:6, nrow = 2),
#   obs = data.frame(group = c("a", "b"), row.names = c("s1", "s2")),
#   var = data.frame(type = c(1L, 2L, 3L), row.names = c("var1", "var2", "var3")),
#   layers = list(
#     spliced = matrix(4:9, nrow = 2),
#     unspliced = matrix(8:13, nrow = 2)
#   ),
#   obsm = list(
#     ones = matrix(rep(1L, 10), nrow = 2),
#     rand = matrix(rnorm(6), nrow = 2),
#     zeros = matrix(rep(0L, 10), nrow = 2)
#   ),
#   varm = list(
#     ones = matrix(rep(1L, 12), nrow = 3),
#     rand = matrix(rnorm(6), nrow = 3),
#     zeros = matrix(rep(0L, 12), nrow = 3)
#   ),
#   uns = list(
#     a = 1,
#     b = data.frame(i = 1:3, j = 4:6, value = runif(3)),
#     c = list(c.a = 3, c.b = 4)
#   )
# )
# ad
## AnnData object with n_obs × n_vars = 2 × 3
## obs: 'group'
## var: 'type'
## uns: 'a', 'b', 'c'
## obsm: 'ones', 'rand', 'zeros'
## varm: 'ones', 'rand', 'zeros'
## layers: 'spliced', 'unspliced'
### GENERATE unified vilas_data_table ###
########################################
# not actually a good idea... better to just keep the annotation separate. as in ANNDATA format
basedir <- "/Users/ahenrie/Projects/NDCN_dev/Vilas/shinydotplot/"
DATA_DIR = paste0(basedir,"starting_data")
#file_path <- file.path(DATA_DIR,csv_name)

norm_csv <- "41467_2020_19737_MOESM15_ESM.csv"
file_path <- file.path(DATA_DIR,norm_csv)

ad <- anndata::read_csv(file_path)
ad
# AnnData object with n_obs × n_vars = 33660 × 16245


#norm_data <- read.csv(file=file_path, header=TRUE, sep=",", row.names=1)
#norm_data2 <- t(norm_data)
#norm_data <- as.matrix(norm_data)
#save(norm_data,file=file.path(DATA_DIR,"norm_tmp.rda"))


annnot_csv <- "41467_2020_19737_MOESM17_ESM.csv"
file_path <- file.path(DATA_DIR,annnot_csv)
clust_data4 <- read.csv(file=file_path, header=TRUE, sep=",", row.names=1)
# out_file_path <- gsub(".csv",".rda",file_path)
# load(out_file_path)  ###genes in rows, cells in columns

ad2 <- anndata::read_csv(file_path)



# functionalize



csv_name <- "41467_2020_19737_MOESM17_ESM.csv"
# cell annotations

ingest_clust_data <- function(csv_name,basedir="") {
  #"norm_data" name DATA = "norm_data"
  DATA_DIR = paste0(basedir,"starting_data/")

  file_path <- file.path(DATA_DIR,csv_name)
  out_file_path <- gsub(".csv",".rda",file_path)
  #16245 rows ---  4 columns - sample_id,cluster_label,batch,color
  if (file.exists(out_file_path)){
    print("loading cluster identity .rda")
    load(out_file_path)  ###genes in rows, cells in columns
    print("loaded cluster identity .rda")
  } else {
    print("converting cluster identity: .csv to .rda")
    clust_data <- read.csv(file=file_path, header=TRUE, sep=",")
    save(clust_data,file=out_file_path)
    print("converted cluster identity  to .rda")
  }
  return (clust_data)
}

#16245 rows ---  4 columns - sample_id,cluster_label,batch,color
clust_data <- ingest_clust_data(csv_name,basedir)





#
# # conversion scripts... # test exporting as SingleCellExperiment / SummarizedExperiment
# # anndata conversion ------------
# require(anndata)
#
# ad <- AnnData(
#   X = X,
#   obs = obs,
#   var = vars
# )
#
# # SummarizedExperiment conversion ------------
#
# require(SummarizedExperiment)
#
# sse <- SummarizedExperiment(assays=list(counts=counts),
#                             colData=obs)
#
# require(SingleCellExperiment)
#
# sce <- SingleCellExperiment(list(counts=counts),
#                             colData=obs,
#                             rowData=vars,
#                             metadata=list(study="test"))
#
#
#
# # out_file_path <- gsub(".csv",".rda",file_path)
# # load(out_file_path)  ###genes in rows, cells in columns
#
#
# clust_data <- as.matrix(clust_data)
# #save(clust_data,file=file.path(DATA_DIR,"clust_tmp.rda"))
#
# clust_data2 <- t(clust_data)
#
# vilas_data_table <- cbind(norm_data,clust_data)
# vilas_data_table <- as(vilas_data_table, "dgeMatrix")
# vilas_data_table <- as(vilas_data_table, "dgTMatrix")
#
# ########################################
# ### GENERATE BROWSABLE FILES ###
# ########################################
#
# annot_mat <- data.frame(sample_name = clust_data$sample_id,
#                         cluster_id = clust_data$cluster_label,
#                         cluster_label = paste0("Cluster_",clust_data$cluster_label),
#                         stringsAsFactors = F)
#
# ###save generated files - note that the expression matrix is re-saved, in case any changes are made to it###
# save(annot_mat,norm_data,file=paste0(basedir,"generated_files/anno_expr_mat.rda"))
#
# ###calculate fractions of cells expressing each gene in each cluster###
# #AH: this could be vectorized?
# # load("generated_files/anno_exp_mat.rda")
# all_clusters <- unique(annot_mat$cluster_label[order(annot_mat$cluster_id)])
# frac_mat <- matrix(0, nrow = nrow(norm_data), ncol = length(all_clusters))
# rownames(frac_mat) <- rownames(norm_data)
# colnames(frac_mat) <- all_clusters
# mean_mat <- frac_mat
#
# for (clust_i in all_clusters) {
#   frac_mat[, clust_i] <- apply(norm_data[, annot_mat$sample_name[annot_mat$cluster_label == clust_i]] > 0, 1, mean)
#   mean_mat[, clust_i] <- apply(norm_data[, annot_mat$sample_name[annot_mat$cluster_label == clust_i]]    , 1, mean)
#   print(c(clust_i, " done"))
# }
# save(frac_mat, file = paste0(basedir,"generated_files/frac_mat.rda"))
# # is this actually used??
# save(mean_mat, file = paste0(basedir,"generated_files/mean_mat.rda"))
#
#
#
# ###  norm_data
#
#
#
#


## additional dataset...

# #
#         Here’s a link to our single-cell microglia data set as a
#         Seurat object, with counts tables, normalized counts tables,
#         UMAP coordinates, and cell metadata in the table. The UMAP
#         coordinates and clusters were generated with a previous v
#         ersion of Seurat, with obsolete normalization and clustering
#         routines. However, for visualization and to kick the tires,
#         I hope this is a good starting data set. This is also what
#         Chris Sifuentes has been playing with in cellxgene, so down
#         the road it could also be a good data set if we want to explore
#         cross-connectivity between visualization and cellxgene options.
# #
#
#setwd("~/Projects/NDCN_dev/Vilas/shinydotplot")
library("Seurat")


basedir <- "/Users/ahenrie/Projects/NDCN_dev/Vilas"


file_name <- "microglia_data_with_counts_RNA_SCT.rda"
file_path <- file.path(basedir,file_name)

load(file_path) # microglia_data


#UPDATE OBJECT
new_microglia_data <- UpdateSeuratObject(object = microglia_data)

# convert to ANNDATA
#
# if (!requireNamespace("remotes", quietly = TRUE)) {
#   install.packages("remotes")
# }
# remotes::install_github("mojaveazure/seurat-disk")
# devtools::install_github('satijalab/seurat-data')

library(SeuratData)

library(SeuratDisk)

out_file_name <- "new_microglia_data.h5Seurat"
out_file_path <- file.path(basedir,out_file_name)


SaveH5Seurat( new_microglia_data, filename = out_file_path, overwrite = TRUE)
# Creating h5Seurat file for version 3.1.5.9900
# Adding counts for RNA
# Adding data for RNA
# No variable features found for RNA
# No feature-level metadata found for RNA
# Adding counts for SCT
# Adding data for SCT
# Adding scale.data for SCT
# Adding variable features for SCT
# Adding feature-level metadata for SCT
# Adding counts for counts
# Adding data for counts
# No variable features found for counts
# No feature-level metadata found for counts
# Adding cell embeddings for pca
# Adding loadings for pca
# No projected loadings for pca
# Adding standard deviations for pca
# No JackStraw data for pca
# Adding cell embeddings for tsne
# No loadings for tsne
# No projected loadings for tsne
# No standard deviations for tsne
# No JackStraw data for tsne
#
#

Convert(out_file_path, dest = "h5ad")
# Validating h5Seurat file
# Adding scale.data from SCT as X
# Transfering meta.features to var
# Adding data from SCT as raw
# Transfering meta.features to raw/var
# Transfering meta.data to obs
# Adding dimensional reduction information for pca
# Adding feature loadings for pca
# Adding dimensional reduction information for tsne
#
new_file_path <- gsub(".h5Seurat",".h5ad",out_file_path)
ad <- read_h5ad(new_file_path)





X <- counts
obs <- meta_data
vars <- data_table

oscar_toy1_X <- X
oscar_toy1_obs <- obs
oscar_toy1_vars <- vars


# ingest script -----------------------
DATA_DIR <- "/Users/ahenrie/Projects/NDCN_dev/dbrowse/ingest"
DATA_DIR <- "ingest"
#basedir <- "/Users/ahenrie/Projects/NDCN_dev/Oscar"
DB_NAME = "Vilas_B"


