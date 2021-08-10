#######################################################################
#######################################################################
##
##  DomenicO "A" dataset.
##  - DIA DATASET developed from "wRapper_example" script
##
#######################################################################
#######################################################################
######################################################################### code to prepare `domenico_A` dataset goes here

# # Set options here
# options(golem.app.prod = FALSE) # TRUE = production mode, FALSE = development mode
# # Detach all loaded packages and clean your environment
# golem::detach_all_attached()
# # rm(list=ls(all.names = TRUE))
# # Document and reload your package
# golem::document_and_reload()

############
##
##  Part 1.  .xls ->  .txt
##    - instrument produced reports to blank text files...
############
################################################
###
### SCRIPT 1
###
#################################################
# from server
databases <- ("uniprot_for_annotation.RData")

DATA_DIR <- "ingest"
DB_NAME = "Domenico_A"

# organism -----------------
organism <- "mmusculus"
# max number of highlighted proteins in volcano/heatmap
top_n <- 100


# = process_sn_prot_quant_report ---------------------------
# compatible with SN >v.9, please use AO proteins scheme to export report from sn
#' here the basic data reading/munging is in this function (process_sn_prot_quant)
#' and the "report" is wrapped in a separate function (make_sn_prot_quant_report)
#' started as a Local copy of ori labs "process.sn.prot.quant.report"
#' compatible with SN >v.9, please use AO proteins scheme to export report from sn
#' not sure why they have so many versions. (e.g. process.sn.pep.quant.report etc.)
#'
#'
#' @param in_file
#' @param exp_name
#'
#' @return
#' @export
#'
#' @examples
#'
#'
#'
process_sn_prot_quant <- function(in_file, exp_name=NULL) {
  mp_all <- read.delim(in_file, sep = "\t", header = TRUE, as.is = TRUE)
  # remove contaminant
  if (length(grep("CON", mp_all$PG.ProteinAccessions)) > 0) {
    mp_all <- mp_all[-grep("CON", mp_all$PG.ProteinAccessions), ]
  } else {
    mp_all <- mp_all
  }
  # remove NaN and pack into matrix
  mp_all <- mp_all[-which(mp_all$PG.Quantity == "NaN"), ]
  data <- matrix(ncol = length(unique(mp_all$R.FileName)), nrow = length(unique(mp_all$PG.ProteinAccessions)))
  index <- unique(mp_all$PG.ProteinAccessions)
  rownames(data) <- index

  for (i in 1:length(unique(mp_all$R.FileName))) {
    t <- mp_all[mp_all$R.FileName == unique(mp_all$R.FileName)[i], ]
    rownames(t) <- t$PG.ProteinAccessions
    data[, i] <- t[index, "PG.Quantity"]
  }
  colnames(data) <- unique(mp_all$R.FileName)

  # write out expression data file
  if ( !is.null(exp_name) ) {
    write.table(data, paste(exp_name, "data_table.txt", sep = "_"), sep = "\t", quote = FALSE)
  }
  return(data)
}
# End process_sn_prot_quant_report ---------------------------

# report table
data_matrix_file <- "20210524_093609_170805_aging_against_SC_merged_all_lib_2_Report.xls"
# candidate table without filters
contrast_matrix_file <- "170805_aging_against_SC_merged_all_lib_2_candidates.xls"
# condition setup
conditions_table_file <- "170805_aging_against_SC_merged_all_lib_2_ConditionSetup.xls"

# make the table and generate report plots
#

file_path <- file.path(DATA_DIR,DB_NAME,data_matrix_file)
raw_data <- process_sn_prot_quant(file_path)


# TODO: turn these into tests
nrow(raw_data)
sum(complete.cases(raw_data))
features <- row.names(raw_data)

# load condition setup
file_path <- file.path(DATA_DIR,DB_NAME,conditions_table_file)
conditions_table <- read.delim(file_path, as.is = TRUE)
# these file name names suck... lets call them "Sample_1" "Sample_2" etc...
obs_names <- paste0("Sample_",conditions_table$X)

row.names(conditions_table) <- obs_names

# differential expression matrix
# TODO: "candidates "is not a good name...  "contrasts"?  "diff_data"?
# resave the "contrast_matrix" as candidates table
file_path <- file.path(DATA_DIR,DB_NAME,contrast_matrix_file)
diff_data <- read.delim(file_path, as.is = TRUE)  #contrasts <- diff_data

# re-order to conditions_table
raw_data <- raw_data[,conditions_table$File.Name]

rownames(raw_data)<- features




############
##
##  Part 2.  create feature/variables annotation var
##
############


length((features))

feat_annots <- unique(diff_data[, c("UniProtIds", "Genes","ProteinDescriptions","ProteinNames",
              "GO.Cellular.Component","GO.Molecular.Function","GO.Biological.Process")])
all_proteins <- unique(diff_data$UniProtIds)

row.names(feat_annots) <- feat_annots$UniProtIds

keep_vars <- which(features %in% feat_annots$UniProtIds)
del_vars <- which(!(features %in% feat_annots$UniProtIds))
# data <- raw_data[keep_vars,]

# force the data matrix to match our obs(conditions_table) and var (annots)
data <- raw_data[feat_annots$UniProtIds,conditions_table$File.Name]
var_names <- row.names(data)

require(tidyverse)

d_data <- diff_data %>% mutate(Comp_Group = gsub(" / ", "_vs_", Comparison..group1.group2.) )

wide_data <- d_data %>% pivot_wider(
                            id_cols = UniProtIds,
                            names_from = c(Comp_Group),
                            values_from = c(AVG.Log2.Ratio, Absolute.AVG.Log2.Ratio,Pvalue,Qvalue,
                                            X..of.Ratios,X..Unique.Total.Peptides, X..Change,Ratio),
                            names_repair = "check_unique"
                          )
to_join <- d_data %>% select(UniProtIds,
                             ProteinDescriptions,
                             ProteinNames, Genes,
                             GO.Cellular.Component,
                             GO.Molecular.Function,
                             GO.Biological.Process) %>%
                      distinct()

comp_data <- wide_data %>% inner_join(to_join)
comp_data <- as.data.frame(comp_data)
rownames(comp_data) <- comp_data$UniProtIds

comp_data <- comp_data[feat_annots$UniProtIds,]



log_data <- log(data)


var_ <- comp_data
var_$geomean <- rowMeans(log_data,na.rm = TRUE)
var_$mean <- rowMeans(data,na.rm = TRUE)

tmp_mat <- data
tmp_mat[which(is.na(tmp_mat), arr.ind = TRUE)] <- 0
var_$var  <- apply(tmp_mat,1,var)

obs <- conditions_table
obs$var_conc <- apply(tmp_mat,2,var)
obs$mean_conc <-apply(tmp_mat,2,mean)


############
##
##  Part 3.  ut into anndata/ omicser formated .rda
##
############

# #
#     X
#     vars
#     obs
#
#     database_name = "Domenico_A",
#     omics_type = "prote", # Transcript-omics, Prote-omics, Lipid-omics, Metabol-omics, misc-
#
#
#     var_names = NULL, # eg.  var_names: Genes, Proteins, Lipids
#     var_annotations = NULL, # colnames of variable annotations.  e.g. gene/protein-families, ontologies, lipid class
#
#     obs_names = NULL, # name of observations e.g. cell ID
#     obs_annotations = NULL,
#     # do i need this here?
#     uns_meta = NULL, # NULL if not used, otherwise list of unstructured annotations...
# #




# pack into X, obs, vars, meta?

Domenico_A_X <- t(data)
Domenico_A_obs <- obs
Domenico_A_var <- var_


X <- Domenico_A_X
obs <- Domenico_A_obs
var_ <- Domenico_A_var
# var_names <- row.names(var)
# obs_names <- row.names(obs)


# usethis::use_data(Domenico_A_X,Domenico_A_obs,Domenico_A_var,overwrite = TRUE)
#
# X <- omicser::Domenico_A_X
# var_ <- omicser::Domenico_A_var
# obs <- omicser::Domenico_A_obs

# pack differential/ aggregated into the obsm lists
# breal tje diff_data up into comparasons?
#
#


# TODO:  change this from top_n to ALL the proteins in features  and then we'll do the actual comparisons...

# THIS IS ALL NONSENSE


add_top_hits_label <- function( diff_data, comp, top_n = 100, pval_cutoff = 0.01, fc_cutoff = 0.58) {
  require(tidyverse)

  all_prots <- unique(diff_data$UniProtIds)

  # # # data table with
  # hits_all <- diff_data[diff_data$Comparison..group1.group2. == comp, ]
  # hits_all <- hits_all[order(hits_all$Qvalue, decreasing = FALSE), ]
  # top_hits <- hits_all[hits_all$Qvalue < pval_cutoff & hits_all$Absolute.AVG.Log2.Ratio > fc_cutoff, ]
  # top_n_hits <- top_hits[c(1:top_n), ]
  # hits_all$label <- ifelse(hits_all$Group %in% top_n_hits$Group, as.character(hits_all$Genes), NA)

    # data table with
    #hits_all <- diff_data[diff_data$Comparison..group1.group2. == comp, ]
    hits_all<- diff_data %>% filter(Comparison..group1.group2. == comp)  %>%
      arrange(Qvalue)  # don't actually need this since we are using top_n


    top_hits <- hits_all %>% filter(Qvalue<pval_cutoff & Absolute.AVG.Log2.Ratio>fc_cutoff ) %>%
      slice_head(n=top_n)

    # top_hits <- diff_data %>% filter(Comparison..group1.group2. == comp) %>%
    #   filter(Qvalue<pval_cutoff & Absolute.AVG.Log2.Ratio>fc_cutoff ) %>%
    #   slice_min(n=top_n,order_by=Qvalue)

    hits_all <- hits_all %>% mutate(label=ifelse(hits_all$Group %in% top_hits$Group, as.character(hits_all$Genes), NA))

    hits_all <- as.data.frame(hits_all)
    rownames(hits_all) <- hits_all$UniProtIds
    return(hits_all[all_prots,])
#
#   ############## Extract quantity from data table and annotate the table###################
#   data <- as.data.frame(data_in)
#   data[is.na(data)] <- 0
#   data <- cbind.data.frame(data, ID = rownames(data)) # mutate?
#   data$Genes <- NA
#
#   ############# Extract quantity from Datatable based on topN table and format table#########
#   data_annotated_top_n <- subset(data, ID %in% top_n_hits$Group)
#   data_annotated_top_n$ID <- factor(data_annotated_top_n$ID)
#   for (i in 1:nrow(data_annotated_top_n)) {
#     if (data_annotated_top_n$ID[i] %in% top_n_hits$ProteinGroups) {
#       data_annotated_top_n$Genes[i] <- top_n_hits[top_n_hits$ProteinGroups == data_annotated_top_n$ID[i], "Genes"]
#     }
#   }
#
#   row.names(data_annotated_top_n) <- make.unique(data_annotated_top_n$Genes)
#   data_annotated_top_n <- dplyr::select(data_annotated_top_n, -ID, -Genes)


}



comps <- unique(diff_data$Comparison..group1.group2.)
dat_top_n <- list()
for (i in 1:length(comps)) {
  comp <- comps[i]
  comp_name <- gsub(" / ", "_vs_", comp)
  # p_val_cutoff<-pval_cutoff  # TODO: parameterize this?
  # fc.cutoff<-0.58  #?
  print(comp_name)
  # get_volcano_selection(candidates,comp,top_n,pval_cutoff=0.01,fc_cutoff=0.58)
  hits <- add_top_hits_label(diff_data, comp, top_n = 50, pval_cutoff = 0.01, fc_cutoff = 0.58)
  dat_top_n[[comp_name]] <- as.matrix(hits)
  #diff_i <- diff_data[diff_data$Comparison..group1.group2. == comp, ]

  print(length(which(diff_data$Comparison..group1.group2. == comp)))
}

varm = dat_top_n



uns <- list(o_vs_y=colnames(varm$o_vs_y),
            g_vs_y = colnames(varm$g_vs_y),
            g_vs_o = colnames(varm$g_vs_o))


# list(
#   old_v_young = old_v_young,
#   geriatric_v_young = geriatric_v_young,
#   geriatric_v_old = geriatric_v_old
# )

obsm = NULL
layers = NULL

Domenico_A_obsm <- obsm
Domenico_A_varm <- varm
Domenico_A_uns <- uns
Domenico_A_layers <- layers


X <- Domenico_A_X
obs <- Domenico_A_obs
var_ <- Domenico_A_var
obsm <- Domenico_A_obsm
varm <- Domenico_A_varm
uns <- Domenico_A_uns
layers <- Domenico_A_layers



# OBSERVABLES
observables <- list(obs = c("var_conc","mean_conc"),
                    var = c("geomean","mean","var"),
                    layers = NA,
                    raw = c("X"),
                    obsm = NA)

# COMPARABLES
comparables <- list(varm = names(varm),
                    obsm = NA)

# Dimred
dimreds <- list(varm = NA,
                obsm = NA)

# usethis::use_data(Domenico_A_obsm,Domenico_A_varm,Domenico_A_uns,Domenico_A_layers,overwrite = TRUE)
#
#
# usethis::use_data_table()
#require("data.table")
## TODO:  check
##      - do we need default_1, 2 and multi?
##

#helper_functions<-('data-raw/data_prep_functions.R')
helper_function<-('data-raw/make_ingest_file_primitives.R')

source(helper_function)

db_dir = "data-raw"
# ui_config <- make_ingest_file_primitives(X,obs,var_,obsm=NA,varm=NA,uns=NA,
#                                           gex.assay = NA, gex.slot = c("data", "scale.data", "counts"),
#                                           gene_mapping = FALSE, db_prefix = "", db_dir = db_dir,
#                                           default_omics1 = NA, default_omics2 = NA, default_multi = NA,
#                                           default_dimred = NA, chunk_size = 500, meta_to_include = NA, legend_cols = 4,
#                                           max_levels_ui = 50)

db_prefix <- "Domenico_A_"
make_ingest_file_primitives(X,obs,var_,obsm,varm,uns, layers,
                            observables, comparables, dimreds,
                            default_omic = NA, default_dimred = NA, meta_to_include = NA,
                            gex.assay = NA, gex.slot = c("data", "scale.data", "counts"),
                            gene_mapping = FALSE, db_prefix = db_prefix, db_dir = db_dir,
                            chunk_size = 500,  legend_cols = 4,
                            max_levels_ui = 50)
#
#vilas_A_conf = readRDS(file.path(db_dir,"test1conf.rds"))
domenico_A_conf = readRDS( paste0(db_dir,"/",db_prefix,"conf.rds") )

# defaults:  list of meta1, meta2, omics1, omics2, omics (list of 10). dimred, grp1, grp2
domenico_A_def = readRDS( paste0(db_dir,"/",db_prefix,"def.rds") )

# list of vars )e/g/ 3000 genes with counts?
domenico_A_omics = readRDS( paste0(db_dir,"/",db_prefix,"omics.rds") )
# use this sorted one to resort everything before packing into anndata

domenico_A_meta = readRDS( paste0(db_dir,"/",db_prefix,"meta.rds") )
domenico_A_X = readRDS( paste0(db_dir,"/",db_prefix,"X.rds") )
domenico_A_obs = readRDS( paste0(db_dir,"/",db_prefix,"obs.rds") )
domenico_A_obsm = readRDS( paste0(db_dir,"/",db_prefix,"obsm.rds") )
domenico_A_var = readRDS( paste0(db_dir,"/",db_prefix,"var.rds") )
domenico_A_varm = readRDS( paste0(db_dir,"/",db_prefix,"varm.rds") )
domenico_A_uns = readRDS( paste0(db_dir,"/",db_prefix,"uns.rds") )
domenico_A_layers = readRDS( paste0(db_dir,"/",db_prefix,"layers.rds") )


#
# usethis::use_data(domenico_A_X,domenico_A_var,domenico_A_obs, overwrite = TRUE)
# usethis::use_data(domenico_A_obsm,domenico_A_varm,domenico_A_layers,domenico_A_uns, overwrite = TRUE)


usethis::use_data(domenico_A_conf, domenico_A_def, domenico_A_omics, domenico_A_meta, overwrite = TRUE)

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

require(anndata)

X <-  domenico_A_X
obs <-  domenico_A_obs
var_ <-  domenico_A_var
obsm <-  domenico_A_obsm
varm <-  domenico_A_varm
uns <-  domenico_A_uns
layers <-  domenico_A_layers



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

#write_h5ad(anndata = ad, filename = file.path(db_dir,"data-raw/domenico_A.h5ad"))
# anndata R wrapper is broken.. .invoke python
#
ad$write_h5ad(filename="data-raw/domenico_A.h5ad")
#> AnnData object with n_obs × n_vars = 2 × 3
#>     obs: 'group'

#source("data-raw/yassene_example.R",echo = FALSE)

#######################################################################
#######################################################################
##
##  ANNDATA example
##
#######################################################################
#######################################################################
#
# library(anndata)
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
#
#
#
#
#
#
#
#
# ad <- AnnData(
#   X = X,
#   obs = obs,
#   var = var_,
#   layers = list(),
#   obsm =  list(),
#   varm =  list(),
#   uns = uns
# )
#
# ad
# #> AnnData object with n_obs × n_vars = 2 × 3
#>     obs: 'group'
#>     var: 'type'
#>     uns: 'a', 'b', 'c'
#>     obsm: 'ones', 'rand', 'zeros'
#>     varm: 'ones', 'rand', 'zeros'
#>     layers: 'spliced', 'unspliced'

