#######################################################################
#######################################################################
##
##  Muscle Stem Cell Proteomics
##  - DIA DATASET developed from "wRapper_example" script
##
#######################################################################
#######################################################################
#  formely `domenico_A`
#
#
# Golem options (Inactive) ---------------------------
# Set options here
# options(golem.app.prod = FALSE) # TRUE = production mode, FALSE = development mode
# # Detach all loaded packages and clean your environment
# golem::detach_all_attached()
# # rm(list=ls(all.names = TRUE))
# # Document and reload your package
# golem::document_and_reload()


# ------------------------------------------
# 0. preamble/setup -------------------------
# ------------------------------------------
require(tidyverse)
require(reticulate)
reticulate::use_condaenv(required = TRUE, condaenv = 'sc39')
require(anndata)

# create the folder to contain the raw data
DB_NAME = "domenico_stem_cell"
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


organism <- "mmusculus"
lab <- "Ori/Ward"
annotation_database <- "uniprot_for_annotation.RData"




# ------------------------------------------
# 2. helper functions ----------------------
# ------------------------------------------

# could also source from elsewhere.

# process_sn_prot_quant_report ---------------------------
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
#' @return data
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



# ------------------------------------------
# 2. load data -----------------------------
# ------------------------------------------

RAW_DIR <- "ingest/Domenico_A"

# report table
data_matrix_file <- "20210524_093609_170805_aging_against_SC_merged_all_lib_2_Report.xls"
raw_data <- process_sn_prot_quant( file.path(RAW_DIR,data_matrix_file) )

# candidate table without filters
contrast_matrix_file <- "170805_aging_against_SC_merged_all_lib_2_candidates.xls"
diff_data <- read.delim( file.path(RAW_DIR,contrast_matrix_file), as.is = TRUE)

# condition setup
conditions_table_file <- "170805_aging_against_SC_merged_all_lib_2_ConditionSetup.xls"
conditions_table <- read.delim( file.path(RAW_DIR,conditions_table_file), as.is = TRUE)  #contrasts <- diff_data


nrow(raw_data)
sum(complete.cases(raw_data))

features <- row.names(raw_data) #protein names

# the file name names suck as obs_names... lets call them "Sample_1" "Sample_2" etc...
obs_names <- paste0("Sample_",conditions_table$X)
row.names(conditions_table) <- obs_names

# re-order to conditions_table
raw_data <- raw_data[,conditions_table$File.Name]

rownames(raw_data) <- features





# ------------------------------------------
# 3. create feature/variables annotation var----------------
# ------------------------------------------

# max number of highlighted proteins in volcano/heatmap
top_n <- 100
length(features)

feat_annots <- unique(diff_data[, c("UniProtIds",
                                    "Genes",
                                    "ProteinDescriptions",
                                    "ProteinNames",
                                    "GO.Cellular.Component",
                                    "GO.Molecular.Function",
                                    "GO.Biological.Process")])
all_proteins <- unique(diff_data$UniProtIds)

row.names(feat_annots) <- feat_annots$UniProtIds

keep_vars <- which(features %in% feat_annots$UniProtIds)
del_vars <- which(!(features %in% feat_annots$UniProtIds))
# data <- raw_data[keep_vars,]

# force the data matrix to match our obs(conditions_table) and var (annots)
data <- raw_data[feat_annots$UniProtIds,conditions_table$File.Name]
var_names <- row.names(data)


## i think that DIFF_DATA in the wide form is the optimal pattern...

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



# ------------------------------------------
# 4. pack into anndata                    --
# ------------------------------------------


var_ <- comp_data
obs <- conditions_table
X <- t(data)


db_prefix = "core_data"
saveRDS(X, file = paste0(DB_DIR, "/", db_prefix, "_X.rds"))
saveRDS(obs, file = paste0(DB_DIR, "/", db_prefix, "_obs.rds"))
saveRDS(var_, file = paste0(DB_DIR, "/", db_prefix, "_var.rds"))


ad <- AnnData(
  X = X,
  obs = obs,
  var = var_,
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

DB_NAME = "domenico_stem_cell"
DB_DIR = file.path("data-raw",DB_NAME)
RAW_DIR <- "ingest/Domenico_A"

db_prefix = "core_data"
X = readRDS( file = paste0(DB_DIR, "/", db_prefix, "_X.rds"))
obs = readRDS(  file = paste0(DB_DIR, "/", db_prefix, "_obs.rds"))
var_ = readRDS( file = paste0(DB_DIR, "/", db_prefix, "_var.rds"))


ad <- read_h5ad(file.path(DB_DIR,"core_data.h5ad"))
ad



sc <- import("scanpy")


test_types <- c('wilcoxon','t-test_overestim_var')


log_adata <- sc$pp$log1p(ad,copy=TRUE)
diff_exp <- data.frame()
for (test_type in test_types) {
  comp_type <- "allVrest"
  key <- paste0(test_type,"_", comp_type)
  sc$tl$rank_genes_groups(log_adata, 'Condition', method='t-test_overestim_var', key_added = key)
  #diff_exp <- list()
  de_table <- sc$get$rank_genes_groups_df(log_adata, NULL, key=key)
  de_table$condition <- comp_type
  de_table$test_type <- test_type
  diff_exp <- dplyr::bind_rows(diff_exp, de_table)
  # for this dataset its straightforward to do all comparisons...

  comp_type <- "oVy"
  key <- paste0(test_type,"_", comp_type)
  sc$tl$rank_genes_groups(log_adata, 'Condition',group = 'o', reference='y', method='t-test_overestim_var', key_added = key)
  de_type <- sc$get$rank_genes_groups_df(log_adata, group = c('o'), key=key)
  de_table$condition <- comp_type
  de_table$test_type <- test_type
  diff_exp <- dplyr::bind_rows(diff_exp, de_table)

  comp_type <- "gVy"
  key <- paste0(test_type,"_", comp_type)
  sc$tl$rank_genes_groups(log_adata, 'Condition',group = 'g', reference='y', method='t-test_overestim_var', key_added = key)
  de_type <- sc$get$rank_genes_groups_df(log_adata, group = c('g'), key=key)
  de_table$condition <- comp_type
  de_table$test_type <- test_type
  diff_exp <- dplyr::bind_rows(diff_exp, de_table)

  comp_type <- "gVo"
  key <- paste0(test_type,"_", comp_type)
  sc$tl$rank_genes_groups(log_adata, 'Condition',group = 'g', reference='o', method='t-test_overestim_var', key_added = key)
  de_type <- sc$get$rank_genes_groups_df(log_adata, group = c('g'), key=key)
  de_table$condition <- comp_type
  de_table$test_type <- test_type
  diff_exp <- dplyr::bind_rows(diff_exp, de_table)

}

new_adata <- log_adata$copy()

new_adata$layers  <- list(log1p = log_adata$X)

new_adata$X <- ad$X

# TODO:  add de_table into uns
#new_adata$uns$de_table <- diff_exp
new_adata$write_h5ad(filename=file.path(DB_DIR,"core_data_plus_de.h5ad"))

db_prefix = "de"
saveRDS(diff_exp, file = paste0(DB_DIR, "/", db_prefix, "_table.rds"))


# also need to pack the diff_exp1 and diff_exp2 into easy to deal wiht tables for volcanos...



# ------------------------------------------
# 5. post processing                      --
# ------------------------------------------



# ------------------------------------------
# 6. create "meta" config + defaults      --
# ------------------------------------------


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

