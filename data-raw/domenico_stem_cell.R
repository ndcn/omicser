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
reticulate::use_condaenv(required = TRUE, condaenv = 'omxr')
require(anndata)
require(matrixStats)
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

# transpose so its in X format
data <- t(data)
log_data <- log(data)

tmp_mat <- data
tmp_mat[is.na(tmp_mat)] <- 0
# tmp_mat <- data
# tmp_mat[which(is.na(tmp_mat), arr.ind = TRUE)] <- 0

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
var_$geomean <- colMeans(log_data,na.rm = TRUE) #exp minus 1?
var_$mean <- colMeans(data,na.rm = TRUE)


var_$var  <- colVars(data,na.rm = TRUE)
var_$frac <- colMeans(tmp_mat>0)

obs <- conditions_table
obs$var_conc <- rowVars(data,na.rm=TRUE)
obs$mean_conc <- rowMeans(data,na.rm=TRUE)
obs$frac <- rowMeans(tmp_mat>0)

X <- data


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

raw_ad <- ad$copy()

sc <- import("scanpy")



# the scanpy tools for differential expression calculations expect log transformed data...
#  and the values are MUCH too big and cause overflow since they are only 32bit values if not
sc$pp$log1p(ad)




test_types <- c('wilcoxon','t-test_overestim_var')
comp_types <- c("grpVrest", "oVy","gVy","gVo")


helper_function<-('data-raw/compute_de_table.R')
source(helper_function)

diff_exp <- compute_de_table(ad,comp_types, test_types, obs_names = "Condition")

#  "log1p" layer and
#  "raw" layer for easy access / clarity
#



norm_adata <- sc$pp$normalize_total(ad,copy=TRUE)
normpc_adata <- sc$pp$normalize_per_cell(ad,copy=TRUE)
z_adata <- sc$pp$scale(ad,copy=TRUE)


rawX <- AnnData(
  X = raw_ad$X,
  var = raw_ad$var)

add_layers = list(
  log1p = ad$X,
  raw = raw_ad$X
)

ad$layers <- add_layers
ad$raw <- rawX
ad$X <- raw_ad$X  # set back to "raw"





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

DB_NAME = "domenico_stem_cell"
DB_DIR = file.path("data-raw",DB_NAME)
#RAW_DIR <- "ingest/Vilas_B"

ad <- read_h5ad(file.path(DB_DIR,"core_data_plus_de.h5ad"))




# ------------------------------------------
# 7 . create config and default files                   --
# ------------------------------------------
require(anndata)
require(reticulate)

DB_NAME = "domenico_stem_cell"
DB_DIR = file.path("data-raw",DB_NAME)
#RAW_DIR <- "ingest/Vilas_B"

ad <- read_h5ad(file.path(DB_DIR,"core_data_plus_de.h5ad"))
db_prefix = "de"
diff_exp = readRDS( file = paste0(DB_DIR, "/", db_prefix, "_table.rds"))



# measures
#  This ordering is the "default"
measures <- list(obs = ad$obs_keys()[11:13],
                 var = ad$var_keys()[32:35])
# [1] "sct.detection_rate"    "sct.gmean"             "sct.variance"
# [4] "sct.residual_mean"     "sct.residual_variance" "sct.variable"

# differentials  #if we care we need to explicitly state. defaults will be the order...
diffs <- list(diff_exp_groups =  levels(factor(diff_exp$group)),
              diff_exp_comp_type =  levels(factor(diff_exp$comp_type)),
              diff_exp_obs_name =  levels(factor(diff_exp$obs_name)),
              diff_exp_tests =  levels(factor(diff_exp$test_type)))

# Dimred
dimreds <- list(varm = NA,
                obsm = NA)


# what ad$obs do we want to make default values for...
# # should just pack according to UI?
default_factors <- c("Condition","Color","Replicate")




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


domenico_stem_cell_conf = readRDS( paste0(DB_DIR,"/",db_prefix,"_conf.rds") )
domenico_stem_cell_def = readRDS( paste0(DB_DIR,"/",db_prefix,"_def.rds") )
domenico_stem_cell_omics = readRDS( paste0(DB_DIR,"/",db_prefix,"_omics.rds") )
domenico_stem_cell_meta = readRDS( paste0(DB_DIR,"/",db_prefix,"_meta.rds") )


usethis::use_data(domenico_stem_cell_conf, domenico_stem_cell_def, domenico_stem_cell_omics, domenico_stem_cell_meta, overwrite = TRUE)


