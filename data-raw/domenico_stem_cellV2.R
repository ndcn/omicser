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
process_sn_prot_quant <- function(matrix_data_path, exp_name=NULL) {
  mp_all <- read.delim(matrix_data_path, sep = "\t", header = TRUE, as.is = TRUE)
  # remove contaminant
  if (length(grep("CON", mp_all$PG.ProteinAccessions)) > 0) {
    mp_all <- mp_all[-grep("CON", mp_all$PG.ProteinAccessions), ]
  }

  # remove NaN and pack into matrix
  mp_all <- mp_all[-which(mp_all$PG.Quantity == "NaN"), ]
  data <- matrix(ncol = length(unique(mp_all$R.FileName)), nrow = length(unique(mp_all$PG.ProteinAccessions)))
  index <- unique(mp_all$PG.ProteinAccessions)
  rownames(data) <- index

  all_files <- unique(mp_all$R.FileName)
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
  return(t(data))
}
# End process_sn_prot_quant_report ---------------------------



process_DIA_annot_de <- function( annot_de_file_path ){
  # full plath
  de_annot_data <- read.delim( annot_de_file_path, as.is = TRUE)

  # converte d_data inot diff_exp
  # TODO:  add the "test_type" currently it is just read from teh file which i think is
  # output from proprietary software.   check with Domenico...
  #   rather... --> pass off "ownership" of DIA prep functions to Domenico.
  diff_exp <- de_annot_data %>% transmute(group = gsub(" / ", "V", Comparison..group1.group2.),
                                      names = Group,
                                      obs_name = "Condition",
                                      test_type = "unknown",
                                      reference = Condition.Denominator,
                                      comp_type = 'grpVref',
                                      logfoldchange = AVG.Log2.Ratio,
                                      scores = NA,
                                      pvals = Pvalue,
                                      pvals_adj = Qvalue )

  feat_annots <- unique(de_annot_data[, c("UniProtIds",
                                      "Genes",
                                      "ProteinDescriptions",
                                      "ProteinNames",
                                      "GO.Cellular.Component",
                                      "GO.Molecular.Function",
                                      "GO.Biological.Process")])

  all_proteins <- unique(de_annot_data$UniProtIds)
  row.names(feat_annots) <- feat_annots$UniProtIds


  return ( list(de=diff_exp,
                annot=feat_annots))

}


get_DIA_conditions <- function(conditions_table_path){
  conditions_table <- read.delim( conditions_table_path, as.is = TRUE)  #contrasts <- diff_data
}


prep_DIA_files <- function(matrix_data_file,annot_de_file,conditions_table_file,path_root){

  raw_data <- process_sn_prot_quant(file.path(RAW_DIR, matrix_data_file))
  de_annot <- process_DIA_annot_de (file.path(RAW_DIR, annot_de_file ))
  conditions_table <- get_DIA_conditions( file.path(RAW_DIR,conditions_table_file) )


  features <- row.names(raw_data) #protein names
  # the file name names suck as obs_names... lets call them "Sample_1" "Sample_2" etc...
  obs_names <- paste0("Sample_",conditions_table$X)
  row.names(conditions_table) <- obs_names

  # re-order to conditions_table
  raw_data <- raw_data[conditions_table$File.Name,]


  # force the data matrix to match our obs(conditions_table) and var (annots)
  data <- raw_data[conditions_table$File.Name,feat_annots$UniProtIds]
  var_names <- row.names(data)

  # add some marginal statistics
  # transpose so its in X format
  tmp_mat <- data
  tmp_mat[is.na(tmp_mat)] <- 0

  de_annot$annot$expr_geomean <- colMeans( log(data),na.rm = TRUE) #exp minus 1?
  de_annot$annot$expr_mean <- colMeans(data,na.rm = TRUE)
  de_annot$annot$expr_var  <- colVars(data,na.rm = TRUE)
  de_annot$annot$expr_frac <- colMeans(tmp_mat>0)

  conditions_table$expr_var <- rowVars(data,na.rm=TRUE)
  conditions_table$expr_mean <- rowMeans(data,na.rm=TRUE)
  conditions_table$expr_frac <- rowMeans(tmp_mat>0)


  # scemaa
  #c("object", "data_mat","obs_meta","var_annot","omics","sample_ID","etc")
  data_list <- list(data_mat = data,
                    obs_meta = conditions_table,
                    var_annot = feat_annots,
                    omics = rownames(feat_annots),
                    sample_ID = rownames(data),
                    etc = NULL,
                    de = de_annot$de )

  return(data_list)
}

# ------------------------------------------
# 2. load data -----------------------------
# ------------------------------------------

RAW_DIR <- "ingest/Domenico_A"

# report table
matrix_data_file <- "20210524_093609_170805_aging_against_SC_merged_all_lib_2_Report.xls"

# candidate table without filters
annot_de_file <- "170805_aging_against_SC_merged_all_lib_2_candidates.xls"

# condition setup
conditions_table_file <- "170805_aging_against_SC_merged_all_lib_2_ConditionSetup.xls"

data_list <- prep_DIA_files(matrix_data_file,annot_de_file,conditions_table_file,RAW_DIR)



diff_exp <- data_list$de

db_prefix = "de"
saveRDS(diff_exp, file = paste0(DB_DIR, "/", db_prefix, "_table.rds"))



# ------------------------------------------
# 4. pack into anndata                    --
# ------------------------------------------


DB_NAME = "domenico_stem_cell"

helper_function<-('R/fct_ingestor.R')
source(helper_function)

ad <- setup_database(database_name=DB_NAME, data_in=data_list, db_meta=NULL , re_pack=TRUE)
ad$write_h5ad(filename=file.path("data-raw",DB_NAME,"core_data.h5ad"))



 ----------------------------
  # 5. post processing                      --
  # ------------------------------------------
require(anndata)
require(reticulate)
DB_NAME = "domenico_stem_cell"


ad <- read_h5ad(file.path("data-raw",DB_NAME,"core_data.h5ad"))
ad
raw <- ad$copy()

sc <- import("scanpy")

#outvals <- sc$pp$filter_cells(ad, min_genes=200, inplace=FALSE)

ad_tmp <- ad$copy()


sc$pp$scale(ad)

# ## Step 3: Do some basic preprocessing to run PCA and compute the neighbor graph
# sc$pp$pca(ad_tmp)
# sc$pp$neighbors(ad)
#
# ## Step 4: Infer clusters with the Louvain algorithm
# #sc$tl$louvain(ad_tmp)
# sc$tl$leiden(ad)
# ## Step 5: Compute tsne and umap embeddings
# sc$tl$umap(ad)

ad$raw <- raw
#ad_tmp$layers <- list(non_regressed=ad$X) #list('count'=layers)


ad$write_h5ad(filename=file.path("data-raw",DB_NAME,"core_data_plus_de.h5ad"))







# ------------------------------------------
# 6. dimension reduction - PCA / umap    --
# ------------------------------------------
require(anndata)
require(reticulate)

DB_NAME = "domenico_stem_cell"
DB_DIR = file.path("data-raw",DB_NAME)
#RAW_DIR <- "ingest/Vilas_B"

#ad <- read_h5ad(file.path(DB_DIR,"core_data_plus_de.h5ad"))
ad <- read_h5ad(file.path(DB_DIR,"core_data.h5ad"))


# ------------------------------------------
# 7 . create config and default files                   --
# ------------------------------------------
require(anndata)
require(reticulate)

DB_NAME = "domenico_stem_cell"
DB_DIR = file.path("data-raw",DB_NAME)
#RAW_DIR <- "ingest/Vilas_B"

#ad <- read_h5ad(file.path(DB_DIR,"core_data_plus_de.h5ad"))
ad <- read_h5ad(file.path(DB_DIR,"core_data.h5ad"))
db_prefix = "de"
diff_exp = readRDS( file = paste0(DB_DIR, "/", db_prefix, "_table.rds"))

# this all needs to go in YML files... so we cna just create the config/defaults at ingest

# measures
#  This ordering is the "default"
measures <- list(obs = ad$obs_keys()[11:14],
                 var = NA)
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


#usethis::use_data(domenico_stem_cell_conf, domenico_stem_cell_def, domenico_stem_cell_omics, domenico_stem_cell_meta, overwrite = TRUE)


