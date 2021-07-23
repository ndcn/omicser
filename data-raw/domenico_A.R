#######################################################################
#######################################################################
##
##  DomenicO "A" dataset.
##  - DIA DATASET developed from "wRapper_example" script
##
#######################################################################
#######################################################################
######################################################################### code to prepare `domenico_A` dataset goes here



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
  # TODO?
  # replace "-" with "_" in colnames because they are replaced automatically by "."in readMSn
  # colnames(data)<-sub("-", "_", colnames(data))

  # correlation #
  # thin wrapper to plot correlation table

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
obs_names <- conditions_table$File.Name
row.names(conditions_table) <- obs_names

# differential expression matrix
# TODO: "candidates "is not a good name...  "contrasts"?  "diff_data"?
# resave the "contrast_matrix" as candidates table
file_path <- file.path(DATA_DIR,DB_NAME,contrast_matrix_file)
diff_data <- read.delim(file_path, as.is = TRUE)  #contrasts <- diff_data

# re-order to conditions_table
raw_data <- raw_data[,conditions_table$File.Name]
log2_data <- log2(raw_data)

############
##
##  Part 2.  create feature/variables annotation var
##
############


length((features))

feat_annots <- unique(diff_data[, c("UniProtIds", "Genes","ProteinDescriptions","ProteinNames",
              "GO.Cellular.Component","GO.Molecular.Function","GO.Biological.Process")])

row.names(feat_annots) <- feat_annots$UniProtIds

keep_vars <- which(features %in% feat_annots$UniProtIds)
del_vars <- which(!(features %in% feat_annots$UniProtIds))
# data <- raw_data[keep_vars,]

# force the data matrix to match our obs(conditions_table) and var (annots)
data <- raw_data[feat_annots$UniProtIds,conditions_table$File.Name]
var_names <- row.names(data)


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
Domenico_A_obs <- conditions_table
Domenico_A_var <- feat_annots

# var_names <- row.names(var)
# obs_names <- row.names(obs)


usethis::use_data(Domenico_A_X,Domenico_A_obs,Domenico_A_var,overwrite = TRUE)



# pack differential/ aggregated into the obsm lists
# breal tje diff_data up into comparasons?
#
#


# TODO:  change this from top_n to ALL the proteins in features  and then we'll do the actual comparisons...


get_comparisons_selection <- function(data_in, diff_data, comp, top_n = 100, pval_cutoff = 0.01, fc_cutoff = 0.58) {

  # data table with
  hits_all <- diff_data[diff_data$Comparison..group1.group2. == comp, ]
  hits_all <- hits_all[order(hits_all$Qvalue, decreasing = FALSE), ]
  top_hits <- hits_all[hits_all$Qvalue < pval_cutoff & hits_all$Absolute.AVG.Log2.Ratio > fc_cutoff, ]
  top_n_hits <- top_hits[c(1:top_n), ]
  hits_all$label <- ifelse(hits_all$Group %in% top_n_hits$Group, as.character(hits_all$Genes), NA)


  ############## Extract quantity from data table and annotate the table###################
  data <- as.data.frame(data_in)
  data[is.na(data)] <- 0
  data <- cbind.data.frame(data, ID = rownames(data)) # mutate?
  data$Genes <- NA

  ############# Extract quantity from Datatable based on topN table and format table#########
  data_annotated_top_n <- subset(data, ID %in% top_n_hits$Group)
  data_annotated_top_n$ID <- factor(data_annotated_top_n$ID)
  for (i in 1:nrow(data_annotated_top_n)) {
    if (data_annotated_top_n$ID[i] %in% top_n_hits$ProteinGroups) {
      data_annotated_top_n$Genes[i] <- top_n_hits[top_n_hits$ProteinGroups == data_annotated_top_n$ID[i], "Genes"]
    }
  }

  row.names(data_annotated_top_n) <- make.unique(data_annotated_top_n$Genes)
  data_annotated_top_n <- dplyr::select(data_annotated_top_n, -ID, -Genes)

  return(t(data_annotated_top_n))

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
  dat_top_n[[comp_name]] <- get_comparisons_selection(data, diff_data, comp, top_n = 50, pval_cutoff = 0.01, fc_cutoff = 0.58)

  #diff_i <- diff_data[diff_data$Comparison..group1.group2. == comp, ]

  print(length(which(diff_data$Comparison..group1.group2. == comp)))
}

uns = dat_top_n
# list(
#   old_v_young = old_v_young,
#   geriatric_v_young = geriatric_v_young,
#   geriatric_v_old = geriatric_v_old
# )

varm = NULL
bosm = NULL

Domenico_A_obsm <- obsm
Domenico_A_varm <- varm
Domenico_A_uns <- uns

obsm<- omicser::domenico_A_obsm
varm<- omicser::domenico_A_varm
varm<- omicser::domenico_A_varm


usethis::use_data(Domenico_A_obsm,Domenico_A_varm,Domenico_A_uns,overwrite = TRUE)




#######################################################################
#######################################################################
##
##  ANNDATA example
##
#######################################################################
#######################################################################

library(anndata)


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
  obsm =  list(),
  varm =  list(),
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

