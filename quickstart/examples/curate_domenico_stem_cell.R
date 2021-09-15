#==== INFO ========================================================
##
##  Muscle Stem Cell Proteomics
##  - DIA DATASET developed from "wRapper_example" script

#==== 0. preamble/setup ==================================================
# assume we are in the [omicser_path]
# getwd()
# pkgload::load_all('.')
require(golem)
golem::document_and_reload()


# BOOTSTRAP the options we have already set up...
# NOTE: we are looking in the "quickstart" folder.  the default is to look for the config in with default getwd()
omxr_options <- omicser::get_config()


CONDA_ENV <- omxr_options$conda_environment
DB_NAME <- omxr_options$database_names[1]
DB_ROOT_PATH <- omxr_options$db_root_path



#DB_NAME = "domenico_stem_cell"
DB_DIR = file.path(DB_ROOT_PATH,DB_NAME)
if (!dir.exists(DB_DIR)) {
  dir.create(DB_DIR)
}


#==== 1. documentation / provenance ==============================================================
# TODO:  markdown file or html with some copy about the database
#  - lab, paper link/name
#  summarize results / data origin whatever

db_meta <- list(
  organism = 'mmusculus',
  lab = "Ori/Ward",
  annotation_database =  "uniprot_for_annotation.RData?",
  title = "DIA proteomics",
  omic_type = "Proteomics",
  measurment = "normalized counts",
  pub = "TBD",
  date = format(Sys.time(), "%a %b %d %X %Y")
)

write_db_meta(db_meta,DB_NAME, db_root = DB_ROOT_PATH)


#==== 2. helper functions =================================================================================

# TODO:could also source from elsewhere.
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
  diff_exp <- de_annot_data %>% dplyr::transmute(group = gsub(" / ", "V", Comparison..group1.group2.),
                                                 names = Group,
                                                 obs_name = "Condition",
                                                 test_type = "unknown",
                                                 reference = Condition.Denominator,
                                                 comp_type = 'grpVref',
                                                 logfoldchanges = AVG.Log2.Ratio,
                                                 scores = NA,
                                                 pvals = Pvalue,
                                                 pvals_adj = Qvalue,
                                                 versus = gsub(" / ", " vs. ", Comparison..group1.group2.) )

  ###
   # the differential expression table has these fields:
   # group - the comparison   {names}V{reference}
   # names - what are we comparing?
   # obs_name  - name of the meta data variable
   # test_type - what statistic are we using
   # reference - the denomenator. or the condition we are comparing expressions values to
   # comp_type - grpVref or grpVrest. rest is all other conditions
   # logfoldchanges - log2(name/reference)
   # scores - statistic score
   # pvals - pvalues from the stats test. e.g. t-test
   # pvals_adj - adjusted pvalue (Q)
   # versus - label which we will choose in the browser
  ###

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
  raw_data <- process_sn_prot_quant(file.path(path_root, matrix_data_file))
  de_annot <- process_DIA_annot_de (file.path(path_root, annot_de_file ))
  conditions_table <- get_DIA_conditions( file.path(path_root,conditions_table_file) )


  features <- row.names(raw_data) #protein names
  # the file name names suck as obs_names... lets call them "Sample_1" "Sample_2" etc...
  obs_names <- paste0("Sample_",conditions_table$X)
  row.names(conditions_table) <- obs_names

  # re-order to conditions_table
  raw_data <- raw_data[conditions_table$File.Name,]

  # force the data matrix to match our obs(conditions_table) and var (annots)
  data <- raw_data[conditions_table$File.Name,de_annot$annot$UniProtIds]
  var_names <- row.names(data)

  # add some marginal statistics
  # transpose so its in X format
  tmp_mat <- data
  tmp_mat[is.na(tmp_mat)] <- 0

  de_annot$annot$expr_geomean <- Matrix::colMeans( log(data),na.rm = TRUE) #exp minus 1?
  de_annot$annot$expr_mean <- Matrix::colMeans(data,na.rm = TRUE)
  de_annot$annot$expr_var  <- matrixStats::colVars(data,na.rm = TRUE)
  de_annot$annot$expr_frac <- Matrix::colMeans(tmp_mat>0)

  conditions_table$expr_var <- matrixStats::rowVars(data,na.rm=TRUE)
  conditions_table$expr_mean <- Matrix::rowMeans(data,na.rm=TRUE)
  conditions_table$expr_frac <- Matrix::rowMeans(tmp_mat>0)


  # scemaa
  #c("object", "data_mat","obs_meta","var_annot","omics","sample_ID","etc")
  data_list <- list(data_mat = data,
                    obs_meta = conditions_table,
                    var_annot = de_annot$annot,
                    omics = rownames(de_annot$annot),
                    sample_ID = rownames(data),
                    etc = NULL,
                    de = de_annot$de )

  return(data_list)
}




# TODO:  define "get marker genes" or some such


#==== 3. load data -========================================================================================

RAW_DIR <- "raw_data/Domenico_A"

# report table
matrix_data_file <- "20210524_093609_170805_aging_against_SC_merged_all_lib_2_Report.xls"

# candidate table without filters
annot_de_file <- "170805_aging_against_SC_merged_all_lib_2_candidates.xls"

# condition setup
conditions_table_file <- "170805_aging_against_SC_merged_all_lib_2_ConditionSetup.xls"

data_list <- prep_DIA_files(matrix_data_file,annot_de_file,conditions_table_file,RAW_DIR)


# save diff expression data for later...
diff_exp <- data_list$de
# saveRDS(diff_exp, file = file.path(DB_DIR, "diff_expr_table.rds"))
saveRDS(diff_exp, file.path(DB_ROOT_PATH,DB_NAME, "db_de_table.rds"))

#==== 4. pack into anndata =========================================================================

# helper_function<-('data-raw/ingest_helpers.R')
# source(helper_function)


ad <- omicser::setup_database(database_name=DB_NAME,
                              db_path=DB_ROOT_PATH,
                              data_in=data_list,
                              db_meta=NULL ,
                              re_pack=TRUE)


ad$write_h5ad(filename=file.path(DB_ROOT_PATH,DB_NAME,"core_data.h5ad"))


#==== 5. post processing =========================================================================               --

# Add the raw field
raw <- ad$copy()
# make raw #includes NA
#     X="zeroed_NA" (for PCA)
#     "scaled"  (zeroed)

# use scanpy to do some scaling and calculations...
sc <- reticulate::import("scanpy")

#create layers, and raw
# na_to_0 (raw but zeroed)
# scaled
# X_is_scaled_na_to_0


zro_na <- ad$copy()
zro_na$X[is.na(zro_na$X)]<-0
zro_na <- zro_na$copy()

scaled <- ad$copy() #scaled
sc$pp$scale(scaled)

ad_out <- ad$copy()

ad_out$X[is.na(ad_out$X)]<-0
ad <- ad_out$copy()
sc$pp$scale(ad)

ad_copy <- ad$copy()
ad$layers <- list(zro_na=zro_na$X,
              scaled=scaled$X,
              X_is_scaled_na_to_0=ad_copy$X) #list('count'=layers)


ad$raw <- raw
ad$raw$to_adata()
ad <- ad$copy()



#  don't know how to make this work....
#sc$pp$highly_variable_genes(ad,n_top_genes=40)
ad$var$var_rank <- order(ad$var$expr_var)
# choose top 40 proteins by variance across dataset as our "targets"
target_omics <- ad$var_names[which(ad$var$var_rank <= 40)]

ad$var$decile <- dplyr::ntile(ad$var$expr_var, 10)


# save an intermediate file (incase we want to revert...)
ad$write_h5ad(filename=file.path(DB_ROOT_PATH,DB_NAME,"normalized_data.h5ad"))

#==== 5-a. dimension reduction - PCA / umap =========================================================================

## Step 3: Do some basic preprocessing to run PCA and compute the neighbor graphs
##

zro_na <- ad$copy()
zro_na$X <- ad$layers$get('zro_na')
sc$pp$pca(zro_na)
sc$pp$neighbors(zro_na)
## Step 4: Infer clusters with the leiden algorithm
sc$tl$leiden(zro_na)
## Step 5: Compute tsne and umap embeddings
sc$tl$umap(zro_na)


sc$pp$pca(ad)
sc$pp$neighbors(ad)
sc$tl$leiden(ad)
sc$tl$umap(ad)


ad$obsm$unscaled_X_pca <- zro_na$obsm$X_pca
ad$obsm$unscaled_X_umap <- zro_na$obsm$X_umap
ad$varm$unscaled_PCs <- zro_na$varm$PCs
ad$obsp$unscaled_distances <- zro_na$obsp$distances
ad$obsp$unscaled_connectivities <- zro_na$obsp$connectivities

# save an intermediate file (incase we want to revert...)
ad$write_h5ad(filename=file.path(DB_ROOT_PATH,DB_NAME,"norm_data_plus_dr.h5ad"))


#==== 6. differential expression  ======================================================================
#already have from the DIA helper fucnctions.
file.exists(file.path(DB_ROOT_PATH,DB_NAME, "db_de_table.rds"))

# diff_exp <- readRDS(file = file.path(DB_ROOT_PATH,DB_NAME, "diff_expr_table.rds"))


#==== 7. create configs =========================================================================
# what ad$obs do we want to make default values for...
default_factors <- c("Condition","Color","Replicate")

# differentials  #if we care we need to explicitly state. defaults will be the order...
config_list <- list(
  x_obs = c("Is.Reference","Condition","Replicate", "Label"),
  y_obs =  c("expr_var", "expr_mean", "expr_frac", "sample_ID", "leiden"), #MEASURES
  obs_groupby = c("Is.Reference","Condition","Replicate", "Label"),
  obs_subset = c("Is.Reference","Condition","Replicate", "Label"),

  x_var = c("decile"),
  y_var = c("expr_geomean", "expr_mean", "expr_var", "expr_frac" ),

  var_groupby = c("decile"),
  var_subset = c("decile"),

  diffs = list(diff_exp_comps = levels(factor(diff_exp$versus)),
               diff_exp_comp_type =  levels(factor(diff_exp$comp_type)), #i don't think we need this
               diff_exp_obs_name =  levels(factor(diff_exp$obs_name)),
               diff_exp_tests =  levels(factor(diff_exp$test_type))
  ),

  layers = c("X","raw","X_is_scaled_na_to_0","scaled","zro_na"),

  # Dimred
  dimreds = list(obsm = ad$obsm_keys(),
                 varm = ad$varm_keys()),

  # what ad$obs do we want to make default values for...
  # # should just pack according to UI?
  default_factors = c("Condition","Color","Replicate"),
  target_omics = target_omics,
  omic_details = c("Genes",
                   "ProteinDescriptions",
                   "ProteinNames",
                   "GO.Cellular.Component",
                   "GO.Molecular.Function",
                   "GO.Biological.Process",
                   "decile")

)



omicser::write_db_conf(config_list,DB_NAME, db_root = DB_ROOT_PATH)

#==== 8. write data file to load  =========================================================================
ad$write_h5ad(filename=file.path(DB_ROOT_PATH,DB_NAME,"db_data.h5ad"))

