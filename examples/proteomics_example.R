# Overview --------------
#     A start to finish example of curating a database and configuring the
#     NDCN omics browser for exploring it.   This database is Muscle Stem Cell Proteomics, and can
#     be obtained from {url / contacti ifo}.   The data consist of bulk protein charachterization from DIA
#     *(Data Independent Acquisition) mass spectrometry, yielding estimates of concentrations for each protein feature
##

#  Step 1: Set paths--------------
#  path where the omicser repo is cloned to (see install_scrip.R)
REPO_PATH <- getwd() # e.g. "~/Projects/NDCN_dev/omicser"
OMICSER_RUN_DIR <- REPO_PATH # file.path(REPO_PATH,"quickstart")

# where are the raw data files for this script?
RAW_DATA_DIR <- file.path(OMICSER_RUN_DIR,"quickstart/raw_data/DOMENICO_A")
# set where we will put the database folder
DB_ROOT_PATH <- file.path(OMICSER_RUN_DIR,"quickstart/test_db")

# what will we call the database
DB_NAME <-  list("Domenico DIA" = "domenico_stem_cell")


# name your python envirnoment
OMICSER_PYTHON <-  "pyenv_omicser"


# Step 2:  Load the omicser package
if (CLONED_OMICSER <- TRUE){
  # this should be a full path... e.g. ~/Projects/NDCN_dev/omicser
  # but for github, we will set relative to the repo BASE
  golem::document_and_reload(pkg = REPO_PATH)
} else {
  require(omicser)
}

# Step 3: Assert python back-end ----------------
#  for the curation we need to have scanpy

packages <- c("anndata")

USE_CONDA <- TRUE  #change this if you prefer virtual environments
if (USE_CONDA){
    # TODO: wrap this in a try
    reticulate::install_miniconda() #in case it is not already installed
    reticulate::conda_create(OMICSER_PYTHON, python_version = 3.8)

    reticulate::conda_install(envname=OMICSER_PYTHON,
                            channel = "conda-forge",
                            packages = packages )

    reticulate::use_condaenv(condaenv = OMICSER_PYTHON,
                                    conda = reticulate::conda_binary(),
                                    required = TRUE)

    # check that we have our python on deck
    reticulate::py_discover_config()

    if (!reticulate::py_module_available(module = "scanpy") ){
      reticulate::conda_install(envname=OMICSER_PYTHON,
                                             packages = "scanpy[leiden]",
                                             forge=TRUE,
                                             conda=reticulate::conda_binary())
    }

} else {
    PYTHON_PATH <-"path/to/python" # EDIT THIS WITH YOUR PYTHON3.8
    # so we'll do everything in conda instead with conda...  possibly doing it in stages makes it faster or less prone to failure... scapyy 1.7.2
    cmd <- c("which python3")
    PYTHON_PATH <- base::system( command = paste(cmd, collapse = " && "))

    if (!reticulate::virtualenv_exists( OMICSER_PYTHON )){
        reticulate::virtualenv_create( envname = OMICSER_PYTHON,
                                    python = PYTHON_PATH,
                                    packages =packages )
    }
    reticulate::use_virtualenv(virtualenv = OMICSER_PYTHON,
                            required = TRUE)


    # check that we have our python on deck
    reticulate::py_discover_config()

    if (!reticulate::py_module_available(module = "scanpy[leiden]") ){
     reticulate::py_install(packages = "scanpy[leiden]")
    }


}

# Step 4:  get the data ---------------
# we already have the data so we just need to locat it

# report table
matrix_data_file <- "20210524_093609_170805_aging_against_SC_merged_all_lib_2_Report.xls"
# candidate table without filters
annot_de_file <- "170805_aging_against_SC_merged_all_lib_2_candidates.xls"
# condition setup
conditions_table_file <- "170805_aging_against_SC_merged_all_lib_2_ConditionSetup.xls"


# Step 5:  define for source helper functions ------------
# process_sn_prot_quant_report
#       here the basic data reading/munging is in this function (process_sn_prot_quant)
#       and the "report" is wrapped in a separate function (make_sn_prot_quant_report)
#       started as a Local copy of ori labs "process.sn.prot.quant.report"
#       compatible with SN >v.9, please use AO proteins scheme to export report from sn
#       not sure why they have so many versions. (e.g. process.sn.pep.quant.report etc.)
#
#
#       @param in_file
#       @param exp_name
#
#       @return data
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
                                                    test_type = "unknown-test",
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


    return(
                    list(de=diff_exp,
                        annot=feat_annots)
                )

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
    tmp_mat <- data
    tmp_mat[is.na(tmp_mat)] <- 0

    de_annot$annot$expr_geomean <- Matrix::colMeans( log(data),na.rm = TRUE) #exp minus 1?
    de_annot$annot$expr_mean <- Matrix::colMeans(data,na.rm = TRUE)
    de_annot$annot$expr_var  <- matrixStats::colVars(data,na.rm = TRUE)
    de_annot$annot$expr_frac <- Matrix::colMeans(tmp_mat>0)

    conditions_table$expr_var <- matrixStats::rowVars(data,na.rm=TRUE)
    conditions_table$expr_mean <- Matrix::rowMeans(data,na.rm=TRUE)
    conditions_table$expr_frac <- Matrix::rowMeans(tmp_mat>0)

    # scema       #c("object", "data_mat","obs_meta","var_annot","omics","sample_ID","etc")
    data_list <- list(data_mat = data,
                        obs_meta = conditions_table,
                    var_annot = de_annot$annot,
                    omics = rownames(de_annot$annot),
                    sample_ID = rownames(data),
                    etc = NULL,
                    de = de_annot$de )

  return(data_list)
}

# Step 6:  pack data into AnnData format --------------
data_list <- prep_DIA_files(matrix_data_file,annot_de_file,conditions_table_file,RAW_DATA_DIR)

adata <- omicser::setup_database(database_name=DB_NAME,
                              db_path=DB_ROOT_PATH,
                              data_in=data_list,
                              db_meta=NULL ,
                              re_pack=TRUE)


if (FALSE)
    adata$write_h5ad(filename=file.path(DB_ROOT_PATH,DB_NAME,"core_data.h5ad"))


# Steps 7-9: CURATION
# Step 7: additional data processing ----

# use scanpy to do some scaling and calculations...
sc <- reticulate::import("scanpy")

# Add the raw field
raw <- adata$copy()

#create layers, and raw
# na_to_0 (raw but zeroed)
# scaled
# X_is_scaled_na_to_0

zro_na <- adata$copy()
zro_na$X[is.na(zro_na$X)]<-0
zro_na <- zro_na$copy()

scaled <- adata$copy() #scaled
sc$pp$scale(scaled)

ad_out <- adata$copy()

ad_out$X[is.na(ad_out$X)]<-0
adata <- ad_out$copy()
sc$pp$scale(adata)

ad_copy <- adata$copy()
adata$layers <- list(zro_na=zro_na$X,
              scaled=scaled$X,
              X_is_scaled_na_to_0=ad_copy$X) #list('count'=layers)

adata$raw <- raw
adata$raw$to_adata()
#adata <- adata$copy()

adata$var$vmr <- adata$var$expr_var/adata$var$expr_mean
# set vmr to zero when mean is zero
adata$var$vmr[adata$var$expr_mean==0] <- 0

# computer no-matter what,
tmp_X <- log1p(adata$X)
tmp_mu <- colMeans(tmp_X,na.rm = TRUE)
tmp_vmr <- matrixStats::colVars(tmp_X,na.rm = TRUE)-tmp_mu

adata$var$logvmr <- tmp_vmr

adata$var$vmr_rank <- order(adata$var$logvmr)

# curate a decile
adata$var$vmr_decile <- dplyr::ntile(adata$var$logvmr, 10)

# fix anotation fators
#
facts <- sapply(adata$var, is.factor)
adata$var[facts] <- lapply(adata$var[facts], as.character)

if (FALSE) {
    # save an intermediate file (incase we want to revert...)
    adata$write_h5ad(filename=file.path(DB_ROOT_PATH,DB_NAME,"normalized_data.h5ad"))
}

#7-b. dimension reduction - PCA / umap
zro_na <- adata$copy()
zro_na$X <- adata$layers$get('zro_na')
sc$pp$pca(zro_na)
sc$pp$neighbors(zro_na)
## Step 4: Infer clusters with the leiden algorithm
sc$tl$leiden(zro_na)
## Step 5: Compute tsne and umap embeddings
sc$tl$umap(zro_na)

sc$pp$pca(adata)
sc$pp$neighbors(adata)
sc$tl$leiden(adata)
sc$tl$umap(adata)

adata$obsm$unscaled_X_pca <- zro_na$obsm$X_pca
adata$obsm$unscaled_X_umap <- zro_na$obsm$X_umap
adata$varm$unscaled_PCs <- zro_na$varm$PCs
adata$obsp$unscaled_distances <- zro_na$obsp$distances
adata$obsp$unscaled_connectivities <- zro_na$obsp$connectivities

if (FALSE){
    # save an intermediate file (incase we want to revert...)
    adata$write_h5ad(filename=file.path(DB_ROOT_PATH,DB_NAME,"norm_data_plus_dr.h5ad"))
}

# Step 8: pre-compute differential expression
# save diff expression data for later...
diff_exp <- data_list$de
saveRDS(diff_exp, file.path(DB_ROOT_PATH,DB_NAME, "db_de_table.rds"))

#diff_exp <- readRDS(file = file.path(DB_ROOT_PATH,DB_NAME, "db_de_table.rds"))

if (FALSE){
    #load intermediat data
    adata <- anndata::read_h5ad(filename=file.path(DB_ROOT_PATH,DB_NAME,"norm_data_plus_dr.h5ad"))
}
# Step 9: Write data files to database directory -------
DB_DIR = file.path(DB_ROOT_PATH,DB_NAME)
if (!dir.exists(DB_DIR)) {
  dir.create(DB_DIR)
}

# save diff expression data
diff_exp <- data_list$de

saveRDS(diff_exp, file.path(DB_ROOT_PATH,DB_NAME, "db_de_table.rds"))

# write the anndata object
adata$write_h5ad(filename=file.path(DB_ROOT_PATH,DB_NAME,"db_data.h5ad"))

#reload
if (FALSE) adata <- anndata::read_h5ad(filename=file.path(DB_ROOT_PATH,DB_NAME,"db_data.h5ad"))



# Step 10:  configure browser ----
omic_type <- "prote" #c("transcript","prote","metabol","lipid","other")
aggregate_by_default <- (if (omic_type=="transcript") TRUE else FALSE ) #e.g.  single cell
# choose top 40 proteins by variance across dataset as our "targets"
target_features <- adata$var_names[which(adata$var$var_rank <= 40)]
#if we care we need to explicitly state. defaults will be the order...
config_list <- list(
  # meta-tablel grouping "factors"
  group_obs = c("Condition","leiden","Is.Reference","Replicate"),
  group_var = c("vmr_decile"),

  # LAYERS
  # each layer needs a label/explanation
  layer_values = c("X","raw","X_is_scaled_na_to_0","scaled","zro_na"),
  layer_names = c("Z-scores","intensity (arb)","Z-scores (na=0)","Z-scores","intensity (na=0)" ),

  # ANNOTATIONS / TARGETS
  # what adata$obs do we want to make default values for...
  # # should just pack according to UI?
  default_obs =  c("Condition","leiden"), #subset & ordering

  obs_annots = c("Condition","leiden", "Is.Reference", "expr_var","expr_mean"),

  default_var = character(0), #just use them in order as defined
  var_annots = c("vmr_decile",
                 "expr_geomean",
                 "expr_mean",
                 "expr_frac",
                 "vmr",
                 "vmr_decile"),


  target_features = target_features,
  feature_deets = c("Genes",
                "ProteinDescriptions",
                "ProteinNames",
                "GO.Cellular.Component",
                "GO.Molecular.Function",
                "GO.Biological.Process",
                "expr_geomean",
                "expr_mean",
                "expr_var" ,
                "expr_frac",
                "mean",
                "std",
                "vmr" ,
                "vmr_rank",
                "vmr_decile" ),

  filter_feature = c("vmr"), #if null defaults to "fano_factor"

  # differential expression
  diffs = list( diff_exp_comps = levels(factor(diff_exp$versus)),
                diff_exp_obs_name =  levels(factor(diff_exp$obs_name)),
                diff_exp_tests =  levels(factor(diff_exp$test_type))
              ),

  # Dimension reduction (depricated)
  dimreds = list(obsm = adata$obsm_keys(),
                 varm = adata$varm_keys()),


  #meta info
  annotation_database =  NA,
  publication = "TBD",
  method = "bulk", # c("single-cell","bulk","other")
  omic_type = omic_type, #c("transcript","prote","metabol","lipid","other")
  aggregate_by_default = aggregate_by_default, #e.g.  single cell

  organism = 'mmusculus',
  lab = "Ori/Ward",
  title = "DIA proteomics",
  date = format(Sys.time(), "%a %b %d %X %Y")
)

omicser::write_db_conf(config_list,DB_NAME, db_root = DB_ROOT_PATH)



# BOOTSTRAP the options we have already set up...
# NOTE: we are looking in the "quickstart" folder.  the default is to look for the config in with default getwd()
omicser_options <- omicser::get_config(in_path = OMICSER_RUN_DIR)
DB_ROOT_PATH <- omicser_options$db_root_path

#DB_NAME <- omicser_options$database_names[1]

if (! (DB_NAME %in% omicser_options$database_names)){
  omicser_options$database_names <- c(omicser_options$database_names,DB_NAME)
  omicser::write_config(omicser_options,in_path = OMICSER_RUN_DIR )
}

# Step 11: Run the browser -------------

