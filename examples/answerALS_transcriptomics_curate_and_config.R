#
# Overview --------------
#### Create an app to browse data from Answer ALS database.
####
#### WARNING::: WIP
####
require("reticulate")

DB_NAME <- list("answer ALS transcriptomics" = "answr_als_transcript") # name our database

#-#_#_#_#_#_#_#_#_#_#_#_#_#_#_#__#_#_#_#_#_#_
#  Step 1: Set paths--------------
OMICSER_RUN_DIR <- getwd() # /path/to/cloned/omicser/examples or just where you run from
OMICSER_RUN_DIR <- file.path(getwd(),"examples") # /path/to/cloned/omicser/examples or just where you run from

#RAW_DATA_DIR <- file.path(OMICSER_RUN_DIR,"raw_data") # set the path for where the raw_data lives...
                                                      # here its going to be in our OMCISER_RUN_DIR
RAW_DATA_DIR <- file.path(Sys.getenv("HOME"),"omicser_databases","answerALS")


if (!dir.exists(RAW_DATA_DIR)) {
  dir.create(RAW_DATA_DIR) #fails if the path has multiple levels to generate
}

# set the path for where the databases live... here its going to be in our OMCISER_RUN_DIR
DB_ROOT_PATH <- file.path(OMICSER_RUN_DIR,"databases")
DB_ROOT_PATH <- "/Users/ahenrie/Projects/NDCN_dev/omicser/quickstart/test_db"

if (!dir.exists(DB_ROOT_PATH)) {
  dir.create(DB_ROOT_PATH)
}

OMICSER_PYTHON <-  "pyenv_omicser"
# installation type (see install_script.R)

#PROTEOMICS
# Assay
# Upon receival, the iPSC pellets were lyophilized and stored in -80 °C. A small portion of the pellet was taken in a 1.5ul microfuge tube and urea lysis buffer (8 M urea, 75 mM NaCl, 50 mM Tris-HCl, pH 8.3) was added in a ratio of 1:2 (pellet to buffer). The resuspended pellets were sonicated at 40% amplitude until the pellets were completely homogenized. Following homogenization, the samples were centrifuged for 5 minutes at 5000g at 4 °C and the supernatant was transferred into a new tube. Protein BCA assay was performed as per manufacturer’s instructions (Pierce™ BCA Protein Assay Kit, Catalog number:  23225). The volume required for 100ug protein was taken from the cell lysate. Protein was reduced using 10mM (final concentration) TCEP (tris(2-carboxyethyl)phosphine) at 56°C for 30 minutes followed by alkylation using 12mM (final concentration) Iodoacetamide (IAA) at room temperature for 50 minutes in dark. Urea concentration was brought down to 1.6M by the addition of 50 mM Tris-HCl, pH 8.3 and subsequently, enzyme trypsin (at 1:50 enzyme to the substrate) was used to digest the proteins. The enzyme-substrate mixture was incubated overnight at 37°C while on a shaker at 1400 rpm. The following day, the reaction was quenched by the addition of 10%Formic acid (1% final concentration) and the pH was ensured to be between 2-3. The sample was then centrifuged at 5000g at RT and the supernatant was transferred to a new tube. Sample desalting was performed using an oasis 96 well HLB plate following the manufacturer’s instructions (Waters, SKU: 186000309). Desalted samples were dried down in a speed vacuum and the dried peptides were reconstituted in 0.1% Formic acid in water. Peptide BCA was performed to obtain the accurate peptide concentration in each sample.
# Profile-mode .wiff files from the shotgun data acquisition were converted to mzML format using the AB Sciex Data Converter (in proteinpilot mode) and then re-converted to mzXML format using ProteoWizard for peaklist generation. The MS2 spectra were queried against the reviewed canonical Swiss-Prot Human complete proteome database appended with iRT protein sequence and shuffled sequence decoys. All data were searched using the X!Tandem and Comet. The search parameters included the following criteria: static modifications of Carbamidomethyl (C) and variable modifications of Oxidation (M). The parent mass tolerance was set to be 50 p.p.m, and mono-isotopic fragment mass tolerance was 100 p.p.m (which was further filtered to be < 0.05 Da for building spectral library); tryptic peptides with up to three missed cleavages were allowed. The identified peptides were processed and analyzed through Trans-Proteomic Pipeline and were validated using the PeptideProphet scoring and the PeptideProphet results were statistically refined using iProphet. All the peptides were filtered at a false discovery rate (FDR) of 1% with peptide probability cutoff >=0.99. The raw spectral libraries were generated from all valid peptide spectrum matches and then refined into the non-redundant consensus libraries using SpectraST. For each peptide, the retention time was mapped into the iRT space with reference to a linear calibration constructed for each shotgun run. The MS assays, constructed from Top six most intense transitions (from ion series: b and y and charge states: 1,2) with Q1 range from 400 to 1,200 m/z excluding the precursor SWATH window, were used for targeted data analysis of SWATH maps.
# SWATH-MS.wiff files from the data-independent acquisition were first converted to profile mzML using ProteoWizard v.3.0.6002. The whole process of SWATH-targeted data analysis was carried out using OpenSWATH running on an internal computing cluster. OpenSWATH utilizes a target-decoy scoring system such as mProphet to estimate the identification of FDR. The best scoring classifier that was built from the sample of most protein identifications was utilized in this study. Based on our final spectral library, OpenSWATH firstly identified the peak groups from all individual SWATH maps at a global peptide FDR of 1% and aligned them between SWATH maps based on the clustering behaviors of retention time in each run with a non-linear alignment algorithm. For this analysis, the MS runs were realigned to each other using LOcally WEighted Scatterplot Smoothing method and the peak group clustering was performed using ‘LocalMST’ method. Specifically, only those peptide peak groups that deviate within 3 standard deviations from the retention time were reported and considered for alignment with the max FDR quality of 5% (quality cutoff to still consider a feature for alignment). Next, to obtain high-quality quantitative data at the protein level, we discarded those proteins whose peptides were shared between multiple different proteins (non-proteotypic peptides). Data pre-processing and statistical analysis of MS runs into quantitative data was performed using mapDIA. The fragment-level intensities were normalized based on TIC (Total Ion Current) from each DIA MS run to remove systematic bias between MS runs. This is followed by outlier removal and peptide/fragment selection that preserves the major quantitative patterns across all samples for each protein. The selected fragments and peptides are used in the final model-based statistical significance analysis of protein-level differential expression between specified groups of samples.

#TRANSCRIPTOMICS
# Assay
# The RNA-seq library preparation and sequencing have been done at the Genomics High Throughput Facility at the University of California, Irvine. For each sample, total RNA was isolated using the Qiagen RNeasy mini kit. RNA samples for each AALS subject (control or ALS) were entered into an electronic tracking system and processed.  RNA QC was conducted using an Agilent Bioanalyzer and Nanodrop. Our primary QC metric for RNA quality is based on RIN values (RNA Integrity Number) ranging from 0-10, 10 being the highest quality RNA. Additionally, we collected QC data on total RNA concentration and 260/280 and 260/230 ratios to evaluate any potential contamination. Only samples with RIN >8 were used for library prep and sequencing. rRNAs were removed and libraries generated using TruSeq Stranded Total RNA library prep kit with Ribo-Zero (Qiagen). RNA-Seq libraries were titrated by qPCR (Kapa), normalized according to size (Agilent Bioanalyzer 2100 High Sensitivity chip). Each cDNA library was then subjected to 100 Illumina (Novaseq 6000) paired-end (PE) sequencing cycles to obtain over 50 million PE reads per sample.
# Data Pipeline
# The quality of the sequencing was first assessed using fastQC tool (v.0.11.9). Raw reads were then adapter and quality trimmed and filtered by a length of 20 bases. Trimmed reads were mapped to the GRCh38 reference genome using Hisat2 (v.2.2.1) and resulting BAM files were indexed using samtools (v.1.10), and gene expression was quantified with featureCounts (subread v.2.0.1) as raw read count files.

# Step 2: Assert python back-end --------------------------------
#  for the curation we need to have scanpy
CONDA_INSTALLED <- reticulate:::miniconda_exists()
OMICSER_PYTHON_EXISTS <- any(reticulate::conda_list()["name"]==OMICSER_PYTHON)

if (!CONDA_INSTALLED){  #you should already have installed miniconda and created the env
  reticulate::install_miniconda() #in case it is not already installed
  }


if (!OMICSER_PYTHON_EXISTS){  #you should already have installed miniconda and created the env
  # simpler pip pypi install
  packages <- c("scanpy[leiden]")
  reticulate::conda_create(OMICSER_PYTHON, python_version = 3.8)
  reticulate::conda_install(envname=OMICSER_PYTHON,
                            # channel = "conda-forge",
                            pip = TRUE,
                            packages =  packages )

}

if ( Sys.getenv("RETICULATE_PYTHON") != "OMICSER_PYTHON" ) {
  Sys.setenv("RETICULATE_PYTHON"=reticulate::conda_python(envname = OMICSER_PYTHON))
}


# check that we have our python on deck
reticulate::py_discover_config()



# step2b:  troublshoot conda install--------
# full conda install
# packages1 <- c("seaborn", "scikit-learn", "statsmodels", "numba", "pytables")
# packages2 <- c("python-igraph", "leidenalg")
#
# reticulate::conda_create(OMICSER_PYTHON, python_version = 3.8,packages = packages1)
# reticulate::conda_install(envname=OMICSER_PYTHON,
#                         channel = "conda-forge",
#                         packages = packages2 )
# reticulate::conda_install(envname=OMICSER_PYTHON,
#                           channel = "conda-forge",
#                           pip = TRUE,
#                           packages = "scanpy" )

# if (!reticulate::py_module_available(module = "scanpy") ) {
#
#   reticulate::conda_install(envname=OMICSER_PYTHON,
#                             packages = "scanpy[leiden]",
#                             pip = TRUE,
#                             conda=reticulate::conda_binary())
#
#
# }

# reticulate::use_condaenv(condaenv = OMICSER_PYTHON,
#                          conda = reticulate::conda_binary(),
#                          required = TRUE)

# reticulate::conda_install(envname=OMICSER_PYTHON,
#                           packages = "leidenalg",
#                           pip = TRUE,
#                           conda=reticulate::conda_binary())



# Step 4:  get the data ---------------
# create directory structure for data and databases
DB_DIR = file.path(DB_ROOT_PATH,DB_NAME)
if (!dir.exists(DB_DIR)) {
  dir.create(DB_DIR)
}




# Step 5:  define for source helper functions -------------------------
# N/A
make_meta_table <- function(RAW_DATA_DIR,mtable,clin_dict,meta_to_agg){
  #forms <- unique(clin_dict[,Form_Name])
  forms <- meta_to_agg
  key_list <- list()
  out_table <- mtable
  for (form in forms) {
    key_list[[form]] <- clin_dict[clin_dict[,Form_Name==form] ]
    tbl <- data.table::fread(file.path(RAW_DATA_DIR,"metadata/clinical",paste0(form,".csv")))

    if ("SubjectUID" %in% colnames(tbl)){
      tbl[,SubjectUID:=NULL]
    }

    if ("Form_Name" %in% colnames(tbl)){
      tbl[,Form_Name:=NULL]
    }

    if ("Visit_Name" %in% colnames(tbl)){
      tbl[,Visit_Name:=NULL]
    }
    if ("Visit_Date" %in% colnames(tbl)){
      tbl[,Visit_Date:=NULL]
    }

    if (dim(tbl)[1]<=dim(out_table)[1]){
      out_table <- merge(out_table,tbl,on="Participant_ID",all = TRUE)
    }
  }

  return (out_table)
}



prep_ANSALS_files <- function(proteomics_data_file,transcriptomics_data_file,meta_data_file,clinical_file,path_root,trans_or_prot){
  #trans_or_prot <- c("trans", "prot")[trans_or_prot]  #i.e. trans=1, prot=2
  # Step 6: Format and ingest raw data
  ptable <- data.table::fread(file.path(path_root,proteomics_data_file))
  ttable <- data.table::fread(file.path(path_root,transcriptomics_data_file))
  mtable <- data.table::fread(file.path(path_root,"metadata",meta_data_file))
  clin_dict <- data.table::fread(file.path(path_root,"metadata",clinical_file))

  meta_to_agg <- c( "AALSDXFX", "AALSHXFX",
                    "ALSFRS_R" ,
                    #"ALS_CBS",
                    "ALS_Gene_Mutations",
                    "ANSASFD",
                    #"ANSWER_ALS_Medications_Log",
                    #"ANSWER_ALS_MobileApp",
                    "Ashworth_Spasticity_Scale",
                    #"Auxiliary_Chemistry",
                    #"Auxiliary_Chemistry_Labs",
                    #"CNS_Lability_Scale",
                    #"Cerebrospinal_Fluid",
                    "DNA_Sample_Collection",
                    "Demographics"
                    #"Diaphragm_Pacing_System_Device",
                    #"Family_History_Log" ,"Feeding_Tube_Placement",
                    #"Grip_Strength_Testing","Hand_Held_Dynamometry", "Medical_History",
                    #"Mortality" ,"NEUROLOG","NIV_Log",
                    #"PBMC_Sample_Collection" , "Permanent_Assisted_Ventilation","Plasma_Sample",
                    #"Reflexes", "Serum_Sample", "Tracheostomy",
                    #"Vital_Capacity", "Vital_Signs","subjects"
  )


  meta_tab <- make_meta_table(path_root,mtable,clin_dict,meta_to_agg)

  proteins <- ptable$Protein
  #p_subjID <- colnames(ptable)[-1]
  p_subjID <- substring(colnames(ptable)[-1],1,16)

  qProt <- (meta_tab$Participant_ID %in% p_subjID)
  colnames(ptable)<-c("Proteins",p_subjID)
  pdata_mat <- as.matrix(ptable,rownames=1)


  ensembleIDs <- ttable$V1
  #p_subjID <- colnames(ptable)[-1]
  t_subjID <- substring(colnames(ttable)[-1],1,16)
  qTrans <- (meta_tab$Participant_ID %in% t_subjID)

  colnames(ttable)<-c("ensembleIDs",t_subjID)


  tdata_mat <- as.matrix(ttable,rownames=1)

  #tdata_mat <- as(tdata_mat, "dgTMatrix")
  # data1 <- as.matrix(tdata_mat) #74.8MB
  # data2 <- as(data1, "dgeMatrix") #141MB
  # data3 <- as(data1, "dgTMatrix") #151MB


  p_omics <- rownames(pdata_mat)
  t_omics <- rownames(tdata_mat)


  # add some marginal statistics
  if (trans_or_prot==1) {
    tmp_mat <- tdata_mat
    sample_ID <- t_subjID
    features <- t_omics
  } else {
    tmp_mat <- pdata_mat
    sample_ID <- p_subjID
    features <- p_omics
  }
  tmp_mat[is.na(tmp_mat)] <- 0 #defensive

  obs_meta <-  meta_tab[sample_ID,]# meta_tab[meta_tab$Participant_ID %in% sample_ID,]

  obs_meta$expr_geomean <- Matrix::colMeans( log1p(tmp_mat),na.rm = TRUE) #exp minus 1?
  obs_meta$expr_mean <- Matrix::colMeans(tmp_mat,na.rm = TRUE)
  obs_meta$expr_var  <- matrixStats::colVars(tmp_mat,na.rm = TRUE)
  obs_meta$expr_frac <- Matrix::colMeans(tmp_mat>0)


  var_annot <- data.table::data.table(features=features)

  var_annot$expr_var <- matrixStats::rowVars(tmp_mat,na.rm=TRUE)
  var_annot$expr_mean <- Matrix::rowMeans(tmp_mat,na.rm=TRUE)
  var_annot$expr_frac <- Matrix::rowMeans(tmp_mat>0)


  # TODO:  strategy for mapping t_omics and p_omics to ids...

  data_list <- list(data_mat = data_mat,
                    obs_meta = obs_meta,
                    var_annot =var_annot,
                    omics = features,
                    sample_ID = sample_ID,
                    etc = NULL)

  return(data_list)
}


  # we have 1145 "subjects", but transcriptomics for just <300 adn proteomics for ~200


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



# Step 6: load helper tools via the "omicser" browser package ---------
CLONED_OMICSER <- TRUE
if ( CLONED_OMICSER ) {
  require("golem")
  REPO_DIR <- getwd()
  golem::document_and_reload(pkg = REPO_DIR)
} else {
  require("omicser")
  #see install_script.R if not installed
}

# Step 6: Format and ingest raw data

# change paths to make data manipulations easier
#

#
proteomics_data_file <- "AnswerALS-204-P-v1-release5_protein-level-matrix.csv"
#
transcriptomics_data_file <- "AnswerALS-290-T-v1-release5_rna-counts-matrix.csv"

meta_data_file <- "aals_participants.csv"
#meta_data_file <- "aals_dataportal_datatable.csv"

clinical_file <- "Clinical_Dictionary.csv"
#
# ptable <- data.table::fread(file.path(RAW_DATA_DIR,proteomics_data_file))
# ttable <- data.table::fread(file.path(RAW_DATA_DIR,transcriptomics_data_file))
# mtable <- data.table::fread(file.path(RAW_DATA_DIR,"metadata",meta_data_file))
#
# clin_dict <- data.table::fread(file.path(RAW_DATA_DIR,"metadata",clinical_file))

# process proteomics & transcriptomics together for meta-info, but only use one..


# Steps 7-9: CURATION
SAVE_INTERMEDIATE_FILES <- FALSE
# Step 7:  pack data into AnnData format --------------
# identify location of raw data
trans_or_prot = 1 #transcriptomics
data_list <- prep_ANSALS_files(proteomics_data_file,transcriptomics_data_file,meta_data_file,clinical_file, RAW_DATA_DIR,trans_or_prot)

# create database formatted as AnnData
adata <- omicser::setup_database(database_name = DB_NAME,
                                 db_path = DB_ROOT_PATH,
                                 data_in = data_list,
                                 re_pack = TRUE)



# Step 8: additional data processing ----
adata$var_names_make_unique()

# use scanpy to do some scaling and calculations...
sc <- reticulate::import("scanpy")

# filter data
sc$pp$filter_cells(adata, min_genes=200)
sc$pp$filter_genes(adata, min_cells=10)


sc$pp$normalize_total(adata, target_sum=1e4)
sc$pp$log1p(adata)
sc$pp$highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)

# transform data
#sc$pp$scale(adata, max_value=10)

# choose top 40 genes by variance across dataset as "targets"
adata$var$var_rank <- order(adata$var$dispersions_norm)

# calculate deciles
adata$var$decile <- dplyr::ntile(adata$var$dispersions_norm, 10)
#raw <- ad$raw$to_adata()

# save intermediate database file
if (SAVE_INTERMEDIATE_FILES){
  adata$write_h5ad(filename=file.path(DB_ROOT_PATH,DB_NAME,"normalized_data.h5ad"))
}

#7-b. dimension reduction - PCA / umap
#pca
sc$pp$pca(adata)
# compute neighbor graph
sc$pp$neighbors(adata)
## infer clusters
sc$tl$leiden(adata)
# compute umap
sc$tl$umap(adata)

# save intermediate database file
if (SAVE_INTERMEDIATE_FILES){
  adata$write_h5ad(filename=file.path(DB_ROOT_PATH,DB_NAME,"norm_data_plus_dr.h5ad"))
}

# Step 8: pre-compute differential expression
# identify stats
# see scanpy documentation for possible stat test choices
test_types <- c('wilcoxon')
comp_types <- c("grpVrest")
obs_names <- c('leiden')
# calculate DE
diff_exp <- omicser::compute_de_table(adata,comp_types, test_types, obs_names,sc)

###
# save intermediate database file
if (SAVE_INTERMEDIATE_FILES){
  adata$write_h5ad(filename=file.path(DB_ROOT_PATH,DB_NAME,"norm_data_with_de.h5ad"))
}

# save DE tables
saveRDS(diff_exp, file = file.path(DB_ROOT_PATH, DB_NAME, "db_de_table.rds"))


# Step 9: Write data files to database directory -----------
# write final database
adata$write_h5ad(filename = file.path(DB_ROOT_PATH, DB_NAME, "db_data.h5ad"))

# set to TRUE and restart from here for re-configuring
if (FALSE) {
  adata <- anndata::read_h5ad(filename=file.path(DB_ROOT_PATH,DB_NAME,"db_data.h5ad"))
  diff_exp <- readRDS( file = file.path(DB_ROOT_PATH,DB_NAME, "db_de_table.rds"))
}
if (FALSE) adata <- anndata::read_h5ad(filename=file.path(DB_ROOT_PATH,DB_NAME,"db_data.h5ad"))



# Step 10:  configure browser ----
omic_type <- "transcript" #c("transcript","prote","metabol","lipid","other")
aggregate_by_default <- (if (omic_type=="transcript") TRUE else FALSE ) #e.g.  single cell
# choose top 40 proteins by variance across dataset as our "targets"
target_features <- adata$var_names[which(adata$var$var_rank <= 40)]
#if we care we need to explicitly state. defaults will be the order...
config_list <- list(
  # meta-tablel grouping "factors"
  group_obs = c("leiden"),
  group_var = c("decile","highly_variable"),

  # LAYERS
  # each layer needs a label/explanation
  layer_values = c("X","raw"),
  layer_names = c("norm-count","counts" ),

  # ANNOTATIONS / TARGETS
  # what adata$obs do we want to make default values for...
  # # should just pack according to UI?
  default_obs =  c("Case_Control","leiden"), #subset & ordering

  obs_annots = c( "leiden", "n_genes","n_genes_by_counts","total_counts","Sex","Diagnosis","Enrollment Status"),

  default_var = c("decile"),#just use them in order as defined
  var_annots = c(
    "highly_variable",
    "dispersions_norm",
    "decile"),


  target_features = target_features,
  feature_details = c( "feature_name",
                     "n_cells",
                     "total_counts",
                     "highly_variable",
                     "means",
                     "dispersions",
                     "dispersions_norm",
                     "var_rank",
                     "decile" ),

  filter_feature = c("dispersions_norm"), #if null defaults to "fano_factor"
  # differential expression
  diffs = list( diff_exp_comps = levels(factor(diff_exp$versus)),
                diff_exp_obs_name =  levels(factor(diff_exp$obs_name)),
                diff_exp_tests =  levels(factor(diff_exp$test_type))
  ),

  # Dimension reduction (depricated)
  dimreds = list(obsm = adata$obsm_keys(),
                 varm = adata$varm_keys()),

  omic_type = omic_type, #c("transcript","prote","metabol","lipid","other")
  aggregate_by_default = aggregate_by_default, #e.g.  single cell

  #meta info
  meta_info = list(
    annotation_database =  NA,
    publication = "TBD",
    method = "single-cell", # c("single-cell","bulk","other")
    organism = "human",
    lab = "?",
    source = "answer ALS transcriptomic",
    title = "answer als transcriptomic",
    measurment = "normalized counts",
    pub = "",
    url = "https://answerALS",
    date = format(Sys.time(), "%a %b %d %X %Y")
  )
)

omicser::write_db_conf(config_list,DB_NAME, db_root = DB_ROOT_PATH)


# BOOTSTRAP the options we have already set up...
# NOTE: we are looking in the "quickstart" folder.  the default is to look for the config in with default getwd()
omicser_options <- omicser::get_config(in_path = OMICSER_RUN_DIR)
#omicser_options <- omicser::get_config()
DB_ROOT_PATH_ <- omicser_options$db_root_path
if (DB_ROOT_PATH_==DB_ROOT_PATH){
  # add the database if we need it...
  if (! (DB_NAME %in% omicser_options$database_names)){
    omicser_options$database_names <- c(omicser_options$database_names,DB_NAME)
  }

} else {
  omicser_options$db_root_path <- DB_ROOT_PATH
  if (any(omicser_options$database_names == "UNDEFINED")) {
    omicser_options$database_names <- DB_NAME
  } else {
    omicser_options$database_names <- c(omicser_options$database_names,DB_NAME)
  }
}


# write the configuration file
omicser::write_config(omicser_options,in_path = OMICSER_RUN_DIR )


# Step 11: Run the browser -------------
omicser::run_defaults()
