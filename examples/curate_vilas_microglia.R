# INFO ========================================================
##  Microglia scRNAseq Transcriptomics
##  - from Vilas:
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
##
#  formely `vilas_B` which was derived from a Seurat file
#
#



DB_NAME <-  list("Vilas Microglia (sceasy)" = "vilas_microglia_sceasy")



require(omicser)
OMICSER_RUN_DIR <- file.path(getwd(),"examples")
setwd(OMICSER_RUN_DIR) # incase we have relative dir
# BOOTSTRAP the options we have already set up...
# NOTE: we are looking in the "quickstart" folder.  the default is to look for the config in with default getwd()

# if you already have an omicser_options.yml for other databases run the following lines
if (FALSE) {
  omicser_options <- omicser::get_config(in_path = OMICSER_RUN_DIR)
  CONDA_ENV <- omicser_options$conda_environment
  DB_ROOT_PATH <- omicser_options$db_root_path
  if (! (DB_NAME %in% omicser_options$database_names)){
    omicser_options$database_names <- c(omicser_options$database_names,DB_NAME)
    omicser::write_config(omicser_options,in_path = OMICSER_RUN_DIR )
  }


} else {
  CONDA_ENV <- "omxr"
  DB_ROOT_PATH <- file.path(OMICSER_RUN_DIR,"test_db") #or as `relative path` = "test_db"

  omicser_options <- list(
    conda_environment = CONDA_ENV,
    db_root_path = DB_ROOT_PATH,
    database_names = DB_NAME)

}



DB_DIR = file.path(DB_ROOT_PATH,DB_NAME)
if (!dir.exists(DB_DIR)) {
  dir.create(DB_DIR)
}





#==== 1. documentation / provenance =========================================================================

db_meta <- list(
  organism = '',
  lab = "Menon",
  annotation_database = NA,
  title = "single-cell microglia data",
  omic_type = "Transcriptomics",
  measurment = "normalized counts",
  pub = "TBD",
  date = format(Sys.time(), "%a %b %d %X %Y")
)


write_db_meta(db_meta,DB_NAME, db_root = DB_ROOT_PATH)



#==== 2. helper functions =========================================================================
# "helper" functions to prep your "raw" data go here or should be sourced here...

#==== 3. load data =========================================================================
RAW_DIR <- file.path(OMICSER_RUN_DIR,"raw_data", "Vilas_B")

# load the og. Seurat Data object
raw_data_file <- "microglia_data_with_counts_RNA_SCT.rda"
load(file.path(RAW_DIR,raw_data_file) )

# pre-pre processing
#  update seurat object -----------------------------
microglia_data_updated <- UpdateSeuratObject(object = microglia_data)
#> OUTPUT:
#   Validating object structure
#   Updating object slots
#   Ensuring keys are in the proper strucutre
#   Ensuring feature names don't have underscores or pipes
#   Object representation is consistent with the most current Seurat version


saveRDS(microglia_data_updated,
        file = file.path(RAW_DIR, "microglia_data_seu_new.rds"))



#==== 4. pack into anndata =========================================================================

data_list <- list(object=file.path(RAW_DIR, "microglia_data_seu_new.rds"))


ad <- omicser::setup_database(database_name=DB_NAME,
                              db_path=DB_ROOT_PATH,
                              data_in=data_list,
                              db_meta=NULL ,
                              re_pack=TRUE)

ad$write_h5ad(filename=file.path(DB_ROOT_PATH,DB_NAME,"core_data.h5ad"))


#==== 5. post processing =========================================================================
sc <- reticulate::import("scanpy")
# re-read the anndata file to instantiate the anndata in R wrappers
ad <- anndata::read_h5ad(filename=file.path(DB_ROOT_PATH,DB_NAME,"core_data.h5ad"))

raw <- ad$raw$copy()

#outvals <- sc$pp$filter_cells(ad, min_genes=200, inplace=FALSE)

sc$pp$filter_cells(ad, min_genes=200)
sc$pp$filter_genes(ad, min_cells=20)


rawX <- raw$X

# FIX RAW
colnames(rawX) <- raw$var_names

rawX <- rawX[which(raw$obs_names %in% ad$obs_names),(raw$var_names %in% ad$var_names)]
raw$obs_names <- raw$obs_names[which(raw$obs_names %in% ad$obs_names)]



newraw <- anndata::AnnData(
  X=rawX,
  var = data.frame(names=colnames(rawX)),
  obs = data.frame(sampleID=rownames(ad$X))

)
ad$raw <- newraw

# to conform with cellXgene scheme gene filtering is tricky.
#
# we need to retain the list / shape of the original genes so we must
# replace zeros/nulls with NA

# outvals <- sc$pp$filter_genes(ad, min_cells=2, inplace=FALSE)

sc$pp$normalize_per_cell(ad, counts_per_cell_after=1e4)
sc$pp$log1p(ad)  # already logg-ed?
sc$pp$highly_variable_genes(ad, min_mean=0.0125, max_mean=3, min_disp=0.5  )
#sc$pp$regress_out(ad_tmp, 'n_counts' )

#sc$pp$highly_variable_genes(ad,n_top_genes=40)


#  don't know how to make this work....
#sc$pp$highly_variable_genes(ad,n_top_genes=40)
ad$var$var_rank <- order(ad$var$sct.residual_variance)
# choose top 40 genes by variance across dataset as our "targets"
target_omics <- ad$var_names[which(ad$var$var_rank <= 40)]


ad$var$decile <- dplyr::ntile(ad$var$sct.residual_variance, 10)



if (FALSE){ # save an intermediate file (incase we want to revert...)
  ad$write_h5ad(filename=file.path(DB_ROOT_PATH,DB_NAME,"normalized_data.h5ad"))
}
#==== 5-a. dimension reduction - PCA / umap  =========================================================================
# ## Step 2: Normalize with a very vanilla recipe
# normalized_data = sc$pp$recipe_seurat(raw, copy=TRUE)

## Step 3: Do some basic preprocessing to run PCA and compute the neighbor graph
sc$pp$pca(ad)
sc$pp$neighbors(ad)

## Step 4: Infer clusters with the Louvain algorithm
#sc$tl$louvain(ad_tmp)
sc$tl$leiden(ad)
## Step 5: Compute tsne and umap embeddings
sc$tl$umap(ad)
#ad_tmp$layers <- list(non_regressed=ad$X) #list('count'=layers)

if (FALSE){ # save an intermediate file (incase we want to revert...)
  ad$write_h5ad(filename=file.path(DB_ROOT_PATH,DB_NAME,"norm_data_plus_dr.h5ad"))
}
#==== 6. differential expression =========================================================================

sc <- import("scanpy")

test_types <- c('wilcoxon','t-test_overestim_var')


comp_types <- c("grpVrest")
obs_names <- c('disease','cell_type')
diff_exp <- omicser::compute_de_table(ad,comp_types, test_types, obs_names,sc)


comp_types <- c("Mild_Cognitive_ImpairmentVAlzheimer_Disease",
                "Mild_Cognitive_ImpairmentVTemporal_Lobe_Epilepsy",
                "Alzheimer_DiseaseVTemporal_Lobe_Epilepsy")
obs_names <- c('disease')
diff_exp2 <- omicser::compute_de_table(ad,comp_types, test_types, obs_names,sc)

diff_exp <- dplyr::bind_rows(diff_exp, diff_exp2)
#diff_exp <- rbind(diff_exp,diff_exp2)
if (FALSE){ # save an intermediate file (incase we want to revert...)
 ad$write_h5ad(filename=file.path(DB_ROOT_PATH,DB_NAME,"norm_data_with_de.h5ad"))
}
saveRDS(diff_exp, file = file.path(DB_ROOT_PATH,DB_NAME, "db_de_table.rds"))

#==== 7. create configs =========================================================================
DB_NAME = "vilas_microglia_sceasy"
if (FALSE) { # load intermediate files if available
  ad <- read_h5ad(filename=file.path("data-raw",DB_NAME,"normalized_data.h5ad"))
  ad <- read_h5ad(filename=file.path("data-raw",DB_NAME,"norm_data_with_de.h5ad"))
  diff_exp <- readRDS( file = file.path("data-raw",DB_NAME, "db_de_table.rds"))
}

# differentials  #if we care we need to explicitly state. defaults will be the order...
config_list <- list(
  x_obs = c("tissue", "disease", "cell_type", "sex", "leiden"),
  y_obs = c("nCount_RNA","nFeature_RNA","nCount_SCT","nFeature_SCT", "n_genes","n_counts"), #MEASURES
  obs_groupby = c("tissue", "disease", "cell_type", "sex", "leiden"),
  obs_subset = c("tissue", "disease", "cell_type", "sex", "leiden"),

  x_var = c("highly_variable","decile"),
  y_var = c("sct.detection_rate", "sct.gmean", "sct.variance","sct.residual_mean","sct.residual_variance", "sct.variable",
            "n_cells",  "means" , "dispersions", "dispersions_norm" ),
  var_groupby = c("highly_variable","decile"),
  var_subset = c("highly_variable","decile"),  # NOTE:  <omic selector> is NOT in the data object so its not actually going to load

  layers = c("X","raw","counts"),

  diffs = list(diff_exp_comps = levels(factor(diff_exp$versus)),
               diff_exp_comp_type =  levels(factor(diff_exp$comp_type)), #i don't think we need this
               diff_exp_obs_name =  levels(factor(diff_exp$obs_name)),
               diff_exp_tests =  levels(factor(diff_exp$test_type))
  ),

  # Dimred
  dimreds = list(obsm = c('X_pca', 'X_tsne'),
                 varm = c('PCs')),

  # what ad$obs do we want to make default values for...
  # # should just pack according to UI?
  default_factors = c("cell_type","disease","tissue"),
  target_omics = target_omics,
  omic_details = c("highly_variable",
                   "decile",
                   "sct.detection_rate",
                   "sct.gmean",
                   "sct.variance",
                   "sct.residual_variance",
                   "sct.variable",
                   "n_cells")

)

omicser::write_db_conf(config_list,DB_NAME, db_root = DB_ROOT_PATH)

#==== 8. write data file to load  =========================================================================

ad$write_h5ad(filename=file.path(DB_ROOT_PATH,DB_NAME,"db_data.h5ad"))



#==== 9. create meta-data narrative for rendering in INGEST  =========================================================================


# ALSO create an `additional_info.Rmd` as in `inst/app/www/additional_info.Rmd`



