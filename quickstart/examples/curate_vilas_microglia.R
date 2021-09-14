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

#==== 0. preamble/setup ==================================================
# assume we are in the [omicser_path]
# getwd()
# pkgload::load_all('.')
require(golem)
golem::document_and_reload()


# BOOTSTRAP the options we have already set up...
# NOTE: we are looking in the "quickstart" folder.  the default is to look for the config in with default getwd()
omxr_options <- omicser::get_config(in_path="quickstart")


CONDA_ENV <- omxr_options$conda_environment
DB_NAME <- omxr_options$database_names[3]
DB_ROOT_PATH <- omxr_options$db_root_path



#DB_NAME = "vilas_microglia_sceasy"
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


write_db_meta(DB_NAME, db_root = DB_ROOT_PATH)



#==== 2. helper functions =========================================================================
# "helper" functions to prep your "raw" data go here or should be sourced here...

#==== 3. load data =========================================================================
RAW_DIR <- "raw_data/Vilas_B"

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
microglia_data_updated <- readRDS(file = file.path(RAW_DIR, "microglia_data_seu_new.rds"))

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
sc$pp$filter_genes(ad, min_cells=3)


# FIX RAW
colnames(raw$X) <- raw$var_names
rownames(raw$X) <- raw$obs_names

raw <- raw[(raw$obs_names %in% ad$obs_names),(raw$var_names %in% ad$var_names)]

#


# to conform with cellXgene scheme gene filtering is tricky.
#
# we need to retain the list / shape of the original genes so we must
# replace zeros/nulls with NA

# outvals <- sc$pp$filter_genes(ad, min_cells=2, inplace=FALSE)

sc$pp$normalize_per_cell(ad, counts_per_cell_after=1e4)
sc$pp$log1p(ad)  # already logg-ed?
sc$pp$highly_variable_genes(ad, min_mean=0.0125, max_mean=3, min_disp=0.5  )
#sc$pp$regress_out(ad_tmp, 'n_counts' )

sc$pp$highly_variable_genes(ad,n_top_genes=40)

ad$raw <- raw

#  don't know how to make this work....
#sc$pp$highly_variable_genes(ad,n_top_genes=40)
ad$var$var_rank <- order(ad$var$expr_var)
# choose top 40 genes by variance across dataset as our "targets"
target_omics <- ad$var_names[which(ad$var$var_rank <= 40)]

# save an intermediate file (incase we want to revert...)
ad$write_h5ad(filename=file.path(DB_ROOT_PATH,DB_NAME,"normalized_data.h5ad"))

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

ad$write_h5ad(filename=file.path("data-raw",DB_NAME,"norm_data_plus_dr.h5ad"))


#==== 6. differential expression =========================================================================

sc <- import("scanpy")
helper_function<-('data-raw/compute_de_table.R')
source(helper_function)

test_types <- c('wilcoxon','t-test_overestim_var')


comp_types <- c("grpVrest")
obs_names <- c('disease','cell_type')
diff_exp <- compute_de_table(ad,comp_types, test_types, obs_names)


comp_types <- c("Mild_Cognitive_ImpairmentVAlzheimer_Disease",
                "Mild_Cognitive_ImpairmentVTemporal_Lobe_Epilepsy",
                "Alzheimer_DiseaseVTemporal_Lobe_Epilepsy")
obs_names <- c('disease')
diff_exp2 <- compute_de_table(ad,comp_types, test_types, obs_names)

diff_exp <- dplyr::bind_rows(diff_exp, diff_exp2)


ad$write_h5ad(filename=file.path("data-raw",DB_NAME,"norm_data_with_de.h5ad"))
saveRDS(diff_exp, file = file.path("data-raw",DB_NAME, "diff_expr_table.rds"))

#==== 7. create configs =========================================================================
DB_NAME = "vilas_microglia_sceasy"

ad <- read_h5ad(filename=file.path("data-raw",DB_NAME,"normalized_data.h5ad"))

ad <- read_h5ad(filename=file.path("data-raw",DB_NAME,"norm_data_with_de.h5ad"))
diff_exp <- readRDS( file = file.path("data-raw",DB_NAME, "diff_expr_table.rds"))


# differentials  #if we care we need to explicitly state. defaults will be the order...
conf_list <- list(
  x_obs = c("tissue", "disease", "cell_type", "sex", "leiden"),
  y_obs = c("nCount_RNA","nFeature_RNA","nCount_SCT","nFeature_SCT", "n_genes","n_counts"), #MEASURES
  obs_groupby = c("tissue", "disease", "cell_type", "sex", "leiden"),
  obs_subset = c("tissue", "disease", "cell_type", "sex", "leiden"),

  x_var = c("highly_variable"),
  y_var = c("sct.detection_rate", "sct.gmean", "sct.variance","sct.residual_mean","sct.residual_variance", "sct.variable",
            "n_cells",  "means" , "dispersions", "dispersions_norm" ),
  var_groupby = c("highly_variable"),
  var_subset = c("highly_variable"),  # NOTE:  <omic selector> is NOT in the data object so its not actually going to load

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
  target_omics = ad$var_names[1:40]
)

configr::write.config(config.dat = conf_list, file.path = file.path("data-raw",DB_NAME,"config.yml" ),
                      write.type = "yaml", indent = 4)


#==== 8. write data file to load  =========================================================================
ad$write_h5ad(filename=file.path("data-raw",DB_NAME,"omxr_data.h5ad"))




#==== DONE- =========================================================================


ad <- read_h5ad(filename=file.path("data-raw",DB_NAME,"omxr_data.h5ad"))
conf_list_out <- configr::read.config( file.path("data-raw",DB_NAME,"config.yml" ) )
