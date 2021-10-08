# curate pbmnc3k
#
#
#==== 0. preamble/setup ==================================================

require(omicser)
OMICSER_RUN_DIR <- getwd()

omicser_options <- omicser::get_config(in_path = OMICSER_RUN_DIR)
CONDA_ENV <- omicser_options$conda_environment
DB_ROOT_PATH <- omicser_options$db_root_path


DB_NAME <-  list("caindb" = "caindb")
DB_NAME <- omicser_options$database_names[1]



#==== 1. documentation / provenance ==============================================================
# ANDY: DID YOU DO SOMETHING LIKE THIS?
#
db_meta <- list(
  organism = 'RAT',
  lab = "",
  source = "",
  annotation_database =  "",
  title = "vilas test",
  omic_type = "Transcriptomics",
  measurment = " counts",
  pub = "TBD",
  date = format(Sys.time(), "%a %b %d %X %Y"),
  url = "k"
)

omicser::write_db_meta(db_meta,DB_NAME, db_root = DB_ROOT_PATH)


#==== 4. pack into anndata =========================================================================
# ANDY: I think you did this:
#

TEST_DIR <- file.path(OMICSER_RUN_DIR,"test_db", "caindb")
adata <- anndata::read_h5ad(file.path(TEST_DIR,'db_data.h5ad'))
diff_exp <- readRDS(file = file.path(TEST_DIR, "db_de_table.rds"))



data_list <- list(object=file.path(DB_ROOT_PATH,DB_NAME,'db_data.h5ad'))

adata <- omicser::setup_database(database_name=DB_NAME,
                                 db_path=DB_ROOT_PATH,
                                 data_in=data_list,
                                 db_meta=NULL ,
                                 re_pack=TRUE)


# ANDY: I think you did this:
#
if (FALSE) {
  test_types <- c('wilcoxon')

  comp_types <- c("grpVrest")
  obs_names <- c('cell.type.manual')

  sc <- reticulate::import("scanpy")
  diff_exp <- omicser::compute_de_table(adata,comp_types, test_types, obs_names,sc)

  saveRDS(diff_exp, file = file.path(DB_ROOT_PATH,DB_NAME, "db_de_table.rds"))
}


# ANDY: Did you do this?
#==== 7. create configs =========================================================================
adata <- anndata::read_h5ad(filename=file.path(DB_ROOT_PATH,DB_NAME,"norm_data_with_de.h5ad"))
diff_exp <- readRDS( file = file.path(DB_ROOT_PATH,DB_NAME, "db_de_table.rds"))


# what ad$obs do we want to make default values for...
default_factors <- c("cell.type.manual")

target_omics <- adata$var_names[1:40] # just take the first ones
# differentials  #if we care we need to explicitly state. defaults will be the order...
config_list <- list(
  x_obs = c("cell.type.manual"),
  y_obs =  c('nCount_RNA', 'nFeature_RNA'), #MEASURES
  obs_groupby = c("cell.type.manual"),
  obs_subset = c("cell.type.manual"),

  x_var = c("vartype"),
  # is there a "measurement" in the meta-data?
  y_var = c(""),

  var_groupby = c("vartype"),
  var_subset = c("vartype"),

  diffs = list(diff_exp_comps = levels(factor(diff_exp$versus)),
               diff_exp_comp_type =  levels(factor(diff_exp$comp_type)), #i don't think we need this
               diff_exp_obs_name =  levels(factor(diff_exp$obs_name)),
               diff_exp_tests =  levels(factor(diff_exp$test_type))
  ),

  layers = c("X","raw"),

  # Dimred
  dimreds = list(obsm = adata$obsm_keys(),
                 varm = adata$varm_keys()),

  # what ad$obs do we want to make default values for...
  # # should just pack according to UI?
  default_factors = default_factors,
  target_omics = target_omics,
  omic_details = c('sample_ID', 'nCount_RNA', 'nFeature_RNA', 'projid', 'orig.ident', 'batch', 'sex', 'pAD', 'Cdx', 'age_death', 'apoe_genotype', 'Experiment', 'cogn_global_random_slope', 'nft_sqrt', 'gpath_sqrt', 'tangles_sqrt', 'amyloid_sqrt', 'diabetes_sr_rx_ever', 'region', 'latent_cell_probability', 'RNA_snn_res.5', 'tree.ident', 'new.cell', 'prev.cell.type', 'prev.cell.sub.type', 'prev.annotated.doublet', 'enriched.cells.de', 'doublet.score', 'high.cell.de.cluster', 'is.doublet', 'is.doublet.text', 'is.doublet.expanded', 'is.doublet.expanded.text', 'RNA_snn_res.1.5', 'RNA_snn_res.1', 'RNA_snn_res.0.7', 'cell.type', 'high.res.clust.cell.level.raw', 'high.res.clust.cell.level.clean', 'cell.type.manual', 'cell.subtype.manual', 'oligodendrocytes.omega.1', 'oligodendrocytes.omega.2', 'oligodendrocytes.omega.3', 'oligodendrocytes.omega.4')

)



omicser::write_db_conf(config_list,DB_NAME, db_root = DB_ROOT_PATH)

#==== 8. write data file to load  =========================================================================
adata$write_h5ad(filename=file.path(DB_ROOT_PATH,DB_NAME,"db_data.h5ad"))
















# ALTERNATE SEURAT STYLE =================
#'
#' # FROM SEURAT::
#' #' Read 10X hdf5 file
#' #'
#' #' Read count matrix from 10X CellRanger hdf5 file.
#' #' This can be used to read both scATAC-seq and scRNA-seq matrices.
#' #'
#' #' @param filename Path to h5 file
#' #' @param use.names Label row names with feature names rather than ID numbers.
#' #' @param unique.features Make feature names unique (default TRUE)
#' #'
#' #' @return Returns a sparse matrix with rows and columns labeled. If multiple
#' #' genomes are present, returns a list of sparse matrices (one per genome).
#' #'
#' #' @export
#' #' @concept preprocessing
#' #'
#' Read10X_h5 <- function(filename, use.names = TRUE, unique.features = TRUE) {
#'   if (!requireNamespace('hdf5r', quietly = TRUE)) {
#'     stop("Please install hdf5r to read HDF5 files")
#'   }
#'   if (!file.exists(filename)) {
#'     stop("File not found")
#'   }
#'   infile <- hdf5r::H5File$new(filename = filename, mode = 'r')
#'   genomes <- names(x = infile)
#'   output <- list()
#'   if (hdf5r::existsGroup(infile, 'matrix')) {
#'     # cellranger version 3
#'     if (use.names) {
#'       feature_slot <- 'features/name'
#'     } else {
#'       feature_slot <- 'features/id'
#'     }
#'   } else {
#'     if (use.names) {
#'       feature_slot <- 'gene_names'
#'     } else {
#'       feature_slot <- 'genes'
#'     }
#'   }
#'   for (genome in genomes) {
#'     counts <- infile[[paste0(genome, '/data')]]
#'     indices <- infile[[paste0(genome, '/indices')]]
#'     indptr <- infile[[paste0(genome, '/indptr')]]
#'     shp <- infile[[paste0(genome, '/shape')]]
#'     features <- infile[[paste0(genome, '/', feature_slot)]][]
#'     barcodes <- infile[[paste0(genome, '/barcodes')]]
#'     sparse.mat <- sparseMatrix(
#'       i = indices[] + 1,
#'       p = indptr[],
#'       x = as.numeric(x = counts[]),
#'       dims = shp[],
#'       giveCsparse = FALSE
#'     )
#'     if (unique.features) {
#'       features <- make.unique(names = features)
#'     }
#'     rownames(x = sparse.mat) <- features
#'     colnames(x = sparse.mat) <- barcodes[]
#'     sparse.mat <- as(object = sparse.mat, Class = 'dgCMatrix')
#'     # Split v3 multimodal
#'     if (infile$exists(name = paste0(genome, '/features'))) {
#'       types <- infile[[paste0(genome, '/features/feature_type')]][]
#'       types.unique <- unique(x = types)
#'       if (length(x = types.unique) > 1) {
#'         message("Genome ", genome, " has multiple modalities, returning a list of matrices for this genome")
#'         sparse.mat <- sapply(
#'           X = types.unique,
#'           FUN = function(x) {
#'             return(sparse.mat[which(x = types == x), ])
#'           },
#'           simplify = FALSE,
#'           USE.NAMES = TRUE
#'         )
#'       }
#'     }
#'     output[[genome]] <- sparse.mat
#'   }
#'   infile$close_all()
#'   if (length(x = output) == 1) {
#'     return(output[[genome]])
#'   } else{
#'     return(output)
#'   }
#' }
# pbmc_data <- Read10X(data.dir = "raw_data/pbmc3k/filtered_gene_bc_matrices/hg19/")
