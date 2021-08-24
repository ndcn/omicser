#######################################################################
#######################################################################
##
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
#######################################################################
#######################################################################
#  formely `vilas_B` which was derived from a Seurat file
#
#

# ------------------------------------------
# 0. preamble/setup -------------------------
# ------------------------------------------


db_meta <- list(
  organizm = '',
  lab = "Menon",
  annotation_database = NA
)

RAW_DIR <- "ingest/Vilas_B"
# Seurat Data object
raw_data_file <- "microglia_data_with_counts_RNA_SCT.rda"

# pre-pre processing
#update seurat object -----------------------------
microglia_data_updated <- UpdateSeuratObject(object = microglia_data)
#> OUTPUT:
#   Validating object structure
#   Updating object slots
#   Ensuring keys are in the proper strucutre
#   Ensuring feature names don't have underscores or pipes
#   Object representation is consistent with the most current Seurat version

saveRDS(microglia_data_updated,
        file = file.path(RAW_DIR, "microglia_data_seu_new.rds"))


data_list <- list(object="microglia_data_seu_new.rds")

# scemaa
data_list <- c("object", "data_mat","obs_meta","var_annot","omics","sample_ID","etc")




omxr_pack_anndata <- function (data_in){

  tools::file_path_sans_ext(data_in)
  if (class(data_in)[1] == "list") {
    # multple files in c("object", "data_mat","obs_meta","var_annot","omics","sample_ID",etc")
    # data_mat - data matrix.  need to assert that this matrix object has omics and sample+ID names



    # obs_meta - ensure that we are factored and sample_ID is first column


    # var_annot - ensure that "omics" is first column


    # etc goes into an uns entry

  } else if (tolower(tools::file_ext(data_in)) == "rds") {
    data_in = readRDS( file = paste0(DB_DIR, "/", db_prefix, "_X.rds"))


    if(class(data_in)[1] == "Seurat"){

      # how stereotyped is this pattern?  check for Oscar...
      ad <- sceasy::convertFormat(data_in, from="seurat", to="anndata",
                                  outFile = NULL,
                                  assay = 'SCT',
                                  main_layer = 'data',
                                  transfer_layers = c('data', 'counts', 'scale.data')
      )

      raw <- sceasy::convertFormat(data_in, from="seurat", to="anndata",
                                   outFile = NULL,
                                   assay = 'RNA',
                                   main_layer = 'counts',
                                   transfer_layers = NULL)


      ad$raw <- raw


    } else if (class(data_in)[1] == "SingleCellExperiment") {

    } else if ("data.frame" in class(data_in)) {
      # could _everything be in a dataframe???
      # yes... lipidomic... strip off first two columns?
    }

  } else if (tolower(tools::file_ext(obj)) == "h5ad"){


  } else if (tolower(tools::file_ext(obj)) == "loom"){

  }





} else {

}

}


setup_database <- function(database_name, data_in, db_meta , re_pack=TRUE){
  #LOAD & PACK into ANNDATA
  ##if data_in contains filenames they must be the full path (i.e. RAW_DIR inlcuded)

  DB_NAME <- database_name

  if tolower(data_type) == "seurat" {
    require("Seurat")
    require("sceasy")
  }

  require(reticulate)
  reticulate::use_condaenv(required = TRUE, condaenv = 'omxr')
  require(anndata)

  DB_DIR = file.path("data-raw",DB_NAME)
  if (!dir.exists(DB_DIR)) {
    dir.create(DB_DIR)
  }

  # # db_meta e.g.
  #   organism <- ""
  #   lab <- "Menon"
  #   annotation_database <- "NA"
  # check to see what level of data we were given in our data_list
  #
  if (length(data_in) & names(data_in)=="object") {
    data_in <- unlist(data_in)
  }

  # check if we have it or are forcing
  if ( !file.exists(past0(DB_DIR,"/core_data.h5ad"))
       | re_pack ) {
    #create it
    # sub-functions to deal with what kind of data we have...
    ad <- omxr_pack_anndata(data_in)

  } else {
    #load it
    ad <- read_h5ad(file.path(DB_DIR,"core_data.h5ad"))
  }

  return(ad)
}



# pre-pre processing
#update seurat object -----------------------------
microglia_data_updated <- UpdateSeuratObject(object = microglia_data)
#> OUTPUT:
#   Validating object structure
#   Updating object slots
#   Ensuring keys are in the proper strucutre
#   Ensuring feature names don't have underscores or pipes
#   Object representation is consistent with the most current Seurat version

saveRDS(microglia_data_updated,
        file = file.path(RAW_DIR, "microglia_data_seu_new.rds"))


data_list <- list(object=file.path(RAW_DIR, "microglia_data_seu_new.rds"))

# scemaa
#c("object", "data_mat","obs_meta","var_annot","omics","sample_ID","etc")
data_list <- list(data_mat = NA,
                  obs_meta = NULL,
                  var_annot = NULL,
                  omics = NULL,
                  sample_ID = NULL,
                  etc = NULL)

ad <- setup_database(database_name, data_in, data_type, db_meta , re_pack=TRUE){


  sc <- import("scanpy")


  test_types <- c('wilcoxon','t-test_overestim_var')
  comp_types <- c("allVrest")


  helper_function<-('data-raw/compute_de_table.R')
  source(helper_function)

  diff_exp <- compute_de_table(ad_,comp_types, test_types, obs_name = c('disease','cell_type'))


  # put the logvals in layers of ad


  # measures
  #  This ordering is the "default"
  measures <- list(obs = c("nCount_SCT","nFeature_SCT","nCount_RNA","nFeature_RNA"),
                   var = NA)
  # [1] "sct.detection_rate"    "sct.gmean"             "sct.variance"
  # [4] "sct.residual_mean"     "sct.residual_variance" "sct.variable"

  # differentials  #if we care we need to explicitly state. defaults will be the order...
  diffs <- list(diff_exp_groups =  levels(factor(diff_exp$group)),
                diff_exp_comp_type =  levels(factor(diff_exp$comp_type)),
                diff_exp_obs_name =  levels(factor(diff_exp$obs_name)),
                diff_exp_tests =  levels(factor(diff_exp$test_type)))

  # Dimred
  dimreds <- list(obsm = c('X_pca', 'X_tsne'),
                  varm = c('PCs'))

  # what ad$obs do we want to make default values for...
  # # should just pack according to UI?
  default_factors <- c("tissue","disease","cell_type")




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




  require("Seurat")

  # if (!requireNamespace("BiocManager", quietly = TRUE))
  #   install.packages("BiocManager")
  # BiocManager::install(c("LoomExperiment", "SingleCellExperiment"))
  #devtools::install_github("cellgeni/sceasy")
  require("sceasy")

  require(reticulate)
  reticulate::use_condaenv(required = TRUE, condaenv = 'omxr')
  require(anndata)

  # create the folder to contain the raw data
  DB_NAME = "vilas_microglia_sceasy"
  DB_DIR = file.path("data-raw",DB_NAME)
  if (!dir.exists(DB_DIR)) {
    dir.create(DB_DIR)
  }


  # ------------------------------------------
  # # 1. documentation / provenance ------------
  # # ------------------------------------------
  #
  # # TODO:  markdown file or html with some copy about the database
  # #  - lab, paper link/name
  # #  summarize results / data origin whatever
  #
  #
  # organism <- ""
  # lab <- "Menon"
  # annotation_database <- "NA"


  # ------------------------------------------
  # 2. helper functions ----------------------
  # ------------------------------------------

  #install.packages("Seurat")
  library("Seurat")

  # better script -------------------------

  DATA_DIR <- "/Users/ahenrie/Projects/NDCN_dev/dbrowse/ingest"
  DATA_DIR <- "ingest"
  #basedir <- "/Users/ahenrie/Projects/NDCN_dev/Oscar"
  DB_NAME = "Oscar_toy1"


  # update data to new seurat object.

  file_name <- file.path(DB_NAME,"toyGeneNames.RData")
  file_path <- file.path(DATA_DIR,file_name)

  load(file_path)
  gene_names = toyGeneNames
  rm( toyGeneNames)


  #Load the differentially epxressed Genes (APOE and 100 other random picked genes)
  file_name <- file.path(DB_NAME,"toyDE.RData")
  file_path <- file.path(DATA_DIR,file_name)

  load(file_path)
  differential_expre = toyDifferentialExpre
  rm( toyDifferentialExpre)


  # give the elements of the list (clusters) more informative names
  names( differential_expre) = c('Neuron', 'Oligodendrocytes', 'Astrocytes', 'Microglia', 'OPC', 'Endothelial')

  file_name <- file.path(DB_NAME,"toyCells.RData")
  file_path <- file.path(DATA_DIR,file_name)
  load(file_path)
  old_data_marker_gene <- toy
  rm( toy)

  #UPDATE OBJECT BY REMAKING IT
  label_data_marker_gene <- CreateSeuratObject(counts = old_data_marker_gene@data, meta.data = old_data_marker_gene@meta.data)



  # Warning: Feature names cannot have underscores ('_'), replacing with dashes ('-')
  label_data_marker_gene <- UpdateSeuratObject(object = label_data_marker_gene)

  #label_data_marker_gene <- SetAllIdent(object = label_data_marker_gene, id = "res.0.6")
  Idents(object = label_data_marker_gene) <- label_data_marker_gene@meta.data$res.0.6


  #ged the dimension reductions
  dim_red <- old_data_marker_gene@dr

  label_data_marker_gene[["pca"]] <- CreateDimReducObject(
    embeddings = dim_red$pca@cell.embeddings,
    loadings = dim_red$pca@gene.loadings,
    key = dim_red$pca@key,
    assay = "RNA")


  label_data_marker_gene[["tsne"]] <- CreateDimReducObject(
    embeddings = dim_red$tsne@cell.embeddings,
    loadings = dim_red$tsne@gene.loadings,
    key = dim_red$tsne@key,
    assay = "RNA")



  #######  save data in "new" Seurat
  # save the data blobs into something easier to work with...
  #
  file_name <- file.path(DB_NAME,"label_data_marker_gene.Rds")
  file_path <- file.path(DATA_DIR,file_name)
  #
  save(differential_expre, file=file_path)

  #cnv_data_list
  file_name <- file.path(DB_NAME,"gene_names.Rds")
  file_path <- file.path(DATA_DIR,file_name)

  save(gene_names, file=file_path)


  differential_expre$Neuron$cell_type<- "Neuron"
  differential_expre$Oligodendrocytes$cell_type <- "Oligodendrocytes"
  differential_expre$Astrocytes$cell_type <- "Astrocytes"
  differential_expre$Microglia$cell_type<- "Microglia"
  differential_expre$OPC$cell_type<- "OPC"
  differential_expre$Endothelial$cell_type<- "Endothelial"

  data_table <- dplyr::bind_cols(list(differential_expre$Neuron,
                                      differential_expre$Oligodendrocytes,
                                      differential_expre$Astrocytes,
                                      differential_expre$Microglia,
                                      differential_expre$OPC,
                                      differential_expre$Endothelial) ,repair = "unique")

  #JDA
  tsnePlot = TSNEPlot(object = label_data_marker_gene, pt.size=0.5, label.size = 2.8, label = TRUE, label.box=TRUE)
  DimPlot(object = label_data_marker_gene, reduction = "tsne", pt.size=0.5, label.size = 2.8, label = TRUE, label.box=TRUE)

  ids <- Idents(object = label_data_marker_gene)


  # make ready for anndata / SCE
  counts <- label_data_marker_gene@assays[1]$RNA@counts
  meta_data <- label_data_marker_gene@meta.data
  dimred_pca <- label_data_marker_gene@reductions$pca@cell.embeddings
  dimred_tsne <- label_data_marker_gene@reductions$tsne@cell.embeddings

  X <- counts
  obs <- meta_data
  vars <- data_table

  oscar_toy1_X <- X
  oscar_toy1_obs <- obs
  oscar_toy1_vars <- vars

  # ------------------------------------------
  # 3. load data -----------------------------
  # ------------------------------------------
  RAW_DIR <- "ingest/Vilas_B"
  # Seurat Data object
  data_file <- "microglia_data_with_counts_RNA_SCT.rda"
  load(file.path(RAW_DIR,data_file) )



  # 3a. update seurat object -----------------------------
  microglia_data_updated <- UpdateSeuratObject(object = microglia_data)



  #> OUTPUT:
  #   Validating object structure
  #   Updating object slots
  #   Ensuring keys are in the proper strucutre
  #   Ensuring feature names don't have underscores or pipes
  #   Object representation is consistent with the most current Seurat version

  saveRDS(microglia_data_updated,
          file = file.path(RAW_DIR, "microglia_data_seu_new.rds"))



  microglia_data_updated <- readRDS(file = file.path(RAW_DIR, "microglia_data_seu_new.rds"))




  ad <- sceasy::convertFormat(microglia_data_updated, from="seurat", to="anndata",
                              outFile = NULL,
                              assay = 'SCT',
                              main_layer = 'data',
                              transfer_layers = c('data', 'counts', 'scale.data')
  )

  raw <- sceasy::convertFormat(microglia_data_updated, from="seurat", to="anndata",
                               outFile = NULL,
                               assay = 'RNA',
                               main_layer = 'counts',
                               transfer_layers = NULL)


  ad$raw <- raw


  ad$write_h5ad(filename=file.path(DB_DIR,"core_data.h5ad"))

  # confirm that casting to anndatainR wrapper doesn't change the file....
  ad_R <- py_to_r(ad)
  ad_R$write_h5ad(filename=file.path(DB_DIR,"core_dataR.h5ad"))

  ad_py <- read_h5ad(file.path(DB_DIR,"core_data.h5ad"))
  ad_R <- read_h5ad(file.path(DB_DIR,"core_dataR.h5ad"))






  # ------------------------------------------
  # 4. pack into anndata   (redundant?)      --
  # ------------------------------------------0

  # R TOOLS
  ad_r <- reticulate::py_to_r(ad)

  var_ <- ad_r$var
  obs <- ad_r$obs
  X <- ad_r$X
  ad_raw <- ad_r$raw
  obsm <- ad_r$obsm
  varm <- ad_r$varm

  layers <- ad_r$layers$get('counts')

  # db_prefix = "core_data"
  # saveRDS(X, file = paste0(DB_DIR, "/", db_prefix, "_X.rds"))
  # saveRDS(obs, file = paste0(DB_DIR, "/", db_prefix, "_obs.rds"))
  # saveRDS(var_, file = paste0(DB_DIR, "/", db_prefix, "_var.rds"))

  ad_ <- AnnData(
    X = X,
    obs = obs,
    var = var_,
    raw = ad_raw,
    obsm = obsm,
    varm=varm,
    layers =  list(counts=ad$layers$get('counts')) #list('count'=layers)
  )

  ad_

  ad_$write_h5ad(filename=file.path(DB_DIR,"recastR.h5ad"))


  detach("package:anndata",unload=TRUE)
  sc <- import("scanpy")
  anndat <- import("anndata")

  var_ <- ad$var
  obs <- ad$obs
  X <- ad$X
  ad_raw <- ad$raw
  obsm <- ad$obsm
  varm <- ad$varm

  layers <- ad$layers


  ad_py <- sc$AnnData(
    X = X,
    obs = obs,
    var = var_,
    raw = ad_raw,
    obsm = obsm,
    varm=varm,
    layers =  layers
  )
  #ad$write_h5ad(filename=file.path(DB_DIR,"core_data.h5ad"))
  # returns AnnDataR6 !?!?!?

  ad_py
  ad_py$write_h5ad(filename=file.path(DB_DIR,"recast_py.h5ad"))


  ad_py2 <- anndat$AnnData(
    X = ad$X,
    obs = ad$obs,
    var = ad$var,
    raw = ad$raw,
    obsm = ad$obsm,
    varm= ad$varm,
    raw = ad$raw,
    layers =  ad$layers
  )
  #ad$write_h5ad(filename=file
  #
  ad_py2$write_h5ad(filename=file.path(DB_DIR,"recast_py2.h5ad"))

  # ------------------------------------------
  # 5. post processing                      --
  # ------------------------------------------
  require(anndata)
  require(reticulate)
  DB_NAME = "vilas_microglia_seu"
  DB_DIR = file.path("data-raw",DB_NAME)
  RAW_DIR <- "ingest/Vilas_B"

  db_prefix = "core_data"
  X = readRDS( file = paste0(DB_DIR, "/", db_prefix, "_X.rds"))
  obs = readRDS(  file = paste0(DB_DIR, "/", db_prefix, "_obs.rds"))
  var_ = readRDS( file = paste0(DB_DIR, "/", db_prefix, "_var.rds"))


  ad2 <- read_h5ad(file.path(DB_DIR,"core_data.h5ad"))
  ad


  sc <- import("scanpy")


  test_types <- c('wilcoxon','t-test_overestim_var')
  comp_types <- c("allVrest")


  helper_function<-('data-raw/compute_de_table.R')
  source(helper_function)

  diff_exp <- compute_de_table(ad_,comp_types, test_types, obs_name = c('disease','cell_type'))


  # put the logvals in layers of ad
  # copy the whole thing and replace X to copy the uns to ad

  ad$write_h5ad(filename=file.path(DB_DIR,"core_data_plus_de.h5ad"))

  # also need to pack the diff_exp1 and diff_exp2 into easy to deal wiht tables for volcanos...
  table(ad$obs$disease)


  db_prefix = "de"
  saveRDS(diff_exp, file = paste0(DB_DIR, "/", db_prefix, "_table.rds"))


  # ------------------------------------------
  # 6. dimension reduction - PCA / umap    --
  # ------------------------------------------
  require(anndata)
  require(reticulate)

  DB_NAME = "vilas_microglia_seu"
  DB_DIR = file.path("data-raw",DB_NAME)
  #RAW_DIR <- "ingest/Vilas_B"

  ad <- read_h5ad(file.path(DB_DIR,"core_data_plus_de.h5ad"))




  # ------------------------------------------
  # 7 . create config and default files                   --
  # ------------------------------------------
  require(anndata)
  require(reticulate)

  DB_NAME = "vilas_microglia_seu"
  DB_DIR = file.path("data-raw",DB_NAME)
  #RAW_DIR <- "ingest/Vilas_B"

  ad <- read_h5ad(file.path(DB_DIR,"core_data_plus_de.h5ad"))
  db_prefix = "de"
  diff_exp = readRDS( file = paste0(DB_DIR, "/", db_prefix, "_table.rds"))




  ###

  # measures
  #  This ordering is the "default"
  measures <- list(obs = c("nCount_SCT","nFeature_SCT","nCount_RNA","nFeature_RNA"),
                   var = NA)
  # [1] "sct.detection_rate"    "sct.gmean"             "sct.variance"
  # [4] "sct.residual_mean"     "sct.residual_variance" "sct.variable"

  # differentials  #if we care we need to explicitly state. defaults will be the order...
  diffs <- list(diff_exp_groups =  levels(factor(diff_exp$group)),
                diff_exp_comp_type =  levels(factor(diff_exp$comp_type)),
                diff_exp_obs_name =  levels(factor(diff_exp$obs_name)),
                diff_exp_tests =  levels(factor(diff_exp$test_type)))

  # Dimred
  dimreds <- list(obsm = c('X_pca', 'X_tsne'),
                  varm = c('PCs'))

  # what ad$obs do we want to make default values for...
  # # should just pack according to UI?
  default_factors <- c("tissue","disease","cell_type")




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


  vilas_microglia_sceasy_conf = readRDS( paste0(DB_DIR,"/",db_prefix,"_conf.rds") )
  vilas_microglia_sceasy_def = readRDS( paste0(DB_DIR,"/",db_prefix,"_def.rds") )
  vilas_microglia_sceasy_omics = readRDS( paste0(DB_DIR,"/",db_prefix,"_omics.rds") )
  vilas_microglia_sceasy_meta = readRDS( paste0(DB_DIR,"/",db_prefix,"_meta.rds") )


  #usethis::use_data(vilas_microglia_sceasy_conf, vilas_microglia_sceasy_def, vilas_microglia_sceasy_omics, vilas_microglia_sceasy_meta, overwrite = TRUE)


