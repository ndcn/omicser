# this is a collection of functions to help with pre-processing/curating data
#


#' Title
#'
#' @param adata the anndata object
#' @param comp_types what kind of comparisons?  "allVrest" or "{a}V{b}"
#' @param test_types statistical tests - t-test, etc
#' @param obs_names name of the adata$obs column defining the comparision groups
#' @param sc scanpy
#'
#' @return
#' @export compute_de_table
#'
#' @importFrom stringr str_match
#'
#' @examples TODO
compute_de_table <- function(adata,comp_types, test_types, obs_names,sc) {
  # this should update adata in place with the diff_exp data...

  diff_exp <- data.frame()
  for (obs_name in obs_names){
    for (test_type in test_types) {
      for (comp_type in comp_types) {
        parts <- str_match(pattern = "^\\{?([a-zA-Z\\+\\._]+)\\}?V\\{?([a-zA-Z\\+\\._]+)\\}?$",
                           string = comp_type)
        reference <- parts[3]
        group <- parts[2]
        print(test_type)
        print(group)
        print(reference)
        key <- paste0(test_type,"_", comp_type)


        if (reference == "rest") { #grpVrest
          sc$tl$rank_genes_groups(adata,
                                  obs_name,
                                  groups = "all",
                                  reference = reference,
                                  method=test_type,
                                  use_raw = FALSE,
                                  key_added = key)
          de_table <- sc$get$rank_genes_groups_df(adata,
                                                  group=NULL,
                                                  key=key)
          de_table$comp_type <- comp_type
        } else { #compare group vs reference
          sc$tl$rank_genes_groups(adata,
                                  obs_name,
                                  groups = list(group),
                                  reference = reference,
                                  method=test_type,
                                  use_raw = FALSE,
                                  key_added = key)
          de_table <- sc$get$rank_genes_groups_df(adata,
                                                  group=group,
                                                  key=key)
          de_table$group <- group
          de_table$comp_type <- 'grpVref'
        }

        de_table$reference <- reference
        de_table$test_type <- test_type
        de_table$obs_name <- obs_name
        de_table$versus <- paste0(de_table$group," vs. ", reference)

        diff_exp <- dplyr::bind_rows(diff_exp, de_table)
      }
    }
  }

  return(diff_exp)
}



#' pack_anndata
#'
#' @param data_in - list of matrices, tables, and lists containing the data to pack into the db
#'
#' @return adata the anndata object we are browsing
#' @export pack_anndata_from_csv
#'
#' @examples  TODO
pack_anndata_from_csv <- function(data_in){

  #tools::file_path_sans_ext(data_in)
  if ( class(data_in)[1] == "list" ) {
    # multple files in c("object", "data_mat","obs_meta","var_annot","omics","sample_ID",etc")
    # data_mat - data matrix.  need to assert that this matrix object has omics and sample+ID names

    if (dim(data_in$data_mat)[2] != dim(data_in$var_annot)[1]) {
      X <-Matrix::t(data_in$data_mat)
    } else {
      X <- data_in$data_mat
    }
    if (is.null(dimnames(X)[1])) { rownames(X) <- data_in$sample_ID }
    if (is.null(dimnames(X)[2])) { colnames(X) <- data_in$omics  }


    # obs_meta - ensure that we are factored and sample_ID is first column
    obs <- data_in$obs_meta
    id_col <- colnames(obs)[1]
    if (all(obs[[id_col]] == data_in$sample_ID)) {
      obs <- obs %>% dplyr::rename(sample_ID=all_of(id_col))
    } else {
      obs <- obs %>% dplyr::mutate(sample_ID=data_in$sample_ID)
    }
    if(is.null(rownames(obs)) ){
      rownames(obs) <- data_in$sample_ID
    }

    # var_annot - ensure that "omics" is first column
    var_ <- data_in$var_annot
    omics_col <- colnames(var_)[1]

    if ( all(var_[[omics_col]] == data_in$omics) ) {
      var_ <- var_ %>% dplyr::rename(feature_name=all_of(omics_col))
    } else {
      var_ <- var_ %>% dplyr::mutate(feature_name=data_in$omics)
    }
    if(is.null(rownames(var_)) ){
      rownames(var_) <- data_in$omics
    }

    # etc goes into an uns entry
    #
    if ( class(data_in$uns)=="list" ){
      uns <- data_in$uns
    } else {
      uns <- list(etc=data_in$uns)
    }

    adata <- anndata::AnnData(
      X = X,
      obs = obs,
      var = var_,
      uns = uns
    )
  } else {
      print("WARNING:  did not receive a list of data as expected")
      adata <- NULL

    }
  return(adata)
}




# sceasy functions rewritten here for simplicity (Bioconductor deps are killing me)
# COPIED FROM https://github.com/cellgeni/sceasy/blob/master/R/functions.R

.regularise_df <- function(df, drop_single_values = TRUE) {
  if (ncol(df) == 0) df[['name']] <- rownames(df)
  if (drop_single_values) {
    k_singular <- sapply(df, function(x) length(unique(x)) == 1)
    if (sum(k_singular) > 0)
      warning(paste('Dropping single category variables:'),
              paste(colnames(df)[k_singular], collapse=', '))
    df <- df[, !k_singular, drop=F]
    if (ncol(df) == 0) df[['name']] <- rownames(df)
  }
  return(df)
}


# modified from sceasy original to take advantage of the anndata package
seurat2anndata <- function(obj,
                           outFile = NULL,
                           assay = 'RNA',
                           main_layer = 'data',
                           transfer_layers = NULL,
                           drop_single_values = TRUE) {
  main_layer <- match.arg(main_layer, c('data', 'counts', 'scale.data'))
  transfer_layers <- transfer_layers[ transfer_layers %in% c('data', 'counts', 'scale.data') ]
  transfer_layers <- transfer_layers[ transfer_layers != main_layer ]

  if (compareVersion(as.character(obj@version), '3.0.0') < 0)
    obj <- Seurat::UpdateSeuratObject(object = obj)

  X <- Seurat::GetAssayData(object = obj, assay = assay, slot = main_layer)

  obs <- .regularise_df(obj@meta.data, drop_single_values = drop_single_values)

  var <- .regularise_df(Seurat::GetAssay(obj, assay = assay)@meta.features, drop_single_values = drop_single_values)

  obsm <- NULL
  reductions <- names(obj@reductions)
  if (length(reductions) > 0) {
    obsm <- sapply(
      reductions,
      function(name) as.matrix(Seurat::Embeddings(obj, reduction=name)),
      simplify = FALSE
    )
    names(obsm) <- paste0('X_', tolower(names(obj@reductions)))
  }

  layers <- list()
  for (layer in transfer_layers) {
    mat <- Seurat::GetAssayData(object = obj, assay = assay, slot = layer)
    if (all(dim(mat) == dim(X))) layers[[layer]] <- Matrix::t(mat)
  }

  anndata <- reticulate::import('anndata', convert = FALSE)
  adata <- anndata::AnnData(
  #adata <- anndata$AnnData(
    X = Matrix::t(X),
    obs = obs,
    var = var,
    obsm = obsm,
    layers = layers
  )

  if (!is.null(outFile))
    adata$write(outFile, compression = 'gzip')

  adata
}



#' pack_anndata_from_seurat
#'
#' @param seurat_obj
#'
#' @return adata the anndata object we are browsing
#' @export pack_anndata_from_seurat
#'
#' @examples  TODO
pack_anndata_from_seurat <- function(seurat_obj_name){


  data_in = readRDS( file = seurat_obj_name )

  if(class(data_in)[1] == "Seurat"){
    # how stereotyped is this pattern?  check for Oscar...
    adata <- seurat2anndata(data_in, outFile = NULL,
                                assay = 'SCT',
                                main_layer = 'data',
                                transfer_layers = c('data', 'counts', 'scale.data')
                                )


    raw <- seurat2anndata(data_in, outFile = NULL,
                              assay = 'RNA',
                              main_layer = 'counts',
                              transfer_layers = NULL
                              )

    if ( !("sample_ID" %in% adata$obs_keys()) ){
      obs <- adata$obs
      obs <- obs %>% dplyr::mutate(sample_ID=adata$obs_names) %>%
                     dplyr::relocate(sample_ID)
      adata$obs <- obs
    }

    #enforce sample_ID
    # TODO:
    #      replace dplyr with data.table
    if ( !("feature_name" %in% adata$var_keys()) ){
      var_ <- adata$var
      var_ <- var_ %>% dplyr::mutate(feature_name=adata$var_names) %>%
                       dplyr::relocate(feature_name)
      adata$var <- var_
    }

    # force RAW?
    if (is.null(adata$raw)) {
      adata$raw <- adata$copy()
    } else {
      obs <- adata$raw$obs
      obs <- obs %>% dplyr::mutate(sample_ID=adata$raw$obs_names) %>%
        dplyr::relocate(sample_ID)
      adata$raw$obs <- obs
      var_ <- adata$raw$var
      var_ <- var_ %>% dplyr::mutate(feature_name=adata$raw$var_names) %>%
                       dplyr::relocate(feature_name)
      adata$raw$var <- var_
    }

    #  DISABLED >> getting bio-conductor dependencies is a pain...
    #     } else if (class(data_in)[1] == "SingleCellExperiment") {
    #       print("SingleCellExperiment not enabled")
    #       adata <- NULL


  } else {
    print("WARNING:  this is not a seurate object")
    adata <- NULL

  }


  return(adata)
}





#' setup_database
#'
#' @param database_name - name of the database being setup
#' @param db_path - system path to the databases
#' @param data_in - a list of data, or a filename
#' @param db_meta - "meta" information about the database origin/ context / etc
#' @param re_pack - force re-packing the anndata object
#'
#' @description A a function to create the anndata database
#'
#' @return The return value, if any, from executing the function.
#' @export setup_database
#'
#' @examples  TODO
setup_database <- function(database_name, db_path, data_in, db_meta , re_pack=TRUE){
  #LOAD & PACK into ANNDATA
  ##if data_in contains filenames they must be the full path (i.e. RAW_DIR inlcuded)

  # require(reticulate)
  # reticulate::use_condaenv(required = TRUE, condaenv = 'omxr')
  # require(anndata)

  DB_DIR = file.path(db_path,database_name)
  if (!dir.exists(DB_DIR)) {
    dir.create(DB_DIR)
  }

  # # db_meta e.g.
  #   organism <- ""
  #   lab <- "Menon"
  #   annotation_database <- "NA"
  # check to see what level of data we were given in our data_list
  #
  if (length(data_in) & names(data_in[1])=="object") {
    data_in <- unlist(data_in)
  }

  # check if we have it or are forcing
  if ( !file.exists(paste0(DB_DIR,"/core_data.h5ad"))
       | re_pack ) {

        # sub-functions to deal with what kind of data we have...
    #tools::file_path_sans_ext(data_in)
    if ( class(data_in)[1] == "list" ) {
      adata <- pack_anndata_from_csv(data_in)

    } else if (tolower(tools::file_ext(data_in)) == "rds") {
      adata <- pack_anndata_from_seurat(data_in)
      #TODO: logic for singleCellExperiment goes here

    } else if (tolower(tools::file_ext(data_in)) == "h5ad") {
      adata <- anndata::read_h5ad(data_in)
      # obs_meta - ensure that we are factored and sample_ID is first column

      if ( !("sample_ID" %in% adata$obs_keys()) ){
         obs <- adata$obs
         obs <- obs %>% dplyr::mutate(sample_ID=adata$obs_names) %>%
                        dplyr::relocate(sample_ID)
         adata$obs <- obs
       }
      #enforce sample_ID
      # TODO:
      #      replace dplyr with data.table
      if ( !("feature_name" %in% adata$var_keys()) ){
        var_ <- adata$var
        var_ <- var_ %>% dplyr::mutate(feature_name=adata$var_names) %>%
                        dplyr::relocate(feature_name)
        adata$var <- var_

      }

      # force RAW?
      if (is.null(adata$raw)) {
        adata$raw <- adata$copy()
      } else {
        obs <- adata$raw$obs
        obs <- obs %>% dplyr::mutate(sample_ID=adata$raw$obs_names) %>%
          dplyr::relocate(sample_ID)
        adata$raw$obs <- obs
        var_ <- adata$raw$var
        var_ <- var_ %>% dplyr::mutate(feature_name=adata$raw$var_names) %>%
                         dplyr::relocate(feature_name)
        adata$raw$var <- var_
      }

    } else if (tolower(tools::file_ext(data_in)) == "loom"){
      print("loom loading not enabled")
      adata <- NULL
    }

  } else {
    #load it
    adata <- anndata::read_h5ad(file.path(DB_DIR,"core_data.h5ad"))
  }

  return(adata)
}
