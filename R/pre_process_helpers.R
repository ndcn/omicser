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
#' @examples TODO
compute_de_table <- function(adata,comp_types, test_types, obs_names,sc) {
  # this should update adata in place with the diff_exp data...

  diff_exp <- data.frame()
  for (obs_name in obs_names){
    for (test_type in test_types) {
      for (comp_type in comp_types) {
        parts <- strsplit(comp_type,"V")[[1]]
        reference <- parts[2]
        group <- parts[1]
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
      var_ <- var_ %>% dplyr::rename(omics_name=all_of(omics_col))
    } else {
      var_ <- var_ %>% dplyr::mutate(omics_name=data_in$omics)
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
      # could _everything be in a dataframe???
      # yes... lipidomic... strip off first two columns?
      print("WARNING:  did not receive a list of data as expected")
      adata <- NULL

    }
  return(adata)
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
    adata <- anndata::AnnDataR6$new(
              sceasy::convertFormat(data_in, from="seurat", to="anndata",
                                outFile = NULL,
                                assay = 'SCT',
                                main_layer = 'data',
                                transfer_layers = c('data', 'counts', 'scale.data')
                                )
    )

    raw <- anndata::AnnDataR6$new(
            sceasy::convertFormat(data_in, from="seurat", to="anndata",
                                 outFile = NULL,
                                 assay = 'RNA',
                                 main_layer = 'counts',
                                 transfer_layers = NULL)
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
    if ( !("omics_name" %in% adata$var_keys()) ){
      var_ <- adata$var
      var_ <- var_ %>% dplyr::mutate(omics_name=adata$var_names) %>%
        dplyr::relocate(omics_name)

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
      var_ <- var_ %>% dplyr::mutate(omics_name=adata$raw$var_names) %>%
        dplyr::relocate(omics_name)

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
      if ( !("omics_name" %in% adata$var_keys()) ){
        var_ <- adata$var
        var_ <- var_ %>% dplyr::mutate(omics_name=adata$var_names) %>%
                        dplyr::relocate(omics_name)

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
        var_ <- var_ %>% dplyr::mutate(omics_name=adata$raw$var_names) %>%
          dplyr::relocate(omics_name)

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
