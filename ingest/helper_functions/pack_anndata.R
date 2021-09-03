#' Title
#'
#' @param data_in
#'
#' @return
#' @export
#'
#' @examples
pack_anndata <- function (data_in){

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

    ad <- anndata::AnnData(
      X = X,
      obs = obs,
      var = var_,
      uns = uns
    )

  } else if (tolower(tools::file_ext(data_in)) == "rds") {

    data_in = readRDS( file = data_in )


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



      #enforce sample_ID
      # TODO: use anndata:: instead of py_to_r wrappers?
      if ( !("sample_ID" %in% py_to_r(ad$obs_keys()) ) ){
        tmp <- dplyr::mutate( py_to_r(ad$obs),
                              sample_ID=py_to_r(ad$obs_names))
        ad$obs <- r_to_py(tmp)

        tmp <- dplyr::mutate( py_to_r(raw$obs),
                              sample_ID=py_to_r(raw$obs_names))
        raw$obs <- r_to_py(tmp)

      }

      ad$raw <- raw

    } else if (class(data_in)[1] == "SingleCellExperiment") {
      print("SingleCellExperiment not enabled")
      ad <- NULL

    } else if ("data.frame" %in% class(data_in)) {
      # could _everything be in a dataframe???
      # yes... lipidomic... strip off first two columns?
      print("enable this for lipidomics? not enabled")
      ad <- NULL

    }

  } else if (tolower(tools::file_ext(data_in)) == "h5ad") {
    ad <- anndata::read_h5ad(data_in)


  } else if (tolower(tools::file_ext(data_in)) == "loom"){
    print("loom loading not enabled")
    ad <- NULL
  }


  return(ad)
}

