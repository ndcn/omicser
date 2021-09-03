
#' setup_database
#'
#' @description A a function to create the anndata database
#'
#' @return The return value, if any, from executing the function.
#' @export setup_database
#'
#' @examples  TODO
setup_database <- function(database_name, data_in, db_meta , re_pack=TRUE){
  #LOAD & PACK into ANNDATA
  ##if data_in contains filenames they must be the full path (i.e. RAW_DIR inlcuded)

  DB_NAME <- database_name


  # require(reticulate)
  # reticulate::use_condaenv(required = TRUE, condaenv = 'omxr')
  # require(anndata)

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
  if (length(data_in) & names(data_in[1])=="object") {
    data_in <- unlist(data_in)
  }

  # check if we have it or are forcing
  if ( !file.exists(paste0(DB_DIR,"/core_data.h5ad"))
       | re_pack ) {
    #create it
    # sub-functions to deal with what kind of data we have...
    ad <- pack_anndata(data_in)

  } else {
    #load it
    ad <- anndata::read_h5ad(file.path(DB_DIR,"core_data.h5ad"))
  }

  return(ad)
}


