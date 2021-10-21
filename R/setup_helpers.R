
#' write_config: helper for Writing the Config file for setting up the omicser databases
#'
#' @param in_options list of database names, conda environment, and database path
#' @param in_path where the config file lives. default cwd
#'
#' @return
#' @export write_config
#' @importFrom configr write.config
#' @examples TODO
write_config <- function(in_options, in_path=NULL) {
  if ( is.null(in_path) ) {
    in_path <- getwd()
  }

  #TODO: check_config_list(in_options)

  OPTIONS_FILE <- file.path(in_path,'omicser_options.yml')
  write.config(config.dat = in_options, file.path = OPTIONS_FILE,
               write.type = "yaml", indent = 4)

}



#TODO:  make a check_config_list <- function(in_options)
#         to ensure that we have everything or we are setting defaults
#


#' get_config: helper for reading the Config file containing information about the databases
#'
#' @param in_path where the config file lives. default getwd() test
#'
#' @return the list of options contained in `omicser_options.yml`
#' @export get_config
#' @importFrom configr read.config
#' @examples TODO
get_config <- function(in_path = NULL, is_running = FALSE) {

  # TODO:  can we make this more dynamic??
  # OPTIONS_FILE <- system.file('omicser_options.yml',
  #                             package = 'omicser')

  if (is_running) {
    OPTIONS_FILE <- system.file('omicser_options.yml', package = 'omicser')
  } else {
    if ( is.null(in_path) ) {
      in_path <- getwd()
    }
    OPTIONS_FILE <- file.path(in_path,'omicser_options.yml')
  }

  options_list <- read.config( file = OPTIONS_FILE )
  #TODO: check_config_list(in_options)
  #
  return(options_list)
}




#' write_db_meta: helper for writing the database meta info
#'
#' @param db_meta data to write
#' @param db_name what is out database is called
#' @param db_root where the folder lives. default cwd
#' @return
#' @export write_db_meta
#' @importFrom configr write.config
#' @examples TODO
write_db_meta <- function(db_meta,db_name, db_root = NULL) {

  if ( is.null(db_root) ) {
    db_root <- getwd()
  }

  out_fn <-  file.path(db_root,db_name,"db_meta.yml")
  write.config(config.dat = db_meta, file.path = out_fn,
               write.type = "yaml", indent = 4)

}


#' get_db_meta: helper for reading the database meta data
#'
#' @param db_name what is out db called
#' @param db_root where the folder lives. default cwd
#'
#' @return meta_info: database meta deta
#' @export get_db_meta
#' @importFrom configr read.config
#' @examples TODO
get_db_meta <- function(db_name, db_root = NULL)  {

  if ( is.null(db_root) ) {
    db_root <- getwd()
  }

  DB_META_FILE <-  file.path(db_root,db_name,"db_meta.yml")

  meta_info <- read.config( file = DB_META_FILE )
  return(meta_info)
}


#TODO:  make a check_db_config_list <- function(config_list)
#         to ensure that we have everything or we are setting defaults
#


#' write_db_conf: helper for writing the ui config options for the database
#'
#' @param config_list the database configuration list
#' @param db_name what is out database is called
#' @param db_root where the folder lives. default cwd
#' @return
#' @export write_db_conf
#' @importFrom configr write.config
#' @examples TODO
write_db_conf <- function(config_list,db_name, db_root = NULL) {

  if ( is.null(db_root) ) {
    db_root <- getwd()
  }

  # TODO:  check_db_config_list(config_list)
  out_fn <-  file.path(db_root,db_name,"db_config.yml")
  write.config(config.dat = config_list, file.path = out_fn,
               write.type = "yaml", indent = 4)

}



#' get_db_conf: helper for reading the database config
#'
#' @param db_name what is out db called
#' @param db_root where the folder lives. default cwd
#'
#' @return config_list: list of UI configurations for accessing data
#' @export get_db_conf
#' @importFrom configr read.config
#' @examples TODO
get_db_conf <- function(db_name, db_root = NULL)  {

  if ( is.null(db_root) ) {
    db_root <- getwd()
  }
  DB_CONF_FILE <-  file.path(db_root,db_name,"db_config.yml")

  config_list <- read.config( file = DB_CONF_FILE )

  # TODO:  check_db_config_list(config_list)

  return(config_list)
}



