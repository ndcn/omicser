
#' write_config: helper for Writing the Config file for setting up the omicser databases
#'
#' @param in_options list of database names, conda environment, and database path
#' @param in_path where the config file lives. default cwd
#'
#' @return
#' @export write_config
#' @importFrom configr write.config
#' @examples TODO
write_config <- function(in_config, in_path=NULL, set_default=FALSE) {

  if ( is.null(in_path) ) {
    if (golem::app_dev()) {
      in_path <- golem::get_golem_wd()
    } else {
      in_path <- getwd()
    }
  }
  #TODO: check_config_list(in_options)
  CONFIG_FILE <- file.path(in_path,'app_config.yml')
  write.config(config.dat = in_config, file.path = CONFIG_FILE,
               write.type = "yaml", indent = 4)

  if (set_default){
    #write a new fallback config in omicser/inst
    CONFIG_FILE <- system.file('app_config.yml', package = 'omicser')
    write.config(config.dat = in_config, file.path = CONFIG_FILE,
                 write.type = "yaml", indent = 4)
  }
}



#TODO:  make a check_config_list <- function(in_options)
#         to ensure that we have everything or we are setting defaults
#
#TODO: check paths in 'app_config.yml' and spawn error if the paths are undefined

#' get_config: helper for reading the Config file containing information about the databases
#'
#' @param in_path where the config file lives. default getwd() test
#'
#' @return the list of options contained in `omicser_options.yml`
#' @export get_config
#' @importFrom configr read.config
#' @examples TODO
get_config <- function(in_path = NULL) {

  if ( is.null(in_path) ) {
    if (golem::app_dev()) {
      in_path <- golem::get_golem_wd()
      install_type <- "dev"
    } else {
      in_path <- getwd()
      install_type <- "configured"
    }
  }


  CONFIG_FILE <- file.path(in_path,'app_config.yml')
  if (file.exists(CONFIG_FILE)) {
    config_list <- read.config( file = CONFIG_FILE )
    install_type <- "configured"
  } else {
    #fallback to default.
    CONFIG_FILE <- system.file('inst/app/app_config.yml', package = 'omicser')
    message(paste("fallback to default:",CONFIG_FILE))
    if (CONFIG_FILE == "") {
      message("can't find config file: hacking defaults")
      config_list <- list(
        database_names=list(UNDEFINED="UNDEFINED"),
        db_root_path = "UNDEFINED",
        install = "hack"
        )
    } else {
      config_list <- read.config( file = CONFIG_FILE )
      install_type <- "default"
    }
  }
  config_list$install <- install_type
  #TODO: check_config_list(in_options)
  #
  return(config_list)
}
#
# #
#
# if {options()$golem.app.prod {
#
#
# } else {
#   golem::get_golem_wd()
#
# }
#
#
#
#   DB_NAMES <- get_golem_config("database_names")
#   DB_ROOT_PATH <- get_golem_config("db_root_path")
#
#   # get_golem_config <- function(
#   # value,
#   # config = Sys.getenv(
#   #   "GOLEM_CONFIG_ACTIVE",
#   #   Sys.getenv(
#   #     "R_CONFIG_ACTIVE",
#   #     "default"
#   #   )
#   # ),
#   # use_parent = TRUE
#   #
#   #   CONFIG <- omicser::get_config()
#   #   DB_NAMES <- CONFIG$database_names
#   #   DB_ROOT_PATH <- CONFIG$db_root_path
#   # PYTHON_ENV <- CONFIG$python_environment
#   # reticulate::use_virtualenv(
#   #     virtualenv = PYTHON_ENV,
#   #     required = TRUE)
#


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



