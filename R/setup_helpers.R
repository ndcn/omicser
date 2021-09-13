
#' helper for Writing the Config file for setting up the omicser databases
#'
#' @param in_options
#' @param in_path where the config file lives. default cwd
#'
#' @return
#' @export write_config
#' @importFrom configr write.config
#' @examples TODO
write_config <- function(in_options, in_path=NULL) {
  # TODO: add a path arg?
  # this should update ad in place with the diff_exp data...
  # OPTIONS_FILE <- system.file('omicser_options.yml',
  #                             package = 'omicser')
  if ( is.null(in_path) ) {
    in_path <- getwd()
  }

  OPTIONS_FILE <- file.path(in_path,'omicser_options.yml')
  write.config(config.dat = in_options, file.path = OPTIONS_FILE,
                        write.type = "yaml", indent = 4)

}



#' helper for reading the Config file... currently from
#'
#' @param in_path where the config file lives. default cwd
#'
#' @return
#' @export get_config
#' @importFrom configr read.config
#' @examples TODO
get_config <- function(in_path = NULL) {

  if ( is.null(in_path) ) {
    in_path <- getwd()
  }
  # TODO:  can we make this more dynamic??
  # OPTIONS_FILE <- system.file('omicser_options.yml',
  #                             package = 'omicser')

  OPTIONS_FILE <- file.path(in_path,'omicser_options.yml')
  options_list <- read.config( file = OPTIONS_FILE )
  return(options_list)
}

