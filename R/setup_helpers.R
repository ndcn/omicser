
#' helper for Writing the Config file for setting up the omicser databases
#'
#' @param
#' @param in_options
#'
#' @return
#' @export write_config
#' @importFrom configr write.config
#' @examples TODO
write_config <- function(in_options) {
  # this should update ad in place with the diff_exp data...
  OPTIONS_FILE <- system.file('omicser_options.yml',
                              package = 'omicser')

  configr::write.config(config.dat = in_options, file.path = OPTIONS_FILE,
                        write.type = "yaml", indent = 4)
}
