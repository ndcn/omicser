#' Run the Shiny Application
#'
#' @param ... arguments to pass to golem_opts.
#' See `?golem::get_golem_options` for more details.
#' @inheritParams shiny::shinyApp
#'
#' @export
#' @importFrom shiny shinyApp
#' @importFrom golem with_golem_options
run_app <- function(
  onStart = NULL,
  options = list(),
  enableBookmarking = NULL,
  uiPattern = "/",
  ...
) {
  with_golem_options(
    app = shinyApp(
      ui = app_ui,
      server = app_server,
      onStart = onStart,
      options = options,
      enableBookmarking = enableBookmarking,
      uiPattern = uiPattern
    ),
    golem_opts = list(...)
  )
}

# runApp <- function(appDir=getwd(),
#                    port=getOption('shiny.port'),
#                    launch.browser = getOption('shiny.launch.browser', interactive()),
#                    host=getOption('shiny.host', '127.0.0.1'),
#                    workerId="", quiet=FALSE,
#                    display.mode=c("auto", "normal", "showcase"),
#                    test.mode=getOption('shiny.testmode', FALSE)) {

#' run_in_browser - shortcut to launch in_browser
#'
#' @param db_root path of the database root
#' @param database_names database names
#' @param install_type installation type
#' @param ... parameters passes on to run_app()
#'
#' @return
#' @export
#'
#' @examples TODO
run_in_browser <- function(db_root=NULL,database_names=NULL,install_type = "empty",...){

  if (is.null(db_root)) db_root <- "UNDEFINED"
  if (is.null(database_names)) database_names <- list(UNDEFINED="UNDEFINED")

  run_app(options = list(launch.browser = TRUE),
          db_root=db_root,
          database_names=database_names,
          install_type=install_type,
          ...) #with = "shinyAppDir"
}



#' run_defai;ts - shortcut to launch with default settings.
#'
#' @return
#' @export
#'
#' @examples TODO
run_defaults <- function(){
  # load yaml and add launch.browser = TRUE to the list
  #
  #
  #
  CONFIG <- omicser::get_config()
  DB_NAMES <- CONFIG$database_names
  DB_ROOT <- CONFIG$db_root_path
  INSTALL_TYPE <- CONFIG$install #production, dev,

   # TODO:  assert ptython settings
  #reticulate::py_module_available("anndata")

  run_app(options = list(launch.browser = TRUE),
          db_root = DB_ROOT,
          database_names = DB_NAMES,
          install_type = INSTALL_TYPE
          ) #with = "shinyAppDir"
}

#'  run_dockers - shortcut to launch with default settings for the docker container
#'
#' @return
#' @export
#'
#' @examples TODO
run_docker <- function(){
  # load yaml and add launch.browser = TRUE to the list
  #
  #
  #
  CONFIG <- omicser::get_config()
  DB_NAMES <- CONFIG$database_names
  DB_ROOT <- CONFIG$db_root_path
  INSTALL_TYPE <- CONFIG$install #production, dev,

  # TODO:  assert ptython settings
  #reticulate::py_module_available("anndata")

  run_app(#options = list(launch.browser = FALSE),
    db_root = DB_ROOT,
    database_names = DB_NAMES,
    install_type = INSTALL_TYPE
  ) #with = "shinyAppDir"
}
