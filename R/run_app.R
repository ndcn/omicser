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
#' @return
#' @export
#'
#' @examples TODO
run_in_browser <- function(){
  run_app(options = list(launch.browser = TRUE)) #with = "shinyAppDir"
}

#'
#'
#' #' run_defai;ts - shortcut to launch with default settings.
#' #'
#' #' @return
#' #' @export
#' #'
#' #' @examples TODO
#' run_defaults <- function(){
#'   # load yaml and add launch.browser = TRUE to the list
#'   #
#'   #
#'   #
#'
#'
#'   if
#'
#'     default_options <-  list(launch.browser = TRUE,
#'                            # the rest are options to pass too the app.
#'                            )
#'
#'   run_app(options =) #with = "shinyAppDir"
#' }
