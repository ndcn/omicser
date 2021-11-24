
# Set options here
options(golem.app.prod = FALSE) # TRUE = production mode, FALSE = development mode
# Detach all loaded packages and clean your environment
#golem::detach_all_attached()
#  DETACHING RETICULATE SEEMS TO CAUSE PROBLEMS...
#  so just detach omicser
#detach("package:omicser")
#rm(list=ls(all.names = TRUE))

# Document and reload your package

OMICSER_RUN_DIR <- getwd()
golem::document_and_reload( OMICSER_RUN_DIR )

# Run the application
OMICSER_PYTHON <-  "pyenv_omicser"
reticulate::use_condaenv(condaenv = OMICSER_PYTHON,
                         conda = reticulate::conda_binary(),
                         required = TRUE)


run_in_browser()
#run_app(options = list(launch.browser = TRUE)) #with = "shinyAppDir"

