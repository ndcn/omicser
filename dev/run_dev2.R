
# Set options here
options(golem.app.prod = FALSE) # TRUE = production mode, FALSE = development mode
# Detach all loaded packages and clean your environment
#golem::detach_all_attached()
#  DETACHING RETICULATE SEEMS TO CAUSE PROBLEMS...
#  so just detach omicser
#detach("package:omicser")
#rm(list=ls(all.names = TRUE))

# Document and reload your package
golem::document_and_reload( "/Users/ahenrie/Projects/NDCN_dev/omicser" )

setwd("quickstart")
# Run the application
run_in_browser()
#run_app(options = list(launch.browser = TRUE)) #with = "shinyAppDir"

