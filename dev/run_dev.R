# Set options here
options(golem.app.prod = FALSE) # TRUE = production mode, FALSE = development mode

# Detach all loaded packages and clean your environment
#golem::detach_all_attached()
#  DETACHING RETICULATE SEEMS TO CAUSE PROBLEMS...

#
# rm(list=ls(all.names = TRUE))
detach("package:omicser")

# Document and reload your package
golem::document_and_reload()

# Run the application
run_app(options = list(launch.browser = TRUE))
          # ,
          # port = "4243",
          # host = "127.0.0.1"))

