#### Software installation ####

## These steps follow vignettes 1-Environment Setup and 2-Installation

# install R dependencies
install.packages("devtools")
install.packages("reticulate")


# find out what python is available... make sure <= 3.6



setwd(OMICSER_RUN_DIR)
anndata <- anndata::read_h5ad( filename=file.path(DB_ROOT_PATH,DB_NAME,"db_data.h5ad" ) )

PYTHON_ENV <-  "./pyenv"
reticulate::virtualenv_create( envname = PYTHON_ENV,
                               python = "/usr/local/bin/python3",
                               packages = c("anndata")    ) # "scanpy"

reticulate::virtualenv_exists( PYTHON_ENV )

reticulate::use_virtualenv(virtualenv = PYTHON_ENV,
                           required = TRUE)

# Document and reload your package
golem::document_and_reload( "/Users/ahenrie/Projects/NDCN_dev/omicser" )


setwd("quickstart")



# install and configure environment
reticulate::install_miniconda()
reticulate::conda_create("omxr", python_version = 3.9)
reticulate::conda_install(envname="omxr",
                          channel = "conda-forge",
                          packages = c("scanpy","leidenalg") )

# install omicser package
devtools::install_github("ndcn/omicser")




