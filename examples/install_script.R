#### Software installation ####

## These steps follow vignettes 1-Environment Setup and 2-Installation

# install R dependencies
install.packages("devtools")
install.packages("reticulate")

# install and configure environment
reticulate::install_miniconda()
reticulate::conda_create("omxr", python_version = 3.9)
reticulate::conda_install(envname="omxr",
                          channel = "conda-forge",
                          packages = c("scanpy","leidenalg") )

# install omicser package
devtools::install_github("ndcn/omicser")
