### install the python virtual enivronment #################################################
# Taken from the 'omicser' package install_script.R


### settings ########################################
# python packages to be installed
# only needs packages for running the app in the browser, no curation needs to be done
packages <- c("anndata")

# set the name of the python environment
OMICSER_PYTHON <-  "pyenv_omicser"

### installation #####################################
# install miniconda
reticulate::install_miniconda()
# create conda envrionment
reticulate::conda_create(OMICSER_PYTHON, python_version = 3.8)
# install python packages
reticulate::conda_install(envname=OMICSER_PYTHON,
                          channel = "conda-forge",
                          packages = packages )
# set which python environment to run
reticulate::use_condaenv(condaenv = OMICSER_PYTHON,
                         conda = reticulate::conda_binary(),
                         required = TRUE)

