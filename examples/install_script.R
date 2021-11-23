#### Software installation ####

## These steps follow vignettes 1-Environment Setup and 2-Installation

# Step 1: install R dependencies -----
install.packages("devtools")
install.packages("reticulate")

# Step 2: make sure we have python available
# find out what python is available... make sure >= 3.6 (3.8 recommended)
#  we recommend using reticulate's installation of mini-conda as below
USE_CONDA <- TRUE

if (FALSE) {
  # which python
  # install : download link?
  # conda or virtual environment...

}

OMICSER_PYTHON <-  "pyenv_omicser"

BROWSER_ONLY <- FALSE
if (BROWSER_ONLY){
    packages <- c("anndata") # needed for browser
} else {
    packages <- c("scanpy[leiden]")  #need scanpy for curation helpers, also installs anndata
}

USE_CONDA <- TRUE  #change this if you prefer virtual environments
if (USE_CONDA){
    reticulate::install_miniconda() #in case it is not already installed
    reticulate::conda_create(OMICSER_PYTHON, python_version = 3.8)
    reticulate::conda_install(envname=OMICSER_PYTHON,
                            channel = "conda-forge",
                            packages = packages )

    reticulate::use_condaenv(condaenv = OMICSER_PYTHON,
                                    conda = reticulate::conda_binary(),
                                    required = TRUE)


} else {
    PYTHON_PATH <-"path/to/python" # EDIT THIS WITH YOUR PYTHON3.8
    # so we'll do everything in conda instead with conda...  possibly doing it in stages makes it faster or less prone to failure... scapyy 1.7.2
    cmd <- c("xxx","yyy")
    base::system( command = paste(cmd, collapse = " && ") )
    reticulate::virtualenv_create( envname = OMICSER_PYTHON,
                               python = PYTHON_PATH,
                               packages = packages   )

    reticulate::use_virtualenv(virtualenv = OMICSER_PYTHON,
                            required = TRUE)
}




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
reticulate::conda_create("conda38", python_version = 3.8)
reticulate::conda_install(envname="omxr",
                          channel = "conda-forge",
                          packages = c("scanpy[leiden]") )

# install omicser package
devtools::install_github("ndcn/omicser")




