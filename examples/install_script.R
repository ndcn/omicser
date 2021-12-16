#### Software installation ####

## These steps follow vignettes 1-Environment Setup and 2-Installation

# brew install cairo


# Step 1: install R dependencies -----
install.packages("devtools")
install.packages("reticulate")

# Step 2: make sure we have python available
# find out what python is available... make sure >= 3.6 (3.8 recommended)
#  we recommend using reticulate's installation of mini-conda as below

OMICSER_PYTHON <-  "pyenv_omicser"

BROWSER_ONLY <- FALSE
if (BROWSER_ONLY){
    packages <- c("anndata") # needed for browser
} else {
    packages <- c("scanpy","leidenalg")  #need scanpy for curation helpers, also installs anndata
}


reticulate::install_miniconda() #in case it is not already installed
reticulate::conda_create(OMICSER_PYTHON, python_version = 3.8)
reticulate::conda_install(envname = OMICSER_PYTHON,
                          pip = TRUE,
                         packages = packages )


if ( Sys.getenv("RETICULATE_PYTHON")!="OMICSER_PYTHON" ) {
    Sys.setenv("RETICULATE_PYTHON"=reticulate::conda_python(envname = OMICSER_PYTHON))
}

# install omicser package
devtools::install_github("ndcn/omicser")

library(omicser)



# # # TROUBLESHOOTING devel install / github clone
# OMICSER_RUN_DIR <- "/path/to/cloned/omicser" #or getcwd() if you are there
#
# devtools::build(OMICSER_RUN_DIR)
# devtools::install_local(OMICSER_RUN_DIR)
#
#
# require("golem")
# golem::document_and_reload(pkg = OMICSER_RUN_DIR)


#
# PYTHON_PATH <-"path/to/python" # EDIT THIS WITH YOUR PYTHON3.8
# # so we'll do everything in conda instead with conda...  possibly doing it in stages makes it faster or less prone to failure... scapyy 1.7.2
# cmd <- c("xxx","yyy")
# base::system( command = paste(cmd, collapse = " && ") )
# reticulate::virtualenv_create( envname = OMICSER_PYTHON,
#                                python = PYTHON_PATH,
#                                packages = packages   )
#
# reticulate::use_virtualenv(virtualenv = OMICSER_PYTHON,
#                            required = TRUE)
#
#
#



# # TROUBLESHOOTING virtual environments
#
# PYTHON_PATH <-"path/to/python" # EDIT THIS WITH YOUR PYTHON3.8
# # so we'll do everything in conda instead with conda...  possibly doing it in stages makes it faster or less prone to failure... scapyy 1.7.2
# cmd <- c("xxx","yyy")
# base::system( command = paste(cmd, collapse = " && ") )
# reticulate::virtualenv_create( envname = OMICSER_PYTHON,
#                                python = PYTHON_PATH,
#                                packages = packages   )
#
# reticulate::use_virtualenv(virtualenv = OMICSER_PYTHON,
#                            required = TRUE)
#
#
#
