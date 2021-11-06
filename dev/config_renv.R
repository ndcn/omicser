

# make sure we are in the repo root
setwd("/Users/ahenrie/Projects/NDCN_dev/omicser")

# makesure we have renv
install.packages("renv")

renv::init()
# install R dependencies
#renv::install("devtools")
#renv::install("reticulate")
#renv::install("BiocManager")
#renv::install("rematch")



# VENV TROUBLE INSTALLING SCANPY
# renv::use_python(type = "virtualenv")
# # choose option 1:  /usr/local/Cellar/python@3.9/3.9.7_1/bin/python3.9
# # with VIRTUAL ENVIRONMENTS..replacing CONDA
# VENV <- "renv/python/virtualenvs/renv-python-3.9"
#
# reticulate::virtualenv_exists(VENV)
# reticulate::virtualenv_python(envname = VENV)
# #FAILS HERE.... BACK TO CONDA
# reticulate::virtualenv_install( envname = VENV,
#                                 packages = "scanpy")




# # install and configure environment
# reticulate::conda_create(envname=CONDA_ENV, conda=CONDA_EXE,python_version = 3.9)
# reticulate::conda_install(envname=CONDA_ENV,
#                           channel = "conda-forge",
#                           packages = c("scanpy","leidenalg") )


CONDA_ENV <- "omxr" #"/Users/ahenrie/Library/r-miniconda/envs/omxr/bin/python"
CONDA_ENV <-"/Users/ahenrie/Library/r-miniconda/envs/omxr" # avoide sicne we have two omxr
CONDA_EXE <- "/Users/ahenrie/Library/r-miniconda/bin/conda" #reticulate::conda_binary()

#CONDA_PY_PATH: "/Users/ahenrie/Library/r-miniconda/envs/omxr/bin/python"
renv::use_python(python = "/Users/ahenrie/Library/r-miniconda/envs/omxr/bin/python", type = "conda")
# can i encapsulate a miniconda in the environment???
#reticulate::install_miniconda()



# conda environment is now encapsulated in teh renv...
CONDA_ENV <-"renv/python/condaenvs/renv-python" # avoide sicne we have two omxr

reticulate::conda_install(envname=CONDA_ENV,
                          channel = "conda-forge",
                          packages = c("scanpy","leidenalg") )

#configure python
# "/opt/anaconda3/envs/base39/bin/python"
py <- reticulate::conda_list()$python[7]


#set up the python environment this way...
# this is the right way to do it, but reenv can't handle the mix of python and conda...
cmd <- c("/Users/ahenrie/Library/r-miniconda/bin/conda create -y -p /Users/ahenrie/Projects/NDCN_dev/omicser/.pyenv python=3.9 pip",
          "/Users/ahenrie/Library/r-miniconda/bin/conda install -p /Users/ahenrie/Projects/NDCN_dev/omicser/.pyenv -y seaborn scikit-learn statsmodels numba pytables ",
          "/Users/ahenrie/Library/r-miniconda/bin/conda install  -p /Users/ahenrie/Projects/NDCN_dev/omicser/.pyenv -y -c conda-forge python-igraph leidenalg",
          "/Users/ahenrie/Projects/NDCN_dev/omicser/.pyenv/bin/python -m pip install scanpy")


# so we'll do everything in conda instead with conda...  possibly doing it in stages makes it faster or less prone to failure... scapyy 1.7.2
cmd <- c("/Users/ahenrie/Library/r-miniconda/bin/conda create -y -p /Users/ahenrie/Projects/NDCN_dev/omicser/.pyenv python=3.9 pip",
         "/Users/ahenrie/Library/r-miniconda/bin/conda install -p /Users/ahenrie/Projects/NDCN_dev/omicser/.pyenv -y seaborn scikit-learn statsmodels numba pytables ",
         "/Users/ahenrie/Library/r-miniconda/bin/conda install  -p /Users/ahenrie/Projects/NDCN_dev/omicser/.pyenv -y -c conda-forge python-igraph leidenalg",
         "/Users/ahenrie/Library/r-miniconda/bin/conda install  -p /Users/ahenrie/Projects/NDCN_dev/omicser/.pyenv -y -c bioconda scanpy")


base::system(command = paste(cmd, collapse = " && "))


#

CONDA_ENV <- "" #"/Users/ahenrie/Library/r-miniconda/envs/omxr/bin/python"
CONDA_ENV <-"/Users/ahenrie/Projects/NDCN_dev/omicser/.pyenv" # avoide sicne we have two omxr
CONDA_EXE <- reticulate::conda_binary()
reticulate::use_condaenv( required = TRUE,
                          condaenv = CONDA_ENV,
                          conda = CONDA_EXE)
#conirm its active.
#

# virtualevn is NOT working because the pip installations on Mac need special compiler flags for scanpy...
reticulate::conda_create(envname = "./.pyenv",python_version = 3.9,
                              packages = c("pip","wheel","setuptools") )

reticulate::use_condaenv("./.pyenv")

reticulate::virtualenv_install(envname = "./.pyenv", packages = c("'scanpy[leiden]'") )


# virtualevn is NOT working because the pip installations on Mac need special compiler flags for scanpy...
reticulate::virtualenv_create(envname = "./.pyenv",
                              python = "python3.9",
                              packages = c("pip","wheel","setuptools"))

reticulate::use_virtualenv("./.pyenv")

reticulate::virtualenv_install(envname = "./.pyenv", packages = c("'scanpy[leiden]'") )


reticulate::py_config()
reticulate::py_exe()

renv::snapshot()
