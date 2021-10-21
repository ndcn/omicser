

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






CONDA_ENV <- "omxr" #"/Users/ahenrie/Library/r-miniconda/envs/omxr/bin/python"
CONDA_ENV <-"/Users/ahenrie/Library/r-miniconda/envs/omxr" # avoide sicne we have two omxr








reticulate::use_condaenv( required = TRUE,
                          condaenv = CONDA_ENV,
                          conda = CONDA_EXE)
#conirm its active.
reticulate::py_config()
reticulate::py_exe()

renv::snapshot()
