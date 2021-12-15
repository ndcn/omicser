
OMICSER_PYTHON <-  "pyenv_omicser" #set your python_env_name

# assuming mini-conda is installed
OMICSER_PYTHON_EXISTS <- any(reticulate::conda_list()["name"]==OMICSER_PYTHON)

if (!OMICSER_PYTHON_EXISTS){  #you should already have installed miniconda and created the env
  # simpler pip pypi install
  packages <- c("anndata")
  reticulate::conda_create(OMICSER_PYTHON, python_version = 3.8)
  reticulate::conda_install(envname=OMICSER_PYTHON,
                            # channel = "conda-forge",
                            pip = TRUE,
                            packages =  packages )
}

# set the `RETICULATE_PYTHON` env variable
if ( !Sys.getenv("RETICULATE_PYTHON")=="OMICSER_PYTHON" ) {
  Sys.setenv("RETICULATE_PYTHON"=reticulate::conda_python(envname = OMICSER_PYTHON))
}

