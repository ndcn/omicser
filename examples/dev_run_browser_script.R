
OMICSER_PYTHON <-  "pyenv_omicser"

# assume CONDA is installed
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

if ( !Sys.getenv("RETICULATE_PYTHON")=="OMICSER_PYTHON" ) {
  Sys.setenv("RETICULATE_PYTHON"=reticulate::conda_python(envname = OMICSER_PYTHON))
}


# check that we have our python on deck
reticulate::py_discover_config()



#

OMICSER_RUN_DIR <- getwd()
golem::document_and_reload(pkg = OMICSER_RUN_DIR)

# read from the app_config.yml found in the current working directory
omicser::run_defaults()



DB_ROOT_PATH <- "/Users/ergonyc/Projects/NDCN_dev/testing/omxr/databases"
# run in system default browser (chrome recommended) with databases found in DB_ROOT_PATH
omicser::run_in_browser(
  db_root = DB_ROOT_PATH,
  database_names = list(UNDEFINED="UNDEFINED"),
  install_type = "arguments"
)





OMICSER_RUN_DIR <- getwd()
#Run browser script
DB_ROOT_PATH <- file.path(OMICSER_RUN_DIR,"examples/test_db") #for example



# run in system default browser with databases found in DB_ROOT_PATH
omicser::run_in_browser(
  db_root = DB_ROOT_PATH,
  database_names = list(`Domenico DIA`="domenico_stem_cell"),
  install_type = "arguments"
)


# run in the rstudio browser with databases found in DB_ROOT_PATH, and modal database choice
omicser::run_app(db_root = DB_ROOT_PATH,
                   database_names = list(UNDEFINED="UNDEFINED"),
                   install_type = "INSTALL_TYPE"
                   )



# run in system default browser with modal choice for defining paths to db_root and databases
omicser::run_in_browser(
  db_root = "UNDEFINED",
  database_names = list(UNDEFINED="UNDEFINED")
)

