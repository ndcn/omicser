
#
library("omicser")


OMICSER_RUN_DIR <- "/path/to/folder_with_app_config"
setwd(OMICSER_RUN_DIR)


# 1. read from the app_config.yml found in the current working directory
# run in system default browser (chrome recommended) with databases found in DB_ROOT_PATH
omicser::run_defaults()


# 2. call `run_in_browser` with arguments specifying where to find databases
DB_ROOT_PATH <- OMICSER_RUN_DIR <- file.path(OMICSER_RUN_DIR,"databases") #/path/to/databases



omicser::run_in_browser(
  db_root = DB_ROOT_PATH,
  database_names = list(UNDEFINED="UNDEFINED"),
  install_type = "arguments"
)

#using UNDEFINED="UNDEFINED" as the database_names will force a modal to load to select them



# run in system default browser with databases found in DB_ROOT_PATH
omicser::run_in_browser(
  db_root = DB_ROOT_PATH,
  database_names = list(`Domenico DIA`="domenico_stem_cell"),
  install_type = "arguments"
)


# 3. call `run_app` from R Studio will
#         run in the rstudio browser with databases found in DB_ROOT_PATH, and modal database choice
omicser::run_app(db_root = DB_ROOT_PATH,
                 database_names = list(UNDEFINED="UNDEFINED"),
                 install_type = "INSTALL_TYPE"
)


# 4. call `run_in_browser`with arguments set to "UNDEFINED" will fall back to modal
omicser::run_in_browser(
  db_root = "UNDEFINED",
  database_names = list(UNDEFINED="UNDEFINED")
)
