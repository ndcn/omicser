

# read from the app_config.yml found in the current working directory
omicser::run_defaults()



OMICSER_RUN_DIR <- getwd()
#Run browser script
DB_ROOT_PATH <- file.path(OMICSER_RUN_DIR,"examples/test_db") #for example

# run in system default browser (chrome recommended) with databases found in DB_ROOT_PATH
omicser::run_in_browser(
           db_root = DB_ROOT_PATH,
           database_names = list(UNDEFINED="UNDEFINED"),
           install_type = "arguments"
   )


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

