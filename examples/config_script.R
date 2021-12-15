
# Now move to the directory where you want to execute the omicser
# and make the omicser_options.yml

DEV_OMICSER <- TRUE


if (DEV_OMICSER){
  REPO_PATH <- getwd()  #/path/to/cloned/repo
  golem::document_and_reload(pkg = REPO_PATH)
} else {

  library(omicser)

}

OMICSER_RUN_DIR <- file.path(REPO_PATH,"examples") # "/path/to/folder_with_app_config"

DB_ROOT_PATH <- OMICSER_RUN_DIR <- file.path(OMICSER_RUN_DIR,"databases") #/path/to/databases

# for example the databases curated with `pbmc3k_curate_and_config.R` & `proteomics_curate_and_config.R`
database_names <- list(
  "Domenico DIA" = "domenico_stem_cell",
  "mypbmc" = "pbmc3k",
)


omicser_options <- list(database_names=database_names,
                        db_root_path=DB_ROOT_PATH,
                        install="configured")

#write the `app_config.yml` into our OMICSER_RUN_DIR (e.g. getwd())
omicser::write_config(omicser_options,in_path = OMICSER_RUN_DIR )



