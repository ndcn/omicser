
# Now move to the directory where you want to execute the omicser
# and make the omicser_options.yml

DEV_OMICSER <- TRUE


if (DEV_OMICSER){
  # this should be a full path... e.g. ~/Projects/NDCN_dev/omicser
  # but for github, we will set relative to the repo BASE
  REPO_PATH <- "/Users/ahenrie/Projects/NDCN_dev/omicser"
  OMICSER_RUN_DIR <- file.path(REPO_PATH,"quickstart")
  golem::document_and_reload(pkg = REPO_PATH)
} else {

  require(omicser)
  OMICSER_RUN_DIR <- file.path(REPO_PATH,"quickstart")

}

db_root_path <- "test_db"

database_names <- list(
  "Domenico DIA" = "domenico_stem_cell",
  "Vilas Microglia (sceasy)" = "vilas_microglia_sceasy",
  "Yassene Lipid concentraions & compositions" ="yassene_lipid"
)


conda_environment = 'omxr'
omicser_options <- list(database_names=database_names,
                        db_root_path=db_root_path,
                        conda_environment=conda_environment)

omicser::write_config(omicser_options,in_path = OMICSER_RUN_DIR )



