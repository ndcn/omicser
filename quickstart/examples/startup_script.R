

install.packages("devtools")
install.packages("BiocManager")
BiocManager::install("SingleCellExperiment")
install.packages("reticulate")

reticulate::install_miniconda()
reticulate::conda_create("omxr",python_version = 3.9)
reticulate::conda_install(envname = "omxr", packages = "scanpy")
reticulate::conda_install(envname="omxr",channel = "conda-forge",packages = c("leidenalg") )


# execute THIS:
devtools::install_github("ergonyc/omicser")

# or if we are doing a DEV installation.  make sure you are in the Omicser based directory and either
# golem::document_and_reload()
# assume we are in the or
getwd()

pkgload::load_all(getwd())
# shortcuts in-case things need to be tweaked
#detach("package:omicser")
#remove.packages("omicser")



# place the configureation file `omicser_optiions.yml` in the quickstart folder of the repo
setwd("quickstart" )
getwd()
#=======================================================
database_names <- list(
  "Domenico DIA" = "domenico_stem_cell",
  "Vilas Microglia" = "vilas_microglia",
  #"Vilas Microglia (seu)" = "vilas_microglia_seu",
  "Vilas Microglia (sceasy)" = "vilas_microglia_sceasy",
  "Yassene Lipid concentraions & compositions" ="yassene_lipid"
  #"Yassene Lipid Concentrations" ="yassene_A_conc",
  #"Yassene Lipid Compositions" ="yassene_A_compos",
  #"Oscar Microglia" ="oscar_microglia"
)

# where do our data files live. # WARNING do not use the ~ alias for Home
# MUST BE FULL or RELATIVE PATH to where this is executed...... ~ will cause loading to
db_root_path = 'test_db'


# python environment
conda_environment = 'omxr'

omicser_options <- list(database_names=database_names,
                        db_root_path=db_root_path,
                        conda_environment=conda_environment)

omicser::write_config(omicser_options)



# CURATE DATA HERE =======================================================
#source("examples/curate_domenico_stem_cell.R")



# CURATE DATA HERE =======================================================


# NOTE this shoudl be run in the quickstart directory...
omicser::run_app(options = list(launch.browser = TRUE))
