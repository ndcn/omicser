

install.packages("devtools")
install.packages("BiocManager")
BiocManager::install("SingleCellExperiment")
install.packages("reticulate")

reticulate::install_miniconda()
reticulate::conda_create("omxr",python_version = 3.9)
reticulate::conda_install(envname = "omxr", packages = "scanpy")
reticulate::conda_install(envname="omxr",channel = "conda-forge",packages = c("leidenalg") )



devtools::install_github("ergonyc/omicser")

# shortcuts incase things need to be tweaked
#detach("package:omicser")
#remove.packages("omicser")



#=======================================================
# example with domenicao datast
## this script sets up the list of datasets and data directory locatiosn for our data into
## the OMXR_databases.yml

# > files to curate a dataset have been put into a "domenico_stem_cell" folder which is a
#   subdirectory of "omicsdata"  in the my current working directory (.e.g. ~/spin-up )
#


dataset_names <- list(
  "Domenico DIA" = "domenico_stem_cell"
)


# where do our data files live. # WARNING do not use the ~ alias for Home
# MUST BE FULL or RELATIVE PATH... ~ will cause loading to

# create the omxr_options.yml ==============================
# could also edit the .yml directly
ds_root_path = 'omicsdata'

# python environment
conda_environment = 'omxr'

omicser_options <- list(dataset_names=dataset_names,
                        ds_root_path=ds_root_path,
                        conda_environment=conda_environment)

require(configr)
configr::write.config(config.dat = omicser_options, file.path = "omxr_options.yml",
                      write.type = "yaml", indent = 4)


#=======================================================
# curate data
# look at this config fileh



#source("examples/curate_domenico_stem_cell.R")


library(omicser)
omicser::run_app(options = list(launch.browser = TRUE))
