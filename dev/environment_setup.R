
install.packages("devtools")
install.packages("BiocManager")
BiocManager::install("SingleCellExperiment")


install.packages("reticulate")

require(reticulate)
reticulate::install_miniconda()
reticulate::conda_create("omxr",python_version = 3.9)
reticulate::conda_install(envname = "omxr", packages = "scanpy")
reticulate::conda_install(envname="omxr",channel = "conda-forge",packages = c("leidenalg") )

# indicate that we want to use a specific condaenv
use_condaenv("omxr")
