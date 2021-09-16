

install.packages("devtools")
# install.packages("BiocManager")
# BiocManager::install("SingleCellExperiment")
install.packages("reticulate")

reticulate::install_miniconda()
reticulate::conda_create("omxr",python_version = 3.9)
reticulate::conda_install(envname = "omxr", packages = "scanpy")
reticulate::conda_install(envname="omxr",channel = "conda-forge",packages = c("leidenalg") )

# execute THIS:
devtools::install_github("ergonyc/omicser")




