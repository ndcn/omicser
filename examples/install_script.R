

install.packages("devtools")
install.packages("reticulate")

reticulate::install_miniconda()
reticulate::conda_create("omxr",python_version = 3.9)
reticulate::conda_install(envname="omxr",
                          channel = "conda-forge",
                          packages = c("scanpy","leidenalg") )


# execute THIS:
devtools::install_github("ndcn/omicser")




