

install.packages("devtools")
# install.packages("BiocManager")
# BiocManager::install("SingleCellExperiment")
install.packages("reticulate")

reticulate::install_miniconda()
reticulate::conda_create("omxr",python_version = 3.9)
reticulate::conda_install(envname="omxr",channel = "conda-forge",packages = c("scanpy","leidenalg") )

# execute THIS:

command <- "git clone git@github.com:ndcn/omicser.git"
system(command)
# move to repo
CWD <- getwd()
pkgload::load_all("omicser")

REPO_PATH <- fie.path(CWD,"omicser")


