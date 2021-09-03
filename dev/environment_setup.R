

install.packages("BiocManager")
BiocManager::install("SingleCellExperiment")


install.packages("reticulate")

require(reticulate)

install_miniconda(path = miniconda_path(), update = TRUE, force = FALSE)

# create a new environment
conda_create("omicser",python_version = 3.9)

# install SciPy
conda_install("omicser", "scanpy")


# indicate that we want to use a specific condaenv
use_condaenv("omicser")

# import SciPy (will use "r-reticulate" as per call to use_condaenv)
scipy <- import("scanpy")

# import SciPy (it will be automatically discovered in "r-reticulate")
scipy <- import("scipy")
