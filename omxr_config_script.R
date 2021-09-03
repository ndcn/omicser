

## this script sets up the list of datasets and data directory locatiosn for our data into
## the OMXR_databases.yml

dataset_names <- list(
  "Domenico DIA" = "domenico_stem_cell",
  "Vilas Microglia" = "vilas_microglia",
  #"Vilas Microglia (seu)" = "vilas_microglia_seu",
  "Vilas Microglia (sceasy)" = "vilas_microglia_sceasy",
  "Yassene Lipid concentraions & compositions" ="yassene_lipid",
  #"Yassene Lipid Concentrations" ="yassene_A_conc",
  #"Yassene Lipid Compositions" ="yassene_A_compos",
  "Oscar Microglia" ="oscar_microglia"
)


# where do our data files live. # WARNING do not use the ~ alias for Home
# MUST BE FULL or RELATIVE PATH... ~ will cause loading to
ds_root_path = 'data-raw'


# python environment
conda_environment = 'omxr'

omicser_options <- list(dataset_names=dataset_names,
                        ds_root_path=ds_root_path,
                        conda_environment=conda_environment)


configr::write.config(config.dat = omicser_options, file.path = "omxr_options.yml",
                      write.type = "yaml", indent = 4)



