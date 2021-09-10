

## this script sets up the list of datasets and data directory locatiosn for our data into
## the OMXR_databases.yml

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



