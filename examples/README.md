# Example files for `omicser` development

If you would like to test the NDCN Omics Browser with example data,
we recommend starting with our example script for 
[PBMC data from 10 Genomics](`pbmc3k_curate_and_config.R`).

Additional scripts include:

- `install_script.R`: installs R and Python packages, creates environment 
- `assert_python_script.R`: short script for asserting that we are using the configured python environment via `reticulate`
- example curation scripts for various data types:
  - `proteomics_curate_and_config.R`: proteomics
  - `pbmc3k_curate_and_config.R`: single-cell transcriptomics
  - COMING SOON: `lipids_curate_and_config.R`: lipidomics
- `run_browser_script.R`: runs the NDCN omics browser via various incantations 

Example files include:
- `raw_data/` and `test_db/`: placeholders for locations to put data and database files, respectively


See also:
- `app_config.yml`: This is an example configuration file for a browser with multiple databases
- `dev_run_browser_script.R`: shows how to load the repo locally and run the browser, as in a _development_ installation via a cloned repo
- `download_examples_script.R`:  short script which will download and extract just the `examples/` directory from the `omicser` directory.

A short [installation & configuration tutorial video](https://www.youtube.com/watch?v=lwJmsxk0vTU) is also available which roughly follows this sequence of vignettes: 

Files not listed here are used by project developers and are not intended for other users.
