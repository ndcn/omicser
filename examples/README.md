# Example files for `omicser` development

If you would like to test the NDCN Omics Browser with example data,
we recommend starting with our example script for 
[PBMC data from 10 Genomics](`pbmc3k_curate_and_config.R`).

Additional scripts include:

- `install_script.R`: installs R and Python packages, creates environment 
- example curation scripts for various data types:
  - `proteomics_curate_and_config.R`: proteomics
  - `pbmc3k_curate_and_config.R`: single-cell transcriptomics
  - COMING SOON: `lipids_curate_and_config.R`: lipidomics
- `run_browser_script.R`: runs the NDCN omics browser via various incantations 

Example files include:
- `raw_data/` and `test_db/`: placeholders for locations to put data and database files, respectively


See also:
- `app_config.yml`: in the repo root.  This is an example configuration file for a browser with multiple databases

Files not listed here are used by project developers and are not intended for other users.
