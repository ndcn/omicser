###### Lipidomics master table #########
#
# Create the lipidomics master table.
#
#
#
#
library(tidyverse)
library(readxl)
library(lipidtranslator)

# Read the lipidnames
# read method 1
method1_lipids <- read_xlsx(path = "./examples/Lipidyzer_lipids.xlsx",
                        sheet = "Method 1",
                        col_names = TRUE)

# read method 2
method2_lipids <- read_xlsx(path = "./examples/Lipidyzer_lipids.xlsx",
                        sheet = "Method 2",
                        col_names = TRUE)

# combine both methods
all_lipids <- rbind(method1_lipids, method2_lipids)


# Cleanup, remove :
#   - some empty lines
#   - lipids starting with 'd'
#   - lipids ending with '_2'
lipids <- all_lipids %>%
  filter(!is.na(ID),
         !str_detect(string = ID,
                     pattern = "^d.*$"),
         !str_detect(string = ID,
                     pattern = "^.*_2$"))

# Translate all lipids
master_table_lipids <- lipidyzer2lm(lipid_name = lipids$ID,
                                    get_ids = TRUE)


# Store as internal data
usethis::use_data(master_table_lipids,
                  internal = TRUE)
