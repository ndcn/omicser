

compute_de_table <- function(ad,comp_types, test_types, obs_names) {
  # this should update ad in place with the diff_exp data...

  diff_exp <- data.frame()
  for (obs_name in obs_names){
    for (test_type in test_types) {
      for (comp_type in comp_types) {
        parts <- strsplit(comp_type,"V")[[1]]
        reference <- parts[2]
        group <- parts[1]
        print(test_type)
        print(group)
        print(reference)

        if (reference == "rest") { #grpVrest
          key <- paste0(test_type,"_", comp_type)
          sc$tl$rank_genes_groups(ad,
                                  obs_name,
                                  method=test_type,
                                  key_added = key)
          #diff_exp <- list()
          de_table <- sc$get$rank_genes_groups_df(ad,
                                                  group=NULL,
                                                  key=key)
          de_table$comp_type <- comp_type
        } else { #compare group vs reference
          key <- paste0(test_type,"_", comp_type)
          sc$tl$rank_genes_groups(ad,
                                  obs_name,
                                  group = group,
                                  reference=reference,
                                  method=test_type,
                                  key_added = key)
          de_table <- sc$get$rank_genes_groups_df(ad,
                                                  group = group,
                                                  key=key)
          de_table$group <- group
          de_table$comp_type <- 'grpVref'
        }
        de_table$reference <- reference
        de_table$test_type <- test_type
        de_table$obs_name <- obs_name

        diff_exp <- dplyr::bind_rows(diff_exp, de_table)
      }
    }
  }

  return(diff_exp)
}

