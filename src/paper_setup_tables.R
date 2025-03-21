# Setup tables code

tab_row_names <- c("Width", "Length", "Density", "Mean distance", 
                   "Count per segment", "Fraction")
tab_col_names <- c("All sensilla", "Olfactory sensilla", "Mechanosensory sensilla")
tab_signal_df <- data.frame("all" = c("n/a", # width
                                      "n/a", # length
                                      paste0(test_all_density_signal$all_phylo$Pf, " (",
                                             test_all_density_signal$all_phylo$n, ")"),  # density
                                      paste0(test_all_mean_dists_signal$all_phylo$Pf," (",
                                             test_all_mean_dists_signal$all_phylo$n, ")"), # mean_dists
                                      paste0(test_all_counts_signal$all_phylo$Pf, " (",
                                             test_all_counts_signal$all_phylo$n, ")"),  # counts
                                      "n/a"  # fraction
), 
"olfactory" = c(paste0(test_widths_signal$olfactory_phylo$Pf, " (",
                       test_widths_signal$olfactory_phylo$n, ")"), # width
                paste0(test_lengths_signal$olfactory_phylo$Pf, " (",
                       test_lengths_signal$olfactory_phylo$n, ")"), # length
                paste0(test_all_density_signal$olfactory_phylo$Pf, " (",
                       test_all_density_signal$olfactory_phylo$n, ")"),  # density
                paste0(test_all_mean_dists_signal$olfactory_phylo$Pf," (",
                       test_all_mean_dists_signal$olfactory_phylo$n, ")"), # mean_dists
                paste0(test_all_counts_signal$olfactory_phylo$Pf, " (",
                       test_all_counts_signal$olfactory_phylo$n, ")"),  # counts
                paste0(test_all_fractions_signal$olfactory_phylo$Pf, " (",
                       test_all_fractions_signal$olfactory_phylo$n, ")")  # fraction
),
"mechano" = c(paste0(test_widths_signal$mechano_phylo$Pf, " (",
                     test_widths_signal$mechano_phylo$n, ")"), # width
              paste0(test_lengths_signal$mechano_phylo$Pf, " (",
                     test_lengths_signal$mechano_phylo$n, ")"), # length
              paste0(test_all_density_signal$mechano_phylo$Pf, " (",
                     test_all_density_signal$mechano_phylo$n, ")"),  # density
              paste0(test_all_mean_dists_signal$mechano_phylo$Pf," (",
                     test_all_mean_dists_signal$mechano_phylo$n, ")"), # mean_dists
              paste0(test_all_counts_signal$mechano_phylo$Pf, " (",
                     test_all_counts_signal$mechano_phylo$n, ")"),  # counts
              paste0(test_all_fractions_signal$mechano_phylo$Pf, " (",
                     test_all_fractions_signal$mechano_phylo$n, ")")  # fraction
)
)

tab_sex_df <- data.frame("all" = c("n/a", # width
                                   "n/a", # length
                                   paste0(test_all_density_sex$all_phylo$Pf, " (",
                                          test_all_density_sex$all_phylo$n, ")"),  # density
                                   paste0(test_all_mean_dists_sex$all_phylo$Pf," (",
                                          test_all_mean_dists_sex$all_phylo$n, ")"), # mean_dists
                                   paste0(test_all_counts_sex$all_phylo$Pf, " (",
                                          test_all_counts_sex$all_phylo$n, ")"),  # counts
                                   "n/a"  # fraction
), 
"olfactory" = c(paste0(test_widths_sex$olfactory_phylo$Pf, " (",
                       test_widths_sex$olfactory_phylo$n, ")"), # width
                paste0(test_lengths_sex$olfactory_phylo$Pf, " (",
                       test_lengths_sex$olfactory_phylo$n, ")"), # length
                paste0(test_all_density_sex$olfactory_phylo$Pf, " (",
                       test_all_density_sex$olfactory_phylo$n, ")"),  # density
                paste0(test_all_mean_dists_sex$olfactory_phylo$Pf," (",
                       test_all_mean_dists_sex$olfactory_phylo$n, ")"), # mean_dists
                paste0(test_all_counts_sex$olfactory_phylo$Pf, " (",
                       test_all_counts_sex$olfactory_phylo$n, ")"),  # counts
                paste0(test_all_fractions_sex$olfactory_phylo$Pf, " (",
                       test_all_fractions_sex$olfactory_phylo$n, ")")  # fraction
),
"mechano" = c(paste0(test_widths_sex$mechano_phylo$Pf, " (",
                     test_widths_sex$mechano_phylo$n, ")"), # width
              paste0(test_lengths_sex$mechano_phylo$Pf, " (",
                     test_lengths_sex$mechano_phylo$n, ")"), # length
              paste0(test_all_density_sex$mechano_phylo$Pf, " (",
                     test_all_density_sex$mechano_phylo$n, ")"),  # density
              paste0(test_all_mean_dists_sex$mechano_phylo$Pf," (",
                     test_all_mean_dists_sex$mechano_phylo$n, ")"), # mean_dists
              paste0(test_all_counts_sex$mechano_phylo$Pf, " (",
                     test_all_counts_sex$mechano_phylo$n, ")"),  # counts
              paste0(test_all_fractions_sex$mechano_phylo$Pf, " (",
                     test_all_fractions_sex$mechano_phylo$n, ")")  # fraction
)
)

rownames(tab_signal_df) <- tab_row_names
rownames(tab_sex_df) <- tab_row_names

