# Setup tables code


tab_row_names <- c("Width", "Length", "Density", "Mean distance", 
                   "Count per segment", "Fraction")
tab_col_signal_names <- c("Mean value (visual)", "Mean value (chemical)",
                          "Sample size","F value", "p value")
tab_col_sex_names <- c("Mean value (female)", "Mean value (male)",
                          "Sample size","F value", "p value")

#### Signal Testing ####

tab_signal_all_df <- data.frame(
  "visual_mean_value" = c("n/a", # width
                        "n/a", # length
                        paste(signif(mean(test_all_density_signal$bar_means$means[ 
                          test_all_density_signal$bar_means$signal == "visual" &
                            test_all_density_signal$bar_means$sensillum_type == "all"], 
                          na.rm = T), digits = 3), "±", signif(sd(test_all_density_signal$bar_means$means[
                            test_all_density_signal$bar_means$signal == "visual" &
                              test_all_density_signal$bar_means$sensillum_type == "all"], 
                            na.rm = T), digits = 1)), # density
                        paste(signif(mean(test_all_mean_dists_signal$bar_means$means[ 
                          test_all_mean_dists_signal$bar_means$signal == "visual" &
                            test_all_mean_dists_signal$bar_means$sensillum_type == "all"], 
                          na.rm = T), digits = 3), "±", signif(sd(test_all_mean_dists_signal$bar_means$means[
                            test_all_mean_dists_signal$bar_means$signal == "visual" &
                              test_all_mean_dists_signal$bar_means$sensillum_type == "all"], 
                            na.rm = T), digits = 1)), # mean_dists
                        paste(signif(mean(test_all_counts_signal$bar_means$means[ 
                          test_all_counts_signal$bar_means$signal == "visual" &
                            test_all_counts_signal$bar_means$sensillum_type == "all"], 
                          na.rm = T), digits = 3), "±", signif(sd(test_all_counts_signal$bar_means$means[
                            test_all_counts_signal$bar_means$signal == "visual" &
                              test_all_counts_signal$bar_means$sensillum_type == "all"], 
                            na.rm = T), digits = 1)), # counts
                        "n/a"  # fraction
  ),
  "chemical_mean_value" = c("n/a", # width
                          "n/a", # length
                          paste(signif(mean(test_all_density_signal$bar_means$means[ 
                            test_all_density_signal$bar_means$signal != "visual" &
                              test_all_density_signal$bar_means$sensillum_type == "all"], 
                            na.rm = T), digits = 3), "±", signif(sd(test_all_density_signal$bar_means$means[
                              test_all_density_signal$bar_means$signal != "visual" &
                                test_all_density_signal$bar_means$sensillum_type == "all"], 
                              na.rm = T), digits = 1)), # density
                          paste(signif(mean(test_all_mean_dists_signal$bar_means$means[ 
                            test_all_mean_dists_signal$bar_means$signal != "visual" &
                              test_all_mean_dists_signal$bar_means$sensillum_type == "all"], 
                            na.rm = T), digits = 3), "±", signif(sd(test_all_mean_dists_signal$bar_means$means[
                              test_all_mean_dists_signal$bar_means$signal != "visual" &
                                test_all_mean_dists_signal$bar_means$sensillum_type == "all"], 
                              na.rm = T), digits = 1)), # mean_dists
                          paste(signif(mean(test_all_counts_signal$bar_means$means[ 
                            test_all_counts_signal$bar_means$signal != "visual" &
                              test_all_counts_signal$bar_means$sensillum_type == "all"], 
                            na.rm = T), digits = 3), "±", signif(sd(test_all_counts_signal$bar_means$means[
                              test_all_counts_signal$bar_means$signal != "visual" &
                                test_all_counts_signal$bar_means$sensillum_type == "all"], 
                              na.rm = T), digits = 1)), # counts
                          "n/a"  # fraction
  ),
  "sample_size" = c("n/a", # width
                "n/a", # length
                test_all_density_signal$all_phylo$n,  # density
                test_all_mean_dists_signal$all_phylo$n, # mean_dists
                test_all_counts_signal$all_phylo$n,  # counts
                "n/a"  # fraction
  ), 
  "f_value" = c("n/a", # width
                "n/a", # length
                signif(test_all_density_signal$olfactory_phylo$F, digits = 3),  # density
                signif(test_all_mean_dists_signal$olfactory_phylo$F, digits = 3), # mean_dists
                signif(test_all_counts_signal$olfactory_phylo$F, digits = 3),  # counts
                "n/a"  # fraction
  ),
  "p_value" = c("n/a", # width
                "n/a", # length
                signif(test_all_density_signal$mechano_phylo$Pf, digits = 2),  # density
                signif(test_all_mean_dists_signal$mechano_phylo$Pf, digits = 2), # mean_dists
                signif(test_all_counts_signal$mechano_phylo$Pf, digits = 2),  # counts
                "n/a"  # fraction
  )
)

rownames(tab_signal_all_df) <- tab_row_names
colnames(tab_signal_all_df) <- tab_col_signal_names

tab_signal_mechano_df <- data.frame(
  "visual_mean_value" = c(paste(signif(mean(test_widths_signal$bar_means$mean[ 
                        test_widths_signal$bar_means$signal == "visual" &
                        test_widths_signal$bar_means$sensillum_type == "mechano"], 
                        na.rm = T), digits = 3), "±", signif(sd(test_widths_signal$bar_means$mean[
                          test_widths_signal$bar_means$signal == "visual" &
                            test_widths_signal$bar_means$sensillum_type == "mechano"], 
                          na.rm = T), digits = 1)), # width
                        paste(signif(mean(test_lengths_signal$bar_means$mean[ 
                          test_lengths_signal$bar_means$signal == "visual" &
                            test_lengths_signal$bar_means$sensillum_type == "mechano"], 
                          na.rm = T), digits = 3), "±", signif(sd(test_lengths_signal$bar_means$mean[
                            test_lengths_signal$bar_means$signal == "visual" &
                              test_lengths_signal$bar_means$sensillum_type == "mechano"], 
                            na.rm = T), digits = 1)), # length
                          paste(signif(mean(test_all_density_signal$bar_means$means[ 
                            test_all_density_signal$bar_means$signal == "visual" &
                              test_all_density_signal$bar_means$sensillum_type == "mechano"], 
                            na.rm = T), digits = 3), "±", signif(sd(test_all_density_signal$bar_means$means[
                              test_all_density_signal$bar_means$signal == "visual" &
                                test_all_density_signal$bar_means$sensillum_type == "mechano"], 
                              na.rm = T), digits = 1)), # density
                          paste(signif(mean(test_all_mean_dists_signal$bar_means$means[ 
                            test_all_mean_dists_signal$bar_means$signal == "visual" &
                              test_all_mean_dists_signal$bar_means$sensillum_type == "mechano"], 
                            na.rm = T), digits = 3), "±", signif(sd(test_all_mean_dists_signal$bar_means$means[
                              test_all_mean_dists_signal$bar_means$signal == "visual" &
                                test_all_mean_dists_signal$bar_means$sensillum_type == "mechano"], 
                              na.rm = T), digits = 1)), # mean_dists
                          paste(signif(mean(test_all_counts_signal$bar_means$means[ 
                            test_all_counts_signal$bar_means$signal == "visual" &
                              test_all_counts_signal$bar_means$sensillum_type == "mechano"], 
                            na.rm = T), digits = 3), "±", signif(sd(test_all_counts_signal$bar_means$means[
                              test_all_counts_signal$bar_means$signal == "visual" &
                                test_all_counts_signal$bar_means$sensillum_type == "mechano"], 
                              na.rm = T), digits = 1)), # counts
                        paste(signif(mean(test_all_fractions_signal$bar_means$means[ 
                          test_all_fractions_signal$bar_means$signal == "visual" &
                            test_all_fractions_signal$bar_means$sensillum_type == "mechano"], 
                          na.rm = T), digits = 3), "±", signif(sd(test_all_fractions_signal$bar_means$means[
                            test_all_fractions_signal$bar_means$signal == "visual" &
                              test_all_fractions_signal$bar_means$sensillum_type == "mechano"], 
                            na.rm = T), digits = 1))  # fraction
  ),
  "chemical_mean_value" = c(paste(signif(mean(test_widths_signal$bar_means$mean[ 
                          test_widths_signal$bar_means$signal != "visual" &
                            test_widths_signal$bar_means$sensillum_type == "mechano"], 
                          na.rm = T), digits = 3), "±", signif(sd(test_widths_signal$bar_means$mean[
                            test_widths_signal$bar_means$signal != "visual" &
                              test_widths_signal$bar_means$sensillum_type == "mechano"], 
                            na.rm = T), digits = 1)), # width
                          paste(signif(mean(test_lengths_signal$bar_means$mean[ 
                            test_lengths_signal$bar_means$signal != "visual" &
                              test_lengths_signal$bar_means$sensillum_type == "mechano"], 
                            na.rm = T), digits = 3), "±", signif(sd(test_lengths_signal$bar_means$mean[
                              test_lengths_signal$bar_means$signal != "visual" &
                                test_lengths_signal$bar_means$sensillum_type == "mechano"], 
                              na.rm = T), digits = 1)), # length, # length
                            paste(signif(mean(test_all_density_signal$bar_means$means[ 
                              test_all_density_signal$bar_means$signal != "visual" &
                                test_all_density_signal$bar_means$sensillum_type == "mechano"], 
                              na.rm = T), digits = 3), "±", signif(sd(test_all_density_signal$bar_means$means[
                                test_all_density_signal$bar_means$signal != "visual" &
                                  test_all_density_signal$bar_means$sensillum_type == "mechano"], 
                                na.rm = T), digits = 1)), # density
                            paste(signif(mean(test_all_mean_dists_signal$bar_means$means[ 
                              test_all_mean_dists_signal$bar_means$signal != "visual" &
                                test_all_mean_dists_signal$bar_means$sensillum_type == "mechano"], 
                              na.rm = T), digits = 3), "±", signif(sd(test_all_mean_dists_signal$bar_means$means[
                                test_all_mean_dists_signal$bar_means$signal != "visual" &
                                  test_all_mean_dists_signal$bar_means$sensillum_type == "mechano"], 
                                na.rm = T), digits = 1)), # mean_dists
                            paste(signif(mean(test_all_counts_signal$bar_means$means[ 
                              test_all_counts_signal$bar_means$signal != "visual" &
                                test_all_counts_signal$bar_means$sensillum_type == "mechano"], 
                              na.rm = T), digits = 3), "±", signif(sd(test_all_counts_signal$bar_means$means[
                                test_all_counts_signal$bar_means$signal != "visual" &
                                  test_all_counts_signal$bar_means$sensillum_type == "mechano"], 
                                na.rm = T), digits = 1)), # counts
                          paste(signif(mean(test_all_fractions_signal$bar_means$means[ 
                            test_all_fractions_signal$bar_means$signal != "visual" &
                              test_all_fractions_signal$bar_means$sensillum_type == "mechano"], 
                            na.rm = T), digits = 3), "±", signif(sd(test_all_fractions_signal$bar_means$means[
                              test_all_fractions_signal$bar_means$signal != "visual" &
                                test_all_fractions_signal$bar_means$sensillum_type == "mechano"], 
                              na.rm = T), digits = 1))  # fraction
  ),
  "sample_size" = c(test_widths_signal$mechano_phylo$n, # width
                    test_lengths_signal$mechano_phylo$n, # length
                    test_all_density_signal$mechano_phylo$n,  # density
                    test_all_mean_dists_signal$mechano_phylo$n, # mean_dists
                    test_all_counts_signal$mechano_phylo$n,  # counts
                    test_all_fractions_signal$mechano_phylo$n  # fraction
  ), 
  "f_value" = c(signif(test_widths_signal$mechano_phylo$F, digits = 3), # width
                signif(test_lengths_signal$mechano_phylo$F, digits = 3), # length
                signif(test_all_density_signal$mechano_phylo$F, digits = 3),  # density
                signif(test_all_mean_dists_signal$mechano_phylo$F, digits = 3), # mean_dists
                signif(test_all_counts_signal$mechano_phylo$F, digits = 3),  # counts
                signif(test_all_fractions_signal$mechano_phylo$F, digits = 3)  # fraction
  ),
  "p_value" = c(signif(test_widths_signal$mechano_phylo$Pf, digits = 3), # width
                signif(test_lengths_signal$mechano_phylo$Pf, digits = 3), # length
                signif(test_all_density_signal$mechano_phylo$Pf, digits = 2),  # density
                signif(test_all_mean_dists_signal$mechano_phylo$Pf, digits = 2), # mean_dists
                signif(test_all_counts_signal$mechano_phylo$Pf, digits = 2),  # counts
                signif(test_all_fractions_signal$mechano_phylo$Pf, digits = 2)  # fraction
  )
)
rownames(tab_signal_mechano_df) <- tab_row_names
colnames(tab_signal_mechano_df) <- tab_col_signal_names

tab_signal_olfactory_df <- data.frame(
  "visual_mean_value" = c(paste(signif(mean(test_widths_signal$bar_means$mean[ 
    test_widths_signal$bar_means$signal == "visual" &
      test_widths_signal$bar_means$sensillum_type == "olfactory"], 
    na.rm = T), digits = 3), "±", signif(sd(test_widths_signal$bar_means$mean[
      test_widths_signal$bar_means$signal == "visual" &
        test_widths_signal$bar_means$sensillum_type == "olfactory"], 
      na.rm = T), digits = 1)), # width
    paste(signif(mean(test_lengths_signal$bar_means$mean[ 
      test_lengths_signal$bar_means$signal == "visual" &
        test_lengths_signal$bar_means$sensillum_type == "olfactory"], 
      na.rm = T), digits = 3), "±", signif(sd(test_lengths_signal$bar_means$mean[
        test_lengths_signal$bar_means$signal == "visual" &
          test_lengths_signal$bar_means$sensillum_type == "olfactory"], 
        na.rm = T), digits = 1)), # length
    paste(signif(mean(test_all_density_signal$bar_means$means[ 
      test_all_density_signal$bar_means$signal == "visual" &
        test_all_density_signal$bar_means$sensillum_type == "olfactory"], 
      na.rm = T), digits = 3), "±", signif(sd(test_all_density_signal$bar_means$means[
        test_all_density_signal$bar_means$signal == "visual" &
          test_all_density_signal$bar_means$sensillum_type == "olfactory"], 
        na.rm = T), digits = 1)), # density
    paste(signif(mean(test_all_mean_dists_signal$bar_means$means[ 
      test_all_mean_dists_signal$bar_means$signal == "visual" &
        test_all_mean_dists_signal$bar_means$sensillum_type == "olfactory"], 
      na.rm = T), digits = 3), "±", signif(sd(test_all_mean_dists_signal$bar_means$means[
        test_all_mean_dists_signal$bar_means$signal == "visual" &
          test_all_mean_dists_signal$bar_means$sensillum_type == "olfactory"], 
        na.rm = T), digits = 1)), # mean_dists
    paste(signif(mean(test_all_counts_signal$bar_means$means[ 
      test_all_counts_signal$bar_means$signal == "visual" &
        test_all_counts_signal$bar_means$sensillum_type == "olfactory"], 
      na.rm = T), digits = 3), "±", signif(sd(test_all_counts_signal$bar_means$means[
        test_all_counts_signal$bar_means$signal == "visual" &
          test_all_counts_signal$bar_means$sensillum_type == "olfactory"], 
        na.rm = T), digits = 1)), # counts
    paste(signif(mean(test_all_fractions_signal$bar_means$means[ 
      test_all_fractions_signal$bar_means$signal == "visual" &
        test_all_fractions_signal$bar_means$sensillum_type == "olfactory"], 
      na.rm = T), digits = 3), "±", signif(sd(test_all_fractions_signal$bar_means$means[
        test_all_fractions_signal$bar_means$signal == "visual" &
          test_all_fractions_signal$bar_means$sensillum_type == "olfactory"], 
        na.rm = T), digits = 1))  # fraction
  ),
  "chemical_mean_value" = c(paste(signif(mean(test_widths_signal$bar_means$mean[ 
    test_widths_signal$bar_means$signal != "visual" &
      test_widths_signal$bar_means$sensillum_type == "olfactory"], 
    na.rm = T), digits = 3), "±", signif(sd(test_widths_signal$bar_means$mean[
      test_widths_signal$bar_means$signal != "visual" &
        test_widths_signal$bar_means$sensillum_type == "olfactory"], 
      na.rm = T), digits = 1)), # width
    paste(signif(mean(test_lengths_signal$bar_means$mean[ 
      test_lengths_signal$bar_means$signal != "visual" &
        test_lengths_signal$bar_means$sensillum_type == "olfactory"], 
      na.rm = T), digits = 3), "±", signif(sd(test_lengths_signal$bar_means$mean[
        test_lengths_signal$bar_means$signal != "visual" &
          test_lengths_signal$bar_means$sensillum_type == "olfactory"], 
        na.rm = T), digits = 1)), # length, # length
    paste(signif(mean(test_all_density_signal$bar_means$means[ 
      test_all_density_signal$bar_means$signal != "visual" &
        test_all_density_signal$bar_means$sensillum_type == "olfactory"], 
      na.rm = T), digits = 3), "±", signif(sd(test_all_density_signal$bar_means$means[
        test_all_density_signal$bar_means$signal != "visual" &
          test_all_density_signal$bar_means$sensillum_type == "olfactory"], 
        na.rm = T), digits = 1)), # density
    paste(signif(mean(test_all_mean_dists_signal$bar_means$means[ 
      test_all_mean_dists_signal$bar_means$signal != "visual" &
        test_all_mean_dists_signal$bar_means$sensillum_type == "olfactory"], 
      na.rm = T), digits = 3), "±", signif(sd(test_all_mean_dists_signal$bar_means$means[
        test_all_mean_dists_signal$bar_means$signal != "visual" &
          test_all_mean_dists_signal$bar_means$sensillum_type == "olfactory"], 
        na.rm = T), digits = 1)), # mean_dists
    paste(signif(mean(test_all_counts_signal$bar_means$means[ 
      test_all_counts_signal$bar_means$signal != "visual" &
        test_all_counts_signal$bar_means$sensillum_type == "olfactory"], 
      na.rm = T), digits = 3), "±", signif(sd(test_all_counts_signal$bar_means$means[
        test_all_counts_signal$bar_means$signal != "visual" &
          test_all_counts_signal$bar_means$sensillum_type == "olfactory"], 
        na.rm = T), digits = 1)), # counts
    paste(signif(mean(test_all_fractions_signal$bar_means$means[ 
      test_all_fractions_signal$bar_means$signal != "visual" &
        test_all_fractions_signal$bar_means$sensillum_type == "olfactory"], 
      na.rm = T), digits = 3), "±", signif(sd(test_all_fractions_signal$bar_means$means[
        test_all_fractions_signal$bar_means$signal != "visual" &
          test_all_fractions_signal$bar_means$sensillum_type == "olfactory"], 
        na.rm = T), digits = 1))  # fraction
  ),
  "sample_size" = c(test_widths_signal$olfactory_phylo$n, # width
                    test_lengths_signal$olfactory_phylo$n, # length
                    test_all_density_signal$olfactory_phylo$n,  # density
                    test_all_mean_dists_signal$olfactory_phylo$n, # mean_dists
                    test_all_counts_signal$olfactory_phylo$n,  # counts
                    test_all_fractions_signal$olfactory_phylo$n  # fraction
  ), 
  "f_value" = c(signif(test_widths_signal$olfactory_phylo$F, digits = 3), # width
                signif(test_lengths_signal$olfactory_phylo$F, digits = 3), # length
                signif(test_all_density_signal$olfactory_phylo$F, digits = 3),  # density
                signif(test_all_mean_dists_signal$olfactory_phylo$F, digits = 3), # mean_dists
                signif(test_all_counts_signal$olfactory_phylo$F, digits = 3),  # counts
                signif(test_all_fractions_signal$olfactory_phylo$F, digits = 3)  # fraction
  ),
  "p_value" = c(signif(test_widths_signal$olfactory_phylo$Pf, digits = 3), # width
                signif(test_lengths_signal$olfactory_phylo$Pf, digits = 3), # length
                signif(test_all_density_signal$olfactory_phylo$Pf, digits = 2),  # density
                signif(test_all_mean_dists_signal$olfactory_phylo$Pf, digits = 2), # mean_dists
                signif(test_all_counts_signal$olfactory_phylo$Pf, digits = 2),  # counts
                signif(test_all_fractions_signal$olfactory_phylo$Pf, digits = 2)  # fraction
  )
)
rownames(tab_signal_olfactory_df) <- tab_row_names
colnames(tab_signal_olfactory_df) <- tab_col_signal_names

#### Sex testing ####

tab_sex_all_df <- data.frame(
  "female_mean_value" = c("n/a", # width
                          "n/a", # length
                          paste(signif(mean(test_all_density_sex$bar_means$means[ 
                            test_all_density_sex$bar_means$sex == "female" &
                              test_all_density_sex$bar_means$sensillum_type == "all"], 
                            na.rm = T), digits = 3), "±", signif(sd(test_all_density_sex$bar_means$means[
                              test_all_density_sex$bar_means$sex == "female" &
                                test_all_density_sex$bar_means$sensillum_type == "all"], 
                              na.rm = T), digits = 1)), # density
                          paste(signif(mean(test_all_mean_dists_sex$bar_means$means[ 
                            test_all_mean_dists_sex$bar_means$sex == "female" &
                              test_all_mean_dists_sex$bar_means$sensillum_type == "all"], 
                            na.rm = T), digits = 3), "±", signif(sd(test_all_mean_dists_sex$bar_means$means[
                              test_all_mean_dists_sex$bar_means$sex == "female" &
                                test_all_mean_dists_sex$bar_means$sensillum_type == "all"], 
                              na.rm = T), digits = 1)), # mean_dists
                          paste(signif(mean(test_all_counts_sex$bar_means$means[ 
                            test_all_counts_sex$bar_means$sex == "female" &
                              test_all_counts_sex$bar_means$sensillum_type == "all"], 
                            na.rm = T), digits = 3), "±", signif(sd(test_all_counts_sex$bar_means$means[
                              test_all_counts_sex$bar_means$sex == "female" &
                                test_all_counts_sex$bar_means$sensillum_type == "all"], 
                              na.rm = T), digits = 1)), # counts
                          "n/a"  # fraction
  ),
  "male_mean_value" = c("n/a", # width
                            "n/a", # length
                            paste(signif(mean(test_all_density_sex$bar_means$means[ 
                              test_all_density_sex$bar_means$sex == "male" &
                                test_all_density_sex$bar_means$sensillum_type == "all"], 
                              na.rm = T), digits = 3), "±", signif(sd(test_all_density_sex$bar_means$means[
                                test_all_density_sex$bar_means$sex == "male" &
                                  test_all_density_sex$bar_means$sensillum_type == "all"], 
                                na.rm = T), digits = 1)), # density
                            paste(signif(mean(test_all_mean_dists_sex$bar_means$means[ 
                              test_all_mean_dists_sex$bar_means$sex == "male" &
                                test_all_mean_dists_sex$bar_means$sensillum_type == "all"], 
                              na.rm = T), digits = 3), "±", signif(sd(test_all_mean_dists_sex$bar_means$means[
                                test_all_mean_dists_sex$bar_means$sex == "male" &
                                  test_all_mean_dists_sex$bar_means$sensillum_type == "all"], 
                                na.rm = T), digits = 1)), # mean_dists
                            paste(signif(mean(test_all_counts_sex$bar_means$means[ 
                              test_all_counts_sex$bar_means$sex == "male" &
                                test_all_counts_sex$bar_means$sensillum_type == "all"], 
                              na.rm = T), digits = 3), "±", signif(sd(test_all_counts_sex$bar_means$means[
                                test_all_counts_sex$bar_means$sex == "male" &
                                  test_all_counts_sex$bar_means$sensillum_type == "all"], 
                                na.rm = T), digits = 1)), # counts
                            "n/a"  # fraction
  ),
  "sample_size" = c("n/a", # width
                    "n/a", # length
                    test_all_density_sex$all_phylo$n,  # density
                    test_all_mean_dists_sex$all_phylo$n, # mean_dists
                    test_all_counts_sex$all_phylo$n,  # counts
                    "n/a"  # fraction
  ), 
  "f_value" = c("n/a", # width
                "n/a", # length
                signif(test_all_density_sex$olfactory_phylo$F, digits = 3),  # density
                signif(test_all_mean_dists_sex$olfactory_phylo$F, digits = 3), # mean_dists
                signif(test_all_counts_sex$olfactory_phylo$F, digits = 3),  # counts
                "n/a"  # fraction
  ),
  "p_value" = c("n/a", # width
                "n/a", # length
                signif(test_all_density_sex$mechano_phylo$Pf, digits = 2),  # density
                signif(test_all_mean_dists_sex$mechano_phylo$Pf, digits = 2), # mean_dists
                signif(test_all_counts_sex$mechano_phylo$Pf, digits = 2),  # counts
                "n/a"  # fraction
  )
)

rownames(tab_sex_all_df) <- tab_row_names
colnames(tab_sex_all_df) <- tab_col_sex_names

tab_sex_mechano_df <- data.frame(
  "female_mean_value" = c(paste(signif(mean(test_widths_sex$bar_means$mean[ 
    test_widths_sex$bar_means$sex == "female" &
      test_widths_sex$bar_means$sensillum_type == "mechano"], 
    na.rm = T), digits = 3), "±", signif(sd(test_widths_sex$bar_means$mean[
      test_widths_sex$bar_means$sex == "female" &
        test_widths_sex$bar_means$sensillum_type == "mechano"], 
      na.rm = T), digits = 1)), # width
    paste(signif(mean(test_lengths_sex$bar_means$mean[ 
      test_lengths_sex$bar_means$sex == "female" &
        test_lengths_sex$bar_means$sensillum_type == "mechano"], 
      na.rm = T), digits = 3), "±", signif(sd(test_lengths_sex$bar_means$mean[
        test_lengths_sex$bar_means$sex == "female" &
          test_lengths_sex$bar_means$sensillum_type == "mechano"], 
        na.rm = T), digits = 1)), # length
    paste(signif(mean(test_all_density_sex$bar_means$means[ 
      test_all_density_sex$bar_means$sex == "female" &
        test_all_density_sex$bar_means$sensillum_type == "mechano"], 
      na.rm = T), digits = 3), "±", signif(sd(test_all_density_sex$bar_means$means[
        test_all_density_sex$bar_means$sex == "female" &
          test_all_density_sex$bar_means$sensillum_type == "mechano"], 
        na.rm = T), digits = 1)), # density
    paste(signif(mean(test_all_mean_dists_sex$bar_means$means[ 
      test_all_mean_dists_sex$bar_means$sex == "female" &
        test_all_mean_dists_sex$bar_means$sensillum_type == "mechano"], 
      na.rm = T), digits = 3), "±", signif(sd(test_all_mean_dists_sex$bar_means$means[
        test_all_mean_dists_sex$bar_means$sex == "female" &
          test_all_mean_dists_sex$bar_means$sensillum_type == "mechano"], 
        na.rm = T), digits = 1)), # mean_dists
    paste(signif(mean(test_all_counts_sex$bar_means$means[ 
      test_all_counts_sex$bar_means$sex == "female" &
        test_all_counts_sex$bar_means$sensillum_type == "mechano"], 
      na.rm = T), digits = 3), "±", signif(sd(test_all_counts_sex$bar_means$means[
        test_all_counts_sex$bar_means$sex == "female" &
          test_all_counts_sex$bar_means$sensillum_type == "mechano"], 
        na.rm = T), digits = 1)), # counts
    paste(signif(mean(test_all_fractions_sex$bar_means$means[ 
      test_all_fractions_sex$bar_means$sex == "female" &
        test_all_fractions_sex$bar_means$sensillum_type == "mechano"], 
      na.rm = T), digits = 3), "±", signif(sd(test_all_fractions_sex$bar_means$means[
        test_all_fractions_sex$bar_means$sex == "female" &
          test_all_fractions_sex$bar_means$sensillum_type == "mechano"], 
        na.rm = T), digits = 1))  # fraction
  ),
  "male_mean_value" = c(paste(signif(mean(test_widths_sex$bar_means$mean[ 
    test_widths_sex$bar_means$sex != "visual" &
      test_widths_sex$bar_means$sensillum_type == "mechano"], 
    na.rm = T), digits = 3), "±", signif(sd(test_widths_sex$bar_means$mean[
      test_widths_sex$bar_means$sex == "male" &
        test_widths_sex$bar_means$sensillum_type == "mechano"], 
      na.rm = T), digits = 1)), # width
    paste(signif(mean(test_lengths_sex$bar_means$mean[ 
      test_lengths_sex$bar_means$sex == "male" &
        test_lengths_sex$bar_means$sensillum_type == "mechano"], 
      na.rm = T), digits = 3), "±", signif(sd(test_lengths_sex$bar_means$mean[
        test_lengths_sex$bar_means$sex == "male" &
          test_lengths_sex$bar_means$sensillum_type == "mechano"], 
        na.rm = T), digits = 1)), # length, # length
    paste(signif(mean(test_all_density_sex$bar_means$means[ 
      test_all_density_sex$bar_means$sex == "male" &
        test_all_density_sex$bar_means$sensillum_type == "mechano"], 
      na.rm = T), digits = 3), "±", signif(sd(test_all_density_sex$bar_means$means[
        test_all_density_sex$bar_means$sex == "male" &
          test_all_density_sex$bar_means$sensillum_type == "mechano"], 
        na.rm = T), digits = 1)), # density
    paste(signif(mean(test_all_mean_dists_sex$bar_means$means[ 
      test_all_mean_dists_sex$bar_means$sex == "male" &
        test_all_mean_dists_sex$bar_means$sensillum_type == "mechano"], 
      na.rm = T), digits = 3), "±", signif(sd(test_all_mean_dists_sex$bar_means$means[
        test_all_mean_dists_sex$bar_means$sex == "male" &
          test_all_mean_dists_sex$bar_means$sensillum_type == "mechano"], 
        na.rm = T), digits = 1)), # mean_dists
    paste(signif(mean(test_all_counts_sex$bar_means$means[ 
      test_all_counts_sex$bar_means$sex == "male" &
        test_all_counts_sex$bar_means$sensillum_type == "mechano"], 
      na.rm = T), digits = 3), "±", signif(sd(test_all_counts_sex$bar_means$means[
        test_all_counts_sex$bar_means$sex == "male" &
          test_all_counts_sex$bar_means$sensillum_type == "mechano"], 
        na.rm = T), digits = 1)), # counts
    paste(signif(mean(test_all_fractions_sex$bar_means$means[ 
      test_all_fractions_sex$bar_means$sex == "male" &
        test_all_fractions_sex$bar_means$sensillum_type == "mechano"], 
      na.rm = T), digits = 3), "±", signif(sd(test_all_fractions_sex$bar_means$means[
        test_all_fractions_sex$bar_means$sex == "male" &
          test_all_fractions_sex$bar_means$sensillum_type == "mechano"], 
        na.rm = T), digits = 1))  # fraction
  ),
  "sample_size" = c(test_widths_sex$mechano_phylo$n, # width
                    test_lengths_sex$mechano_phylo$n, # length
                    test_all_density_sex$mechano_phylo$n,  # density
                    test_all_mean_dists_sex$mechano_phylo$n, # mean_dists
                    test_all_counts_sex$mechano_phylo$n,  # counts
                    test_all_fractions_sex$mechano_phylo$n  # fraction
  ), 
  "f_value" = c(signif(test_widths_sex$mechano_phylo$F, digits = 3), # width
                signif(test_lengths_sex$mechano_phylo$F, digits = 3), # length
                signif(test_all_density_sex$mechano_phylo$F, digits = 3),  # density
                signif(test_all_mean_dists_sex$mechano_phylo$F, digits = 3), # mean_dists
                signif(test_all_counts_sex$mechano_phylo$F, digits = 3),  # counts
                signif(test_all_fractions_sex$mechano_phylo$F, digits = 3)  # fraction
  ),
  "p_value" = c(signif(test_widths_sex$mechano_phylo$Pf, digits = 3), # width
                signif(test_lengths_sex$mechano_phylo$Pf, digits = 3), # length
                signif(test_all_density_sex$mechano_phylo$Pf, digits = 2),  # density
                signif(test_all_mean_dists_sex$mechano_phylo$Pf, digits = 2), # mean_dists
                signif(test_all_counts_sex$mechano_phylo$Pf, digits = 2),  # counts
                signif(test_all_fractions_sex$mechano_phylo$Pf, digits = 2)  # fraction
  )
)
rownames(tab_sex_mechano_df) <- tab_row_names
colnames(tab_sex_mechano_df) <- tab_col_sex_names


### !!!!! need to add means to the following table
tab_sex_olfactory_df <- data.frame(
  "female_mean_value" = c(paste(signif(mean(test_widths_sex$bar_means$mean[ 
    test_widths_sex$bar_means$sex == "female" &
      test_widths_sex$bar_means$sensillum_type == "olfactory"], 
    na.rm = T), digits = 3), "±", signif(sd(test_widths_sex$bar_means$mean[
      test_widths_sex$bar_means$sex == "female" &
        test_widths_sex$bar_means$sensillum_type == "olfactory"], 
      na.rm = T), digits = 1)), # width
    paste(signif(mean(test_lengths_sex$bar_means$mean[ 
      test_lengths_sex$bar_means$sex == "female" &
        test_lengths_sex$bar_means$sensillum_type == "olfactory"], 
      na.rm = T), digits = 3), "±", signif(sd(test_lengths_sex$bar_means$mean[
        test_lengths_sex$bar_means$sex == "female" &
          test_lengths_sex$bar_means$sensillum_type == "olfactory"], 
        na.rm = T), digits = 1)), # length
    paste(signif(mean(test_all_density_sex$bar_means$means[ 
      test_all_density_sex$bar_means$sex == "female" &
        test_all_density_sex$bar_means$sensillum_type == "olfactory"], 
      na.rm = T), digits = 3), "±", signif(sd(test_all_density_sex$bar_means$means[
        test_all_density_sex$bar_means$sex == "female" &
          test_all_density_sex$bar_means$sensillum_type == "olfactory"], 
        na.rm = T), digits = 1)), # density
    paste(signif(mean(test_all_mean_dists_sex$bar_means$means[ 
      test_all_mean_dists_sex$bar_means$sex == "female" &
        test_all_mean_dists_sex$bar_means$sensillum_type == "olfactory"], 
      na.rm = T), digits = 3), "±", signif(sd(test_all_mean_dists_sex$bar_means$means[
        test_all_mean_dists_sex$bar_means$sex == "female" &
          test_all_mean_dists_sex$bar_means$sensillum_type == "olfactory"], 
        na.rm = T), digits = 1)), # mean_dists
    paste(signif(mean(test_all_counts_sex$bar_means$means[ 
      test_all_counts_sex$bar_means$sex == "female" &
        test_all_counts_sex$bar_means$sensillum_type == "olfactory"], 
      na.rm = T), digits = 3), "±", signif(sd(test_all_counts_sex$bar_means$means[
        test_all_counts_sex$bar_means$sex == "female" &
          test_all_counts_sex$bar_means$sensillum_type == "olfactory"], 
        na.rm = T), digits = 1)), # counts
    paste(signif(mean(test_all_fractions_sex$bar_means$means[ 
      test_all_fractions_sex$bar_means$sex == "female" &
        test_all_fractions_sex$bar_means$sensillum_type == "olfactory"], 
      na.rm = T), digits = 3), "±", signif(sd(test_all_fractions_sex$bar_means$means[
        test_all_fractions_sex$bar_means$sex == "female" &
          test_all_fractions_sex$bar_means$sensillum_type == "olfactory"], 
        na.rm = T), digits = 1))  # fraction
  ),
  "male_mean_value" = c(paste(signif(mean(test_widths_sex$bar_means$mean[ 
    test_widths_sex$bar_means$sex != "visual" &
      test_widths_sex$bar_means$sensillum_type == "olfactory"], 
    na.rm = T), digits = 3), "±", signif(sd(test_widths_sex$bar_means$mean[
      test_widths_sex$bar_means$sex == "male" &
        test_widths_sex$bar_means$sensillum_type == "olfactory"], 
      na.rm = T), digits = 1)), # width
    paste(signif(mean(test_lengths_sex$bar_means$mean[ 
      test_lengths_sex$bar_means$sex == "male" &
        test_lengths_sex$bar_means$sensillum_type == "olfactory"], 
      na.rm = T), digits = 3), "±", signif(sd(test_lengths_sex$bar_means$mean[
        test_lengths_sex$bar_means$sex == "male" &
          test_lengths_sex$bar_means$sensillum_type == "olfactory"], 
        na.rm = T), digits = 1)), # length, # length
    paste(signif(mean(test_all_density_sex$bar_means$means[ 
      test_all_density_sex$bar_means$sex == "male" &
        test_all_density_sex$bar_means$sensillum_type == "olfactory"], 
      na.rm = T), digits = 3), "±", signif(sd(test_all_density_sex$bar_means$means[
        test_all_density_sex$bar_means$sex == "male" &
          test_all_density_sex$bar_means$sensillum_type == "olfactory"], 
        na.rm = T), digits = 1)), # density
    paste(signif(mean(test_all_mean_dists_sex$bar_means$means[ 
      test_all_mean_dists_sex$bar_means$sex == "male" &
        test_all_mean_dists_sex$bar_means$sensillum_type == "olfactory"], 
      na.rm = T), digits = 3), "±", signif(sd(test_all_mean_dists_sex$bar_means$means[
        test_all_mean_dists_sex$bar_means$sex == "male" &
          test_all_mean_dists_sex$bar_means$sensillum_type == "olfactory"], 
        na.rm = T), digits = 1)), # mean_dists
    paste(signif(mean(test_all_counts_sex$bar_means$means[ 
      test_all_counts_sex$bar_means$sex == "male" &
        test_all_counts_sex$bar_means$sensillum_type == "olfactory"], 
      na.rm = T), digits = 3), "±", signif(sd(test_all_counts_sex$bar_means$means[
        test_all_counts_sex$bar_means$sex == "male" &
          test_all_counts_sex$bar_means$sensillum_type == "olfactory"], 
        na.rm = T), digits = 1)), # counts
    paste(signif(mean(test_all_fractions_sex$bar_means$means[ 
      test_all_fractions_sex$bar_means$sex == "male" &
        test_all_fractions_sex$bar_means$sensillum_type == "olfactory"], 
      na.rm = T), digits = 3), "±", signif(sd(test_all_fractions_sex$bar_means$means[
        test_all_fractions_sex$bar_means$sex == "male" &
          test_all_fractions_sex$bar_means$sensillum_type == "olfactory"], 
        na.rm = T), digits = 1))  # fraction
  ),
  "sample_size" = c(test_widths_sex$olfactory_phylo$n, # width
                    test_lengths_sex$olfactory_phylo$n, # length
                    test_all_density_sex$olfactory_phylo$n,  # density
                    test_all_mean_dists_sex$olfactory_phylo$n, # mean_dists
                    test_all_counts_sex$olfactory_phylo$n,  # counts
                    test_all_fractions_sex$olfactory_phylo$n  # fraction
  ), 
  "f_value" = c(signif(test_widths_sex$olfactory_phylo$F, digits = 3), # width
                signif(test_lengths_sex$olfactory_phylo$F, digits = 3), # length
                signif(test_all_density_sex$olfactory_phylo$F, digits = 3),  # density
                signif(test_all_mean_dists_sex$olfactory_phylo$F, digits = 3), # mean_dists
                signif(test_all_counts_sex$olfactory_phylo$F, digits = 3),  # counts
                signif(test_all_fractions_sex$olfactory_phylo$F, digits = 3)  # fraction
  ),
  "p_value" = c(signif(test_widths_sex$olfactory_phylo$Pf, digits = 3), # width
                signif(test_lengths_sex$olfactory_phylo$Pf, digits = 3), # length
                signif(test_all_density_sex$olfactory_phylo$Pf, digits = 2),  # density
                signif(test_all_mean_dists_sex$olfactory_phylo$Pf, digits = 2), # mean_dists
                signif(test_all_counts_sex$olfactory_phylo$Pf, digits = 2),  # counts
                signif(test_all_fractions_sex$olfactory_phylo$Pf, digits = 2)  # fraction
  )
)
rownames(tab_sex_olfactory_df) <- tab_row_names
rownames(tab_sex_olfactory_df) <- tab_row_names
