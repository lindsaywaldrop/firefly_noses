# Run Phylogenetically corrected statistics

test_all_mean_dists_sex <- phylo_test_hairs(bug_tree, all_xy, 
                                            "mean_dists", "sex",
                                            "all", T, 0.15)
test_all_mean_dists_signal <- phylo_test_hairs(bug_tree, all_xy, 
                                               "mean_dists", "signal",
                                               "all", T, 0.15)

#test_olf_mean_dists_sex <- phylo_test_hairs(bug_tree, all_xy, 
#                                        "mean_dists", "sex",
#                                        "olfactory", T)
#test_olf_mean_dists_signal <- phylo_test_hairs(bug_tree, all_xy, 
#                                        "mean_dists", "signal",
#                                        "olfactory", T)

test_all_density_sex <- phylo_test_hairs(bug_tree, all_xy, 
                                         "density", "sex",
                                         "all", T)
test_all_density_signal <- phylo_test_hairs(bug_tree, all_xy, 
                                            "density", "signal",
                                            "all", T)

test_olf_density_sex <- phylo_test_hairs(bug_tree, all_xy, 
                                         "density", "sex",
                                         "olfactory", T)
test_olf_density_signal <- phylo_test_hairs(bug_tree, all_xy, 
                                            "density", "signal",
                                            "olfactory", T)

test_all_counts_sex <- phylo_test_hairs(bug_tree, all_xy, 
                                        "est_seg_count", "sex",
                                        "all", T)
test_all_counts_signal <- phylo_test_hairs(bug_tree, all_xy, 
                                           "est_seg_count", "signal",
                                           "all", T)

test_olf_counts_signal <- phylo_test_hairs(bug_tree, all_xy, 
                                           "est_seg_count", "signal",
                                           "olfactory", T)
test_olf_counts_sex <- phylo_test_hairs(bug_tree, all_xy, 
                                        "est_seg_count", "sex",
                                        "olfactory", T)

test_all_fractions_sex <- phylo_test_hairs(bug_tree, all_xy, 
                                           "fraction", "sex",
                                           "all", T)
test_all_fractions_signal <- phylo_test_hairs(bug_tree, all_xy, 
                                              "fraction", "signal",
                                              "all", T)

test_olf_fractions_sex <- phylo_test_hairs(bug_tree, all_xy, 
                                        "fraction", "sex",
                                        "olfactory", T)
test_olf_fractions_signal <- phylo_test_hairs(bug_tree, all_xy, 
                                        "fraction", "signal",
                                        "olfactory", T)

# Uncomment the line below to run the analysis corrected for body length (not presented in the main manuscript)
#morpho_data$value <- morpho_data$value/morpho_data$body_length
test_widths_signal <- phy_test_morphometric(bug_tree, morpho_data, 
                                            "sensillum width med", "signal", 0.25)
test_widths_sex <- phy_test_morphometric(bug_tree, morpho_data, 
                                         "sensillum width med", "sex", 0.25)
test_lengths_signal <- phy_test_morphometric(bug_tree, morpho_data, 
                                             "sensillum length", "signal", 0.25)
test_lengths_sex <- phy_test_morphometric(bug_tree, morpho_data, 
                                          "sensillum length", "sex", 0.25)
