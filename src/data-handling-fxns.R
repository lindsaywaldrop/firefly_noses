#### Data handling functions ####

## Loads specimen data
load_specimen_data <- function(underscore = FALSE, include_unmeasured = TRUE){
  require(readr)
  require(tidyr)
  require(janitor)
  specimen_data <- read_csv("./data/Specimen_list.csv") %>%
    remove_empty(which = c("rows", "cols")) %>%
    clean_names() 
  specimen_data$requested <- factor(specimen_data$requested, levels = c("yes", "no"))
  specimen_data$in_custody <- factor(specimen_data$in_custody, levels = c("yes", "no", "lost", "not available"))
  specimen_data$scanned <- factor(specimen_data$scanned, levels = c("yes", "no"))
  specimen_data$measured <- factor(specimen_data$measured, levels = c("yes", "no"))
  specimen_data$antenna_complete <- factor(specimen_data$antenna_complete, levels = c("yes", "no"))
  if(underscore) specimen_data$species <- sub(" ", "_", specimen_data$species)
  if(include_unmeasured == FALSE) {
    specimen_data <- specimen_data[specimen_data$measured == "yes", ]
  } else{
    drop_these <- which(specimen_data$in_custody == "not available" | 
                          specimen_data$in_custody == "no" | 
                          specimen_data$in_custody == "lost")
    specimen_data <- specimen_data[-drop_these,]
  }
  specimen_data$specimen <- factor(specimen_data$specimen)
  specimen_data$species <- factor(specimen_data$species)
  specimen_data$signal <- factor(specimen_data$signal, ordered = T,
                                 levels = c("visual","chemical","both"))
  specimen_data$sex <- factor(specimen_data$sex, levels = c("male", "female"))
  specimen_data$type_of_antennae <- factor(specimen_data$type_of_antennae)
  return(specimen_data)
}

## Loads morphometric measured data 
load_morphometric_data <- function(specimen_data, analysis_date = NULL){
  require(tidyr)
  require(janitor)
  if(is.null(analysis_date)){
    files_list <- list.files(path = "./data/Morphometric_files", full.names = TRUE, 
                             pattern = "*.csv")
    combined_dat <- data.frame()
    for(i in files_list){
      dat <- read.csv(i, skip = 1)
      if(is.character(dat$Value)) warning(print(i))
      combined_dat <- rbind(combined_dat, dat)
    }
    morpho_data <- combined_dat %>%
      remove_empty(which = c("rows", "cols")) %>%
      clean_names()
    # cleaning out specimen data not in species_data list
    morpho_data$value <- as.numeric(morpho_data$value)
    morpho_data$specimen <- factor(morpho_data$specimen, levels = levels(specimen_data$specimen))
    morpho_data <- morpho_data[!is.na(morpho_data$specimen),]
    morpho_data$species <- factor(morpho_data$species, levels = levels(specimen_data$species))
    # creating and assigning signal data from species_data list
    morpho_data$signal <- rep(NA, nrow(morpho_data))
    morpho_data$sex <- rep(NA, nrow(morpho_data))
    for (u in 1:nrow(specimen_data)){
      for (k in 1:nrow(morpho_data)){
        if(morpho_data$species[k] == specimen_data$species[u]){
          morpho_data$signal[k] <- as.character(specimen_data$signal[u])
          morpho_data$sex[k] <- as.character(specimen_data$sex[u])
        }
      }
    }
    for (k in 1:nrow(morpho_data)){
      for (u in 1:nrow(specimen_data)){
        if(morpho_data$species[k] == specimen_data$species[u]){
          morpho_data$body_length[k] <- as.numeric(specimen_data$body_length_mm[u])
        }
      }
    }
  }else{
    morpho_data <- read.csv(paste0("./results/csv-files/morpho_data_", analysis_date, ".csv"))
    morpho_data$specimen <- factor(morpho_data$specimen)
    morpho_data$species <- factor(morpho_data$species)
  }
  
  if(nlevels(factor(morpho_data$signal)) != 3){
    warning("Warning: incorrect number of levels for signal. ")
  }
  morpho_data$signal <- factor(morpho_data$signal, ordered = T,
                               levels = c("visual","chemical","both"))
  if(nlevels(factor(morpho_data$sex)) != 2){
    warning("Warning: incorrect number of levels for sex. ")
  }
  morpho_data$sex <- factor(morpho_data$sex, levels = c("male", "female"))
  if(nlevels(factor(morpho_data$measurement_type)) != 7){ 
    warning("Warning: incorrect number of levels for measurement type. ")
  }
  morpho_data$measurement_type <- factor(morpho_data$measurement_type)
  if(nlevels(factor(morpho_data$location)) != 10){ 
    warning("Warning: incorrect number of levels for location. ")
  }
  morpho_data$location <- factor(morpho_data$location, ordered = T, 
                                 levels = paste0("F", 1:10))
  if(nlevels(factor(morpho_data$sensillum_type)) != 2){ 
    warning("Warning: incorrect number of levels for sensillum type. ")
  }
  morpho_data$sensillum_type <- factor(morpho_data$sensillum_type, ordered = T,
                                       levels = c("mechano", "olfactory", "hygro-thermo", "all"))
  if(nlevels(factor(morpho_data$olfactory_type)) != 2){ 
    warning("Warning: incorrect number of levels for sensillum type. ")
  }
  morpho_data$olfactory_type <- factor(morpho_data$olfactory_type, ordered = T,
                                       levels = c("basiconica", "trichoid", "capitular"))
  return(morpho_data)
}

load_xy_data <- function(analysis_date){
  all_xy <- read.csv(paste0("./results/csv-files/all_xy_", analysis_date,".csv"))
  all_xy$specimen <- factor(all_xy$specimen)
  all_xy$species <- factor(all_xy$species)
  all_xy$sex <- factor(all_xy$sex)
  all_xy$location <- factor(all_xy$location, ordered = T, 
                            levels = paste0("F", seq(1:10)))
  all_xy$signal <- factor(all_xy$signal, ordered = T, 
                          levels = c("visual", "chemical", "both"))
  all_xy$sensillum_type <- factor(all_xy$sensillum_type, ordered = T, 
                                  levels = c("mechano", "olfactory", "hygro-thermo", "all"))
  all_xy$olfactory_type <- ifelse(all_xy$sensillum_type == "olfactory" &          
                                    is.na(all_xy$olfactory_type),
                                  "all", as.character(all_xy$olfactory_type))
  #all_xy$olfactory_type <- ifelse(all_xy$olfactory_type == "basiconica" | 
  #                                         all_xy$olfactory_type == "capitular", 
  #                                       "other", all_xy$olfactory_type)
  all_xy$olfactory_type <- factor(all_xy$olfactory_type, ordered = T,
                                  levels = c("basiconica", "trichoid", "capitular", "all"))
  return(all_xy)
}

## Loads x,y coordinate data
load_coord_data <- function(specimen_data){
  require(tidyr)
  require(janitor)
  require(purrr)
  coord_data <- list.files(path = "./data/Coordinate_files", 
                                       full.names = TRUE, pattern = "*.csv") %>%
    map_df(~read.csv(.)) %>%
    remove_empty(which = c("rows", "cols")) %>%
    clean_names()
  # Keeping only data in specimen list
  coord_data$specimen <- factor(coord_data$specimen, levels = levels(specimen_data$specimen))
  coord_data <- coord_data[!is.na(coord_data$specimen),]
  coord_data$species <- factor(coord_data$species, levels = levels(specimen_data$species))
  coord_data$signal <- rep(NA, nrow(coord_data))
  coord_data$sex <- rep(NA, nrow(coord_data))
  coord_data$body_length <- rep(NA, nrow(coord_data))
  for (k in 1:nrow(coord_data)){
    for (u in 1:nrow(specimen_data)){
      if(coord_data$species[k] == specimen_data$species[u]){
        coord_data$signal[k] <- as.character(specimen_data$signal[u])
        coord_data$sex[k] <- as.character(specimen_data$sex[u])
        coord_data$body_length[k] <- as.character(specimen_data$body_length_mm[u])
      }
    }
  }
  
  if(nlevels(factor(coord_data$signal)) != 3){ 
    warning("Warning: incorrect number of levels for signal. ")
  }
  coord_data$signal <- factor(coord_data$signal, ordered = T,
                              levels = c("visual", "chemical", "both"))
  if(nlevels(factor(coord_data$sex)) != 2){ 
    warning("Warning: incorrect number of levels for sex. ")
  }
  coord_data$sex <- factor(coord_data$sex, levels = c("male", "female"))
  if(nlevels(factor(coord_data$sensillum_type)) != 3){ 
    warning("Warning: incorrect number of levels for sensillum type. ")
  }
  coord_data$sensillum_type <- factor(coord_data$sensillum_type, ordered = T,
                                      levels = c("mechano", "olfactory", "hygro-thermo"))
  if(nlevels(factor(coord_data$olfactory_type)) != 3){ 
    warning("Warning: incorrect number of levels for olfactory type. ")
  }
  coord_data$olfactory_type <- factor(coord_data$olfactory_type, ordered = T,
                                      levels = c("basiconica", "trichoid", "capitular"))
  if(nlevels(factor(coord_data$location)) != 10){ 
    warning("Warning: incorrect number of levels for location. ")
  }
  coord_data$location <- factor(coord_data$location,  ordered = T, 
                                levels = paste0("F", 1:10))
  return(coord_data)
}

## Calculates distances between hairs using nearest neighbor 
calc_nn_dist <- function(dat, specimen, loc, s_type = NULL, o_type = NULL){
  require(RTriangle)
  if(!is.character(specimen)) stop("species must be a character")
  if(!is.character(loc)) stop("loc must be a character")
  if(is.null(s_type)){
    coords <- as.matrix(dat[dat$location == loc & dat$specimen == specimen, 7:8])
    coords <- na.omit(coords)
  } else {
    if(!is.null(o_type)){
      coords <- as.matrix(dat[dat$location == loc & dat$specimen == specimen & 
                                dat$sensillum_type == s_type &
                                dat$olfactory_type == o_type, 7:8])
    }else {
      coords <- as.matrix(dat[dat$location == loc & dat$specimen == specimen & 
                                dat$sensillum_type == s_type, 7:8])
    }
    
  }
  if(dim(na.omit(coords))[1] < 3) {
    dist <- NA
    x <- NA
  } else {
    if(any(duplicated(coords))){
      coords <- coords[!duplicated(coords),]
    }
    F.pslg <- pslg(na.omit(coords))
    x <- triangulate(F.pslg)
    dist <- sqrt((x$P[x$E[, 2], 1] - x$P[x$E[, 1], 1])^2 + 
                 (x$P[x$E[, 2], 2] - x$P[x$E[, 1], 2])^2
    )
  }
  return(list(dist, x))
}

## This function finds the densities and fractions of each type of sensillum. 

find_hair_densities <- function(dat, specimen, loc) {
  require(pracma)
  mean_dists <- rep(NA, 7)
  density <- rep(NA, 7)
  frac <- rep(NA, 7)
  counts <- rep(NA, 7)
  if(nrow(dat[dat$specimen == specimen & dat$location == loc,]) == 0 | 
     nrow(dat[dat$specimen == specimen & dat$location == loc,]) == 1){
    
  } else {
    all_dots <- calc_nn_dist(dat, specimen, loc)
    total_area <- polyarea(all_dots[[2]]$P[all_dots[[2]]$S[, 1], 1],
                           all_dots[[2]]$P[all_dots[[2]]$S[, 1], 2])
    total_counts <- length(all_dots[[2]]$P[, 1]) 
    counts[length(density)] <- total_counts
    density[length(density)] <- length(all_dots[[2]]$P[, 1])/total_area
    mean_dists[length(mean_dists)] <- mean(all_dots[[1]], na.rm = TRUE)
    for (j in 1:(length(density) - 1)){
      if(j == 1) { # set mechano
        sensillum_type <- levels(dat$sensillum_type)[1]
        olfactory_type <- NULL
      }else if(j == 2){ # set hygro-thermo
        sensillum_type <- levels(dat$sensillum_type)[3]
        olfactory_type <- NULL
      }else if(j == 3) {
        sensillum_type <- levels(dat$sensillum_type)[2]
        olfactory_type <- "trichoid"
      }else if(j == 4) {
        sensillum_type <- levels(dat$sensillum_type)[2]
        olfactory_type <- "basiconica"
      }else if(j == 5) {
        sensillum_type <- levels(dat$sensillum_type)[2]
        olfactory_type <- "capitular"
      }else{
        sensillum_type <- levels(dat$sensillum_type)[2]
        olfactory_type <- NULL
      }
      dots <- calc_nn_dist(dat, specimen, loc, sensillum_type, olfactory_type)
      if(is.na(dots[[1]][1])){
        mean_dists[j] <- NA   
        density[j] <- NA
        frac[j] <- NA
        counts[j] <- NA
      } else{
        mean_dists[j] <- mean(dots[[1]], na.rm = TRUE)
        density[j] <- length(dots[[2]]$P[, 1])/total_area
        
        counts[j] <- length(dots[[2]]$P[, 1])
        if(sensillum_type=="olfactory" & is.null(olfactory_type)){
          frac[j] <- length(dots[[2]]$P[, 1])/length(all_dots[[2]]$P[, 1])
        } else{
          frac[j] <- length(dots[[2]]$P[, 1])/length(all_dots[[2]]$P[, 1])
        }
      }
      
    }
    counts[is.na(counts)] <- 0
    counts[6] <- sum(counts[3:5], na.rm = TRUE)
    frac[length(frac)] <- sum(frac[c(1,2,6)], na.rm = TRUE)
    frac[3] <- counts[3]/counts[6]
    frac[4] <- counts[4]/counts[6]
    frac[5] <- counts[5]/counts[6]
    if(isTRUE(all.equal(frac[length(frac)], 1.0)) & 
       isTRUE(all.equal(frac[1], 1.0)) & 
       is.na(frac[2])){
      frac[2] <- 0.0
    }
  }
  return(data.frame("specimen" = factor(rep(specimen, length(density)), 
                                       levels = levels(dat$specimen)), 
                    "species" = factor(rep(dat$species[dat$specimen == specimen][1])),
                    "signal" = factor(rep(dat$signal[dat$specimen == specimen][1])),
                    "sex" = factor(rep(dat$sex[dat$specimen == specimen][1])),
                    "body_length" = rep(dat$body_length[dat$specimen == specimen][1]),
                    "location" = factor(rep(loc, length(density)), 
                                        levels = levels(dat$location)),
                    "sensillum_type" = factor(c("mechano", "hygro-thermo",
                                                "olfactory", "olfactory", "olfactory", 
                                                "olfactory", "all"),
                                              levels = 
                                                c("mechano", "hygro-thermo",
                                                  "olfactory", "all")),
                    "olfactory_type" = factor(c(NA, NA, "trichoid", "basiconica", 
                                                "capitular", NA, NA), 
                                              levels = c("trichoid", "basiconica", "capitular")),
                    "mean_dists" = mean_dists,
                    "density" = density, 
                    "fraction" = frac,
                    "counts" = counts))
}

## Use both functions to calculate across the whole data set:
calc_hair_densities <- function(dat){
  ns <- nlevels(dat$specimen)
  nf <- nlevels(dat$location)
  df <- data.frame()
  for (j in 1:ns){ # specimen
    for (k in 1:nf){ # location
      #print(paste(j, k))
      df2 <- find_hair_densities(dat, levels(dat$specimen)[j], 
                            levels(dat$location)[k])
      df <- rbind(df, df2) 
      rm(df2)
    }
  }
  return(df)
}

## Calculate the estimated number of sensilla per segement: 
estimate_counts <- function(specimen_data, morpho_data, all_xy){
  require(tidyr)
  require(janitor)
  require(dplyr)
  morph_segments <- morpho_data %>%
    select(- c("sensillum_type", "olfactory_type", "label_on_image", 
               "filename", "comments")) %>%
    filter(measurement_type == "segment width" | measurement_type == "segment length") %>%
    pivot_wider(id_cols = c(species, specimen, signal, sex, body_length, location), 
                names_from = measurement_type) %>%
    clean_names()
  
  morph_segments$est_area <- 2 * morph_segments$segment_length*
    morph_segments$segment_width
  
  all_xy <- left_join(all_xy, morph_segments)
  all_xy$est_seg_count <- floor(all_xy$density * all_xy$est_area)
  all_xy$est_seg_count <- ifelse(all_xy$counts == 0, 0, all_xy$est_seg_count)
  return(all_xy)
}

## Reconcile species lists
reconcile_species <- function(bug_tree, notyet){
  require(stringr)
  specieswehave <- load_specimen_data(T, notyet)
  specieswehave$io_species <- rep(NA, length = dim(specieswehave)[1])
  for(i in 1:length(specieswehave$species)){
    a <- which(specieswehave$species[i] == bug_tree$tip.label)
    if(length(a) == 0){
      specieswehave$io_species[i] <- 0
    } else {
      specieswehave$io_species[i] <- a
    }
    
  }
  
  list_species <- specieswehave[, c("species", "io_species")]
  list_species2 <- list_species[which(list_species$io_species == 0),1]
  list_species2 <- unique(list_species2)
  list_species2 <- list_species2[order(list_species2$species),]
  list_species2$genus <- pull_genus_names(list_species2$species)
  return(list_species2)
}

## Pulls out the genus name from each species
pull_genus_names <- function(species_list){
  require(stringr)
  list_split <- str_split(species_list, "_", n = 2)
  df_split <- data.frame(list_split)
  genus_list <- as.character(t(df_split[1,]))
  return(genus_list)
}

## Drops species on the Martin tree that we don't have in our analysis 
final_drop <- function(bug_tree, specieswehave=NULL){
  require(ape)
  if(is.null(specieswehave)) specieswehave <- load_specimen_data(T,F)
  # Final drop of all tips not in data set: 
  drop <- rep(NA, length(bug_tree$tip.label))
  for(i in 1:length(drop)){
    b <- which(specieswehave$species == bug_tree$tip.label[i])
    if(length(b) == 0){
      drop[i] <- 0
    } else {
      drop[i] <- 1
    }
  }
  bug_tree_final <- drop.tip(bug_tree, bug_tree$tip.label[drop == 0])
  return(bug_tree_final)
}


## Test hairs against with phylogenetically corrected anova
phylo_test_hairs <- function(bug_tree, all_xy, mean_type, test_var,
                             sensillum_type = "all",
                             plotit = T, pad = 0.15){
  require(ggplot2)
  require(phytools)
  require(dplyr)
  
  test_results <- list()
  
  #if(mean_type == "counts") all_xy <- all_xy[all_xy$location == "F4", ]
  
  if(sensillum_type == "all"){
    #not_all <- all_xy[which(all_xy$sensillum_type != "all"), ]
    not_all <- all_xy
    not_all <- not_all[-which(not_all$olfactory_type == "trichoid" | 
                                not_all$olfactory_type == "basiconica" | 
                                not_all$olfactory_type == "capitular"), ]
    not_all$sensillum_type <- as.character(not_all$sensillum_type)
    column_type <- "sensillum_type"
    
    bar_means <- not_all %>% 
      group_by(species, sensillum_type, .data[[test_var]]) %>%
      #na.omit() %>%
      summarise_at(vars(mean_type), list(means = mean, sds = sd), na.rm = T)
    
    bar_means$sds[is.na(bar_means$sds)] <- 0
    bar_means$high_sd <- bar_means$means + bar_means$sds
    bar_means$low_sd <- bar_means$means - bar_means$sds
    bar_means$low_sd[bar_means$low_sd <= 0] <- 0
    
    bar_means <- bar_means[-which(bar_means$sensillum_type == "hygro-thermo"),]
    bar_means$sensillum_type <- factor(bar_means$sensillum_type, 
                                       ordered = TRUE,  
                                       levels = c("mechano", "olfactory", "all"))
    
  }else if(sensillum_type == "olfactory"){
    not_all <- all_xy[which(all_xy$sensillum_type == "olfactory"), ]
    not_all <- not_all[not_all$olfactory_type != "all", ]
    not_all$sensillum_type <- as.character(not_all$sensillum_type)
    column_type <- "olfactory_type"
    
    bar_means <- not_all %>% 
      group_by(species, olfactory_type, .data[[test_var]]) %>%
      #na.omit() %>%
      summarise_at(vars(mean_type), list(means = mean, sds = sd), 
                   na.rm = T)
    bar_means$sds[is.na(bar_means$sds)] <- 0
    bar_means$high_sd <- bar_means$means + bar_means$sds
    bar_means$low_sd <- bar_means$means - bar_means$sds
    bar_means$low_sd[bar_means$low_sd <= 0] <- 0
    #bar_means <- bar_means[-which(bar_means$olfactory_type == "capitular"),]
    bar_means$olfactory_type <- factor(bar_means$olfactory_type)
  }
  
  
  bar_means$species <- factor(bar_means$species, 
                              ordered = T, levels = bug_tree$tip.label)
  test_results[["bar_means"]] <- bar_means
  
  
  if(sensillum_type == "all"){
    for(i in 1:nlevels(bar_means$sensillum_type)){
      test_set <- bar_means[bar_means$sensillum_type == 
                              levels(bar_means$sensillum_type)[i],]
      if(test_var == "signal") test_set$signal[test_set$signal == "both"] <- "chemical"
      if(mean_type == "fraction") {
        test_set$means <- log(test_set$means / (1 - test_set$means))
        test_set$means[is.infinite(test_set$means)] <- NA
      }
      test_set <- droplevels(test_set)
      test_set <- test_set[!is.na(test_set$means), ]
      if(nrow(test_set) < 2) next
      if(test_var == "signal") {
        test_vars <- test_set$signal
      }else if(test_var == "sex"){
        test_set <- test_set[!is.na(test_set$sex), ]
        test_vars <- test_set$sex
      }
      names(test_vars) <- test_set$species
      test_means <- test_set$means
      names(test_means) <- test_set$species
      if(length(bug_tree$tip.label) != length(test_means)){
        bug_tree2 <- final_drop(bug_tree, specieswehave = test_set)
      }else{
        bug_tree2 <- bug_tree
      }
      print(paste("Sample size n = ", length(test_means)))
      print(paste("phylANOVA results for", test_var, "and",
                  levels(bar_means$sensillum_type)[i]))
      test_results[[paste(levels(bar_means$sensillum_type)[i], 
                          "phylo", sep = "_")]] <-  
        phylANOVA(bug_tree2, test_vars, test_means)
      test_results[[paste(levels(bar_means$sensillum_type)[i], 
                         "phylo", sep = "_")]]$n <- 
        length(test_means)
      print(test_results[[paste(levels(bar_means$sensillum_type)[i], 
                                "phylo", sep = "_")]])
      print(paste("Regular ANOVA results for", test_var, "and", 
                  levels(bar_means$sensillum_type)[i]))
      if(test_var == "signal") {
        test_results[[paste(levels(bar_means$sensillum_type)[i], 
                          "reg", sep = "_")]] <- 
        anova(lm(means ~ signal, data = test_set))
      test_results[[paste(levels(bar_means$sensillum_type)[i], 
                          "reg", sep = "_")]]$n <- 
        length(test_means)
      print(test_results[[paste(levels(bar_means$sensillum_type)[i], 
                                "reg", sep = "_")]])
      }else{
        test_results[[paste(levels(bar_means$sensillum_type)[i], 
                            "reg", sep = "_")]] <- 
          anova(lm(means ~ sex, data = test_set))
        test_results[[paste(levels(bar_means$sensillum_type)[i], 
                            "reg", sep = "_")]]$n <- 
          length(test_means)
        print(test_results[[paste(levels(bar_means$sensillum_type)[i], 
                                  "reg", sep = "_")]])
      }
      
      
    }
    
  } else if (sensillum_type=="olfactory"){
    for(i in 1:nlevels(bar_means$olfactory_type)){
      test_set <- bar_means[bar_means$olfactory_type ==
                              levels(bar_means$olfactory_type)[i],]
      if(test_var == "signal") test_set$signal[test_set$signal == "both"] <- "chemical"
      if(mean_type == "fraction") {
        test_set$means <- log(test_set$means / (1 - test_set$means))
        test_set$means[is.infinite(test_set$means)] <- NA
      }
      test_set <- droplevels(test_set)
      test_set <- test_set[!is.na(test_set$means),]
      if(nrow(test_set) < 2) next
      if(test_var == "signal") {
        test_vars <- test_set$signal
      }else if(test_var == "sex"){
        test_set <- test_set[!is.na(test_set$sex), ]
        test_vars <- test_set$sex
      }
      names(test_vars) <- test_set$species
      test_means <- test_set$means
      names(test_means) <- test_set$species
      if(length(bug_tree$tip.label) != length(test_means)){
        bug_tree2 <- final_drop(bug_tree, specieswehave = test_set)
      }else{
        bug_tree2 <- bug_tree
      }
      print(paste("Sample size n = ", length(test_means)))
      print(paste("phylANOVA results for", test_var, "and",  
                  levels(bar_means$olfactory_type)[i]))
      test_results[[paste(levels(bar_means$olfactory_type)[i], 
                          "phylo", sep = "_")]] <- 
        phylANOVA(bug_tree2, test_vars, test_means)
      test_results[[paste(levels(bar_means$olfactory_type)[i], 
                            "phylo", sep = "_")]]$n <- 
        length(test_means)
      print(test_results[[paste(levels(bar_means$olfactory_type)[i], 
                                "phylo", sep = "_")]])
      if(test_var == "signal") {
        print(paste("Regular ANOVA results for", test_var, "and",  
                  levels(bar_means$olfactory_type)[i]))
      test_results[[paste(levels(bar_means$olfactory_type)[i], 
                          "reg", sep = "_")]] <- 
        anova(lm(means ~ signal, data = test_set))
      test_results[[paste(levels(bar_means$olfactory_type)[i], 
                          "reg", sep = "_")]]$n <- 
        length(test_means)
      print(test_results[[paste(levels(bar_means$olfactory_type)[i], 
                                "reg", sep = "_")]])
      }else{
        print(paste("Regular ANOVA results for", test_var, "and",  
                    levels(bar_means$olfactory_type)[i]))
        test_results[[paste(levels(bar_means$olfactory_type)[i], 
                            "reg", sep = "_")]] <- 
          anova(lm(means ~ sex, data = test_set))
        test_results[[paste(levels(bar_means$olfactory_type)[i], 
                            "reg", sep = "_")]]$n <- 
          length(test_means)
        print(test_results[[paste(levels(bar_means$olfactory_type)[i], 
                                  "reg", sep = "_")]])
      }
    }
  }
  
  if(plotit){
    if(mean_type == "fraction"){
      sumsqrs <- function(sds) sqrt(sum(na.omit(sds)^2))
      if(sensillum_type == "all") bar_means <- bar_means[bar_means$sensillum_type != "all",]
      bar_means_sum <- bar_means %>% 
        group_by(.data[[column_type]], .data[[test_var]]) %>%
        summarise(means = sum(means, na.rm = T), sds = sumsqrs(sds))
      bar_means_sum <- na.omit(bar_means_sum)
      bar_means_sum$fractional_means <- NA
      bar_means_sum$fractional_sds <- NA  
      if(sensillum_type == "all"){
        for(var_level in levels(bar_means_sum[[test_var]])){
          for(hair_level in levels(bar_means_sum$sensillum_type)[1:2]){
            bar_means_sum$fractional_means[bar_means_sum$sensillum_type ==  hair_level &
                                             bar_means_sum[[test_var]] == var_level] <- 
              bar_means_sum$means[bar_means_sum$sensillum_type == hair_level &
                                    bar_means_sum[[test_var]] == var_level] /
              sum(bar_means_sum$means[bar_means_sum[[test_var]] == var_level])
            bar_means_sum$fractional_sds[bar_means_sum$sensillum_type ==  hair_level &
                                           bar_means_sum[[test_var]] == var_level] <- 
              bar_means_sum$sds[bar_means_sum$sensillum_type == hair_level &
                                  bar_means_sum[[test_var]] == var_level] /
              sum(bar_means_sum$means[bar_means_sum[[test_var]] == var_level])
          }
        }
        bar_means_sum$frac_sds_low <- bar_means_sum$fractional_means - bar_means_sum$fractional_sds
        bar_means_sum$frac_sds_high <- bar_means_sum$fractional_means + bar_means_sum$fractional_sds
        bar_means_sum$frac_sds_low[bar_means_sum$frac_sds_low < 0] <- 0
        bar_means_sum$frac_sds_low[bar_means_sum$sensillum_type == "mechano"] <- NA
        bar_means_sum$frac_sds_high[bar_means_sum$sensillum_type == "mechano"] <- NA
        #bar_means_sum <- bar_means_sum[-which(bar_means_sum[[column_type]] == "all"), ]
      } else if (sensillum_type == "olfactory"){
          for(var_level in levels(bar_means_sum[[test_var]])){
            for(hair_level in levels(bar_means_sum$olfactory_type)){
              bar_means_sum$fractional_means[bar_means_sum$olfactory_type ==  hair_level &
                                               bar_means_sum[[test_var]] == var_level] <- 
              bar_means_sum$means[bar_means_sum$olfactory_type == hair_level &
                                    bar_means_sum[[test_var]] == var_level] /
                sum(bar_means_sum$means[bar_means_sum[[test_var]] == var_level])
              bar_means_sum$fractional_sds[bar_means_sum$olfactory_type ==  hair_level &
                                               bar_means_sum[[test_var]] == var_level] <- 
                bar_means_sum$sds[bar_means_sum$olfactory_type == hair_level &
                                      bar_means_sum[[test_var]] == var_level] /
                sum(bar_means_sum$means[bar_means_sum[[test_var]] == var_level])
            }
          }
        bar_means_sum$frac_sds_high <- bar_means_sum$fractional_means + bar_means_sum$fractional_sds
        bar_means_sum$frac_sds_low <- bar_means_sum$fractional_means - bar_means_sum$fractional_sds
        bar_means_sum$frac_sds_high[bar_means_sum$olfactory_type == "basiconica"] <- NA
        bar_means_sum$frac_sds_low[bar_means_sum$olfactory_type == "basiconica"] <- NA
      }
      
      p1 <- ggplot(bar_means_sum, aes(x = .data[[test_var]], y = fractional_means, 
                                         fill = .data[[column_type]])) +
        geom_bar(stat = "identity", col = "black") +
        geom_errorbar(mapping = aes(ymin = frac_sds_low, 
                                    ymax = frac_sds_high), color = "black", 
                      width = 0.5) +
        coord_cartesian() + 
        ylab("Fraction of sensilla")  + xlab(" ") +
        theme_minimal() + 
        theme(axis.text.x = element_text(angle = 25))
        #theme(axis.text.x = element_text(angle = 25), 
        #                        legend.direction = "horizontal", 
        #                        legend.justification = "center", 
        #                        legend.box.just = "bottom",
        #                        legend.position = "bottom") 
    }else{
      p1 <- ggplot(bar_means, aes(.data[[column_type]], means, 
                                         fill = .data[[test_var]])) +
        geom_boxplot(position = "dodge", col = "black", outlier.shape = NA) +
        stat_boxplot(position = "dodge", geom = "errorbar") +
        #geom_jitter(alpha = 0.4, shape = 21, position = position_jitterdodge()) +
        geom_pointrange(mapping = aes(ymin = low_sd, ymax = high_sd),
                        position = position_jitterdodge(jitter.width = 0.2), 
                        pch = 21, alpha = 0.6) +
        xlab("") +
        theme_minimal() + theme(axis.text.x = element_text(angle = 25)) 
    }
    if(mean_type=="mean_dists") {
      p1 <- p1 + ylab("Distances between\nsensilla (\u00B5m)")
    } else if(mean_type=="counts") {
      p1 <- p1 + ylab("Number of\nsensilla per segment") 
    } else if(mean_type=="density"){
      p1 <- p1 + ylab(bquote(Density~of~sensilla~(mm^2)))
    }
    # if(sensillum_type == "all" & mean_type != "fraction"){
    #   p1 <- p1 + coord_cartesian(ylim = c(layer_scales(p1)$y$range$range[1] - 
    #                                           diff(layer_scales(p1)$y$range$range)*pad,
    #                                         layer_scales(p1)$y$range$range[2])) +
    #     annotate("text", x = 1, y = layer_scales(p1)$y$range$range[1] - 
    #                   diff(layer_scales(p1)$y$range$range)*0.65*pad, 
    #                 label = paste0("n = ", test_results[["mechano_phylo"]]$n, 
    #                                "\np = ", test_results[["mechano_phylo"]]$Pf)) +
    #     annotate("text", x = 2, y = layer_scales(p1)$y$range$range[1] - 
    #                diff(layer_scales(p1)$y$range$range)*0.65*pad, 
    #              label = paste0("n = ", test_results[["olfactory_phylo"]]$n,
    #                             "\np = ", test_results[["olfactory_phylo"]]$Pf))  +
    #     annotate("text", x = 3, y = layer_scales(p1)$y$range$range[1] - 
    #                diff(layer_scales(p1)$y$range$range)*0.65*pad, 
    #              label = paste0("n = ", test_results[["all_phylo"]]$n, 
    #                             "\np = ", test_results[["all_phylo"]]$Pf)) 
    # }else if(sensillum_type == "olfactory" & mean_type != "fraction"){
    #   p1 <- p1 + coord_cartesian(ylim = c(layer_scales(p1)$y$range$range[1] - 
    #                                         diff(layer_scales(p1)$y$range$range)*pad,
    #                                       layer_scales(p1)$y$range$range[2])) +
    #     annotate("text", x = 1, y = layer_scales(p1)$y$range$range[1] - 
    #                diff(layer_scales(p1)$y$range$range)*0.65*pad, 
    #              label = paste0("n = ", test_results[["basiconica_phylo"]]$n,
    #                             "\np = ", test_results[["basiconica_phylo"]]$Pf)) +
    #     annotate("text", x = 2, y = layer_scales(p1)$y$range$range[1] - 
    #                diff(layer_scales(p1)$y$range$range)*0.65*pad, 
    #              label = paste0("n = ", test_results[["trichoid_phylo"]]$n,
    #                             "\np = ", test_results[["trichoid_phylo"]]$Pf))
    # }  
      
    test_results[["plot"]] <- p1
  }
  
  return(test_results)
}

phy_test_morphometric <- function(bug_tree, dat, meas_type, test_var,
                                  plotit = T, pad = 0.15){
  require(ggplot2)
  require(phytools)
  
  test_results<-list()
  
  dat_sub <- dat[which(dat$measurement_type == meas_type), ]
  
  if(test_var == "signal"){
    bar_means <- dat_sub %>% 
      group_by(species, sensillum_type, signal) %>%
      summarise(mean = mean(value, na.rm = T), sds = sd(value, na.rm = T))
  }else{
    bar_means <- dat_sub %>% 
      group_by(species, sensillum_type, sex) %>%
      summarise(mean = mean(value, na.rm = T), sds = sd(value, na.rm = T))
  }
  
  bar_means$sensillum_type <- factor(bar_means$sensillum_type, ordered = T, 
                      levels = c("mechano", "olfactory"))
  bar_means <- bar_means[!is.na(bar_means$sensillum_type),]
  
  test_results[["bar_means"]] <- bar_means

  for(i in 1:nlevels(bar_means$sensillum_type)){
    test_set <- bar_means[bar_means$sensillum_type == levels(bar_means$sensillum_type)[i], ]
    if(test_var == "signal") test_set$signal[test_set$signal == "both"] <- "chemical"
    test_set <- droplevels(test_set)
    test_set <- na.omit(test_set)
    if(nrow(test_set) < 2) next
    if(test_var == "signal") {
      test_vars <- test_set$signal
    }else if(test_var == "sex"){
      test_vars <- test_set$sex
    }
    names(test_vars) <- test_set$species
    test_means <- test_set$mean
    names(test_means) <- test_set$species
    if(length(bug_tree$tip.label) != length(test_means)){
      bug_tree2 <- final_drop(bug_tree, specieswehave = test_set)
    }else{
      bug_tree2 <- bug_tree
    }
    print(paste("Sample size n = ", length(test_means)))
    print(paste("phylANOVA results for", test_var, "and", 
                levels(bar_means$sensillum_type)[i]))
    test_results[[paste(levels(bar_means$sensillum_type)[i], 
                        "phylo", sep = "_")]] <-  
      phylANOVA(bug_tree2, test_vars, test_means)
    test_results[[paste(levels(bar_means$sensillum_type)[i], 
                        "phylo", sep = "_")]]$n <-
      length(test_means)
    print(test_results[[paste(levels(bar_means$sensillum_type)[i], 
                              "phylo", sep = "_")]])
    if(test_var == "signal"){
      print(paste("Regular ANOVA results for ", test_var, "and",
                  levels(bar_means$sensillum_type)[i]))
      test_results[[paste(levels(bar_means$sensillum_type)[i], 
                          "reg", sep = "_")]] <- 
        anova(lm(mean ~ signal, data = test_set))
      test_results[[paste(levels(bar_means$sensillum_type)[i], 
                          "reg", sep = "_")]]$n <- 
        length(test_means)
      print(test_results[[paste(levels(bar_means$sensillum_type)[i], 
                                "reg", sep = "_")]])
    }else{
      print(paste("Regular ANOVA results for ", test_var, "and",
                  levels(bar_means$sensillum_type)[i]))
      test_results[[paste(levels(bar_means$sensillum_type)[i], 
                          "reg", sep = "_")]] <- 
        anova(lm(mean ~ sex, data = test_set))
      test_results[[paste(levels(bar_means$sensillum_type)[i], 
                          "reg", sep = "_")]]$n <- 
        length(test_means)
      print(test_results[[paste(levels(bar_means$sensillum_type)[i], 
                                "reg", sep = "_")]])
    }
    
  }
  
  if(plotit){
    p1 <- ggplot(bar_means, aes(sensillum_type, mean, 
                                  fill = .data[[test_var]])) +
        geom_boxplot(position="dodge", col = "black", 
                     outlier.shape = NA) +
        stat_boxplot(position = "dodge", geom = "errorbar") +
        #geom_jitter(alpha = 0.4, shape = 21, position=position_jitterdodge()) +
        geom_pointrange(mapping = aes(ymin = mean - sds, ymax = mean + sds),
                    position = position_jitterdodge(jitter.width = 0.2), 
                    pch = 21, alpha = 0.6) +
        xlab("") +
        theme_minimal() + theme(axis.text.x = element_text(angle = 25)) 
    if(grepl("width", meas_type)) {
      p1 <- p1 + ylab("Mean sensillum\nwidth (\u00B5m)")
    } else if(grepl("sensillum length", meas_type)) {
      p1 <- p1 + ylab("Mean sensillum\nlength (\u00B5m)")
    } else if(grepl("antenna length", meas_type)){
      p1 <- p1 + ylab("Antenna length (\u00B5m)")
    } else if(grepl("segment length", meas_type)){
      p1 <- p1 + ylab("Antennomere length (\u00B5m)")
    }
    # p1 <- p1 + coord_cartesian(ylim = c(layer_scales(p1)$y$range$range[1] - 
    #                                       diff(layer_scales(p1)$y$range$range)*pad,
    #                                     layer_scales(p1)$y$range$range[2])) +
    #   annotate("text", x = 1, y = layer_scales(p1)$y$range$range[1]- 
    #              diff(layer_scales(p1)$y$range$range)*0.65*pad, 
    #            label = paste0("n = ", test_results[["mechano_phylo"]]$n,
    #                           "\np = ", test_results[["mechano_phylo"]]$Pf)) +
    #   annotate("text", x = 2, y = layer_scales(p1)$y$range$range[1] - 
    #              diff(layer_scales(p1)$y$range$range)*0.65*pad, 
    #            label = paste0("n = ", test_results[["olfactory_phylo"]]$n, 
    #                           "\np = ", test_results[["olfactory_phylo"]]$Pf))  
    
    test_results[["plot"]] <- p1
    #ggsave(paste0("./results/plots/", sensillum_type, "-", 
    #              mean_type, "-1.pdf"), plot = p1, 
    #       width = 4, height = 4)
  }
  
  return(test_results)
}

create_pie_charts <- function(xy_means, colname){
  require(ggplot2)
  require(cowplot)
  a <- list()
  for (k in 1:nlevels(xy_means$species)){
    xy_species1 <- xy_means[xy_means$species == levels(xy_means$species)[k], ]
    a[[k]] <- ggplot(xy_species1, aes(x = "", y = .data[[colname]], 
                                      fill = sensillum_type)) + 
      geom_bar(stat = "identity") + coord_polar("y", start = 0) + 
      scale_fill_viridis(discrete = TRUE, drop = FALSE) + 
      #ggtitle(levels(xy_means2$species)[k])+
      theme_void() + theme(legend.position = "none") 
    #rm(xy_species1)
  }
  
  p.legend <- ggplot(xy_species1, aes(x = "", y = .data[[colname]], 
                                      fill = sensillum_type)) + 
    geom_bar(stat = "identity") + coord_polar("y", start = 0) + 
    scale_fill_viridis(discrete = TRUE, drop = FALSE, name = "Sensillum type") + 
    #ggtitle(levels(xy_means2$species)[k])+
    theme_void() 
  a[[k+1]] <- get_legend(p.legend + theme(legend.direction = "vertical", 
                                                   legend.justification = "center" , 
                                                   legend.box.just = "right"))
  names(a) <- c(levels(xy_means$species), "legend")
  return(a)
}

