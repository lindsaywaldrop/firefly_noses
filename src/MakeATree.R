# Script for producing tree for analysis
library(tidyverse)
library(phytools)

source("./src/data-handling-fxns.R")
# Would you like to plot the trees?
plot_it <- TRUE
# Include species we don't have yet? 
notyet <- FALSE
el <- 0.01 #default edge length for all the inserted species

#Read the tree file in: 
bug_tree <- read.tree("./data/trees/Lampyridae_AHE.tree") # File from Martin et al. 2019 

# Clean the tree of some unknown species and replicate species
drop_set <- c(4, 13, 26, 27, 30, 40, 48, 59, 60, 61, 69, 70, 71, 83, 
              84, 85, 86, 88)
bug_tree <- drop.tip(bug_tree, bug_tree$tip.label[drop_set]) 
bug_tree$tip.label <- sub("_1", "", bug_tree$tip.label)
bug_tree$tip.label <- sub("[.]", "", bug_tree$tip.label)

# Plot the tree to have a look, original tree from reference removing replicates
if(plot_it) plotTree(bug_tree, node.numbers = T)

# Save the tree tip labels
write.csv(bug_tree$tip.label,"./data/trees/bug_tree_tips.csv")

# Read in the second set to drop
drop_set2 <- read_csv("./data/trees/bug-include.csv")
colnames(drop_set2)<-c("number", "species", "keep")
bug_tree1 <- drop.tip(bug_tree, 
                      bug_tree$tip.label[drop_set2$number[drop_set2$keep == 0]]) 

# Cleaned tree dropping genera that aren't in our data set
if(plot_it) plotTree(bug_tree1, node.numbers = T)


list_species2 <- reconcile_species(bug_tree, T)

# Adding species that share a genus name with others in the original tree, 
# with notes on how tips are added.
genera_to_add <- unique(list_species2$genus)


# Adding Bicellonycha wickershamorum as sister to Bicellonycha sp.
Bicellonycha_species <- grep(genera_to_add[1], list_species2$species)
sp_pos <- which(bug_tree1$tip.label == "Bicellonycha_sp")
bug_tree2 <- bind.tip(bug_tree1, as.character(list_species2$species[Bicellonycha_species[1]]), 
                      edge.length = el, where = sp_pos, position = 0.01)



# Adding Ellychnia species in a polytomy

Ellychnia_species <- grep(genera_to_add[2], list_species2$species)
sp_pos <- which(bug_tree2$tip.label == "Ellychnia_sp")
bug_tree2 <- bind.tip(bug_tree2, as.character(list_species2$species[Ellychnia_species[1]]), 
                      edge.length = el, where = sp_pos, position = 0.01)
bug_tree2 <- bind.tip(bug_tree2, as.character(list_species2$species[Ellychnia_species[2]]),
                      edge.length = el, where = sp_pos, position = 0.01)


# Adding Lucidota_sp. and Lucidota punctata as sister to Lucidota atria
Lucidota_species <- grep(genera_to_add[3], list_species2$species)
sp_pos <- which(bug_tree2$tip.label == "Lucidota_atra")
bug_tree2 <- bind.tip(bug_tree2, as.character(list_species2$species[Lucidota_species[1]]),
                      edge.length = el, where = sp_pos, position = 0.01)
bug_tree2 <- bind.tip(bug_tree2, as.character(list_species2$species[Lucidota_species[2]]),
                      edge.length = el, where = sp_pos, position = 0.01)

# Adding Lucio_sp. as sister to Lucio blattinum
Lucio_species <- grep(genera_to_add[4], list_species2$species)
sp_pos <- which(bug_tree2$tip.label == "Lucio_blattinum")
bug_tree2 <- bind.tip(bug_tree2, as.character(list_species2$species[Lucio_species[1]]), 
                      edge.length = el, where = sp_pos, position = 0.01)


# Adding Microphotus species as a polytomy with Microphotus sp.
Microphotus_species <- grep(genera_to_add[5], list_species2$species)
sp_pos <- which(bug_tree2$tip.label == "Microphotus_sp")
bug_tree2 <- bind.tip(bug_tree2, as.character(list_species2$species[Microphotus_species[1]]), 
                      edge.length = el, where = sp_pos, position = 0.01)

# Adding Photinus species
Photinus_species <- grep(genera_to_add[6], list_species2$species)

# Adding Photinus marginellus as sister to Photinus floridanus (Stanger-Hall & Lloyd 2015 Evolution)
sp_pos <- which(bug_tree2$tip.label == "Photinus_floridanus")
bug_tree2 <- bind.tip(bug_tree2, as.character(list_species2$species[Photinus_species[3]]), 
                      edge.length = el, where = sp_pos, position = 0.01)

# Adding Photinus indictus as sister to Photinus macdermotti & Photinus consanguineus (Stanger-Hall & Lloyd 2015 Evolution)
sp_pos <- which(bug_tree2$tip.label=="Photinus_macdermotti")
bug_tree2 <- bind.tip(bug_tree2, as.character(list_species2$species[Photinus_species[2]]), 
                      edge.length = 2*el, where = sp_pos, position = 2*0.01)

# Adding Photinus consanguineus as sister to Photinus macdermotti (Stanger-Hall & Lloyd 2015 Evolution)
sp_pos <- which(bug_tree2$tip.label=="Photinus_macdermotti")
bug_tree2 <- bind.tip(bug_tree2, as.character(list_species2$species[Photinus_species[1]]), 
                      edge.length = el, where = sp_pos, position = 0.01)

# Adding Photinus scintillans as sister to P. carolinus - P. ardens clade (Stanger-Hall & Lloyd 2015 Evolution)
sp_pos <- 75
bug_tree2 <- bind.tip(bug_tree2, as.character(list_species2$species[Photinus_species[4]]), 
                      edge.length = el, where = sp_pos, position = 0.01)

# Adding Photuris species as a polytomy with Photuris frontalis and other Photuris
Photuris_species <- grep(genera_to_add[7], list_species2$species)
sp_pos <- which(bug_tree2$tip.label == "Photuris_frontalis")
bug_tree2 <- bind.tip(bug_tree2, as.character(list_species2$species[Photuris_species[1]]), 
                      edge.length = el, where = sp_pos, position = 0.01)
bug_tree2 <- bind.tip(bug_tree2, as.character(list_species2$species[Photuris_species[2]]), 
                      edge.length = el, where = sp_pos, position = 0.01)

# Adding Pyractomena species as a polytomy with Pyractomena borealis
Pyractomena_species <- grep(genera_to_add[8], list_species2$species)
sp_pos <- which(bug_tree2$tip.label == "Pyractomena_borealis")
bug_tree2 <- bind.tip(bug_tree2, as.character(list_species2$species[Pyractomena_species[1]]), 
                      edge.length = el, where = sp_pos, position = 0.01)
bug_tree2 <- bind.tip(bug_tree2, as.character(list_species2$species[Pyractomena_species[2]]), 
                      edge.length = el, where = sp_pos, position = 0.01)

# Adding Vesta basalis as a polytomy to Vesta impressicollis and Vesta sp/Vesta saturnalis
Vesta_species <- grep(genera_to_add[9], list_species2$species)
sp_pos <- which(bug_tree2$tip.label == "Vesta_sp")
bug_tree2 <- bind.tip(bug_tree2, as.character(list_species2$species[Vesta_species[1]]), 
                      edge.length = el, where = sp_pos, position = 0.01)


# Tree with tips added:
if(plot_it) plotTree(bug_tree2, node.numbers = F)


bug_tree3 <- final_drop(bug_tree2)

if(plot_it) plotTree(bug_tree3, node.numbers = F)

# Save tree for later
if(notyet){
  write.tree(bug_tree3, "./data/trees/Treefor29species.tree")
} else {
  write.tree(bug_tree3, "./data/trees/Treefor26species.tree")
}



