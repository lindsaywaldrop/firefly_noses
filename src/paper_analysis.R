# Morphology paper analysis script

source("./src/data-handling-fxns.R")
dir.create("./results/csv-files/", showWarnings = FALSE)

# Reads in all specimen list data
specimen_data <- load_specimen_data(include_unmeasured = F)

# Reads in all morphometric data from ./data/Morphometric_files:
morpho_data <- load_morphometric_data(specimen_data)

# Saves morphometric data as a summary csv in ./results/csv-files folder:
write.csv(morpho_data, file = paste0("./results/csv-files/morpho_data_", Sys.Date(), ".csv"),
          row.names = F)

# Load coordinate data files from ./data/Coordinate_files/:
coord_data <- load_coord_data(specimen_data)

#Calculate the density of sensillum, count, and fraction of each type of sensillum:
all_xy <- calc_hair_densities(coord_data)  

# Saves coordinate data calculations as a summary csv in ./results/csv-files folder:
write.csv(all_xy, file = paste0("./results/csv-files/all_xy_", Sys.Date(), ".csv"),
          row.names = F)
