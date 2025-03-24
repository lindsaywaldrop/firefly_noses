# Firefly morphology project

Morphometric measurements and calculations of firefly antennae.


## Navigating this repository

This repository is organized in the following way: 

 - The main directory holds: 
    - `README` file with instructions on reproducing analyses and manuscript files
    - `firefly-morphology.Rproj`, and R-project file which can be opened in RStudio
    - `LICENSE` file which outlines acceptable use and reuse of code and other intellectual property held herein

- `data` folder which contains all raw measurement data necessary for reproducing figures and analyses
- `doc` folder which contains all documentation in the project, including the main manuscript and supplemental information RMD files. 
- `results` folder which contains all results produced by the analyses, including figures. 
- `src` folder which contains all code, including an installation script for packages required. 

## Required software and packages

Figures, tables, and statistical values are fully reproducible in the main manuscript as well as the supplemental information using RStudio version 2024.12.1+563 and R version 4.4.1 or higher from the R project file. 

The following R packages are available on CRAN and required to run this code: 

 - `ape`
 - `bookdown`
 - `cowplot`
 - `dplyr`
 - `ggplot2`
 - `ggrepel`
 - `janitor`
 - `knitr`
 - `patchwork`
 - `pracma`
 - `phytools`
 - `purr`
 - `RColorBrewer`
 - `readr`
 - `rmarkdown`
 - `RTriangle`
 - `stringr`
 - `tidyverse`
 - `tidyr`
 - `viridis`
 
With the exception of `ggtree` which must be installed via BiocManager: 

```
install.packages("BiocManager")
BiocManager::install("ggtree")
```

An install packages script is located in `src` to install all packages needed to run the code. Please note that it may require running this script more than once to complete the install. 


## Reproducing the analyses and manuscript files

In order to regenerate the content of both documents:  

 1. Clone the Github repository to your machine. 
 2. Open the Rproj file in the main directory in RStudio. 
 3. Navigate to `src/` folder and run `install_R_packages.R` to install the required packages for this analysis. Please note that it may require running this script more than once to complete the install. 
 4. Navigate to `doc/` from the main directory and open `Morphology_paper.Rmd`. 
 5. Adjust the settings under `Knit` to knit the project directory. 
 6. Knit the document OR select `Run all` to run the code in situ. 
 
Repeat this process for the supplemental information by opening `Morphology_paper_SI.Rmd` in step 3 and following steps 4-6. 

All code to reproduce the statistical tests are located in `src/` in the script `paper_run_stats.R`. This script is meant to be run as a part of the RMD files mentioned above and will not run independently. 

To reproduce the analysis independently from raw measurement data, in line 54 of `Morphology_paper.Rmd` change the line to: 

`rerun_analyses <- TRUE`

You may also regenerate and save each figure as an independent file by changing lines 52 and 56 to read `TRUE`.

