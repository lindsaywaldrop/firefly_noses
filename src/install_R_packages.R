#### Installs required R packages ####

packages <- c("ape", "bookdown", "cowplot", "dplyr", "ggplot2", "ggrepel",
              "janitor", "knitr", "patchwork", "pracma", "phytools", "purr",
              "RColorBrewer", "readr", "rmarkdown", "RTriangle", "stringr", 
              "tidyverse", "tidyr", "viridis")

package.check <- lapply(
  packages,
  FUN <- function(x) {
    if (!require(x, character.only = TRUE)) {
      install.packages(x, dependencies = TRUE, repos = "http://cran.us.r-project.org")
      library(x, character.only = TRUE)
    }
  }
)

# For the ggtree package:
install.packages("BiocManager")
BiocManager::install("ggtree")
