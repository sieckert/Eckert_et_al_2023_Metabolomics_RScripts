# PACKAGES ====================================================================
library(here)
library(zen4R)

# INFO ========================================================================
# This R script is part of the statistical analysis published in
# Eckert et al. (2023) in the scientific journal Metabolomics
# To load the corresponding raw data use the following R code:
zen4R::download_zenodo(doi = "10.5281/zenodo.7687634",
                       path = here(".")) # add target directory to save files
# unzip downloaded files
utils::unzip("Eckert_et_al_2023_RawData_Zenodo.zip")

# SESSION INFO ================================================================
sessionInfo()