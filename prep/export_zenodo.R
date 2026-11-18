################################################################################
#' Public copy of the analysis datasets for Zenodo
#'
#' Needs outputs/data/ (R/01_data.R). Run from the
#' repository root. Writes outputs/zenodo/: the analysis datasets read by
#' R/02_fit_main.R onwards, and README.txt; and outputs/ssa_incidence_analysis_data.zip
#' of those files, for upload.
################################################################################

library(tidyverse)

source("R/functions.R")

path_zenodo <- "outputs/zenodo"
dir.create(path_zenodo, showWarnings = FALSE, recursive = TRUE)

file.copy(file.path(path_derived, "analysis_data.csv"), path_zenodo, overwrite = TRUE)
file.copy(file.path(path_derived, "flowchart_counts.csv"), path_zenodo, overwrite = TRUE)

writeLines(c(
  "Analysis datasets for 'HIV incidence declines and changing age patterns of new infections",
  "in sub-Saharan Africa', as built by R/01_data.R in the paper repository.",
  "",
  "analysis_data.csv has one row per observation. extrapolated_value is the modelled",
  "incidence used as the offset in the main analysis; extrapolated_value_no_phia,",
  "extrapolated_value_spectrum and extrapolated_value_spectrum_national are the offsets",
  "for the sensitivity analyses (Table S6).",
  "",
  "To use them, place the files in outputs/data/ of the repository and run the analysis",
  "from R/02_fit_main.R onwards, skipping R/01_data.R.",
  "flowchart_counts.csv gives the study counts for Figure 1."
), file.path(path_zenodo, "README.txt"))

#' The archive uploaded to Zenodo
zip_file <- normalizePath(file.path(dirname(path_zenodo), "ssa_incidence_analysis_data.zip"), mustWork = FALSE)
unlink(zip_file)
withr::with_dir(path_zenodo, zip(zip_file, list.files(), flags = "-X"))
