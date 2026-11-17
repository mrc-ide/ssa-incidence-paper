################################################################################
#' Public copy of the analysis datasets for Zenodo
#'
#' Needs outputs/data/ (R/01_data.R, run with private/ present). Run from the
#' repository root. Writes outputs/zenodo/: the analysis datasets read by
#' R/02_fit_main.R onwards, and README.txt; and outputs/ssa_incidence_analysis_data.zip
#' of those files, for upload.
#'
#' Release rule: Namibia's modelled incidence (the offset) comes from the
#' non-public Namibia 2026 Spectrum file, so Namibia observations are removed.
#' They are dropped rather than given a blank offset because INLA reads a missing
#' offset as zero, which would keep them in the fits with the wrong offset.
#' Results from this copy therefore differ slightly from the manuscript.
################################################################################

library(tidyverse)

source("R/functions.R")

path_zenodo <- "outputs/zenodo"
dir.create(path_zenodo, showWarnings = FALSE, recursive = TRUE)

private_iso3 <- "NAM"
is_private <- function(area_id) str_detect(area_id, paste0("(^|, )", private_iso3, "(_|$)"))

data <- read_csv(file.path(path_derived, "analysis_data.csv"), show_col_types = FALSE)
public <- filter(data, !is_private(area_id))
message("Removed ", nrow(data) - nrow(public), " Namibia observations")
write_csv(public, file.path(path_zenodo, "analysis_data.csv"), na = "")

#' Figure 1 counts describe the review and are not redacted
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
  "",
  "Namibia observations have been removed: their modelled incidence comes from a",
  "non-public Spectrum file. Model fits, figures and numbers from these files will",
  "therefore not exactly reproduce the manuscript results. flowchart_counts.csv",
  "(Figure 1) gives the counts of the full dataset, including Namibia."
), file.path(path_zenodo, "README.txt"))

#' The archive uploaded to Zenodo
zip_file <- normalizePath(file.path(dirname(path_zenodo), "ssa_incidence_analysis_data.zip"), mustWork = FALSE)
unlink(zip_file)
withr::with_dir(path_zenodo, zip(zip_file, list.files(), flags = "-X"))
