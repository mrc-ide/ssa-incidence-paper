################################################################################
#' Build the analysis datasets
#'
#' Reads the review extracts and additional study data (data/raw/), matches each
#' observation to UNAIDS areas (data/geography/) and to modelled incidence
#' (data/modelled/), and classifies studies by recruitment type.
#'
#' Writes (outputs/data/, not committed):
#'   analysis_data.csv                       one row per observation. extrapolated_value is the
#'                                           main offset (Naomi-extrapolated incidence);
#'                                           extrapolated_value_<variant> the sensitivity
#'                                           offsets (Table S6): no_phia, spectrum,
#'                                           spectrum_national
#'   flowchart_counts.csv                    studies/observations at each step (Figure 1)
#'
#' Zenodo: rather than regenerate these from the original files, they can be
#' downloaded from https://doi.org/10.5281/zenodo.24284889 into outputs/data/,
#' then start at R/02_fit_main.R.
################################################################################

library(tidyverse)

source("R/functions.R")
source("R/data_functions.R")

dir.create(path_derived, showWarnings = FALSE, recursive = TRUE)

pipeline <- align_observations()

modelled <- read_csv(file.path(path_modelled, "matched_incidence.csv.gz"), show_col_types = FALSE)

build_dataset <- function(offset_col) {
  pipeline$aligned %>%
    match_modelled(modelled, offset_col) %>%
    classify_studies(pipeline$inclusion) %>%
    mutate(sex_int = case_when(sex == "female" ~ 0, sex == "male" ~ 1, sex == "both" ~ 0.5),
           csi_dummy = as.integer(study_type == "Cross-sectional incidence")) %>%
    select(author, doi, study_type, idx, area_id, year, sex, age_group, extrapolated_value, ir, pys,
           pys_source, infections, representative, recruitment_type, study_idx, sex_int, csi_dummy) %>%
    arrange(idx)
}

#' Sensitivity offsets as extra columns: the observations are the same, only the
#' matched modelled incidence differs
sensitivity_offsets <- c(no_phia = "naomi_no_phia", spectrum = "spectrum", spectrum_national = "spectrum_national")

analysis_data <- build_dataset("naomi")
for (variant in names(sensitivity_offsets)) {
  offsets <- build_dataset(sensitivity_offsets[[variant]]) %>%
    select(idx, !!paste0("extrapolated_value_", variant) := extrapolated_value)
  stopifnot(setequal(offsets$idx, analysis_data$idx))
  analysis_data <- left_join(analysis_data, offsets, by = "idx")
}
analysis_data <- analysis_data %>%
  relocate(starts_with("extrapolated_value_"), .after = extrapolated_value) %>%
  filter(!is.na(extrapolated_value))
write_csv(analysis_data, file.path(path_derived, "analysis_data.csv"), na = "")

################################################################################
#' ## Figure 1 counts
################################################################################

#' Review studies that pass the exclusions but cannot be matched to an area, a year
#' or person-years (and so to modelled incidence) drop out in build_dataset()
review_final <- filter(analysis_data, idx <= pipeline$review_max_idx)
stopifnot(all(pipeline$primary_report_dois %in% str_trim(analysis_data$doi)))
flowchart_counts <- bind_rows(
  pipeline$tracker %>% mutate(box = "review"),
  tibble(step = "Matched to location, year and person-years", studies = n_distinct(review_final$doi),
         observations = nrow(review_final), box = "review"),
  pipeline$additional_counts %>% rename(step = source) %>% mutate(box = "additional"),
  pipeline$additional_sources %>% rename(step = source) %>% mutate(box = "additional studies identified"),
  tibble(step = "Final data", studies = n_distinct(analysis_data$doi), observations = nrow(analysis_data),
         box = "final")
)
write_csv(flowchart_counts, file.path(path_derived, "flowchart_counts.csv"))

################################################################################
#' ## Check against the submitted analysis dataset (private copy)
################################################################################

submitted_file <- file.path(path_private, "final_data_2026.csv")
if (file.exists(submitted_file)) {
  submitted <- read_csv(submitted_file, show_col_types = FALSE)
  cmp <- full_join(submitted, analysis_data, by = "idx", suffix = c("_sub", "_new"))
  message("Rows: submitted ", nrow(submitted), ", rebuilt ", nrow(analysis_data),
          "; only submitted ", sum(is.na(cmp$doi_new)), ", only rebuilt ", sum(is.na(cmp$doi_sub)))
  both <- filter(cmp, !is.na(doi_sub), !is.na(doi_new))
  for (v in c("ir", "pys", "extrapolated_value")) {
    d <- abs(both[[paste0(v, "_new")]] / both[[paste0(v, "_sub")]] - 1)
    message(v, ": ", sum(d > 1e-8, na.rm = TRUE), " differ (max relative ", signif(max(d, na.rm = TRUE), 3), ")")
  }
  for (v in c("doi", "area_id", "year", "sex", "age_group", "recruitment_type", "study_idx")) {
    message(v, ": ", sum(as.character(both[[paste0(v, "_new")]]) != as.character(both[[paste0(v, "_sub")]]), na.rm = TRUE), " differ")
  }
}
