################################################################################
#' Leave-one-study-out analysis (Supplementary Figure S2)
#'
#' Refits the main model once per study (doi) with that study's observations
#' removed. Writes outputs/fits/loo_coefficients.csv: exp(fixed effects) for every
#' refit, with the left-out study's recruitment type.
#'
#' Zenodo: outputs/data/ can be downloaded rather than regenerated (see
#' R/02_fit_main.R); results then will not exactly reproduce the manuscript.
#'
#' Rscript R/03_fit_loo.R [cores]   (default 8; ~20 s per fit and core)
################################################################################

source("R/02_fit_main.R")

cores <- as.integer(commandArgs(trailingOnly = TRUE)[1])
if (is.na(cores)) cores <- 8L

data <- read_analysis_data()
studies <- data %>% distinct(doi, study_idx, recruitment_type)
stopifnot(!anyDuplicated(studies$doi))

progress_file <- file.path(path_fits, "loo_progress.log")
cat("", file = progress_file)

# Fit the cohort slope prior once here, so the forked workers share it
invisible(slope_prior())

loo <- parallel::mclapply(studies$doi, function(x) {
  df <- model_frame(filter(data, doi != x), years_main)
  # One thread per fit: parallelism is across studies
  fit <- fit_model(df, years_main, num.threads = "1:1")
  cat(x, "\n", file = progress_file, append = TRUE)
  summarise_inla(fit) %>% filter(type == "fixed") %>% mutate(doi = x)
}, mc.cores = cores)

failed <- !map_lgl(loo, is.data.frame)
if (any(failed)) stop("Failed refits: ", paste(studies$doi[failed], collapse = ", "))

bind_rows(loo) %>%
  mutate(var = str_remove(var, "^recruitment_type")) %>%
  left_join(studies, by = "doi") %>%
  select(doi, study_idx, recruitment_type, var, mean, lower, upper, txt) %>%
  write_csv(file.path(path_fits, "loo_coefficients.csv"))
