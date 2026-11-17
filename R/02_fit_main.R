################################################################################
#' Main regression model and its sensitivity analyses
#'
#' Ratio of empirical to modelled incidence, by recruitment type, sex and year
#' (Methods; Table S6). Writes to outputs/fits/:
#'   main_fit.rds                        the main INLA fit
#'   coefficients_<variant>.csv          exp(fixed effects) and exp(internal hyperparameters)
#'   main_predictions.csv                Fig 2: log ratio by recruitment type, sex and year
#'   modelled_correlation.csv            Fig S1: weighted correlation of empirical and modelled incidence
#'   cohort_slope_prior.csv              prior for the SD of study random slopes, from the cohorts
#'
#' Zenodo: to start here without regenerating the data from the original files,
#' download the analysis datasets from https://doi.org/10.5281/zenodo.24284889
#' into outputs/data/ and skip R/01_data.R. That copy has no Namibia observations (non-public
#' modelled incidence), so results will not exactly reproduce the manuscript.
#'
#' Run all variants:   Rscript R/02_fit_main.R
#' Run some variants:  Rscript R/02_fit_main.R main nb
################################################################################

source("R/functions.R")
library(INLA)

dir.create(path_fits, showWarnings = FALSE, recursive = TRUE)

#' Analysis data with the main offset, or with extrapolated_value_<variant> as the offset
read_analysis_data <- function(variant = "") {
  data <- read_csv(file.path(path_derived, "analysis_data.csv"), show_col_types = FALSE)
  if (variant != "") data$extrapolated_value <- data[[paste0("extrapolated_value_", variant)]]
  data %>%
    select(-starts_with("extrapolated_value_")) %>%
    filter(!is.na(extrapolated_value))
}

#' Prediction grid (rows first, no outcome) followed by the observations within `years`.
#' Grid rows have no offset, so their linear predictor is the log ratio of
#' empirical to modelled incidence for each recruitment type, sex and year.
model_frame <- function(data, years, anchor_year = reference_year) {
  grid <- crossing(year = years, sex = c("female", "male"), recruitment_type = recruitment_levels, pys = 1) %>%
    filter(!(sex == "male" & recruitment_type %in% c("ANC", "FSW")),
           !(sex == "female" & recruitment_type == "MSM"))

  bind_rows(grid, data %>% mutate(obs_iid = row_number())) %>%
    filter(year %in% years) %>%
    mutate(
      year2 = year,
      recruitment_type = factor(recruitment_type, levels = recruitment_levels),
      sex_int = case_when(sex == "female" ~ 0, sex == "both" ~ 0.5, sex == "male" ~ 1),
      binary_sex = ifelse(sex == "female", 1, 2),
      std_year = (year - anchor_year) / 5,
      # Measurement method (Table S7): prospective follow-up, cross-sectional incidence or retrospective cohort
      study2 = case_when(
        study_type %in% c("Intervention Trials", "RCT", "Prospective Cohort") ~ "aaprospective",
        study_type %in% c("Cross-sectional incidence", "Retrospective Cohort") ~ study_type
      ) %>%
        factor(levels = c("aaprospective", "Cross-sectional incidence", "Retrospective Cohort")),
      slope_idx = study_idx
    ) %>%
    # Covariate for the study random slope: scaled year centred on the study's
    # pys-weighted mean year. It is 0 for single-year studies, so only studies
    # observed in 2+ years inform (or are affected by) the slopes. NA on grid rows.
    group_by(study_idx) %>%
    mutate(slope_x = std_year - weighted.mean(std_year, pys)) %>%
    ungroup() %>%
    mutate(slope_x = ifelse(is.na(study_idx), NA_real_, slope_x))
}

#' Prior for the SD of study random slopes on scaled year, from the population
#' cohorts (the studies with long within-study follow-up): a cohort-only model with
#' a cohort intercept and random slope (PC priors, 1% probability an SD exceeds 1).
#' Returns the posterior mean and precision of the slopes' log precision, used as
#' a normal prior on the log precision in the main model.
cohort_slope_prior <- function(data = read_analysis_data()) {
  co <- data %>%
    filter(doi %in% cohort_dois) %>%
    mutate(std_year = (year - reference_year) / 5) %>%
    group_by(doi) %>%
    mutate(slope_x = std_year - weighted.mean(std_year, pys)) %>%
    ungroup() %>%
    mutate(cohort_idx = as.integer(factor(doi)), slope_idx = cohort_idx, obs_iid = row_number())
  pc_sd <- list(prec = list(prior = "pc.prec", param = c(1, 0.01)))
  fit <- inla(infections ~ 1 + std_year * sex +
                f(cohort_idx, model = "iid", hyper = pc_sd) +
                f(slope_idx, slope_x, model = "iid", hyper = pc_sd) +
                f(obs_iid, model = "iid", hyper = pc_sd) +
                offset(log(extrapolated_value)),
              family = "xpoisson", data = co, E = pys)
  m <- fit$internal.marginals.hyperpar[["Log precision for slope_idx"]]
  mean <- inla.emarginal(identity, m)
  sd <- sqrt(inla.emarginal(function(x) x^2, m) - mean^2)
  q <- inla.qmarginal(c(0.975, 0.5, 0.025), m)
  list(mean = mean, prec = 1 / sd^2, sd_slope = exp(-q / 2))
}

#' Computed once per session (the same prior for every variant and leave-one-out fit)
slope_prior <- local({
  cache <- NULL
  function() {
    if (is.null(cache)) cache <<- cohort_slope_prior()
    cache
  }
})

#' Constraint identifying the shared RW2 by sex: level and (scaled-year) slope
#' are zero at the reference year, so the linear trend is carried by the fixed
#' effects. (A recruitment type x year RW2 deviation was dropped: with study
#' random slopes its SD was estimated at zero and the fixed effects were unchanged.)
model_constraints <- function(years, anchor_year = reference_year) {
  stopifnot(anchor_year %in% years)
  list(A_year = rbind(as.numeric(years == anchor_year), as.vector(scale(years - anchor_year))))
}

fit_model <- function(df, years, family = "xpoisson", measurement = FALSE, num.threads = NULL,
                      anchor_year = reference_year, slope_hyper = slope_prior()) {
  cons <- model_constraints(years, anchor_year)
  slope_prec <- list(prec = list(prior = "normal", param = c(slope_hyper$mean, slope_hyper$prec)))
  formula <- infections ~ 0 + sex_int * std_year + recruitment_type +
    f(year2, model = "rw2", replicate = binary_sex, scale.model = TRUE,
      extraconstr = list(A = cons$A_year, e = matrix(0, nrow = 2))) +
    f(study_idx, model = "iid", constr = TRUE) +
    # Study random slopes, so a single study's within-study trend (e.g. Seri 2018,
    # 18% of person-years) cannot set the shared year effect; prior from the cohorts
    f(slope_idx, slope_x, model = "iid", hyper = slope_prec) +
    f(obs_iid, model = "iid") +
    offset(log(extrapolated_value))
  if (measurement) formula <- update(formula, . ~ . + study2)

  inla(formula, family = family, data = df, E = pys, control.compute = list(config = TRUE),
       num.threads = num.threads %||% inla.getOption("num.threads"))
}

################################################################################
#' ## Model variants (Table S6, S7, S12)
################################################################################

years_main <- 1990:2022

variants <- list(
  main = list(data = "", years = years_main),
  # Negative binomial likelihood requires integer counts
  nb = list(data = "", years = years_main, family = "nbinomial",
            prepare = function(d) mutate(d, infections = round(infections))),
  # Matched to Spectrum/EPP fits without PHIA incidence data
  no_phia = list(data = "no_phia", years = years_main),
  # Spectrum subnational where available, otherwise national (no Naomi)
  spectrum_subnational = list(data = "spectrum", years = years_main),
  spectrum_national = list(data = "spectrum_national", years = years_main),
  since_2010 = list(data = "", years = 2010:2022),
  # Observations from 2015, so this variant is anchored at 2015 (no data in the
  # main reference year); the grid keeps 2014 as in the submitted analysis
  since_2015 = list(data = "", years = 2014:2022, anchor_year = 2015L, prepare = function(d) filter(d, year > 2014)),
  # Measurement method, total population only (FSW/MSM observations excluded; their grid rows kept)
  measurement = list(data = "", years = years_main, measurement = TRUE,
                     prepare = function(d) filter(d, !recruitment_type %in% c("FSW", "MSM"))),
  # Table S12: drop observations whose person-years were imputed from the
  # regression on number of participants (Text S2 steps 2-3). The submitted
  # "37 observations from 13 studies" was counted on an earlier version of the
  # data (36 / 13 of these observations are in the source's final_data.csv).
  denominators = list(data = "", years = years_main,
                      prepare = function(d) filter(d, !str_detect(pys_source, "participants")))
)

run_variant <- function(name) {
  v <- variants[[name]]
  data <- read_analysis_data(v$data)
  if (!is.null(v$prepare)) data <- v$prepare(data)
  anchor <- v$anchor_year %||% reference_year
  df <- model_frame(data, v$years, anchor)
  t0 <- Sys.time()
  fit <- fit_model(df, v$years, family = v$family %||% "xpoisson", measurement = isTRUE(v$measurement),
                   anchor_year = anchor)
  message(sprintf("%s: %d observations, fitted in %.1f min", name, sum(!is.na(df$infections)),
                  as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  summarise_inla(fit) %>%
    mutate(variant = name, n_obs = sum(!is.na(df$infections)), n_studies = n_distinct(df$study_idx, na.rm = TRUE)) %>%
    write_csv(file.path(path_fits, paste0("coefficients_", name, ".csv")))
  list(fit = fit, df = df)
}

#' Fig 2: posterior log ratio for each recruitment type and sex, restricted to
#' the years spanned by that type's sex-specific observations.
main_predictions <- function(fit, df, data) {
  observed <- data %>%
    filter(sex != "both") %>%
    group_by(sex, recruitment_type) %>%
    summarise(first_year = min(year), last_year = max(year), .groups = "drop")

  sample_predictor(fit, df) %>%
    select(year, sex, recruitment_type, mean, median, lower, upper) %>%
    mutate(recruitment_type = as.character(recruitment_type)) %>%
    left_join(observed, by = c("sex", "recruitment_type")) %>%
    mutate(in_observed_range = !is.na(first_year) & year >= first_year & year <= last_year) %>%
    select(-first_year, -last_year)
}

#' Fig S1: pys-weighted correlation of log empirical and log modelled incidence,
#' sex-specific total-population observations. Zero-incidence observations are
#' omitted (main) or given 0.1 infections (zero-cell correction, Fig S1 caption).
modelled_correlation <- function(data) {
  d <- data %>% filter(sex != "both", !recruitment_type %in% c("FSW", "MSM"))
  wcor <- function(d) {
    est <- weights::wtd.cor(x = log(d$ir), y = log(d$extrapolated_value), weight = d$pys)
    tibble(r = est[1, "correlation"], se = est[1, "std.err"], n_obs = nrow(d))
  }
  bind_rows(
    omit_zero = wcor(filter(d, ir != 0)),
    zero_cell_correction = wcor(mutate(d, ir = ifelse(infections == 0, 0.1, infections) / pys)),
    .id = "zero_incidence"
  ) %>%
    mutate(lower = r - 1.96 * se, upper = r + 1.96 * se)
}

################################################################################

# Fit only when run as a script; source() defines the functions and variants
if (sys.nframe() == 0L) {
  to_run <- commandArgs(trailingOnly = TRUE)
  if (length(to_run) == 0) to_run <- names(variants)
  stopifnot(all(to_run %in% names(variants)))

  with(slope_prior(), tibble(log_prec_mean = mean, log_prec_sd = 1 / sqrt(prec),
                             sd_slope_lower = sd_slope[1], sd_slope_median = sd_slope[2], sd_slope_upper = sd_slope[3])) %>%
    write_csv(file.path(path_fits, "cohort_slope_prior.csv"))

  for (name in to_run) {
    res <- run_variant(name)
    if (name == "main") {
      saveRDS(res$fit, file.path(path_fits, "main_fit.rds"))
      data <- read_analysis_data()
      write_csv(main_predictions(res$fit, res$df, data), file.path(path_fits, "main_predictions.csv"))
      write_csv(modelled_correlation(data), file.path(path_fits, "modelled_correlation.csv"))
    }
  }
}
