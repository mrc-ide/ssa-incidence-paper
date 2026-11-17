################################################################################
#' Main regression model and its sensitivity analyses
#'
#' Ratio of empirical to modelled incidence, by recruitment type, sex and year
#' (Methods; Table S6). Writes to outputs/fits/:
#'   main_fit.rds                        the main INLA fit
#'   coefficients_<variant>.csv          exp(fixed effects) and exp(internal hyperparameters)
#'   main_predictions.csv                Fig 2: log ratio by recruitment type, sex and year
#'   modelled_correlation.csv            Fig S1: weighted correlation of empirical and modelled incidence
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
      # Recruitment type varies fastest within year
      id.year.recruitment = (match(year, years) - 1L) * length(recruitment_levels) + as.integer(recruitment_type)
    )
}

#' Constraints identifying the time trends.
#' Shared RW2 by sex: level and (scaled-year) slope are zero at the reference
#' year, so the linear trend is carried by the fixed effects.
#' Year x recruitment type RW2 deviations: zero level and zero slope (from the
#' previous year) at the reference year for each type, and summing to zero over
#' types in every other year.
model_constraints <- function(years, anchor_year = reference_year) {
  n_years <- length(years)
  n_types <- length(recruitment_levels)
  anchor_i <- match(anchor_year, years)
  previous_i <- match(anchor_year - 1L, years)
  stopifnot(!is.na(anchor_i), !is.na(previous_i))

  index <- function(year_i, type_i) (year_i - 1L) * n_types + type_i

  A_level <- A_slope <- matrix(0, nrow = n_types, ncol = n_years * n_types)
  for (type_i in seq_len(n_types)) {
    A_level[type_i, index(anchor_i, type_i)] <- 1
    A_slope[type_i, index(previous_i, type_i)] <- -1
    A_slope[type_i, index(anchor_i, type_i)] <- 1
  }
  sum_years <- setdiff(seq_len(n_years), c(previous_i, anchor_i))
  A_sum <- matrix(0, nrow = length(sum_years), ncol = n_years * n_types)
  for (row_i in seq_along(sum_years)) A_sum[row_i, index(sum_years[row_i], seq_len(n_types))] <- 1

  #' RW2 structure scaled to unit generalised variance (as scale.model = TRUE), so
  #' the precision is on the same scale as the sex-specific RW2
  Q_rw2 <- inla.scale.model(dfertility::make_rw_structure_matrix(n_years, 2),
                            constr = list(A = rbind(rep(1, n_years), seq_len(n_years)), e = c(0, 0)))

  list(
    A_year = rbind(as.numeric(years == anchor_year), as.vector(scale(years - anchor_year))),
    Q_interaction = kronecker(Q_rw2, as(diag(1, n_types), "sparseMatrix")),
    A_interaction = rbind(A_level, A_slope, A_sum),
    rankdef = n_years + n_types - 1
  )
}

fit_model <- function(df, years, family = "xpoisson", measurement = FALSE, num.threads = NULL,
                      anchor_year = reference_year) {
  cons <- model_constraints(years, anchor_year)
  formula <- infections ~ 0 + sex_int * std_year + recruitment_type +
    f(year2, model = "rw2", replicate = binary_sex, scale.model = TRUE,
      extraconstr = list(A = cons$A_year, e = matrix(0, nrow = 2))) +
    # PC prior: 1% probability that the typical deviation of a recruitment type's
    # trend from the shared trend exceeds 0.5 on the log scale
    f(id.year.recruitment, model = "generic0", Cmatrix = cons$Q_interaction, rankdef = cons$rankdef,
      hyper = list(prec = list(prior = "pc.prec", param = c(0.5, 0.01))),
      extraconstr = list(A = cons$A_interaction, e = matrix(0, nrow = nrow(cons$A_interaction)))) +
    f(study_idx, model = "iid", constr = TRUE) +
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
  # main reference year); the grid keeps 2014 for the 2014-15 slope constraint
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
