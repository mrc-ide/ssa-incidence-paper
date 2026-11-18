################################################################################
#' Population cohorts and PHIA surveys
#'
#' Writes to outputs/fits/:
#'   coefficients_cohort_vs_model.csv   Text S6: cohort incidence relative to matched modelled incidence
#'   cohort_trends.csv                  Table S8: 15-49 incidence by cohort and year, crude and smoothed,
#'                                      relative to 2010 (2011 in Kisesa and Siaya)
#'   cohort_declines.csv                declines since 2010 to the latest year, with 95% CrI (Results)
#'   coefficients_cohort_age.csv        Table S9: cohort age/sex/time model
#'   cohort_age_trends.csv              annual IRR by sex and age group (Results; year + interactions)
#'   cohort_age_predictions.csv         Fig 5A / S3: incidence by cohort, year, sex and age (15-24, 25-49)
#'   cohort_population.csv              Table S11: HIV-negative population in each cohort's district
#'   cohort_infection_shares.csv        Fig 5C / S4: share of new 15-49 infections by sex and age
#'   coefficients_phia_age.csv          Table S10: PHIA age/sex/time model
#'   phia_age_trends.csv                annual IRR by sex and age group (Results)
#'   phia_spectrum.csv                  Fig 4A: PHIA 15-49 incidence and national modelled incidence
#'   phia_spectrum_fit.csv              Fig 4B: fitted relationship of PHIA to modelled incidence
#'   phia_spectrum_correlation.csv      Fig 4: weighted correlation of log incidence, by PHIA phase
#'
#' outputs/fits is not committed.
#'
#' Zenodo: outputs/data/ can be downloaded rather than regenerated (see
#' R/02_fit_main.R); results then will not exactly reproduce the manuscript.
################################################################################

source("R/functions.R")
library(INLA)

dir.create(path_fits, showWarnings = FALSE, recursive = TRUE)

data <- read_csv(file.path(path_derived, "analysis_data.csv"), show_col_types = FALSE)
cohort <- filter(data, doi %in% cohort_dois)

#' Base year for declines: 2010, or 2011 where the cohort has no 2010 round
base_year <- c("Kisesa" = 2011, "Kisumu HDSS" = 2011)
cohort_base_year <- function(doi) ifelse(doi %in% names(base_year), base_year[doi], 2010)

################################################################################
#' ## Text S6: cohort incidence relative to modelled incidence
################################################################################

cohort_vs_model_df <- cohort %>%
  mutate(id.year = year - min(year),
         std_year = (year - reference_year) / 5,
         id.doi = as.integer(factor(doi)),
         obs_iid = row_number())

#' Each cohort is one study, and its AR1 trend (rho ~0.97) already carries the
#' cohort's level, so there is no separate study intercept (its variance went to
#' zero when included). PC priors: 1% probability that an SD exceeds 1 on the log
#' scale; P(rho > 0.7) = 0.7.
pc_sd <- list(prec = list(prior = "pc.prec", param = c(1, 0.01)))

cohort_vs_model <- inla(
  infections ~ 1 + std_year * sex +
    f(id.year, model = "ar1", group = id.doi, control.group = list(model = "iid"),
      hyper = list(prec = pc_sd$prec, rho = list(prior = "pc.cor1", param = c(0.7, 0.7)))) +
    f(obs_iid, model = "iid", hyper = pc_sd) +
    offset(log(extrapolated_value)),
  family = "xpoisson", data = cohort_vs_model_df, E = pys, control.compute = list(config = TRUE))

summarise_inla(cohort_vs_model) %>%
  write_csv(file.path(path_fits, "coefficients_cohort_vs_model.csv"))

################################################################################
#' ## Table S8: 15-49 incidence trends within cohorts
################################################################################

cohort_15to49 <- cohort %>%
  filter(age_group %in% c(age_groups_15to49, "Y045_054")) %>%
  group_by(study_idx, area_id, doi, year) %>%
  summarise(infections = sum(infections), pys = sum(pys), .groups = "drop") %>%
  mutate(std_year = year - reference_year,
         id.year = year - min(year),
         id.doi = as.integer(factor(doi)))

trend_fit <- inla(
  infections ~ 1 + std_year +
    f(id.year, model = "ar1", group = id.doi, control.group = list(model = "iid")) +
    f(study_idx, model = "iid", constr = TRUE),
  family = "xpoisson", data = cohort_15to49, E = pys, control.compute = list(config = TRUE))

trend_samples <- exp(predictor_samples(trend_fit, seq_len(nrow(cohort_15to49))))

cohort_trends <- cohort_15to49 %>%
  mutate(ir = infections / pys,
         smoothed = rowMeans(trend_samples),
         base = cohort_base_year(doi)) %>%
  group_by(doi) %>%
  mutate(crude_irr = ir / ir[year == base],
         smoothed_irr = smoothed / smoothed[year == base]) %>%
  ungroup() %>%
  select(doi, year, infections, pys, ir, crude_irr, smoothed, smoothed_irr)

write_csv(cohort_trends, file.path(path_fits, "cohort_trends.csv"))

#' Decline from the base year to the latest year, per posterior sample
cohort_declines <- cohort_15to49 %>%
  mutate(row = row_number(), base = cohort_base_year(doi)) %>%
  group_by(doi) %>%
  filter(any(year == base[1]), max(year) > base[1]) %>%
  summarise(base_year = base[1], latest_year = max(year),
            decline = list(1 - trend_samples[row[year == max(year)], ] / trend_samples[row[year == base[1]], ]),
            .groups = "drop") %>%
  mutate(median = map_dbl(decline, median),
         lower = map_dbl(decline, quantile, 0.025),
         upper = map_dbl(decline, quantile, 0.975)) %>%
  select(-decline)

write_csv(cohort_declines, file.path(path_fits, "cohort_declines.csv"))

################################################################################
#' ## Table S9: age and sex patterns of cohort incidence (15-24, 25-49)
################################################################################

cohort_age_data <- cohort %>%
  filter(age_group %in% c(age_groups_15to49, "Y035_044", "Y045_054")) %>%
  mutate(age_group = ifelse(age_group %in% c("Y015_019", "Y020_024"), "Y015_024", "Y025_049")) %>%
  group_by(doi, area_id, year, sex, age_group) %>%
  summarise(infections = sum(ir * pys), pys = sum(pys), .groups = "drop") %>%
  mutate(ir = infections / pys)

#' Prediction rows first: every cohort, year 1990-2022, sex and age group
cohort_age_df <- crossing(year = 1990:2022, distinct(cohort_age_data, doi, area_id),
                          sex = c("male", "female"), age_group = c("Y015_024", "Y025_049")) %>%
  mutate(pys = 1) %>%
  bind_rows(cohort_age_data) %>%
  mutate(year_min = year - reference_year,
         year_min2 = year_min,
         year_min_f = ifelse(sex == "female", year_min, NA),
         year_min_m = ifelse(sex == "male", year_min, NA),
         doi_int = as.integer(factor(doi)),
         doi3 = doi,
         age_group2 = factor(age_group, levels = c("Y015_024", "Y025_049")),
         age_int = ifelse(age_group == "Y025_049", 1L, 2L),
         #' Cohort-specific linear change in the 25-49 : 15-24 log ratio, coded +/- 1/2
         #' so it shifts the two age groups symmetrically
         slope_id = doi_int,
         age_slope = year_min * ifelse(age_group == "Y025_049", 0.5, -0.5))

cohort_age_fit <- inla(
  infections ~ year_min * sex + age_group2 * sex + age_group2 * year_min +
    f(year_min2, model = "ar", order = 2, replicate = age_int, constr = TRUE) +
    f(year_min_f, model = "ar", order = 2, group = doi_int, control.group = list(model = "iid")) +
    f(year_min_m, model = "ar", order = 2, group = doi_int, control.group = list(model = "iid")) +
    #' Sum-to-zero within each age group, so the age-group fixed effect is identified
    #' (unconstrained, it trades off against the mean of the cohort intercepts)
    f(doi3, replicate = age_int, constr = TRUE) +
    #' PC prior: 1% probability that the SD of cohort slopes exceeds 0.1 per year
    f(slope_id, age_slope, model = "iid", hyper = list(prec = list(prior = "pc.prec", param = c(0.1, 0.01)))),
  family = "xpoisson", data = cohort_age_df, E = pys, control.compute = list(config = TRUE))

summarise_inla(cohort_age_fit) %>%
  write_csv(file.path(path_fits, "coefficients_cohort_age.csv"))

age_sex_trends(cohort_age_fit, "year_min:age_group2Y025_049") %>%
  write_csv(file.path(path_fits, "cohort_age_trends.csv"))

cohort_age_pred <- sample_predictor(cohort_age_fit, cohort_age_df) %>%
  select(doi, area_id, year, sex, age_group, mean, median, lower, upper) %>%
  mutate(across(mean:upper, exp))

write_csv(cohort_age_pred, file.path(path_fits, "cohort_age_predictions.csv"))

################################################################################
#' ## Table S11, Fig 5C: shares of new infections by sex and age
#' Modelled incidence x HIV-negative population of the cohort district
################################################################################

cohort_range <- cohort %>%
  group_by(doi) %>%
  summarise(first_year = min(year), last_year = max(year))

cohort_hiv_neg <- read_csv(file.path(path_modelled, "cohort_population.csv"), show_col_types = FALSE) %>%
  mutate(hiv_negative = population * (1 - prevalence)) %>%
  inner_join(distinct(cohort, doi, area_id), by = "area_id", relationship = "many-to-many")

cohort_hiv_neg %>%
  left_join(cohort_range, by = "doi") %>%
  filter(year >= first_year, year <= last_year) %>%
  select(doi, area_id, year, sex, age_group, population, prevalence, hiv_negative) %>%
  arrange(doi, sex, age_group, year) %>%
  write_csv(file.path(path_fits, "cohort_population.csv"))

cohort_age_pred %>%
  left_join(cohort_hiv_neg, by = c("doi", "area_id", "year", "sex", "age_group")) %>%
  mutate(infections = mean * hiv_negative) %>%
  group_by(doi, year) %>%
  mutate(share = infections / sum(infections)) %>%
  ungroup() %>%
  left_join(cohort_age_data %>% group_by(doi) %>% summarise(first_year = min(year), last_year = max(year)),
            by = "doi") %>%
  mutate(in_observed_range = year >= first_year & year <= last_year) %>%
  select(doi, year, sex, age_group, incidence = mean, hiv_negative, infections, share, in_observed_range) %>%
  write_csv(file.path(path_fits, "cohort_infection_shares.csv"))

################################################################################
#' ## Table S10: age and sex patterns of PHIA incidence (public extract)
################################################################################

phia_age_data <- read_csv(file.path(path_raw, "phia-incid-extract_2025-11-13b.csv"), show_col_types = FALSE) %>%
  filter(age_group %in% c("Y015_024", "Y025_049"), sex != "both") %>%
  select(iso3, survey_id, year = survey_year, age_group, sex, incid, pys = n_pys_approx)

phia_age_df <- crossing(year = 2016:2022, sex = c("male", "female"), age_group = c("Y015_024", "Y025_049")) %>%
  mutate(pys = 1) %>%
  bind_rows(phia_age_data) %>%
  mutate(infections = incid * pys,
         year_min = year - reference_year,
         iso3 = factor(iso3),
         iso32 = iso3,
         age_group = factor(age_group, levels = c("Y015_024", "Y025_049")))

phia_age_fit <- inla(
  infections ~ year_min * sex + age_group * sex + age_group * year_min +
    f(iso32) +
    f(iso3, year_min, model = "iid", hyper = list(prec = list(prior = "normal", param = c(2, 1)))),
  family = "xpoisson", data = phia_age_df, E = pys, control.compute = list(config = TRUE))

summarise_inla(phia_age_fit) %>%
  write_csv(file.path(path_fits, "coefficients_phia_age.csv"))

age_sex_trends(phia_age_fit, "year_min:age_groupY025_049") %>%
  write_csv(file.path(path_fits, "phia_age_trends.csv"))

################################################################################
#' ## Fig 4: PHIA vs national modelled incidence, adults 15-49 (internal extract)
################################################################################

#' National file only (the extract also holds subnational EPP and Thembisa areas)
spectrum_national <- read_csv(file.path(path_modelled, "spectrum_incidence.csv.gz"), show_col_types = FALSE) %>%
  filter(variant == "published", sex == "both", age_group == "Y015_049", area_id == iso3) %>%
  transmute(iso3, year, modelled = infections / (population - plhiv))

phia_spectrum <- read_csv(file.path(path_raw, "phia-incid-extract_INTERNAL_2024-05-07.csv"), show_col_types = FALSE) %>%
  filter(age_group == "Y015_049", sex == "both") %>%
  select(survey_id, iso3, year = survey_year, phia_phase, incid, lower = ci_l.incid, upper = ci_u.incid,
         pys = n_pys_approx) %>%
  left_join(spectrum_national, by = c("iso3", "year")) %>%
  mutate(modelled_in_ci = modelled >= lower & modelled <= upper)

write_csv(phia_spectrum, file.path(path_fits, "phia_spectrum.csv"))

phia_cor <- function(d) {
  tibble(r = wCorr::weightedCorr(log(d$incid), log(d$modelled), weights = d$pys, method = "Pearson"),
         n_surveys = nrow(d), n_modelled_in_ci = sum(d$modelled_in_ci))
}
bind_rows(
  all = phia_cor(phia_spectrum),
  PHIA1 = phia_cor(filter(phia_spectrum, phia_phase == "PHIA1")),
  PHIA2 = phia_cor(filter(phia_spectrum, phia_phase == "PHIA2")),
  .id = "phase"
) %>%
  write_csv(file.path(path_fits, "phia_spectrum_correlation.csv"))

phia_fit_df <- tibble(modelled = seq(0.0001, 0.023, 0.001), pys = 1) %>%
  bind_rows(mutate(phia_spectrum, infections = incid * pys))

phia_fit <- inla(infections ~ log(modelled) + f(survey_id, model = "iid"),
                 family = "xpoisson", data = phia_fit_df, E = pys, control.compute = list(config = TRUE))

sample_predictor(phia_fit, phia_fit_df) %>%
  select(modelled, mean, median, lower, upper) %>%
  mutate(across(mean:upper, exp)) %>%
  write_csv(file.path(path_fits, "phia_spectrum_fit.csv"))
