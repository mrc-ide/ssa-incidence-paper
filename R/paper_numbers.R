################################################################################
#' # Numbers and tables quoted in the manuscript
#'
#' Every number quoted in "HIV incidence declines and changing age patterns of new
#' infections in sub-Saharan Africa" and its supplement goes through `record()`,
#' with where it is used and its value from this analysis. All are written to
#' outputs/paper_numbers.csv.
#'
#' Tables are written to outputs/tables/, as they appear in the manuscript.
#'
#' Reads only: outputs/data/ (R/01_data.R) and outputs/fits/ (R/02-04). Run from
#' the repository root: Rscript R/paper_numbers.R
#'
#' Zenodo: outputs/data/ can be downloaded rather than regenerated (see
#' R/02_fit_main.R); numbers then will not exactly reproduce the manuscript.
#'
#' Numbers not reproducible from code (Covidence screening counts, numbers from
#' cited publications, the Spectrum file metadata) are recorded with a `note`.
################################################################################

source("R/functions.R")

dir.create(path_tables, showWarnings = FALSE, recursive = TRUE)

data <- read_csv(file.path(path_derived, "analysis_data.csv"), show_col_types = FALSE)
flowchart <- read_csv(file.path(path_derived, "flowchart_counts.csv"), show_col_types = FALSE)

read_fit <- function(file) read_csv(file.path(path_fits, file), show_col_types = FALSE)

variants <- c("main", "nb", "no_phia", "spectrum_subnational", "spectrum_national",
              "since_2010", "since_2015", "measurement", "denominators")
coefs <- map_dfr(variants, ~ read_fit(paste0("coefficients_", .x, ".csv")))
main <- filter(coefs, variant == "main")

################################################################################
#' ## Recording
################################################################################

paper_numbers <- tibble(where = character(), label = character(), value = character(),
                        note = character())

#' Store and print one number quoted in the manuscript
record <- function(where, label, value, note = NA) {
  value <- as.character(value)
  paper_numbers <<- bind_rows(paper_numbers, tibble(where = where, label = label, value = value,
                                                    note = as.character(note)))
  message(sprintf("[%s] %s: %s", where, label, value))
  invisible(value)
}

#' Coefficient of fit `fits` (a coefficients table) as "1.05 (0.71, 1.56)"
coef_txt <- function(fits, term, d = 2) {
  x <- filter(fits, var == term)
  stopifnot(nrow(x) == 1)
  fmt_ci(x$mean, x$lower, x$upper, d)
}

coef_value <- function(fits, term, stat = "mean") filter(fits, var == term)[[stat]]

#' Percentage change per year from a ratio and its interval, e.g. 0.94 (0.92, 0.96) -> "6% (4, 8%)"
#' (a decline is positive)
pct_decline <- function(median, lower, upper) {
  sprintf("%.0f%% (%.0f, %.0f%%)", 100 * (1 - median), 100 * (1 - upper), 100 * (1 - lower))
}

#' Write a manuscript table to outputs/tables/<name>.csv
write_table <- function(table, name) {
  write_csv(table, file.path(path_tables, paste0(name, ".csv")), na = "")
  message(sprintf("%s: %d rows", name, nrow(table)))
  invisible(table)
}

################################################################################
#' ## Data description (Abstract, Results paragraph 1, Table S5)
################################################################################

n_obs <- nrow(data)
n_studies <- n_distinct(data$doi)

#' One row per observation and country it covers (multi-country studies list several areas)
obs_country <- data %>%
  select(doi, idx, area_id, recruitment_type) %>%
  separate_longer_delim(area_id, ", ") %>%
  mutate(iso3 = substr(area_id, 1, 3)) %>%
  distinct(doi, idx, iso3, recruitment_type) %>%
  left_join(region_lookup, by = "iso3")
stopifnot(!anyNA(obs_country$region))

n_countries <- n_distinct(obs_country$iso3)

record("Abstract; Results", "Observations", fmt_int(n_obs))
record("Abstract; Results", "Studies", n_studies)
record("Abstract; Results", "Countries", n_countries)
record("Results", "Infections", fmt_int(sum(data$infections)))
record("Results", "Person-years (million)", sprintf("%.1f", sum(data$pys) / 1e6))

studies_by_country <- obs_country %>%
  distinct(doi, iso3) %>%
  count(iso3, sort = TRUE)
top5 <- head(studies_by_country$iso3, 5)
n_top5 <- n_distinct(obs_country$doi[obs_country$iso3 %in% top5])

record("Results", "Share of studies from five countries",
       sprintf("%s (%d/%d)", fmt_pct(n_top5 / n_studies), n_top5, n_studies))
for (iso in top5) {
  record("Results", paste("Studies in", country_name(iso)),
         studies_by_country$n[studies_by_country$iso3 == iso])
}
#' The text names these five countries
stopifnot(setequal(top5, c("ZAF", "KEN", "UGA", "ZWE", "MWI")))

type_obs <- count(data, recruitment_type)
n_type <- function(types) sum(type_obs$n[type_obs$recruitment_type %in% types])
record("Results", "Observations from population cohorts and UTT trials",
       sprintf("%d; %s", n_type("UTT/Cohort"), fmt_pct(n_type("UTT/Cohort") / n_obs)))
record("Results", "Observations from household surveys",
       sprintf("%d; %s", n_type(c("National HHS", "Other HHS")), fmt_pct(n_type(c("National HHS", "Other HHS")) / n_obs)))

csi_studies <- data %>% filter(csi_dummy == 1) %>% distinct(doi)
record("Results", "Studies using cross-sectional incidence assays", nrow(csi_studies))
record("Results", "of which PHIA surveys", sum(str_detect(csi_studies$doi, "PHIA$")))

wca <- filter(obs_country, region == "WCA")
record("Results", "Studies in Western and Central Africa", n_distinct(wca$doi))
record("Results", "of which total population (non-FSW/MSM)",
       n_distinct(wca$doi[!wca$recruitment_type %in% c("FSW", "MSM")]))

n_by_sex <- count(data, sex)
n_sex <- function(s) n_by_sex$n[n_by_sex$sex == s]
record("Results", "Observations stratified by sex",
       sprintf("%s; %d/%d", fmt_pct(1 - n_sex("both") / n_obs), n_obs - n_sex("both"), n_obs))
record("Results", "Share of observations among women", fmt_pct(n_sex("female") / n_obs))
record("Results", "Share of observations among men", fmt_pct(n_sex("male") / n_obs))

record("Results", "Median year of observation (IQR)", med_iqr(data$year))
record("Methods (Regression analysis); Text S5", "Reference year (median study year)", reference_year,
       note = sprintf("Median year of observation is %d", median(data$year)))

#' Table S5: description of observations (as make_descriptive_table in the source, clean.R)
s5_block <- function(df, var, study_df = df) {
  studies <- study_df %>% distinct(doi, .data[[var]]) %>% count(.data[[var]], name = "n_studies")
  obs <- df %>% distinct(idx, .data[[var]]) %>% count(.data[[var]], name = "n_obs")
  full_join(studies, obs, by = var) %>%
    mutate(pct_studies = round(100 * n_studies / n_studies_all), pct_obs = round(100 * n_obs / n_obs_all)) %>%
    rename(variable = all_of(var))
}
n_studies_all <- n_studies
n_obs_all <- n_obs

s5_sex <- s5_block(data, "sex") %>%
  mutate(variable = recode(variable, female = "Female", male = "Male", both = "Both (sex disaggregates unavailable)")) %>%
  arrange(match(variable, c("Male", "Female", "Both (sex disaggregates unavailable)")))

s5_year <- data %>%
  mutate(period = paste0(5 * floor(year / 5), "-", 5 * floor(year / 5) + 4)) %>%
  s5_block("period") %>%
  arrange(variable)

s5_region <- s5_block(mutate(obs_country, region = region_names[region]), "region")

s5_country <- obs_country %>%
  mutate(country = ifelse(iso3 %in% top5, country_name(iso3), "All other countries")) %>%
  s5_block("country") %>%
  arrange(variable == "All other countries", desc(n_studies))

s5_method_label <- c("UTT/Cohort" = "UTT/Cohort", "National HHS" = "National HHS", "Other HHS" = "Other HHS",
                     "ANC" = "Pregnant/postpartum women", "Intervention Trials" = "Intervention Trials",
                     "Non-representative - Other" = "Other studies", "No inclusion criteria" = "No inclusion criteria",
                     "FSW" = "FSW", "MSM" = "MSM")
s5_method <- data %>%
  mutate(method = s5_method_label[recruitment_type]) %>%
  s5_block("method") %>%
  arrange(match(variable, s5_method_label))

s5_age <- data %>%
  mutate(age = case_when(
    age_group %in% c("Y015_049", "Y015_059", "Y015_064", "Y015_079") ~ "All adults (15-49, 15-59, 15-64, 15-79)",
    age_group_end(age_group) - age_group_start(age_group) == 4 ~ age_group_label(age_group),
    TRUE ~ "Other age groups")) %>%
  s5_block("age") %>%
  mutate(n_studies = NA, pct_studies = NA) %>%
  arrange(!str_detect(variable, "^\\d+-\\d+$"), variable)

header <- function(x) tibble(variable = x)
table_s5 <- bind_rows(
  header("Sex"), s5_sex,
  header("Year of observation"), s5_year,
  header("Study region"), s5_region, s5_country,
  header("Study method"), header("Total population"), filter(s5_method, !variable %in% c("FSW", "MSM")),
  header("Key population"), filter(s5_method, variable %in% c("FSW", "MSM")),
  header("Age group"), s5_age
) %>%
  select(variable, n_studies, pct_studies, n_obs, pct_obs)

write_table(table_s5, "table_s5")

################################################################################
#' ## Figure 1: flowchart
#' Screening counts (titles, abstracts, full texts) come from Covidence. The
#' study counts for the review steps are from the code; the pptx box definitions
#' differ from the code's steps (flowchart_counts.csv), so both are reported.
################################################################################

fc <- function(step, col) flowchart[[col]][flowchart$step == step]

record("Figure 1", "Pooled systematic review dataset (studies)", fc("Pooled systematic review data", "studies"),
       note = "Pptx box counts studies before de-duplication across reviews; code counts after")
review_steps <- flowchart %>%
  filter(box == "review") %>%
  mutate(excluded_studies = lag(studies) - studies, excluded_obs = lag(observations) - observations)
pptx_steps <- c("Risk-based inclusion criteria", "PrEP", "BED assay", "Secondary analyses",
                "Sex/age/year stratified data available", "Unable to calculate incidence rate")
for (i in seq_len(nrow(review_steps))[-1]) {
  s <- review_steps$step[i]
  record("Figure 1", paste("Studies excluded:", s), review_steps$excluded_studies[i],
         note = if (s %in% pptx_steps) "Pptx step definitions differ from the code's" else
           "Step not a separate box in the pptx")
}
record("Figure 1", "Review studies retained", fc("Unable to calculate incidence rate", "studies"))
record("Figure 1", "PHIA surveys", fc("PHIA surveys", "studies"))
record("Figure 1", "PHIA survey years", paste(range(data$year[str_detect(data$doi, "PHIA$")]), collapse = "-"),
       note = "Data years are survey midpoint years")
for (d in cohort_dois) {
  yrs <- range(data$year[data$doi == d])
  record("Figure 1", paste("Cohort years:", cohort_labels[[d]]), paste(yrs, collapse = "-"),
         note = "Midpoint years of the observations")
}
record("Figure 1", "Analysis dataset studies", fc("Final data", "studies"))
record("Figure 1", "Analysis dataset observations", fc("Final data", "observations"))

################################################################################
#' ## Main regression (Abstract, Results, Table 1)
################################################################################

type_term <- function(type) paste0("recruitment_type", type)

record("Results", "Correlation of log empirical and log modelled incidence (zeros omitted)",
       with(read_fit("modelled_correlation.csv") %>% filter(zero_incidence == "omit_zero"),
            sprintf("%.2f (%.2f-%.2f)", r, lower, upper)))

record("Methods", "Prior for study slope SD per 5 years, from population cohorts (median, 95% CrI)",
       with(read_fit("cohort_slope_prior.csv"),
            sprintf("%.2f (%.2f, %.2f)", sd_slope_median, sd_slope_lower, sd_slope_upper)))
record("Abstract; Results", "Calendar time per five years", coef_txt(main, "std_year"))
record("Abstract; Results", "Population cohorts and UTT trials", coef_txt(main, type_term("UTT/Cohort")))
record("Abstract; Results", "Nationally-representative household surveys", coef_txt(main, type_term("National HHS")))
record("Results", "Other household surveys", coef_txt(main, type_term("Other HHS")))
record("Results", "Population-representative sources: % above modelled (text: '20% higher')",
       paste(fmt_pct(map_dbl(type_term(c("UTT/Cohort", "National HHS", "Other HHS")), ~ coef_value(main, .x)) - 1),
             collapse = ", "))
record("Abstract; Results", "Control arms of clinical trials", coef_txt(main, type_term("Intervention Trials")))
record("Abstract; Results", "Pregnant/post-partum women", coef_txt(main, type_term("ANC")))
record("Results", "Other studies with inclusion criteria", coef_txt(main, type_term("Non-representative - Other")))
record("Abstract; Results", "Female sex workers", coef_txt(main, type_term("FSW")))
record("Abstract; Results", "Men who have sex with men", coef_txt(main, type_term("MSM")))
record("Abstract; Results", "Men in 2014", coef_txt(main, "sex_int"))
record("Results", "Men in 2014: % below modelled (text: '15% lower')", fmt_pct(1 - coef_value(main, "sex_int")))
record("Abstract; Results", "Sex:year interaction per five years", coef_txt(main, "sex_int:std_year"))
record("Results", "Men: additional decline per year (text: '2.2%/year')",
       fmt_pct(1 - coef_value(main, "sex_int:std_year")^(1 / 5), 1))
record("Discussion", "Fold higher incidence in trial control arms, pregnant women (text: 'roughly 2.5-fold')",
       sprintf("%.1f, %.1f", coef_value(main, type_term("Intervention Trials")), coef_value(main, type_term("ANC"))))

#' Table 1
table_1_rows <- tribble(
  ~parameter, ~term,
  "Study recruitment type/population (N studies)", NA,
  "Matched modelled estimates", "reference",
  "Studies representative of the total adult population", NA,
  "Population cohort studies and UTT trials [control clusters]", "UTT/Cohort",
  "Nationally-representative household surveys", "National HHS",
  "Other household surveys", "Other HHS",
  "Studies not representative of the total adult population", NA,
  "Clinical trials, excluding UTT trials [control arms]", "Intervention Trials",
  "Pregnant/post-partum women", "ANC",
  "Other studies", "Non-representative - Other",
  "Convenience sample from total population", "No inclusion criteria",
  "Female sex workers", "FSW",
  "Men who have sex with men", "MSM",
  "Sex", NA,
  "Female", "reference",
  "Male", "sex_int",
  "Calendar time: Per 5 years [reference year = 2014]", "std_year",
  "Sex [male] : five-year interaction", "sex_int:std_year",
  "Random effect standard deviation", NA,
  "Study (IID)", "SD for study_idx",
  "Study slope per 5 years (IID)", "SD for slope_idx",
  "Observation (IID)", "SD for obs_iid",
  "Sex-specific time trend (RW2)", "SD for year2"
)

studies_by_type <- data %>% distinct(doi, recruitment_type) %>% count(recruitment_type, name = "n_studies")

#' Coefficient text for a table row: types are recruitment levels, otherwise a term name.
#' A type with no observations in the fit (convenience samples since 2015) has
#' only its prior, shown as "-".
term_estimate <- function(fits, term) {
  map_chr(term, function(t) {
    if (is.na(t)) return(NA_character_)
    if (t == "reference") return("1.00 (Reference)")
    if (t %in% recruitment_levels) t <- type_term(t)
    if (!t %in% fits$var) return("-")
    if (with(filter(fits, var == t), upper / lower > 1e6)) return("-")
    coef_txt(fits, t)
  })
}

table_1 <- table_1_rows %>%
  left_join(studies_by_type, by = c("term" = "recruitment_type")) %>%
  mutate(estimate = term_estimate(main, term),
         estimate = ifelse(parameter == "Random effect standard deviation", "Estimate (95% CrI)", estimate))

write_table(select(table_1, parameter, n_studies, estimate), "table_1")

for (t in c("UTT/Cohort", "National HHS", "Other HHS", "Intervention Trials", "ANC",
            "Non-representative - Other", "No inclusion criteria", "FSW", "MSM")) {
  record("Table 1", paste("N studies:", recruitment_labels[[t]]),
         studies_by_type$n_studies[studies_by_type$recruitment_type == t])
}
for (i in which(!is.na(table_1$term) & table_1$term != "reference")) {
  record("Table 1", table_1$parameter[i], table_1$estimate[i])
}

################################################################################
#' ## Sensitivity analyses (Table S6, S7, S12, Supplementary Figure S2)
################################################################################

table_s6_rows <- tribble(
  ~parameter, ~term,
  "Modelled estimate", "reference",
  "Cohort/UTT trials", "UTT/Cohort",
  "National HHS", "National HHS",
  "Other HHS", "Other HHS",
  "Clinical trials", "Intervention Trials",
  "Preg/pp. women", "ANC",
  "Other studies", "Non-representative - Other",
  "Convenience", "No inclusion criteria",
  "FSW", "FSW",
  "MSM", "MSM",
  "Sex", NA,
  "Female", "reference",
  "Male", "sex_int",
  "Five calendar years", "std_year",
  "Sex:five-year interaction", "sex_int:std_year"
)

s6_variants <- c("main", "no_phia", "spectrum_subnational", "spectrum_national", "since_2010", "since_2015", "nb")
table_s6 <- table_s6_rows
for (v in s6_variants) {
  table_s6[[v]] <- term_estimate(filter(coefs, variant == v), table_s6$term)
}
#' Since 2015 has no convenience-sample observations, shown as "-"
write_table(select(table_s6, -term), "table_s6")

s6_long <- table_s6 %>%
  filter(!is.na(term), term != "reference") %>%
  pivot_longer(all_of(s6_variants), names_to = "variant", values_to = "value")
pwalk(s6_long, function(parameter, variant, value, ...) {
  record("Table S6", paste0(variant, ": ", parameter), value)
})

record("Table S6", "NB size (1/overdispersion)",
       coef_txt(filter(coefs, variant == "nb"), "log size for the nbinomial observations (1/overdispersion)"))

#' Table S7: measurement method
measurement <- filter(coefs, variant == "measurement")
table_s7_rows <- table_1_rows %>%
  filter(!term %in% c("FSW", "MSM"), !str_detect(coalesce(term, ""), "^SD for"),
         parameter != "Random effect standard deviation") %>%
  mutate(parameter = recode(parameter, "Calendar time: Per 5 years [reference year = 2014]" =
                              "Per 5 calendar years [reference year = 2014]"))
table_s7_rows <- bind_rows(
  slice(table_s7_rows, 1:which(table_s7_rows$parameter == "Convenience sample from total population")),
  tribble(~parameter, ~term,
          "HIV incidence measurement method", NA,
          "Prospective cohorts", "reference",
          "Retrospective cohorts", "study2Retrospective Cohort",
          "Cross-sectional assay", "study2Cross-sectional incidence"),
  slice(table_s7_rows, (which(table_s7_rows$parameter == "Convenience sample from total population") + 1):n())
)
table_s7 <- mutate(table_s7_rows, estimate = term_estimate(measurement, term))
write_table(select(table_s7, parameter, estimate), "table_s7")
for (i in which(!is.na(table_s7$term) & table_s7$term != "reference")) {
  record("Table S7", table_s7$parameter[i], table_s7$estimate[i])
}
record("Results", "Cross-sectional vs prospective measurement",
       coef_txt(measurement, "study2Cross-sectional incidence"))
record("Results", "Retrospective vs prospective cohorts",
       coef_txt(measurement, "study2Retrospective Cohort"))

#' Table S12: dropping observations with regression-imputed person-years
denominators <- filter(coefs, variant == "denominators")
imputed <- filter(data, str_detect(pys_source, "^participants"))
record("Supplementary Table S12 title", "Observations dropped (imputed person-years)", nrow(imputed))
record("Supplementary Table S12 title", "Studies with dropped observations", n_distinct(imputed$doi))
record("Table S12", "Studies remaining", denominators$n_studies[1])

table_s12_rows <- table_1_rows %>%
  filter(!str_detect(coalesce(term, ""), "^SD for"), parameter != "Random effect standard deviation") %>%
  mutate(parameter = recode(parameter, "Female sex workers" = "FSW", "Men who have sex with men" = "MSM",
                            "Calendar time: Per 5 years [reference year = 2014]" =
                              "Per 5 calendar years [reference year = 2014]"))
table_s12 <- mutate(table_s12_rows, estimate = term_estimate(denominators, term))
write_table(select(table_s12, parameter, estimate), "table_s12")
for (i in which(!is.na(table_s12$term) & table_s12$term != "reference")) {
  record("Table S12", table_s12$parameter[i], table_s12$estimate[i])
}

#' Supplementary Figure S2: leave-one-out ranges for each recruitment type, over refits leaving out a study of that type
loo <- read_fit("loo_coefficients.csv")
loo_ranges <- loo %>%
  filter(var == recruitment_type) %>%
  group_by(recruitment_type) %>%
  summarise(n = n(), min = min(mean), max = max(mean))
pwalk(loo_ranges, function(recruitment_type, n, min, max) {
  record("Supplementary Figure S2", paste("Leave-one-out range:", recruitment_type),
         sprintf("%.2f-%.2f (%d refits)", min, max, n))
})

################################################################################
#' ## Text S6: cohort incidence relative to modelled incidence
################################################################################

cohort_vs_model <- read_fit("coefficients_cohort_vs_model.csv")
#' The AR1 correlation is reported as rho, from exp(internal) E: rho = (E - 1) / (E + 1)
rho <- filter(cohort_vs_model, var == "Rho_intern for id.year") %>%
  mutate(across(c(mean, lower, upper), ~ (.x - 1) / (.x + 1)))

text_s6 <- tribble(
  ~parameter, ~estimate,
  "Intercept", coef_txt(cohort_vs_model, "(Intercept)"),
  "Five-calendar years (Reference year = 2014)", coef_txt(cohort_vs_model, "std_year"),
  "Sex", NA,
  "Female", "1.00 (Reference)",
  "Male", coef_txt(cohort_vs_model, "sexmale"),
  "Sex:five-year interaction", coef_txt(cohort_vs_model, "std_year:sexmale"),
  "Hyperparameters", "SD (95% CrI)",
  "Year (AR1)", coef_txt(cohort_vs_model, "SD for id.year"),
  NA, "Autocorrelation (95% CrI)",
  "Year", fmt_ci(rho$mean, rho$lower, rho$upper)
)
write_table(text_s6, "text_s6")
for (p in c("Intercept", "Five-calendar years (Reference year = 2014)", "Male", "Sex:five-year interaction",
            "Year (AR1)", "Year")) {
  record("Text S6", p, text_s6$estimate[which(text_s6$parameter == p)])
}
record("Results", "Cohorts only: five calendar years", coef_txt(cohort_vs_model, "std_year"))
record("Results", "Cohorts only: men", coef_txt(cohort_vs_model, "sexmale"))
record("Results", "Cohorts only: sex:five-year interaction", coef_txt(cohort_vs_model, "std_year:sexmale"))

################################################################################
#' ## Cohort incidence declines (Abstract, Results, Discussion, Table S8)
################################################################################

declines <- read_fit("cohort_declines.csv")
#' The text quotes the four cohorts with data to 2019 or later
text_cohorts <- c("AHRI", "RCCS", "Kyamulibwa GPC", "Kisumu HDSS")
for (d in text_cohorts) {
  x <- filter(declines, doi == d)
  record("Results", paste("Decline in 15-49 incidence:", cohort_short[[d]]),
         sprintf("%.0f%% (%.0f-%.0f%%; %d-%d)", 100 * x$median, 100 * x$lower, 100 * x$upper, x$base_year, x$latest_year))
}
record("Abstract; Discussion", "Median decline since 2010 across cohorts",
       fmt_pct(median(declines$median[declines$doi %in% text_cohorts])),
       note = "Median of the four cohorts quoted in the Results (Kisesa and Manicaland end before 2016)")

#' Table S8
cohort_trends <- read_fit("cohort_trends.csv")
table_s8 <- cohort_trends %>%
  transmute(cohort = cohort_short[doi], year,
            seroconversions = round(infections), person_years = round(pys),
            incidence = round(100 * ir, 1), crude_irr = round(crude_irr, 2),
            smoothed_incidence = round(100 * smoothed, 2), smoothed_irr = round(smoothed_irr, 2)) %>%
  mutate(across(-cohort, as.character))
write_table(table_s8, "table_s8")

################################################################################
#' ## Age and sex patterns: cohorts (Table S9, S11) and PHIA (Table S10)
################################################################################

cohort_age <- read_fit("coefficients_cohort_age.csv")
cohort_trend <- read_fit("cohort_age_trends.csv")
phia_age <- read_fit("coefficients_phia_age.csv")
phia_trend <- read_fit("phia_age_trends.csv")

trend_of <- function(trends, s, a) filter(trends, sex == s, age_group == a)

#' Abstract: "declined 4% faster per year (95%CrI 1-6%) among young adults 15-24 years than 25+ years"
#' = 1 - 1/exp(year:age 25-49), with 15-24 as the reference age group
age_int <- filter(cohort_age, var == "year_min:age_group2Y025_049")
record("Abstract", "Faster annual decline among 15-24 than 25+ (cohorts)",
       sprintf("%.0f%% (%.0f-%.0f%%)", 100 * (1 - 1 / age_int$mean), 100 * (1 - 1 / age_int$lower), 100 * (1 - 1 / age_int$upper)))

f1524 <- trend_of(cohort_trend, "female", "Y015_024")
f2549 <- trend_of(cohort_trend, "female", "Y025_049")
record("Results", "Cohorts: annual decline among women 15-24", pct_decline(f1524$median, f1524$lower, f1524$upper))
record("Results", "Cohorts: annual IRR among women 15-24 (year term)", coef_txt(cohort_age, "year_min"))
record("Results", "Cohorts: annual decline among women 25-49", pct_decline(f2549$median, f2549$lower, f2549$upper))
record("Results", "Cohorts: year:age 25-49 interaction", coef_txt(cohort_age, "year_min:age_group2Y025_049"))
record("Results", "Cohorts: year:sex interaction", coef_txt(cohort_age, "year_min:sexmale"))

table_s9_rows <- tribble(
  ~parameter, ~term,
  "Year (Reference year = 2014)", "year_min",
  "Age group", NA,
  "15-24", "reference",
  "25+", "age_group2Y025_049",
  "Gender", NA,
  "Female", "reference",
  "Male", "sexmale",
  "Sex : year interaction", NA,
  "Male : year", "year_min:sexmale",
  "Year : age interaction", NA,
  "Year : 25+", "year_min:age_group2Y025_049",
  "Sex : age interaction", NA,
  "Male: 25+", "sexmale:age_group2Y025_049",
  "Hyperparameters", "header",
  #' Labels follow the model terms ; AR2 SDs are marginal SDs
  "Age-specific cohort intercepts (IID)", "SD for doi3",
  "Age-specific time trend (AR2)", "SD for year_min2",
  "Cohort-specific time trend, women (AR2)", "SD for year_min_f",
  "Cohort-specific time trend, men (AR2)", "SD for year_min_m",
  "Cohort-specific age : time slope (IID)", "SD for slope_id"
)
table_s9 <- table_s9_rows %>%
  mutate(estimate = ifelse(term %in% "header", "SD (95% CrI)", term_estimate(cohort_age, na_if(term, "header"))))
#' The cohort slope SD is per year, so it needs three decimals
slope_sd <- filter(cohort_age, var == "SD for slope_id")
table_s9$estimate[table_s9$term %in% "SD for slope_id"] <- fmt_ci(slope_sd$mean, slope_sd$lower, slope_sd$upper, d = 3)
write_table(select(table_s9, parameter, estimate), "table_s9")
for (i in which(!is.na(table_s9$term) & !table_s9$term %in% c("reference", "header"))) {
  record("Table S9", table_s9$parameter[i], table_s9$estimate[i])
}

#' PHIA (Table S10)
p1524 <- trend_of(phia_trend, "female", "Y015_024")
p2549 <- trend_of(phia_trend, "female", "Y025_049")
record("Results; Discussion", "PHIA: annual decline among women 15-24", pct_decline(p1524$median, p1524$lower, p1524$upper))
record("Results", "PHIA: annual decline among women 25-49", pct_decline(p2549$median, p2549$lower, p2549$upper))
record("Results", "PHIA: year:age 25-49 interaction", coef_txt(phia_age, "year_min:age_groupY025_049"))
record("Results", "PHIA: year:sex interaction", coef_txt(phia_age, "year_min:sexmale"))

table_s10 <- table_s9_rows %>%
  filter(!str_detect(parameter, regex("cohort|age-specific", ignore_case = TRUE))) %>%
  mutate(term = str_replace(term, "age_group2", "age_group")) %>%
  bind_rows(tribble(~parameter, ~term,
                    "Country-specific intercepts (IID)", "SD for iso32",
                    "Country-specific time trend (IID)", "SD for iso3")) %>%
  mutate(estimate = ifelse(term %in% "header", "SD (95% CrI)", term_estimate(phia_age, na_if(term, "header"))))
write_table(select(table_s10, parameter, estimate), "table_s10")
for (i in which(!is.na(table_s10$term) & !table_s10$term %in% c("reference", "header"))) {
  record("Table S10", table_s10$parameter[i], table_s10$estimate[i])
}

#' Shares of new infections in 2022 (cohort district populations x cohort model incidence)
shares <- read_fit("cohort_infection_shares.csv")
shares_2022 <- shares %>%
  filter(year == 2022) %>%
  group_by(doi) %>%
  summarise(u25 = sum(share[age_group == "Y015_024"]), female = sum(share[sex == "female"]))
record("Results", "2022 share of new infections among ages 15-24, median across cohorts (IQR)",
       med_iqr(shares_2022$u25, scale = 100) %>% str_replace("\\)$", "%)") %>% str_replace("^(\\d+)", "\\1%"))
record("Results", "2022 share of new infections among women, median across cohorts (IQR)",
       med_iqr(shares_2022$female, scale = 100) %>% str_replace("\\)$", "%)") %>% str_replace("^(\\d+)", "\\1%"))

latest_shares <- shares %>%
  filter(in_observed_range) %>%
  group_by(doi) %>%
  filter(year == max(year)) %>%
  summarise(year = first(year), female_25 = sum(share[sex == "female" & age_group == "Y025_049"]))
record("Discussion", "Largest share of new infections among women 25-49 in a cohort's latest year (text: 'up to half')",
       fmt_pct(max(latest_shares$female_25)))

age_2022 <- read_fit("cohort_age_predictions.csv") %>%
  filter(year == 2022, sex == "female") %>%
  select(doi, age_group, median) %>%
  pivot_wider(names_from = age_group, values_from = median)
record("Discussion", "Cohorts where 2022 incidence among women 15-24 is below women 25-49",
       sprintf("%d of %d", sum(age_2022$Y015_024 < age_2022$Y025_049), nrow(age_2022)))

#' Crude sex ratio of new infections in the latest data (Results, Supplementary Figure S3B)
latest_sex <- data %>%
  filter(doi %in% c("Kyamulibwa GPC", "Kisumu HDSS"), sex != "both") %>%
  group_by(doi) %>%
  filter(year >= max(year) - 2) %>%
  group_by(doi, sex) %>%
  summarise(infections = sum(infections), years = paste(range(year), collapse = "-"), .groups = "drop") %>%
  pivot_wider(names_from = sex, values_from = infections)
kya <- filter(latest_sex, doi == "Kyamulibwa GPC")
record("Results", "Kyamulibwa: female:male new infections, latest three years of data",
       sprintf("%.1f:1 (%s)", kya$female / kya$male, kya$years),
       note = "Text reads the latest point of Fig S3B; crude ratio over the last three data years here")
siaya <- filter(latest_sex, doi == "Kisumu HDSS")
record("Results", "Siaya: male infections, latest data", sprintf("%.0f (%s)", siaya$male, siaya$years))

#' Table S11: HIV-negative population of the cohort districts
cohort_population <- read_fit("cohort_population.csv")
table_s11 <- cohort_population %>%
  transmute(cohort = cohort_short[doi], year = as.character(year), sex,
            age_group = recode(age_group, Y015_024 = "15-24", Y025_049 = "25+"),
            population = as.character(round(population)), prevalence = sprintf("%.1f%%", 100 * prevalence),
            hiv_negative = as.character(round(hiv_negative)))
write_table(table_s11, "table_s11")

################################################################################
#' ## PHIA vs national modelled incidence (Results, Figure 4)
#' Namibia's modelled value is from a non-public Spectrum file: only the summaries
#' are recorded here.
################################################################################

phia_cor <- read_fit("phia_spectrum_correlation.csv")
pc <- function(p, col) phia_cor[[col]][phia_cor$phase == p]
record("Results", "PHIA vs modelled: correlation of log incidence, all surveys", sprintf("%.2f", pc("all", "r")))
record("Results", "PHIA vs modelled: correlation, first round", sprintf("%.2f [n=%d]", pc("PHIA1", "r"), pc("PHIA1", "n_surveys")))
record("Results", "PHIA vs modelled: correlation, second round", sprintf("%.2f [n=%d]", pc("PHIA2", "r"), pc("PHIA2", "n_surveys")))
record("Results", "PHIA surveys whose 95% CI contains the modelled estimate",
       sprintf("%d/%d", pc("all", "n_modelled_in_ci"), pc("all", "n_surveys")))
phia_spectrum <- read_fit("phia_spectrum.csv")
second_lower <- phia_spectrum %>%
  group_by(iso3) %>%
  filter(n() == 2, year == max(year), !modelled_in_ci, modelled < lower) %>%
  pull(iso3)
record("Results", "Countries with two PHIAs where modelled < second-survey CI",
       if (length(second_lower)) paste(sort(second_lower), collapse = ", ") else "none")
record("Methods", "PHIA surveys; countries",
       sprintf("%d; %d", n_distinct(phia_spectrum$survey_id), n_distinct(phia_spectrum$iso3)))
record("Results", "Countries with two PHIA surveys", sum(table(phia_spectrum$iso3) == 2))

################################################################################
#' ## Numbers from outside this analysis
################################################################################

record("Methods", "Spectrum files with a constant age pattern of new infections", "58/74",
       note = "From the 2026 Spectrum file metadata, not this code")
record("Text S3", "Unique geocoding search strings; with >1 location; same district; resolved manually",
       "125; 17; 11; 6", note = "Geocoding cache summary (private/matched_gps.rds); not recomputed")

#' Discussion: 15-24 incidence decline 2010-2020 in the 2026 Spectrum files, all
#' sub-Saharan African national files pooled (infections / HIV-negative population),
#' with its annual equivalent, to compare with the cohort and PHIA annual declines.
#' Namibia's file is private: it is included when private/ is present, otherwise
#' left out (the pooled figure moves by <0.5 percentage points).
spectrum_file <- file.path(path_private, "spectrum_incidence_full.csv.gz")
if (!file.exists(spectrum_file)) spectrum_file <- file.path(path_modelled, "spectrum_incidence.csv.gz")
ssa_1524 <- read_csv(spectrum_file, show_col_types = FALSE) %>%
  filter(variant == "published", area_id == iso3, sex == "both", age_group == "Y015_024",
         year %in% c(2010, 2020), !is.na(infections)) %>%
  group_by(year) %>%
  summarise(incidence = sum(infections) / sum(population - plhiv), n_iso3 = n_distinct(iso3))
ssa_ratio <- ssa_1524$incidence[ssa_1524$year == 2020] / ssa_1524$incidence[ssa_1524$year == 2010]
record("Discussion", "Spectrum 2026: 15-24 incidence decline 2010-2020, sub-Saharan Africa; per year",
       sprintf("%.0f%%; %.0f%%/year", 100 * (1 - ssa_ratio), 100 * (1 - ssa_ratio^(1 / 10))),
       note = sprintf("%d national Spectrum files (%s)", max(ssa_1524$n_iso3), basename(spectrum_file)))
record("Discussion", "Median modelled decline across countries 2010-2023", "75%",
       note = "Cited UNAIDS estimates (refs 40, 51)")
record("Introduction", "UNAIDS incidence decline 2010-2024", "nearly 70%", note = "Cited (refs 1, 2)")
record("Results", "Share of new infections among 15-24 in modelled estimates", "40-50%",
       note = "Read from Supplementary Figure S3A")
record("Results", "Modelled female:male sex ratio of new infections 1995-2008", "1.5:1",
       note = "Read from Supplementary Figure S3B")

################################################################################
#' ## Write
################################################################################

write_csv(paper_numbers, "outputs/paper_numbers.csv", na = "")
message(sprintf("\n%d numbers recorded", nrow(paper_numbers)))
