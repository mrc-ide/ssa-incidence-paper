################################################################################
#' Shared helpers: age groups, regions, INLA summaries and plot styling
################################################################################

library(tidyverse)

#' Paths, relative to the repository root
path_raw <- "data/raw"
path_geography <- "data/geography"
path_modelled <- "data/modelled"
path_private <- "private"
path_derived <- "outputs/data"
path_fits <- "outputs/fits"
path_figures <- "outputs/figures"
path_tables <- "outputs/tables"

#' Reference year for all regressions: the median year of observation
reference_year <- 2014L

#' Recruitment types, in regression order; the first is the reference... except
#' that the model has no intercept, so every level gets a coefficient relative to
#' the matched modelled estimate.
recruitment_levels <- c("UTT/Cohort", "National HHS", "Other HHS", "ANC", "Intervention Trials",
                        "Non-representative - Other", "No inclusion criteria", "FSW", "MSM")

recruitment_labels <- c(
  "UTT/Cohort" = "Population cohort studies and UTT trials [control clusters]",
  "National HHS" = "Nationally-representative household surveys",
  "Other HHS" = "Other household surveys",
  "Intervention Trials" = "Clinical trials, excluding UTT trials [control arms]",
  "ANC" = "Pregnant/post-partum women",
  "Non-representative - Other" = "Other studies",
  "No inclusion criteria" = "Convenience sample from total population",
  "FSW" = "Female sex workers",
  "MSM" = "Men who have sex with men"
)

#' Population-representative cohorts (study identifiers as in the analysis data)
cohort_dois <- c("AHRI", "Karonga", "Kisesa", "Kisumu HDSS", "Kyamulibwa GPC", "Manicaland", "RCCS")

cohort_labels <- c(
  "AHRI" = "uMkhanyakude, South Africa",
  "Kisumu HDSS" = "Siaya, Kenya",
  "Karonga" = "Karonga, Malawi",
  "Kisesa" = "Kisesa, Tanzania",
  "Kyamulibwa GPC" = "Kyamulibwa, Uganda",
  "Manicaland" = "Manicaland, Zimbabwe",
  "RCCS" = "Rakai, Uganda"
)

cohort_short <- c("AHRI" = "uMkhanyakude", "Kisumu HDSS" = "Siaya", "Karonga" = "Karonga", "Kisesa" = "Kisesa",
                  "Kyamulibwa GPC" = "Kyamulibwa", "Manicaland" = "Manicaland", "RCCS" = "Rakai")

################################################################################
#' ## Age groups
#' Codes follow the Naomi convention: "Y015_019" is ages 15-19, "Y080_999" 80+.
################################################################################

age_group_code <- function(start, end) {
  sprintf("Y%03d_%03d", as.integer(start), as.integer(end))
}

age_group_start <- function(age_group) as.integer(substr(age_group, 2, 4))
age_group_end <- function(age_group) as.integer(substr(age_group, 6, 8))

#' Five-year age groups 0-4 to 75-79, and 80+
five_year_age_groups <- c(age_group_code(seq(0, 75, 5), seq(4, 79, 5)), "Y080_999")

#' Five-year age groups within 15-49
age_groups_15to49 <- age_group_code(seq(15, 45, 5), seq(19, 49, 5))

#' Convert a label such as "15-19" or "80+" to a five-year age group code; other
#' spans (e.g. "18-19", "35-44") return NA, as they did when matched against the
#' Naomi age group table.
age_label_to_group <- function(label) {
  out <- rep(NA_character_, length(label))
  m <- str_match(label, "^(\\d+)-(\\d+)$")
  five <- !is.na(m[, 1]) & as.integer(m[, 3]) - as.integer(m[, 2]) == 4 & as.integer(m[, 2]) %% 5 == 0
  out[five] <- age_group_code(m[five, 2], m[five, 3])
  out[label %in% "80+"] <- "Y080_999"
  out
}

age_group_label <- function(age_group) {
  ifelse(age_group_end(age_group) == 999,
         paste0(age_group_start(age_group), "+"),
         paste0(age_group_start(age_group), "-", age_group_end(age_group)))
}

#' Five-year age groups (15-19 to 75-79) wholly contained in each requested age group
five_year_groups_within <- function(age_groups) {
  five <- tibble(age_group = age_group_code(seq(15, 75, 5), seq(19, 79, 5))) %>%
    mutate(start = age_group_start(age_group), end = age_group_end(age_group))
  tibble(age_group_data = unique(age_groups)) %>%
    mutate(data_start = age_group_start(age_group_data), data_end = age_group_end(age_group_data)) %>%
    cross_join(five) %>%
    filter(start >= data_start, end <= data_end) %>%
    select(age_group, age_group_data)
}

#' Add both-sexes rows by summing count indicators over sex
sex_aggregation <- function(df, indicators) {
  both <- df %>%
    group_by(across(-all_of(c(indicators, "sex")))) %>%
    summarise(across(all_of(indicators), sum), .groups = "drop") %>%
    mutate(sex = "both")
  bind_rows(df, both) %>% ungroup()
}

################################################################################
#' ## Regions and names
################################################################################

region_lookup <- tribble(
  ~iso3, ~region, ~four_region,
  "AGO", "ESA", "Southern", "BDI", "ESA", "Eastern", "BEN", "WCA", "Western", "BFA", "WCA", "Western",
  "CIV", "WCA", "Western", "CMR", "WCA", "Central", "COD", "WCA", "Central", "COG", "WCA", "Central",
  "ETH", "ESA", "Eastern", "ERI", "ESA", "Eastern", "GAB", "WCA", "Central", "GIN", "WCA", "Western",
  "GMB", "WCA", "Western", "KEN", "ESA", "Eastern", "LBR", "WCA", "Western", "LSO", "ESA", "Southern",
  "MLI", "WCA", "Western", "MOZ", "ESA", "Southern", "MWI", "ESA", "Southern", "NAM", "ESA", "Southern",
  "NER", "WCA", "Western", "NGA", "WCA", "Western", "RWA", "ESA", "Eastern", "SWZ", "ESA", "Southern",
  "TCD", "WCA", "Central", "TGO", "WCA", "Western", "TZA", "ESA", "Eastern", "UGA", "ESA", "Eastern",
  "ZMB", "ESA", "Southern", "ZWE", "ESA", "Southern", "GHA", "WCA", "Western", "SEN", "WCA", "Western",
  "SLE", "WCA", "Western", "CAF", "WCA", "Central", "BWA", "ESA", "Southern", "ZAF", "ESA", "Southern",
  "GNB", "WCA", "Western", "GNQ", "WCA", "Central", "SSD", "ESA", "Eastern"
)

region_names <- c(ESA = "Eastern and Southern Africa", WCA = "Western and Central Africa")

country_name <- function(iso3) {
  countrycode::countrycode(iso3, "iso3c", "country.name",
                           custom_match = c(CIV = "Côte d'Ivoire", SWZ = "eSwatini", COD = "DR Congo"))
}

################################################################################
#' ## Formatting and summaries
################################################################################

#' "1.22 (0.65, 2.31)"
fmt_ci <- function(mean, lower, upper, d = 2) {
  sprintf(paste0("%.", d, "f (%.", d, "f, %.", d, "f)"), mean, lower, upper)
}

fmt_pct <- function(x, d = 0) sprintf(paste0("%.", d, "f%%"), 100 * x)

fmt_int <- function(x) format(round(x), big.mark = ",", trim = TRUE)

#' Median and IQR (optionally weighted), as "median (IQR lower-upper)"
med_iqr <- function(x, weights = NULL, d = 0, scale = 1) {
  q <- DescTools::Quantile(x, weights, c(0.5, 0.25, 0.75), names = FALSE, na.rm = TRUE) * scale
  sprintf(paste0("%.", d, "f (IQR %.", d, "f-%.", d, "f)"), q[1], q[2], q[3])
}

#' Exponentiated fixed effects and hyperparameters of an INLA fit. Hyperparameters
#' are reported as in the manuscript tables: exp() of the posterior summary of the
#' internal (log precision / transformed) scale.
summarise_inla <- function(fit) {
  fixed <- fit$summary.fixed %>%
    rownames_to_column("var") %>%
    select(var, mean, lower = `0.025quant`, upper = `0.975quant`) %>%
    mutate(across(-var, exp), type = "fixed")
  hyper <- fit$internal.summary.hyperpar %>%
    rownames_to_column("var") %>%
    select(var, mean, lower = `0.025quant`, upper = `0.975quant`) %>%
    mutate(across(-var, exp), type = "hyper")
  #' Standard deviations from log precisions: SD = exp(-theta / 2), so the
  #' interval limits swap; the mean is taken over the marginal
  log_prec <- grep("^Log precision for ", names(fit$internal.marginals.hyperpar), value = TRUE)
  sds <- map_dfr(log_prec, function(v) {
    m <- fit$internal.marginals.hyperpar[[v]]
    q <- INLA::inla.qmarginal(c(0.025, 0.975), m)
    tibble(var = sub("^Log precision for ", "SD for ", v),
           mean = INLA::inla.emarginal(function(x) exp(-x / 2), m),
           lower = exp(-q[2] / 2), upper = exp(-q[1] / 2), type = "sd")
  })
  bind_rows(fixed, hyper, sds) %>%
    mutate(txt = fmt_ci(mean, lower, upper))
}

#' Posterior samples (rows x samples matrix) of the linear predictor for data rows `rows`
predictor_samples <- function(fit, rows, n = 1000, seed = 1L) {
  samples <- suppressWarnings(INLA::inla.posterior.sample(n, fit, seed = seed))
  contents <- fit$misc$configs$contents
  start <- contents$start[contents$tag == "Predictor"]
  mat <- sapply(samples, function(x) x$latent[start - 1 + rows])
  if (length(rows) == 1) mat <- matrix(mat, nrow = 1)
  mat
}

#' Posterior samples (terms x samples matrix) of fixed effects `terms`
fixed_samples <- function(fit, terms, n = 1000, seed = 1L) {
  samples <- suppressWarnings(INLA::inla.posterior.sample(n, fit, seed = seed))
  mat <- sapply(samples, function(x) x$latent[paste0(terms, ":1"), 1])
  if (length(terms) == 1) mat <- matrix(mat, nrow = 1)
  rownames(mat) <- terms
  mat
}

#' Annual incidence rate ratio by sex and age group from a model with
#' year_min * sex + age * year_min terms, as posterior median and 95% interval.
#' `age_term` is the name of the 25-49 level's interaction with year_min.
age_sex_trends <- function(fit, age_term, n = 1000, seed = 1L) {
  terms <- c("year_min", "year_min:sexmale", age_term)
  b <- fixed_samples(fit, terms, n, seed)
  groups <- list(
    "female Y015_024" = b["year_min", ],
    "female Y025_049" = b["year_min", ] + b[age_term, ],
    "male Y015_024" = b["year_min", ] + b["year_min:sexmale", ],
    "male Y025_049" = b["year_min", ] + b["year_min:sexmale", ] + b[age_term, ]
  )
  imap_dfr(groups, function(x, g) {
    q <- quantile(exp(x), c(0.5, 0.025, 0.975))
    tibble(group = g, median = q[1], lower = q[2], upper = q[3])
  }) %>%
    separate(group, c("sex", "age_group"), sep = " ")
}

#' Posterior samples of the linear predictor for the prediction rows, which are
#' placed first in the data frame (the rows with no outcome).
#' Returns the prediction rows with mean and 95% interval on the log scale, and
#' (optionally) the sample matrix. A non-zero seed makes the samples reproducible
#' (INLA then samples serially).
sample_predictor <- function(fit, df, outcome = "infections", n = 1000, keep_samples = FALSE, seed = 1L) {
  n_pred <- sum(is.na(df[[outcome]]))
  stopifnot(all(is.na(df[[outcome]][seq_len(n_pred)])))
  idx <- seq_len(n_pred)
  mat <- predictor_samples(fit, idx, n, seed)
  qtls <- apply(mat, 1, quantile, c(0.025, 0.5, 0.975))
  out <- df[idx, ] %>%
    mutate(mean = rowMeans(mat), median = qtls[2, ], lower = qtls[1, ], upper = qtls[3, ])
  if (keep_samples) attr(out, "samples") <- mat
  out
}

################################################################################
#' ## Plot styling (as moz.utils::standard_theme)
################################################################################

theme_paper <- function() {
  theme_minimal() +
    theme(legend.position = "bottom",
          plot.title = element_text(size = 16),
          axis.text = element_text(size = 12),
          axis.title = element_text(size = 14, face = "bold"),
          legend.text = element_text(size = 12),
          strip.text = element_text(size = 13, face = "bold"),
          strip.background = element_rect(fill = NA, colour = "white"),
          plot.tag = element_text(size = 16, face = "bold"),
          panel.background = element_rect(fill = NA, color = "black"),
          panel.grid = element_blank())
}

pal_two <- wesanderson::wes_palette("Zissou1")[c(1, 4)]
pal_four <- c(wesanderson::wes_palette("Darjeeling2")[2], wesanderson::wes_palette("Darjeeling1")[c(2, 4, 5)])
pal_blue <- wesanderson::wes_palette("Zissou1")[1]

#' Pixel sizes and font sizes are as in the original figure code (72 dpi)
save_png <- function(plot, file, width, height, res = NA) {
  dir.create(path_figures, showWarnings = FALSE, recursive = TRUE)
  png(file.path(path_figures, file), width = width, height = height, res = res)
  print(plot)
  dev.off()
  invisible(file.path(path_figures, file))
}
