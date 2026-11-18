################################################################################
#' Limited risk-of-bias assessment: participation and attrition (Supplementary Text S9)
#'
#' Reads the participation/attrition extractions in private/reference/ (Supplementary
#' File S2; copied by prep/copy_private_inputs.R):
#'   hiv_study_extraction.xlsx                   systematic review studies and cohorts
#'   kp_participation_attrition_extraction.xlsx  FSW and MSM review studies
#'   PHIA_blood_draw_per_eligible.xlsx           PHIA surveys (blood draws / eligible)
#' and restricts them to studies in the analysis dataset.
#'
#' Writes:
#'   outputs/fits/risk_of_bias.csv               one row per analysis study
#'   outputs/fits/risk_of_bias_correlation.csv   Pearson R of participation/attrition with
#'                                               log10 median IRR and log10 person-years, by group
#'   outputs/figures/figure_text_s9.png          Text S9 figure
#'
#' Rscript R/05_risk_of_bias.R [analysis_data.csv]
################################################################################

source("R/functions.R")

ref <- file.path("private", "reference")
args <- commandArgs(trailingOnly = TRUE)
data_file <- if (length(args)) args[1] else file.path(path_derived, "analysis_data.csv")
out_suffix <- if (length(args) > 1) args[2] else ""

data <- read_csv(data_file, show_col_types = FALSE)

#' Population-based cohorts from demographic surveillance sites: participation and
#' attrition are defined differently from other prospective studies (Text S9)
dss_cohorts <- c("RCCS", "AHRI", "Manicaland", "Kyamulibwa GPC", "Kisesa", "Karonga", "Kisumu HDSS")

as_rate <- function(x) suppressWarnings(as.numeric(x))  # text ("16-48%", "N/A") becomes NA

review_rob <- readxl::read_excel(file.path(ref, "hiv_study_extraction.xlsx")) %>%
  transmute(doi = recode(doi, Siaya = "Kisumu HDSS"),
            participation = as_rate(as.character(participation)),
            attrition = as_rate(as.character(attrition)))

phia_rob <- readxl::read_excel(file.path(ref, "PHIA_blood_draw_per_eligible.xlsx")) %>%
  transmute(doi = survey_id, participation = `Blood/eligible TOTAL`, attrition = NA_real_)

#' Key population extraction: match each row to the analysis study by review and
#' author surname (and year where the review has two studies by that author)
kp_ids <- data %>%
  filter(str_detect(doi, "^(Stannah|Jones & Anderson) \\| ")) %>%
  distinct(doi) %>%
  mutate(review = ifelse(str_detect(doi, "^Stannah"), "MSM", "FSW"),
         name = str_remove(doi, "^.* \\| "))

kp_rob <- readxl::read_excel(file.path(ref, "kp_participation_attrition_extraction.xlsx")) %>%
  filter(!str_detect(`Design / relevance`, "^See ")) %>%
  mutate(surname = str_extract(Study, "^[A-Za-z]+( Simaleko| Damme)?"),
         surname = recode(surname, Naicker = "Naiker"),
         year = str_extract(Study, "\\d{4}")) %>%
  rowwise() %>%
  mutate(doi = {
    cand <- kp_ids %>% filter(review == Review, str_detect(name, fixed(surname)))
    # where both give a year, it must agree (e.g. Wahome 2018 is not Wahome 2020)
    if (!is.na(year)) cand <- filter(cand, !str_detect(name, "\\d{4}") | str_detect(name, year))
    if (nrow(cand) == 1) cand$doi else NA_character_
  }) %>%
  ungroup() %>%
  transmute(doi,
            # cross-sectional surveys and routine programme data have no follow-up
            kp_cross_sectional = str_detect(`Design / relevance`, "(?i)cross-sectional|routine programme") &
              !str_detect(`Design / relevance`, "(?i)longitudinal|nested cohort"),
            participation = `Consented / enrolled (n)` / `Eligible (n)`,
            attrition = as_rate(as.character(`Attrition / LTFU rate`))) %>%
  filter(!is.na(doi)) %>%
  #' Wahome 2020 (doi:10.1016/j.eclinm.2020.100541) does not report a study participation
  #' rate: 179 is the number of registered cohort MSM and 172 the number assessed for PrEP
  #' eligibility, so 179/172 (104%) is not a participation rate. Recorded as not reported.
  mutate(participation = ifelse(doi == "Stannah | Wahome 2020", NA_real_, participation))

stopifnot(!anyDuplicated(kp_rob$doi))

studies <- data %>%
  mutate(irr = ir / extrapolated_value) %>%
  group_by(doi, recruitment_type, study_type) %>%
  summarise(median_irr = DescTools::Quantile(irr, weights = pys, probs = 0.5),
            pys = sum(pys), .groups = "drop") %>%
  group_by(doi) %>%
  slice_max(pys, n = 1, with_ties = FALSE) %>%
  ungroup()
stopifnot(!anyDuplicated(studies$doi))

rob <- bind_rows(review_rob, phia_rob, kp_rob) %>%
  filter(!is.na(participation) | !is.na(attrition) | !is.na(kp_cross_sectional)) %>%
  distinct(doi, .keep_all = TRUE) %>%
  right_join(studies, by = "doi") %>%
  mutate(longitudinal = study_type != "Cross-sectional incidence" & !coalesce(kp_cross_sectional, FALSE),
         # attrition is only defined for studies with follow-up
         attrition = ifelse(longitudinal, attrition, NA),
         group = case_when(recruitment_type %in% c("FSW", "MSM") ~ "Key population",
                           doi %in% dss_cohorts ~ "Population-representative cohorts",
                           TRUE ~ "Remaining studies") %>%
           factor(levels = c("Remaining studies", "Population-representative cohorts", "Key population"))) %>%
  arrange(group, doi)

write_csv(rob, file.path(path_fits, paste0("risk_of_bias", out_suffix, ".csv")))

#' Pearson correlation of a rate with log10(x), as ggpmisc::stat_correlation on a log10 axis
correlation <- rob %>%
  pivot_longer(c(participation, attrition), names_to = "measure", values_to = "rate") %>%
  pivot_longer(c(median_irr, pys), names_to = "against", values_to = "x") %>%
  filter(!is.na(rate), x > 0) %>%
  group_by(measure, against, group) %>%
  filter(n() >= 4) %>%
  summarise(n = n(),
            r = cor(log10(x), rate),
            lower = cor.test(log10(x), rate)$conf.int[1],
            upper = cor.test(log10(x), rate)$conf.int[2],
            .groups = "drop")

write_csv(correlation, file.path(path_fits, paste0("risk_of_bias_correlation", out_suffix, ".csv")))

################################################################################
#' Figure (Text S9)
################################################################################

group_colours <- c("Remaining studies" = "#619CFF", "Population-representative cohorts" = "#00BA38",
                   "Key population" = "#F8766D")

rob_panel <- function(measure, against) {
  xlab <- c(median_irr = "Median incidence rate ratio (log scale)", pys = "Observed person-years (log scale)")[[against]]
  d <- rob %>% filter(!is.na(.data[[measure]]), .data[[against]] > 0)
  labels <- correlation %>%
    filter(measure == !!measure, against == !!against) %>%
    arrange(group) %>%
    mutate(label = sprintf("R = %.2f, 95%% CI [%.2f, %.2f]", r, lower, upper),
           row = row_number())
  y_top <- measure == "attrition"
  ggplot(d, aes(x = .data[[against]], y = .data[[measure]], color = group)) +
    geom_point() +
    geom_text(data = labels, aes(label = label, color = group,
                                 y = if (y_top) 0.99 - 0.08 * (row - 1) else 0.03 + 0.08 * (nrow(labels) - row)),
              x = -Inf, hjust = -0.04, size = 6.2, inherit.aes = FALSE, show.legend = FALSE) +
    scale_color_manual(values = group_colours, drop = FALSE) +
    scale_y_continuous(labels = scales::label_percent(), limits = c(0, 1.05), breaks = seq(0, 1, 0.25)) +
    scale_x_log10(labels = scales::label_number(big.mark = ",", drop0trailing = TRUE)) +
    expand_limits(x = if (against == "pys") c(50, 1e6) else c(0.1, 1500)) +
    labs(x = xlab, y = str_to_sentence(measure), color = NULL) +
    theme_paper() +
    theme_large() +
    theme(aspect.ratio = 1, plot.margin = margin(5, 30, 5, 5))
}

library(patchwork)
fig <- (rob_panel("participation", "median_irr") | rob_panel("participation", "pys")) /
  (rob_panel("attrition", "median_irr") | rob_panel("attrition", "pys")) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

if (out_suffix == "") save_png(fig, "figure_text_s9.png", width = 1100, height = 1023)

################################################################################
#' Text S9 numbers
################################################################################

q <- function(x) quantile(x, c(0.5, 0.25, 0.75), na.rm = TRUE)
message(sprintf("%d of %d studies reported a participation rate and %d of %d longitudinal studies reported an attrition rate",
                sum(!is.na(rob$participation)), nrow(rob), sum(!is.na(rob$attrition)), sum(rob$longitudinal)))
message(sprintf("participation median %.0f%% (IQR %.0f-%.0f%%); attrition %.0f%% (IQR %.0f-%.0f%%)",
                100 * q(rob$participation)[1], 100 * q(rob$participation)[2], 100 * q(rob$participation)[3],
                100 * q(rob$attrition)[1], 100 * q(rob$attrition)[2], 100 * q(rob$attrition)[3]))
print(as.data.frame(correlation), digits = 2)
