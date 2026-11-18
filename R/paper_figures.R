################################################################################
#' # Figures for the manuscript and supplement
#'
#' Reads only: outputs/data/ (R/01_data.R), outputs/fits/ (R/02-04) and the
#' modelled extracts in data/modelled/. Writes PNGs to outputs/figures/, with the
#' pixel sizes and font sizes of the original figure code. Figure 1 (the review
#' flowchart) is flowchart.pptx; its numbers are recorded in paper_numbers.R.
#'
#' Zenodo: outputs/data/ can be downloaded rather than regenerated (see
#' R/02_fit_main.R); figures then will not exactly reproduce the manuscript.
#'
#' Run from the repository root: Rscript R/paper_figures.R
################################################################################

source("R/functions.R")
library(patchwork)

data <- read_csv(file.path(path_derived, "analysis_data.csv"), show_col_types = FALSE)
read_fit <- function(file) read_csv(file.path(path_fits, file), show_col_types = FALSE)

spectrum <- read_csv(file.path(path_modelled, "spectrum_incidence.csv.gz"), show_col_types = FALSE)

#' Cohort labels as factors in the paper's panel order
cohort_factor <- function(doi, short = FALSE) {
  labels <- if (short) cohort_short else cohort_labels
  factor(labels[doi], levels = labels)
}

cohort <- filter(data, doi %in% cohort_dois)
cohort_area <- distinct(cohort, doi, area_id)

#' Cohort data by sex and age group 15-24, 25-49 (as R/04_fit_cohort_phia.R)
cohort_age_data <- cohort %>%
  filter(age_group %in% c(age_groups_15to49, "Y035_044", "Y045_054")) %>%
  mutate(age_group = ifelse(age_group %in% c("Y015_019", "Y020_024"), "Y015_024", "Y025_049")) %>%
  group_by(doi, area_id, year, sex, age_group) %>%
  summarise(infections = sum(ir * pys), pys = sum(pys), .groups = "drop") %>%
  mutate(ir = infections / pys)

larger_text <- function(base = 16) {
  theme(axis.text = element_text(size = base), legend.text = element_text(size = base),
        legend.title = element_text(size = base), axis.title = element_text(size = base + 2),
        strip.text = element_text(size = base + 1, face = "plain"), plot.tag = element_text(size = base + 2))
}

################################################################################
#' ## Figure 2: empirical:modelled incidence ratio by recruitment type, sex and year
################################################################################

type_panel_labels <- c("UTT/Cohort" = "Population cohort/\nUTT trials", "National HHS" = "National\nhousehold survey",
                       "Other HHS" = "Other\nhousehold survey", "ANC" = "Pregnant/\npostpartum women",
                       "Intervention Trials" = "Intervention trials", "Non-representative - Other" = "Other studies")

#' Regression line over the years with data, for each recruitment type and sex
predictions <- read_fit("main_predictions.csv") %>%
  filter(in_observed_range, recruitment_type != "No inclusion criteria") %>%
  mutate(across(c(mean, lower, upper), exp), sex = str_to_sentence(sex))

ratio_points <- data %>%
  filter(sex != "both") %>%
  mutate(irr = ir / extrapolated_value, sex = str_to_sentence(sex))

ratio_panel <- function(types, title, facets, y_breaks = c(1, 10), y_max = 16) {
  #' Zoom with coord_cartesian (not scale limits), so points and ribbons beyond y_max
  #' are cut at the panel edge rather than dropped
  #' Zero-incidence observations cannot be shown on the log scale
  pts <- filter(ratio_points, recruitment_type %in% types, irr > 0)
  y_lim <- if (is.na(y_max)) NULL else c(min(pts$irr), y_max)
  ggplot(pts, aes(x = year)) +
    geom_hline(yintercept = 1, linetype = 2) +
    geom_point(aes(y = irr, size = pys), shape = 16, alpha = 0.4, show.legend = FALSE) +
    geom_ribbon(data = filter(predictions, recruitment_type %in% types), aes(ymin = lower, ymax = upper),
                fill = pal_blue, alpha = 0.5) +
    geom_line(data = filter(predictions, recruitment_type %in% types), aes(y = mean), linewidth = 1.5, color = pal_blue) +
    facets +
    scale_x_continuous(breaks = seq(1990, 2020, 10)) +
    scale_y_continuous(transform = "log", labels = scales::label_number(accuracy = 1, big.mark = ""),
                       breaks = y_breaks) +
    coord_cartesian(xlim = c(1990, 2023), ylim = y_lim) +
    labs(x = NULL, y = "Incidence rate ratio (log scale)", title = title) +
    theme_paper() +
    larger_text() +
    theme(strip.background = element_rect(fill = NA, colour = NA),
          plot.title = element_text(size = 22, face = "bold"),
          plot.margin = margin(0, 2, 1, 0, unit = "cm"))
}

type_facet <- function(df) factor(type_panel_labels[df], levels = type_panel_labels)

#' Panels tagged A-C, with space above the B and C titles to separate them from
#' the x-axis of the panel above
fig2a <- ratio_panel(c("UTT/Cohort", "National HHS", "Other HHS"), "Population representative studies",
                     facet_grid(type_facet(recruitment_type) ~ sex)) +
  labs(tag = "A")
fig2b <- ratio_panel(c("ANC", "Intervention Trials", "Non-representative - Other"), "Non-representative studies",
                     ggh4x::facet_grid2(type_facet(recruitment_type) ~ sex, render_empty = FALSE)) +
  labs(tag = "B") +
  theme(plot.margin = margin(1, 2, 1, 0, unit = "cm"))
fig2c <- ratio_panel(c("FSW", "MSM"), "Key populations", facet_grid(" \n " ~ recruitment_type, labeller = labeller(recruitment_type = recruitment_labels)),
                     y_breaks = c(1, 50, 1000), y_max = NA) +
  labs(y = " ", tag = "C") +
  theme(plot.margin = margin(1, 2, 0, 0, unit = "cm"))

#' Tags drawn over the top-left corner, so they take no width from the panels
save_png(fig2a / fig2b / fig2c + plot_layout(heights = c(3, 3, 1)) &
           theme(plot.tag = element_text(size = 22, face = "bold"), plot.tag.position = c(0.01, 0.99)),
         "figure_2.png", width = 600, height = 1400)

################################################################################
#' ## Supplementary Figure S1: empirical vs matched modelled incidence
################################################################################

fig_s1 <- data %>%
  filter(sex != "both", !recruitment_type %in% c("FSW", "MSM")) %>%
  mutate(sex = str_to_sentence(sex)) %>%
  ggplot(aes(x = extrapolated_value, y = ir)) +
  geom_point(aes(size = pys)) +
  geom_abline(intercept = 0, slope = 1, linetype = 2) +
  scale_x_continuous(transform = "log", breaks = c(0.0001, 0.001, 0.01, 0.1),
                     labels = scales::label_number(scale = 100), limits = c(0.0001, 0.13)) +
  scale_y_continuous(transform = "log", breaks = c(0.0001, 0.001, 0.01, 0.1),
                     labels = scales::label_number(scale = 100), limits = c(0.0001, 0.13)) +
  facet_wrap(~sex) +
  labs(x = "Modelled incidence/100py (log scale)", y = "Empirical incidence/100py (log scale)", size = "Person-years") +
  theme_paper() +
  larger_text(15) +
  theme(legend.box = "vertical", aspect.ratio = 1)

save_png(fig_s1, "figure_s1.png", width = 700, height = 400)

################################################################################
#' ## Supplementary Figure S2: leave-one-study-out estimates
################################################################################

main_coefs <- read_fit("coefficients_main.csv") %>%
  filter(type == "fixed") %>%
  mutate(var = str_remove(var, "^recruitment_type"), doi = "00000")

loo <- read_fit("loo_coefficients.csv") %>%
  filter(var == recruitment_type | var %in% c("sex_int", "std_year", "sex_int:std_year")) %>%
  bind_rows(main_coefs) %>%
  mutate(main = doi == "00000")

loo_panel <- function(vars, labels = vars, breaks = scales::breaks_pretty(n = 3), limits = NULL) {
  loo %>%
    filter(var %in% vars) %>%
    mutate(var = factor(var, levels = vars, labels = labels)) %>%
    ggplot(aes(x = fct_rev(doi), y = mean, color = main)) +
    geom_point(show.legend = FALSE) +
    geom_linerange(aes(ymin = lower, ymax = upper), show.legend = FALSE) +
    geom_hline(yintercept = 1, linetype = 2, color = "red") +
    scale_y_continuous(transform = "log", labels = function(x) format(x, drop0trailing = TRUE, trim = TRUE),
                       breaks = breaks, limits = limits) +
    scale_color_manual(values = c("black", "blue")) +
    facet_wrap(~var, scales = "free", nrow = 1, labeller = label_wrap_gen(16)) +
    coord_flip() +
    labs(x = NULL, y = "IRR") +
    theme_paper() +
    theme_large() +
    # Extra space between panels so the end labels of adjacent free axes don't merge
    theme(axis.text.y = element_blank(), aspect.ratio = 1, panel.spacing.x = unit(2.5, "lines"))
}

fig_s2 <- loo_panel(c("sex_int", "std_year", "sex_int:std_year"), c("Sex", "Year", "Sex:Year")) /
  loo_panel(c("UTT/Cohort", "National HHS", "Other HHS"), breaks = c(0.6, 1, 1.5, 2.5), limits = c(0.5, 3)) /
  loo_panel(c("ANC", "Intervention Trials", "Non-representative - Other", "No inclusion criteria"), breaks = c(0.3, 1, 3, 10)) /
  loo_panel(c("FSW", "MSM"), breaks = c(1, 3, 10, 30, 100))

save_png(fig_s2, "figure_s2.png", width = 1100, height = 1500)

################################################################################
#' ## Figure 3: incidence trends in population cohorts (ages 15-49, both sexes)
################################################################################

cohort_trends <- read_fit("cohort_trends.csv") %>%
  mutate(cohort = cohort_factor(doi),
         #' Kisesa and Siaya have no 2010 round: their 2011 base is drawn at 2010
         year = ifelse(year == 2011 & doi %in% c("Kisesa", "Kisumu HDSS"), 2010, year),
         se = sqrt(ir * (1 - ir) / pys), lower = ir - 1.96 * se, upper = ir + 1.96 * se)

#' Cohort names, placed by hand in free space next to each line so that none
#' overlaps another label, a point or a smoothed line (the lines converge on 2010,
#' so end-of-line labels collide). Karonga (two rounds, no line) is labelled below
#' its 2006 point with a leader. Positions are on the data scale; check them by eye
#' if the data or figure size change.
fig3a_labels <- tribble(
  ~doi,             ~year,  ~crude_irr, ~hjust,
  "Manicaland",     1997.4, 2.25,       1,
  "Kisesa",         2015.8, 1.00,       0,
  "Karonga",        2003.2, 0.50,       0.5,
  "AHRI",           2023.0, 0.45,       0,
  "Kyamulibwa GPC", 2023.3, 0.28,       0,
  "Kisumu HDSS",    2022.4, 0.17,       0,
  "RCCS",           2018.3, 0.20,       1
) %>%
  mutate(cohort = cohort_factor(doi), label = cohort_short[doi])
fig3a_leader <- tibble(x = 2005.8, y = 0.93, xend = 2003.6, yend = 0.56,
                       cohort = cohort_factor("Karonga"))

fig3a <- ggplot(cohort_trends, aes(x = year, y = crude_irr, color = cohort)) +
  annotate("segment", x = 2010, xend = 2010, y = 0, yend = 1, linetype = 2) +
  annotate("segment", x = 1990, xend = 2010, y = 1, yend = 1, linetype = 2) +
  geom_point(aes(size = pys), alpha = 0.5) +
  geom_smooth(aes(weight = pys), se = FALSE, span = 0.9, linewidth = 1.5, show.legend = FALSE,
              method = "loess", formula = y ~ x) +
  scale_y_continuous(transform = "log", breaks = c(0.1, 0.25, 0.5, 1, 2, 5), expand = expansion(add = 0)) +
  scale_x_continuous(breaks = seq(1995, 2020, 5), expand = expansion(add = 0)) +
  scale_size(guide = "none") +
  geom_segment(data = fig3a_leader, aes(x = x, y = y, xend = xend, yend = yend, color = cohort),
               linewidth = 0.6, inherit.aes = FALSE, show.legend = FALSE) +
  geom_text(data = fig3a_labels, aes(label = label, hjust = hjust), size = 7, fontface = "bold",
            show.legend = FALSE) +
  coord_cartesian(xlim = c(1990, 2032), ylim = c(0.1, 3)) +
  labs(x = NULL, y = "Normalised 15-49 HIV incidence", color = NULL, tag = "A") +
  theme_paper() +
  larger_text() +
  theme(legend.position = "none")

cohort_modelled <- read_csv(file.path(path_modelled, "figure_incidence.csv"), show_col_types = FALSE) %>%
  inner_join(cohort_area, by = "area_id") %>%
  filter(year <= 2023) %>%
  mutate(source = "Modelled estimates")

fig3b_order <- c("Kisumu HDSS", "Kisesa", "Kyamulibwa GPC", "RCCS", "Karonga", "Manicaland", "AHRI")
fig3b <- cohort_trends %>%
  mutate(source = "Cohort data") %>%
  ggplot(aes(x = year, color = source)) +
  geom_line(data = cohort_modelled, aes(y = incidence), linewidth = 1.5) +
  geom_pointrange(aes(y = ir, ymin = lower, ymax = upper)) +
  facet_wrap(~factor(cohort_labels[doi], levels = cohort_labels[fig3b_order]), ncol = 3) +
  scale_color_manual(values = pal_two) +
  scale_x_continuous(breaks = seq(1995, 2020, 10), expand = expansion(add = 0)) +
  scale_y_continuous(labels = scales::label_number(scale = 100)) +
  coord_cartesian(ylim = c(0, 0.06)) +
  labs(x = NULL, y = "15-49 HIV incidence/100py", color = NULL, tag = "B") +
  theme_paper() +
  larger_text() +
  theme(strip.text = element_text(size = 18))

save_png(ggpubr::ggarrange(fig3a, fig3b, nrow = 2, heights = c(1, 1.3)), "figure_3.png", width = 800, height = 1200)

################################################################################
#' ## Figure 4: PHIA surveys and national modelled incidence (ages 15-49)
################################################################################

phia_spectrum <- read_fit("phia_spectrum.csv")

fig4a <- phia_spectrum %>%
  select(iso3, year, PHIA = incid, lower, upper, `National modelled estimates` = modelled) %>%
  pivot_longer(c(PHIA, `National modelled estimates`), names_to = "source", values_to = "incid") %>%
  mutate(across(c(lower, upper), ~ ifelse(source == "PHIA", .x, NA)),
         source = factor(source, levels = c("National modelled estimates", "PHIA"),
                         labels = c("National modelled estimates", "Population-based HIV Impact Assessment surveys (PHIA)")),
         lab = paste0(iso3, "\n", year)) %>%
  left_join(region_lookup, by = "iso3") %>%
  mutate(region = factor(ifelse(region == "WCA", "WCA", four_region), levels = c("WCA", "Eastern", "Southern"),
                         labels = c("West & Central\nAfrica", "Eastern Africa", "Southern Africa"))) %>%
  ggplot(aes(x = lab, y = incid, group = source)) +
  geom_col(aes(fill = source), position = position_dodge2(0.9)) +
  geom_linerange(aes(ymin = lower, ymax = upper), position = position_dodge2(0.9), na.rm = TRUE) +
  scale_fill_manual(values = pal_two) +
  scale_y_continuous(labels = scales::label_number(scale = 100)) +
  ggh4x::facet_grid2(~region, space = "free_x", scales = "free_x") +
  labs(x = NULL, y = "Incidence/100py, adults age 15-49 years", fill = NULL, tag = "A") +
  theme_paper() +
  theme(plot.tag = element_text(size = 25), axis.text.x = element_text(size = 16), axis.text.y = element_text(size = 20),
        axis.title = element_text(size = 22), strip.text = element_text(size = 20), legend.text = element_text(size = 20))

phia_fit <- read_fit("phia_spectrum_fit.csv")
fig4b <- ggplot(phia_spectrum, aes(x = modelled, y = incid)) +
  geom_ribbon(data = phia_fit, aes(y = mean, ymin = lower, ymax = upper), fill = pal_blue, alpha = 0.3) +
  geom_line(data = phia_fit, aes(y = mean), color = pal_blue, linewidth = 1.5) +
  geom_point(aes(size = pys), show.legend = FALSE) +
  geom_linerange(aes(ymin = lower, ymax = upper), linewidth = 0.7) +
  geom_abline(intercept = 0, slope = 1, linetype = 2) +
  scale_x_continuous(transform = "log", labels = scales::label_number(scale = 100, accuracy = 0.01)) +
  scale_y_continuous(transform = "log", labels = scales::label_number(scale = 100, accuracy = 0.01)) +
  coord_fixed(xlim = c(0.0002, 0.023), ylim = c(0.0002, 0.023)) +
  labs(x = "Modelled incidence/100py\nadults 15-49 (log scale)", y = "PHIA incidence/100py\nadults 15-49 (log scale)", tag = "B") +
  theme_paper() +
  theme(plot.tag = element_text(size = 25), axis.text = element_text(size = 20), axis.title = element_text(size = 22),
        aspect.ratio = 1)

two_phia <- phia_spectrum %>% count(iso3) %>% filter(n == 2) %>% pull(iso3)
national_15to49 <- spectrum %>%
  filter(variant == "published", sex == "both", age_group == "Y015_049", area_id == iso3,
         iso3 %in% two_phia, year %in% 2015:2022) %>%
  mutate(incidence = infections / (population - plhiv), source = "National modelled estimates")

fig4c <- phia_spectrum %>%
  filter(iso3 %in% two_phia) %>%
  mutate(source = "PHIA survey") %>%
  ggplot(aes(x = year, color = source)) +
  geom_line(data = national_15to49, aes(y = incidence), linewidth = 1.5) +
  geom_point(aes(y = incid), size = 3) +
  geom_linerange(aes(ymin = lower, ymax = upper), linewidth = 1.1) +
  facet_wrap(~country_name(iso3), nrow = 2) +
  scale_color_manual(values = pal_two, guide = "none") +
  scale_x_continuous(breaks = c(2016, 2020)) +
  scale_y_continuous(labels = scales::label_number(scale = 100)) +
  expand_limits(y = 0) +
  labs(x = NULL, y = "Incidence/100py\nadults 15-49", tag = "C") +
  theme_paper() +
  theme(plot.tag = element_text(size = 25), axis.text = element_text(size = 20), strip.text = element_text(size = 18),
        axis.title = element_text(size = 22), aspect.ratio = 1)

save_png(ggpubr::ggarrange(fig4a, ggpubr::ggarrange(fig4b, fig4c, nrow = 1, widths = c(1, 1.3)), nrow = 2),
         "figure_4.png", width = 1200, height = 1200)

################################################################################
#' ## Figure 5 and Supplementary Figure S4: age patterns in cohorts and PHIA surveys
################################################################################

age_ratio_theme <- function() {
  theme_paper() +
    larger_text() +
    theme(plot.title = element_text(size = 18, face = "bold"), strip.text = element_text(size = 16, hjust = 0),
          plot.tag = element_text(size = 22), aspect.ratio = 1)
}

#' A: 25-49 : 15-24 incidence rate ratio, cohort model (lines) and data (points)
cohort_age_ratio <- read_fit("cohort_age_predictions.csv") %>%
  inner_join(cohort_age_data %>% group_by(doi) %>% summarise(first = min(year), last = max(year)), by = "doi") %>%
  filter(year >= first, year <= last) %>%
  group_by(doi, year, sex) %>%
  summarise(ratio = mean[age_group == "Y025_049"] / mean[age_group == "Y015_024"], .groups = "drop")

cohort_age_points <- cohort_age_data %>%
  group_by(doi, year, sex) %>%
  filter(n() == 2) %>%
  summarise(ratio = ir[age_group == "Y025_049"] / ir[age_group == "Y015_024"], .groups = "drop") %>%
  filter(is.finite(ratio), ratio > 0)  # a zero in either age group cannot be shown on the log scale

fig5a <- ggplot(cohort_age_ratio, aes(x = year, y = ratio, color = cohort_factor(doi))) +
  geom_point(data = cohort_age_points, alpha = 0.6) +
  geom_line(linewidth = 1.3) +
  geom_hline(yintercept = 1, linetype = 2) +
  geom_text(data = tibble(sex = c("male", "female"), x = c(2006, 2004), y = c(0.3, 2.8),
                          label = c("Incidence higher among\nages 15-24 years", "Incidence higher among\nages 25-49 years")),
            aes(x = x, y = y, label = label), inherit.aes = FALSE, size = 5, fontface = "bold") +
  geom_segment(data = tibble(sex = "female"), aes(x = 1991, xend = 1991, y = 1.2, yend = 4), inherit.aes = FALSE,
               arrow = arrow(angle = 15, type = "closed")) +
  geom_segment(data = tibble(sex = "male"), aes(x = 2020, xend = 2020, y = 0.8, yend = 0.25), inherit.aes = FALSE,
               arrow = arrow(angle = 15, type = "closed")) +
  facet_wrap(~str_to_title(sex)) +
  scale_y_log10() +
  coord_cartesian(ylim = c(min(cohort_age_points$ratio), 4)) +
  scale_x_continuous(breaks = c(1990, 2000, 2010, 2020)) +
  guides(color = guide_legend(override.aes = list(alpha = 1), ncol = 3, position = "bottom")) +
  labs(x = NULL, y = "Incidence rate ratio,\nages 25-49 : 15-24 years", title = "Population-based cohorts", color = NULL, tag = "A") +
  age_ratio_theme()

#' B: PHIA surveys (public extract); model line from the fixed effects (country
#' effects zero), 25-49 : 15-24 ratio = exp(age + sex:age + year:age * year_min)
phia_age_data <- read_csv(file.path(path_raw, "phia-incid-extract_2025-11-13b.csv"), show_col_types = FALSE) %>%
  filter(age_group %in% c("Y015_024", "Y025_049"), sex != "both") %>%
  group_by(survey_id, iso3, year = survey_year, sex) %>%
  summarise(ratio = incid[age_group == "Y025_049"] / incid[age_group == "Y015_024"], .groups = "drop") %>%
  filter(is.finite(ratio), ratio > 0)

phia_coef <- read_fit("coefficients_phia_age.csv") %>% select(var, mean) %>% deframe() %>% log()
phia_age_line <- crossing(year = 2016:2022, sex = c("female", "male")) %>%
  mutate(ratio = exp(phia_coef["age_groupY025_049"] + (sex == "male") * phia_coef["sexmale:age_groupY025_049"] +
                       (year - reference_year) * phia_coef["year_min:age_groupY025_049"]))

fig5b <- ggplot(phia_age_data, aes(x = year, y = ratio)) +
  geom_line(aes(group = iso3), linewidth = 0.8, color = "grey75") +
  geom_point(size = 2, alpha = 0.5) +
  geom_line(data = phia_age_line, linewidth = 2, color = pal_blue) +
  geom_hline(yintercept = 1, linetype = 2) +
  facet_wrap(~str_to_title(sex)) +
  scale_x_continuous(breaks = c(2017, 2019, 2021)) +
  scale_y_log10() +
  coord_cartesian(ylim = c(min(phia_age_data$ratio), 4)) +
  labs(x = NULL, y = "Incidence rate ratio,\nages 25-49 : 15-24 years", title = "PHIA surveys (national incidence rate)", tag = "B") +
  age_ratio_theme()

#' C / S4: share of new 15-49 infections by sex and age, over the years with data
shares <- read_fit("cohort_infection_shares.csv") %>%
  filter(in_observed_range) %>%
  mutate(group = factor(paste(sex, age_group),
                        levels = c("female Y015_024", "female Y025_049", "male Y015_024", "male Y025_049"),
                        labels = c("Female 15-24", "Female 25-49", "Male 15-24", "Male 25-49"))) %>%
  arrange(doi, year, desc(group)) %>%
  group_by(doi, year) %>%
  mutate(end = cumsum(share), start = lag(end, default = 0), mid = (start + end) / 2) %>%
  ungroup()

share_plot <- function(dois, nrow) {
  d <- filter(shares, doi %in% dois)
  labels <- d %>%
    group_by(doi) %>%
    filter(year %in% range(year)) %>%
    ungroup() %>%
    mutate(lab = paste0(round(100 * share), "%"))
  key <- tibble(doi = intersect(c("Kyamulibwa GPC", dois), dois)[1],
                label = levels(shares$group), x = c(2005, 2005, 2002, 2005), y = c(0.85, 0.55, 0.37, 0.125),
                color = "white")
  ggplot(d, aes(x = year)) +
    geom_rect(aes(xmin = year - 0.5, xmax = year + 0.5, ymin = start, ymax = end, fill = group)) +
    geom_label(data = labels, aes(y = mid, label = lab, color = group), size = 5, fontface = "bold",
               linewidth = 1.5, show.legend = FALSE) +
    geom_text(data = key, aes(x = x, y = y, label = label), color = key$color, size = 5.5, fontface = "bold",
              inherit.aes = FALSE) +
    facet_wrap(~cohort_factor(doi), nrow = nrow) +
    expand_limits(x = c(1995, 2025)) +
    scale_y_continuous(labels = scales::label_percent(), expand = expansion(0, 0)) +
    scale_fill_manual(values = pal_four) +
    scale_color_manual(values = pal_four) +
    labs(x = NULL, y = "Share of new adult HIV infections,\nages 15-49 years", fill = NULL) +
    theme_paper() +
    larger_text() +
    theme(legend.position = "none", strip.text = element_text(size = 16, hjust = 0), plot.tag = element_text(size = 22))
}

fig5c <- share_plot(c("AHRI", "Kyamulibwa GPC", "Manicaland", "RCCS"), nrow = 2) + labs(tag = "C")
save_png(fig5a / fig5b / fig5c + plot_layout(heights = c(1, 0.9, 1.7)), "figure_5.png", width = 900, height = 1400)

save_png(share_plot(cohort_dois, nrow = 3), "figure_s4.png", width = 1100, height = 1300)

################################################################################
#' ## Supplementary Figure S3: share of new infections under 25 and sex ratio,
#' ## cohorts vs national (provincial) modelled estimates
################################################################################

#' Modelled estimates for each cohort: Thembisa KwaZulu-Natal, EPP Nyanza (KEN_1_7)
#' and Manicaland (ZWE_1_12), national Spectrum files for Malawi, Tanzania and Uganda
cohort_spectrum_area <- c("AHRI" = "ZAF_1_KZN", "Kisumu HDSS" = "KEN_1_7", "Karonga" = "MWI", "Kisesa" = "TZA",
                          "Kyamulibwa GPC" = "UGA", "Manicaland" = "ZWE_1_12", "RCCS" = "UGA")
cohort_spectrum <- spectrum %>%
  filter(variant == "published", area_id %in% cohort_spectrum_area, year %in% 1990:2023) %>%
  inner_join(enframe(cohort_spectrum_area, "doi", "area_id"), by = "area_id", relationship = "many-to-many")

s3_cohorts <- c("AHRI", "Kyamulibwa GPC", "Manicaland", "RCCS", "Kisesa")
s3_panel <- function(points, lines, y, ylab, tag) {
  ggplot(points, aes(x = year, y = .data[[y]], color = source)) +
    geom_point(size = 2) +
    geom_smooth(data = filter(points, doi %in% s3_cohorts), se = FALSE, span = 1.5, method = "loess", formula = y ~ x) +
    geom_line(data = lines, linewidth = 1.5) +
    facet_wrap(~cohort_factor(doi), ncol = 4) +
    scale_color_manual(values = pal_two) +
    scale_x_continuous(breaks = c(2000, 2015)) +
    labs(x = NULL, y = ylab, color = NULL, tag = tag) +
    theme_paper() +
    larger_text() +
    theme(strip.text = element_text(size = 16), aspect.ratio = 1)
}

u25_points <- cohort_age_data %>%
  group_by(doi, year) %>%
  summarise(prop = sum(infections[age_group == "Y015_024"]) / sum(infections), .groups = "drop") %>%
  mutate(source = "Cohort data")
u25_lines <- cohort_spectrum %>%
  filter(sex == "both", age_group %in% c("Y015_024", "Y025_049")) %>%
  group_by(doi, year) %>%
  summarise(prop = infections[age_group == "Y015_024"] / sum(infections), .groups = "drop") %>%
  mutate(source = "National modelled estimates")
fig_s3a <- s3_panel(u25_points, u25_lines, "prop", "Proportion of new infections 15-49 among 15-24", "A") +
  expand_limits(y = c(0, 1)) +
  scale_y_continuous(labels = scales::label_percent())

sr_points <- cohort_age_data %>%
  group_by(doi, year, sex) %>%
  summarise(ir = sum(infections) / sum(pys), .groups = "drop") %>%
  group_by(doi, year) %>%
  filter(n() == 2) %>%
  summarise(sr = ir[sex == "female"] / ir[sex == "male"], .groups = "drop") %>%
  mutate(source = "Cohort data")
sr_lines <- cohort_spectrum %>%
  filter(sex != "both", age_group == "Y015_049", year >= 1995) %>%
  group_by(doi, year) %>%
  summarise(sr = incidence[sex == "female"] / incidence[sex == "male"], .groups = "drop") %>%
  mutate(source = "National modelled estimates")
fig_s3b <- s3_panel(sr_points, sr_lines, "sr", "Female : male HIV incidence\nrate ratio (15-49)", "B") +
  geom_hline(yintercept = 1, linetype = 2) +
  scale_y_log10(breaks = c(0.6, 1, 2, 4))

save_png(ggpubr::ggarrange(fig_s3a, fig_s3b, nrow = 2), "figure_s3.png", width = 1000, height = 1200)

################################################################################
#' ## Text S4 figure: KwaZulu-Natal district trajectories vs Thembisa province
################################################################################

kzn_districts <- read_csv(file.path(path_geography, "areas.csv"), show_col_types = FALSE) %>%
  filter(spectrum_area_id == "ZAF_1_KZN", area_level == 2) %>%
  pull(area_id)
kzn <- read_csv(file.path(path_modelled, "figure_incidence.csv"), show_col_types = FALSE) %>%
  filter(area_id %in% c("ZAF_1_KZN", kzn_districts)) %>%
  mutate(source = ifelse(area_id == "ZAF_1_KZN", "Thembisa: KwaZulu-Natal", "Extrapolated district estimates"))

fig_text_s4 <- ggplot(kzn, aes(x = year, y = incidence, group = area_id, color = source)) +
  geom_line(data = filter(kzn, area_id != "ZAF_1_KZN"), alpha = 0.6) +
  geom_line(data = filter(kzn, area_id == "ZAF_1_KZN"), linewidth = 1.5) +
  scale_color_manual(values = c("black", "red")) +
  scale_y_continuous(labels = scales::label_number(scale = 100)) +
  labs(x = NULL, y = "HIV incidence/100py, adults 15-49", color = NULL) +
  theme_paper() +
  larger_text()

save_png(fig_text_s4, "figure_text_s4.png", width = 700, height = 500)

################################################################################
#' ## Text S7 figure: national 15-49 incidence with and without PHIA calibration
################################################################################

phia_countries <- unique(phia_spectrum$iso3)
#' Panels for Western and Central, Eastern and Southern Africa; y axes start at 0 and
#' share an upper limit within each region
fig_text_s7_data <- spectrum %>%
  filter(sex == "both", age_group == "Y015_049", area_id == iso3, iso3 %in% phia_countries, year %in% 2000:2025) %>%
  left_join(select(region_lookup, iso3, four_region), by = "iso3") %>%
  mutate(variant = factor(variant, levels = c("published", "no_phia"),
                          labels = c("Published estimates", "Without PHIA incidence")),
         region = case_when(four_region %in% c("Western", "Central") ~ "Western and Central Africa",
                            TRUE ~ paste(four_region, "Africa")) %>%
           factor(levels = c("Western and Central Africa", "Eastern Africa", "Southern Africa"))) %>%
  arrange(region, country_name(iso3)) %>%
  mutate(country = fct_inorder(country_name(iso3)))

#' One panel per region, titled with the region; countries within a region share the y axis
text_s7_region <- function(region, ncol, y_title = NULL) {
  fig_text_s7_data %>%
    filter(region == !!region) %>%
    ggplot(aes(x = year, y = incidence, color = variant)) +
    geom_line(linewidth = 1) +
    facet_wrap(~country, ncol = ncol) +
    scale_color_manual(values = pal_two) +
    scale_x_continuous(breaks = seq(2000, 2020, 10)) +
    scale_y_continuous(labels = scales::label_number(scale = 100), limits = c(0, NA),
                       expand = expansion(mult = c(0, 0.05))) +
    labs(title = region, x = NULL, y = y_title, color = NULL) +
    theme_paper() +
    theme_large() +
    theme(panel.spacing.x = unit(1, "lines"))
}

#' Four countries per row: Western and Central Africa one row, Eastern and Southern two each
#' (the three-country row keeps the width of three columns of four)
fig_text_s7 <- (text_s7_region("Western and Central Africa", 3) + plot_spacer() + plot_layout(widths = c(3, 1))) /
  text_s7_region("Eastern Africa", 4, y_title = "HIV incidence/100py, adults 15-49") /
  text_s7_region("Southern Africa", 4) +
  plot_layout(heights = c(1, 2, 2), guides = "collect") &
  theme(legend.position = "bottom", legend.text = element_text(size = 12 * 1.6))

save_png(fig_text_s7, "figure_text_s7.png", width = 1100, height = 1350)
