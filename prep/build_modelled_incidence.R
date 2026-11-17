################################################################################
#' Build the modelled-incidence extracts in data/modelled/ (Supplementary Text S4)
#'
#' Needs private/ (prep/copy_private_inputs.R) and data/geography/
#' (prep/build_geography.R). Run from the repository root.
#'
#' Method. Naomi district estimates are available for one year (2025 for most
#' countries). For each district, sex and age group, incidence in that year is
#' multiplied by the ratio of incidence in each year 1990-2025 to incidence in the
#' Naomi year, taken from the Spectrum estimates for the area containing the
#' district (national Spectrum file, provincial EPP subpopulation, or Thembisa
#' province for South Africa). Observations recorded at a Spectrum area use that
#' area's Spectrum estimate directly.
#'
#' Release rule: values derived from the Namibia 2026 Spectrum file (private) are
#' written to private/ only, and blank (suppressed = TRUE) in data/modelled/.
#'
#' Writes:
#'   data/modelled/matched_incidence.csv.gz  one row per area x year x sex x age
#'     group with data: modelled incidence under the main matching (naomi) and the
#'     sensitivity analyses (naomi_no_phia, spectrum, spectrum_national)
#'   data/modelled/spectrum_incidence.csv.gz Spectrum-area estimates (infections,
#'     population, PLHIV, incidence) by sex for ages 15-24, 25-49 and 15-49,
#'     published files and files refit without PHIA incidence
#'   data/modelled/cohort_population.csv     population and HIV prevalence in the
#'     cohort districts, ages 15-24 and 25-49 (Supplementary Table S11)
################################################################################

library(tidyverse)

source("R/functions.R")
source("R/data_functions.R")

dir.create(path_modelled, showWarnings = FALSE, recursive = TRUE)

private_spectrum_iso3 <- "NAM"

areas <- read_csv(file.path(path_geography, "areas.csv"), show_col_types = FALSE)
aligned <- suppressWarnings(align_observations())$aligned

iso3_in_data <- areas %>%
  filter(area_id %in% aligned$area_id) %>%
  distinct(iso3) %>%
  pull(iso3)

data_age_groups <- unique(aligned$age_group)

################################################################################
#' ## Inputs
################################################################################

naomi <- readRDS(file.path(path_private, "naomi_indicators.rds")) %>%
  filter(indicator %in% c("population", "plhiv", "incidence", "prevalence"))

#' Naomi estimate year by country
naomi_year <- distinct(naomi, iso3, naomi_year = year)

spectrum_national <- bind_rows(readRDS(file.path(path_private, "spectrum_national.rds")))
spectrum_national_no_phia <- bind_rows(readRDS(file.path(path_private, "spectrum_national_no_phia.rds")))

#' EPP subpopulation fits; for the no-PHIA variant, the refit files replace the
#' published files of the same name
epp <- readRDS(file.path(path_private, "subnational_epp.rds"))
epp_no_phia <- epp %>%
  filter(!file %in% read_csv(file.path(path_private, "no_phia_files.csv"), show_col_types = FALSE)$file) %>%
  bind_rows(readRDS(file.path(path_private, "subnational_epp_no_phia.rds")))

#' Single-year ages to five-year groups (0-4, ..., 75-79, 80+)
to_five_year <- function(df) {
  df %>%
    mutate(age = as.integer(age),
           age_group = ifelse(age >= 80, "Y080_999", age_group_code(5 * (age %/% 5), 5 * (age %/% 5) + 4))) %>%
    select(-age)
}

ages_15to79 <- age_group_code(seq(15, 75, 5), seq(19, 79, 5))

################################################################################
#' ## Population by district over time: Naomi age/sex proportions of the national
#' ## Spectrum population, held constant (all area levels)
################################################################################

extrapolate_population <- function(spectrum_national) {

  spec_pop <- spectrum_national %>%
    filter(iso3 %in% naomi$iso3, year %in% 1990:2025) %>%
    to_five_year() %>%
    group_by(iso3, year, sex, age_group) %>%
    summarise(totpop = sum(totpop), .groups = "drop")

  naomi_pop <- naomi %>%
    group_by(iso3) %>%
    filter(indicator == "population", age_group %in% five_year_age_groups,
           area_level == max(area_level), sex != "both") %>%
    ungroup() %>%
    select(iso3, area_id, sex, age_group, naomi_year = year, mean)

  prop <- naomi_pop %>%
    left_join(spec_pop, by = c("iso3", "sex", "age_group", "naomi_year" = "year")) %>%
    transmute(iso3, area_id, sex, age_group, prop = mean / totpop)

  lowest <- spec_pop %>%
    inner_join(prop, by = c("iso3", "sex", "age_group"), relationship = "many-to-many") %>%
    mutate(pop = prop * totpop) %>%
    select(iso3, area_id, year, sex, age_group, pop)

  # Aggregate to every ancestor level
  parent <- setNames(areas$parent_area_id, areas$area_id)
  ancestors <- tibble(area_id = unique(lowest$area_id)) %>%
    mutate(ancestor = map(area_id, function(a) {
      out <- a
      while (!is.na(parent[[a]])) { a <- parent[[a]]; out <- c(out, a) }
      out
    })) %>%
    unnest(ancestor)

  lowest %>%
    inner_join(ancestors, by = "area_id", relationship = "many-to-many") %>%
    group_by(iso3, area_id = ancestor, year, sex, age_group) %>%
    summarise(pop = sum(pop), .groups = "drop")
}

################################################################################
#' ## Spectrum-area estimates by five-year age group
################################################################################

epp_iso3 <- c("TZA", "MOZ", "NAM", "SWZ", "NGA", "ZWE", "KEN")
epp_level <- c(TZA = 2, MOZ = 1, NAM = 1, SWZ = 1, NGA = 2, ZWE = 1, KEN = 1)

#' National (or ETH regional) Spectrum files
spectrum_file_estimates <- function(spectrum_national) {
  regions <- areas %>%
    filter(spectrum_region_code != 0, area_level == 1) %>%
    select(iso3, area_id, spectrum_region_code)

  spec <- spectrum_national %>%
    filter(iso3 %in% iso3_in_data, iso3 != "ZAF") %>%
    left_join(regions, by = c("iso3", "spectrum_region_code")) %>%
    mutate(area_id = ifelse(is.na(area_id), iso3, area_id))

  bind_rows(spec,
            spectrum_national %>% filter(iso3 == "ETH") %>% mutate(area_id = iso3)) %>%
    to_five_year() %>%
    filter(age_group %in% ages_15to79) %>%
    group_by(iso3, area_id, year, sex, age_group) %>%
    summarise(across(c(totpop, hivpop, infections), sum), .groups = "drop")
}

#' EPP subpopulations: national age/sex incidence and prevalence scaled by the
#' subpopulation:national ratio (ages 15-49) from EPP; population from Naomi
epp_estimates <- function(epp, spectrum_national, pop) {

  sub <- epp %>%
    filter(iso3 %in% epp_iso3, area_name != "National") %>%
    mutate(area_name = case_when(area_name == "Karas" ~ "!Karas",
                                 area_name == "Zambezia" ~ "Zambézia",
                                 area_name == "CrossRiver" ~ "Cross River",
                                 area_name == "Akwa-Ibom" ~ "Akwa Ibom",
                                 TRUE ~ area_name),
           area_name = str_trim(str_remove(area_name, "Region")),
           year = as.integer(year)) %>%
    left_join(areas %>%
                filter(area_level == epp_level[iso3]) %>%
                select(iso3, area_id, area_name),
              by = c("iso3", "area_name")) %>%
    filter(year < 2026)

  if (any(is.na(sub$area_id))) stop("Unmatched EPP subpopulations: ", toString(unique(sub$area_name[is.na(sub$area_id)])))

  ratios <- sub %>%
    filter(year > 1989) %>%
    mutate(infections = incid * pop * (1 - prev), plhiv = pop * prev, hiv_neg = pop * (1 - prev)) %>%
    group_by(iso3, year) %>%
    mutate(irr = incid / (sum(infections) / sum(hiv_neg)),
           pr = prev / (sum(plhiv) / sum(pop))) %>%
    ungroup() %>%
    select(iso3, area_id, year, irr, pr)

  national <- spectrum_national %>%
    filter(iso3 %in% epp_iso3) %>%
    to_five_year() %>%
    group_by(iso3, year, sex, age_group) %>%
    summarise(across(c(totpop, hivpop, infections), sum), .groups = "drop") %>%
    filter(age_group %in% ages_15to79) %>%
    transmute(iso3, year, sex, age_group, incidence = infections / (totpop - hivpop), prevalence = hivpop / totpop)

  # Sex-specific population weights the five-year groups and the two sexes when
  # aggregated. (The submitted analysis used the both-sex population for each sex.)
  ratios %>%
    left_join(national, by = c("iso3", "year"), relationship = "many-to-many") %>%
    left_join(pop, by = c("iso3", "area_id", "year", "sex", "age_group")) %>%
    mutate(incidence = incidence * irr,
           prevalence = prevalence * pr,
           infections = incidence * pop * (1 - prevalence),
           hivpop = prevalence * pop) %>%
    select(iso3, area_id, year, sex, age_group, totpop = pop, hivpop, infections)
}

#' Thembisa v5.0: national and provincial, single age, 1990-2025
thembisa_estimates <- function() {
  file <- file.path(path_raw, "modelled", "thembisa_5.0_age_specific_outputs.xlsx")
  provinces <- c("South Africa", "Eastern Cape", "Free State", "Gauteng", "KwaZulu-Natal", "Limpopo",
                 "Mpumalanga", "Northern Cape", "North West", "Western Cape")
  map2(2:11, provinces, function(sheet, province) {
    readxl::read_xlsx(file, sheet = sheet, .name_repair = "unique_quiet") %>%
      filter(!is.na(`...1`)) %>%
      mutate(indicator = ifelse(str_detect(`...1`, " "), `...1`, NA)) %>%
      fill(indicator, .direction = "down") %>%
      filter(!is.na(`2021`)) %>%
      pivot_longer(-c(indicator, `...1`), names_to = "year") %>%
      rename(age = `...1`) %>%
      mutate(area_name = province)
  }) %>%
    bind_rows() %>%
    left_join(areas %>% filter(iso3 == "ZAF", area_level == 1) %>% select(area_id, area_name), by = "area_name") %>%
    mutate(area_id = ifelse(area_name == "South Africa", "ZAF", area_id)) %>%
    filter(indicator %in% c("Male population", "Female population", "Male HIV incidence", "Female HIV incidence",
                            "Male HIV prevalence", "Female HIV prevalence")) %>%
    mutate(sex = ifelse(str_detect(indicator, "Female"), "female", "male"),
           indicator = str_extract(indicator, "population|prevalence|incidence")) %>%
    type.convert(as.is = TRUE) %>%
    filter(age %in% 15:79, year %in% 1990:2025) %>%
    select(area_id, year, sex, age, indicator, value) %>%
    pivot_wider(names_from = indicator, values_from = value) %>%
    mutate(hivpop = population * prevalence, infections = incidence * (population - hivpop)) %>%
    to_five_year() %>%
    group_by(iso3 = "ZAF", area_id, year, sex, age_group) %>%
    summarise(totpop = sum(population), hivpop = sum(hivpop), infections = sum(infections), .groups = "drop")
}

#' Aggregate five-year estimates to the age groups in the data (five-year groups
#' wholly within each), add both sexes
to_data_age_groups <- function(five_year, age_groups) {
  five_year %>%
    filter(year %in% 1990:2025) %>%
    inner_join(five_year_groups_within(age_groups), by = "age_group", relationship = "many-to-many") %>%
    select(-age_group) %>%
    sex_aggregation(c("totpop", "hivpop", "infections")) %>%
    group_by(iso3, area_id, year, sex, age_group = age_group_data) %>%
    summarise(across(c(totpop, hivpop, infections), sum), .groups = "drop") %>%
    mutate(incidence = infections / (totpop - hivpop))
}

#' District population is always extrapolated from the published Spectrum files,
#' including for the no-PHIA variant (as in the submitted analysis)
naomi_population <- extrapolate_population(spectrum_national)

spectrum_area_estimates <- function(spectrum_national, epp, age_groups) {
  pop <- naomi_population
  five_year <- bind_rows(spectrum_file_estimates(spectrum_national),
                         epp_estimates(epp, spectrum_national, pop),
                         thembisa_estimates())
  to_data_age_groups(five_year, age_groups)
}

################################################################################
#' ## Naomi estimates in the data age groups (non-Spectrum areas)
################################################################################

#' Naomi five-year groups overlapping each data age group (for groups that Naomi
#' does not report directly)
naomi_age_map <- function(age_groups) {
  direct <- intersect(age_groups, unique(naomi$age_group))
  other <- tibble(reqd = setdiff(age_groups, direct)) %>%
    mutate(start = age_group_start(reqd), end = age_group_end(reqd)) %>%
    filter(start < end)
  five <- tibble(age_group = age_group_code(seq(15, 75, 5), seq(19, 79, 5))) %>%
    mutate(age = map(age_group_start(age_group), ~ .x + 0:4)) %>%
    unnest(age)
  bind_rows(
    other %>%
      cross_join(five) %>%
      filter(age >= start, age <= end) %>%
      distinct(reqd, age_group),
    tibble(reqd = direct, age_group = direct)
  )
}

naomi_in_data_age_groups <- function(age_groups, spectrum_areas) {
  naomi %>%
    filter(!area_id %in% spectrum_areas, indicator != "prevalence") %>%
    inner_join(naomi_age_map(age_groups), by = "age_group", relationship = "many-to-many") %>%
    select(iso3, area_id, sex, reqd, age_group, indicator, mean) %>%
    pivot_wider(names_from = indicator, values_from = mean) %>%
    mutate(infections = incidence * (population - plhiv)) %>%
    group_by(iso3, area_id, sex, age_group = reqd) %>%
    summarise(across(c(population, plhiv, infections), sum), .groups = "drop") %>%
    mutate(incidence = infections / (population - plhiv))
}

################################################################################
#' ## Extrapolated incidence
################################################################################

extrapolate_incidence <- function(spectrum_dat, area_ids = aligned$area_id) {

  ratio <- spectrum_dat %>%
    left_join(naomi_year, by = "iso3") %>%
    mutate(naomi_year = ifelse(iso3 == "ZAF", 2025, naomi_year)) %>%
    group_by(iso3, area_id, sex, age_group) %>%
    mutate(ratio = incidence / incidence[year == naomi_year]) %>%
    ungroup() %>%
    select(iso3, spectrum_area_id = area_id, year, sex, age_group, ratio)

  naomi_dat <- naomi_in_data_age_groups(data_age_groups, unique(spectrum_dat$area_id)) %>%
    filter(area_id %in% area_ids) %>%
    left_join(areas %>% select(area_id, spectrum_area_id), by = "area_id")

  extrapolated <- naomi_dat %>%
    inner_join(ratio, by = c("iso3", "spectrum_area_id", "sex", "age_group"), relationship = "many-to-many") %>%
    transmute(iso3, area_id, year, sex, age_group, value = incidence * ratio)

  bind_rows(extrapolated,
            spectrum_dat %>% transmute(iso3, area_id, year, sex, age_group, value = incidence))
}

spectrum_dat <- spectrum_area_estimates(spectrum_national, epp, data_age_groups)
spectrum_dat_no_phia <- spectrum_area_estimates(spectrum_national_no_phia, epp_no_phia, data_age_groups)

main <- extrapolate_incidence(spectrum_dat)
no_phia <- extrapolate_incidence(spectrum_dat_no_phia)

#' Check against the extrapolated incidence used in the submitted analysis
check <- function(new, file) {
  old <- readRDS(file.path(path_private, file)) %>%
    filter(area_id %in% aligned$area_id) %>%
    inner_join(new, by = c("area_id", "year", "sex", "age_group"))
  message(file, ": ", nrow(old), " rows matched; max relative difference ",
          signif(max(abs(old$value / old$extrapolated_value - 1), na.rm = TRUE), 3))
  invisible(old)
}
check(main, "extrapolated_incidence_2026.rds")
check(no_phia, "extrapolated_incidence_2026_no_phia.rds")

################################################################################
#' ## Matched incidence for every area/year/sex/age group in the data
################################################################################

needed <- aligned %>%
  distinct(area_id, sex, age_group) %>%
  left_join(areas %>% select(area_id, iso3, spectrum_area_id), by = "area_id") %>%
  mutate(iso3 = coalesce(iso3, substr(area_id, 1, 3))) %>%
  cross_join(tibble(year = 1990:2025))

spec_inc <- spectrum_dat %>% select(spectrum_area_id = area_id, year, sex, age_group, spectrum = incidence)

matched_incidence <- needed %>%
  left_join(main %>% select(area_id, year, sex, age_group, naomi = value), by = c("area_id", "year", "sex", "age_group")) %>%
  left_join(no_phia %>% select(area_id, year, sex, age_group, naomi_no_phia = value), by = c("area_id", "year", "sex", "age_group")) %>%
  left_join(spec_inc, by = c("spectrum_area_id", "year", "sex", "age_group")) %>%
  left_join(spec_inc %>% rename(iso3 = spectrum_area_id, spectrum_national = spectrum), by = c("iso3", "year", "sex", "age_group")) %>%
  filter(if_any(c(naomi, naomi_no_phia, spectrum, spectrum_national), ~ !is.na(.x))) %>%
  mutate(suppressed = iso3 %in% private_spectrum_iso3) %>%
  select(area_id, year, sex, age_group, naomi, naomi_no_phia, spectrum, spectrum_national, suppressed) %>%
  arrange(area_id, sex, age_group, year)

#' Full values (private) and released values (Namibia blank)
write_csv(matched_incidence, file.path(path_private, "matched_incidence_full.csv.gz"))
matched_incidence %>%
  mutate(across(c(naomi, naomi_no_phia, spectrum, spectrum_national), ~ ifelse(suppressed, NA, .x))) %>%
  write_csv(file.path(path_modelled, "matched_incidence.csv.gz"), na = "")

################################################################################
#' ## Spectrum-area estimates for figures (ages 15-24, 25-49, 15-49)
################################################################################

fig_age_groups <- c("Y015_024", "Y025_049", "Y015_049")

spectrum_all <- function(spectrum_national, epp) {
  # all countries' national files, not only those with data
  iso3_keep <- iso3_in_data
  iso3_in_data <<- unique(spectrum_national$iso3)
  on.exit(iso3_in_data <<- iso3_keep)
  spectrum_area_estimates(spectrum_national, epp, fig_age_groups)
}

spectrum_incidence <- bind_rows(
  spectrum_all(spectrum_national, epp) %>% mutate(variant = "published"),
  spectrum_all(spectrum_national_no_phia, epp_no_phia) %>% mutate(variant = "no_phia")
) %>%
  rename(population = totpop, plhiv = hivpop) %>%
  mutate(suppressed = iso3 %in% private_spectrum_iso3) %>%
  select(variant, iso3, area_id, year, sex, age_group, population, plhiv, infections, incidence, suppressed)

write_csv(spectrum_incidence, file.path(path_private, "spectrum_incidence_full.csv.gz"))
spectrum_incidence %>%
  mutate(across(c(population, plhiv, infections, incidence), ~ ifelse(suppressed, NA, .x))) %>%
  write_csv(file.path(path_modelled, "spectrum_incidence.csv.gz"), na = "")

################################################################################
#' ## Cohort districts: population and HIV prevalence, ages 15-24 and 25-49
#' Naomi prevalence in its estimate year, extrapolated with the Spectrum-area
#' prevalence trend; population from extrapolate_population()
################################################################################

cohort_areas <- c("ZAF_2_DC27", "MWI_5_02", "TZA_4_136ev", "KEN_2_41", "UGA_3_048df", "ZWE_1_12", "UGA_3_121")

pop <- naomi_population %>%
  filter(area_id %in% cohort_areas, age_group %in% age_groups_15to49) %>%
  mutate(age_group = ifelse(age_group %in% c("Y015_019", "Y020_024"), "Y015_024", "Y025_049")) %>%
  group_by(area_id, year, sex, age_group) %>%
  summarise(population = sum(pop), .groups = "drop")

spec_prev <- spectrum_incidence %>%
  filter(variant == "published", age_group %in% c("Y015_024", "Y025_049"), sex != "both") %>%
  transmute(spectrum_area_id = area_id, iso3, year, sex, age_group, prev = plhiv / population) %>%
  left_join(naomi_year, by = "iso3") %>%
  mutate(naomi_year = ifelse(iso3 == "ZAF", 2025, naomi_year)) %>%
  group_by(spectrum_area_id, sex, age_group) %>%
  mutate(prev_ratio = prev / prev[year == naomi_year]) %>%
  ungroup()

cohort_population <- naomi %>%
  filter(area_id %in% cohort_areas, indicator == "prevalence", age_group %in% c("Y015_024", "Y025_049"), sex != "both") %>%
  select(area_id, sex, age_group, naomi_prev = mean) %>%
  left_join(areas %>% select(area_id, spectrum_area_id), by = "area_id") %>%
  left_join(spec_prev %>% select(spectrum_area_id, year, sex, age_group, prev_ratio),
            by = c("spectrum_area_id", "sex", "age_group"), relationship = "many-to-many") %>%
  transmute(area_id, year, sex, age_group, prevalence = naomi_prev * prev_ratio) %>%
  left_join(pop, by = c("area_id", "year", "sex", "age_group")) %>%
  select(area_id, year, sex, age_group, population, prevalence)

write_csv(cohort_population, file.path(path_modelled, "cohort_population.csv"))

################################################################################
#' ## Extrapolated 15-49 incidence, both sexes, for figures
#' Cohort districts (Figure 3B) and all KwaZulu-Natal districts with the
#' Thembisa province estimate they are extrapolated from (Text S4 figure)
################################################################################

kzn_areas <- c("ZAF_1_KZN", areas$area_id[areas$spectrum_area_id == "ZAF_1_KZN" & areas$area_level == 2])
figure_areas <- c(cohort_areas, kzn_areas)

extrapolate_incidence(spectrum_dat, figure_areas) %>%
  filter(area_id %in% figure_areas, sex == "both", age_group == "Y015_049") %>%
  select(area_id, year, incidence = value) %>%
  arrange(area_id, year) %>%
  write_csv(file.path(path_modelled, "figure_incidence.csv"))
