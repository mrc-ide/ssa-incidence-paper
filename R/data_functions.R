################################################################################
#' Functions to build the analysis dataset from the systematic review extracts
#' and the additional study data
#'
#' Pipeline (run by R/01_data.R):
#'   read_review_data() -> clean_review_data() -> impute_person_years()
#'   read_additional_data()
#'   match_areas()            observation -> UNAIDS subnational area(s)
#'   match_modelled()         observation -> matched modelled incidence (offset)
#'   classify_studies()       study type and recruitment type
################################################################################

#' Track the number of studies and observations through each exclusion step
track <- function(tracker, dat, step) {
  bind_rows(tracker, tibble(step = step, studies = n_distinct(dat$doi), observations = nrow(dat)))
}

################################################################################
#' ## Systematic review data
################################################################################

#' Joshi et al. review extract plus the 2024 update (three extraction sheets)
read_review_data <- function(dir = file.path(path_raw, "review")) {

  # Studies without a DOI (PHIA, SHIMS and THIS survey reports and three papers) are
  # keyed by author and year, so each counts as its own study in Figure 1
  joshi <- read_csv(file.path(dir, "data.combine.csv"), show_col_types = FALSE) %>%
    mutate(`doi/url` = ifelse(is.na(`doi/url`), paste(author, pubyear), `doi/url`))

  #' Rows reporting cumulative incidence over a period, not a rate, are dropped
  read_update <- function(file) {
    x <- read_csv(file.path(dir, file), show_col_types = FALSE) %>%
      mutate(row = row_number())
    period <- filter(x, str_detect(IR_PER_100, "year|period"))
    x %>%
      filter(!row %in% period$row) %>%
      select(-row, -FollowUp) %>%
      type.convert(as.is = TRUE)
  }

  update <- bind_rows(read_update("2024_i.csv"), read_update("2024_ii_fixed.csv"), read_update("2024_iii.csv")) %>%
    select(author = Author, `doi/url`, StudyClass, TimeLow, TimeHigh, `Country (ISO0)`, ISO1, `ISO3 [Name]`,
           Sex, AgeLow, AgeHigh, starts_with("N_"), starts_with("IR"), PersonYears, group, remove) %>%
    rename(N_Seroconvert = N_Seroconversions, N_total = N_Total) %>%
    mutate(StudyClass = recode(StudyClass,
                               "Cross Sectional Incidence" = "CSI", "Clinical Trial" = "CT",
                               "Prospective Cohort" = "PC", "Retrospective Cohort" = "RC"),
           InterventionArm = case_when(
             StudyClass == "CT" & str_detect(group, "ontrol") ~ "control",
             StudyClass == "CT" & !str_detect(group, "ontrol") ~ "intervention/overall",
             TRUE ~ NA),
           idx = row_number() + max(joshi$idx))

  bind_rows(joshi, update)
}

#' 2024 update studies extracted in Covidence that give no incidence rate (Figure 1):
#' those reporting only cumulative incidence over a period (rows dropped in
#' read_review_data) and those with nothing to extract (extracted_no_rate.csv)
read_update_no_rate <- function(dir = file.path(path_raw, "review")) {
  sheets <- bind_rows(lapply(c("2024_i.csv", "2024_ii_fixed.csv", "2024_iii.csv"), function(f)
    read_csv(file.path(dir, f), show_col_types = FALSE, col_types = cols(.default = "c"))))
  rate_rows <- filter(sheets, !str_detect(IR_PER_100, "year|period") | is.na(IR_PER_100))
  cumulative <- setdiff(unique(sheets$`doi/url`), rate_rows$`doi/url`)
  bind_rows(tibble(doi = cumulative, reason = "Cumulative incidence only"),
            read_csv(file.path(dir, "extracted_no_rate.csv"), show_col_types = FALSE) %>%
              select(doi, reason))
}

#' Joshi et al. papers that are secondary analyses of trials whose primary reports
#' are extracted by hand (new_data.csv, source "primary_report"). One row per Joshi
#' paper; `primary_report_doi` lists the primary report(s), separated by "; "
read_primary_reports <- function(dir = file.path(path_raw, "review")) {
  read_csv(file.path(dir, "replaced_by_primary_report.csv"), show_col_types = FALSE)
}

primary_report_dois <- function(primary = read_primary_reports()) {
  unique(str_trim(unlist(str_split(primary$primary_report_doi, ";"))))
}

read_inclusion <- function(dir = file.path(path_raw, "review")) {
  read_csv(file.path(dir, "inclusion.csv"), show_col_types = FALSE) %>%
    rename(doi = `doi/url`)
}

#' Apply the review exclusion criteria (Figure 1). Returns the cleaned
#' observations and a tracker of studies/observations after each step.
clean_review_data <- function(dat, inclusion, match_level, no_rate = read_update_no_rate(),
                              primary = read_primary_reports()) {

  notes_doi <- function(pattern, exclude = NULL) {
    x <- filter(inclusion, str_detect(notes, pattern))
    if (!is.null(exclude)) x <- filter(x, !str_detect(notes, exclude))
    x$doi
  }

  dat <- dat %>%
    rename(iso3 = `Country (ISO0)`, area_name = `ISO3 [Name]`, pys = PersonYears, N_neg = N_Neg, doi = `doi/url`) %>%
    left_join(match_level, by = c("iso3", "ISO1", "ISO2", "area_name"))

  # Extracted update studies with no incidence rate have no observations; they are
  # counted in the pooled data and excluded at "Unable to calculate incidence rate"
  stopifnot(!any(no_rate$doi %in% dat$doi))
  n_no_rate <- n_distinct(no_rate$doi)
  # Joshi et al. secondary analyses replaced by the trials' primary reports, which join
  # after the exclusions (new_data.csv; read_additional_data)
  stopifnot(all(primary$doi %in% dat$doi))
  track <- function(tracker, dat, step) {
    bind_rows(tracker, tibble(step = step, studies = n_distinct(dat$doi) + n_no_rate,
                              observations = nrow(dat)))
  }

  tracker <- track(tibble(), dat, "Pooled systematic review data")

  dat <- filter(dat, !doi %in% primary$doi)
  tracker <- track(tracker, dat, "Replaced by primary trial report")

  # Total population studies with risk-based inclusion criteria, and key populations.
  # Studies of FSW, MSM or transgender people only are counted first (these populations
  # come from the dedicated reviews instead); then all other risk-based studies
  risk_ok <- function(x) {
    (x$HighRiskGroup %in% c("GP", "PW", "PW/GP") | is.na(x$HighRiskGroup)) &
      !x$doi %in% notes_doi("HR", exclude = "AHRI")
  }
  kp_groups <- c("SW", "MSM", "MSM*", "MSM/TG", "TG", "SW/MSM")
  kp_only <- dat %>%
    filter(!doi %in% dat$doi[risk_ok(dat)]) %>%
    group_by(doi) %>%
    filter(all(HighRiskGroup %in% kp_groups)) %>%
    pull(doi) %>%
    unique()
  dat <- filter(dat, !doi %in% kp_only)
  tracker <- track(tracker, dat, "FSW and MSM studies")

  dat <- dat[risk_ok(dat), ]
  tracker <- track(tracker, dat, "Risk-based inclusion criteria")

  dat <- filter(dat, !doi %in% filter(inclusion, remove == "PrEP")$doi)
  tracker <- track(tracker, dat, "PrEP")

  # Keep control arms (and overall ECHO, which was not an HIV prevention trial)
  dat <- filter(dat,
                (is.na(InterventionArm) |
                   InterventionArm %in% c("control", "placebo", "Placebo Gel", "oral placebo") |
                   str_detect(group, "control") |
                   (doi == "10.1016/S0140-6736(19)31288-7" & InterventionArm == "overall")),
                !doi %in% filter(inclusion, notes == "Remove: intervention study")$doi)
  tracker <- track(tracker, dat, "Intervention arms or intervention study")

  dat <- filter(dat, !doi %in% filter(inclusion, notes == "Remove: BED")$doi)
  tracker <- track(tracker, dat, "BED assay")

  dat <- dat %>%
    filter(is.na(remove),
           !doi %in% notes_doi("long term follow up"),
           !doi %in% notes_doi("eplicate"))
  dat <- remove_overlapping_review(dat)
  tracker <- track(tracker, dat, "Overlapping observations and duplicates")

  dat <- filter(dat, !doi %in% filter(inclusion, str_detect(notes, "Remove"), str_detect(notes, "analysis"))$doi)
  tracker <- track(tracker, dat, "Secondary analyses")

  # Replaced by more finely stratified data from study teams (read_additional_data)
  replaced <- function(x) {
    filter(x, !doi %in% c("10.1056/NEJMoa1814556",         # HPTN071/PopART
                          "10.1056/NEJMoa1809866",         # SEARCH
                          "10.1056/NEJMoa1812281",         # Ya Tsie/BCPP
                          "10.1016/S2352-3018(17)30205-9", # TasP
                          notes_doi("agg")),
           !str_detect(author, "PHIA|SHIMS|THIS"))
  }
  tracker <- track(tracker, replaced(dat), "Sex/age/year stratified data available")

  # Person-years are imputed from all studies left after the secondary-analysis step;
  # observations with no incidence data drop out here and are counted at the next step
  dat <- impute_person_years(dat, inclusion) %>%
    replaced() %>%
    filter(!doi %in% filter(inclusion, notes == "Remove: No incidence data")$doi)
  n_no_rate <- 0
  tracker <- track(tracker, dat, "Unable to calculate incidence rate")

  stopifnot(!any(dat$doi %in% filter(inclusion, str_detect(notes, "Remove"))$doi))

  dat <- dat %>%
    clean_age_groups() %>%
    clean_year() %>%
    # Study type as re-reviewed (inclusion.csv)
    left_join(inclusion %>%
                filter(!str_detect(notes, "Remove") | is.na(notes)) %>%
                select(doi, StudyClass, study_type),
              by = c("doi", "StudyClass")) %>%
    select(-StudyClass) %>%
    rename(StudyClass = study_type)

  list(dat = dat, tracker = tracker)
}

#' Overlapping observations not flagged in the review extracts, found when the
#' final dataset was checked for double counting (peer review, November 2026).
#' Each study's person-time should enter the model once.
remove_overlapping_review <- function(dat) {

  huerga <- "10.1186/s12879-017-2740-y"
  grabbe <- "10.1521/aeap.2015.27.4.350"
  machekano <- "10.1371/journal.pone.0209782"
  gray <- "10.1097/QAD.0b013e3283504a3f."
  mehta <- "10.1097/01.aids.0000432444.30308.2d"
  kinuthia_2010 <- "10.2174/157016210793499213"

  # Huerga 2017 reports every stratum three times, in extract order: incidence with no
  # ARV measure, with self-reported ART, and with ARV detection. Keep ARV detection.
  huerga_rows <- filter(dat, doi == huerga)
  stopifnot(nrow(huerga_rows) == 12, all(count(huerga_rows, Sex, AgeLow, AgeHigh)$n == 3))

  dat %>%
    group_by(doi, Sex, AgeLow, AgeHigh) %>%
    filter(!doi %in% huerga | row_number() == 3) %>%
    ungroup() %>%
    filter(
      # Grabbe 2015 reports the same 259 infections overall, by district, by year, by
      # sex/age and by another two-way split: keep the calendar-year rows
      !(doi %in% grabbe & TimeLow %in% "01/01/2006" & TimeHigh %in% "30/06/2010"),
      # Machekano 2018: the postpartum-only estimate overlaps the age-specific rows,
      # which cover pregnancy and postpartum follow-up
      !(doi %in% machekano & HighRiskGroup %in% "GP"),
      # Gray 2012: post-trial surveillance of the Rakai circumcision trial ran within the
      # RCCS, whose men 15-49 come from the study team (read_additional_data)
      !doi %in% gray,
      # Mehta 2013: long-term follow-up of the Kisumu circumcision trial, already
      # included as Bailey 2007 (data/raw/review/new_data.csv)
      !doi %in% mehta) %>%
    # Kinuthia et al. 2010 (Curr HIV Res) carries the
    # 2011-13 dates of Kinuthia 2015. The paper gives no dates; it was published in
    # October 2010 and cites 2008 Kenyan births, so take 2008-09.
    mutate(TimeLow = ifelse(doi %in% kinuthia_2010, "01/01/2008", TimeLow),
           TimeHigh = ifelse(doi %in% kinuthia_2010, "31/12/2009", TimeHigh))
}

#' Person-years where missing (Supplementary Text S2). `pys_source` records how
#' each observation's person-years were obtained.
impute_person_years <- function(dat, inclusion) {

  dat <- dat %>%
    filter(N_total != N_positive_baseline | is.na(N_positive_baseline) | is.na(N_total), # all HIV+ at baseline
           !(is.na(IR_PER_100) & is.na(N_Seroconvert) & is.na(N_RI))) %>%                 # no incidence data
    mutate(ir = IR_PER_100 / 100,
           lower = IR_CI_LOW / 100,
           upper = IR_CI_HIGH / 100,
           ir = ifelse(is.na(ir) & StudyClass != "CSI", N_Seroconvert / pys, ir),
           pys_source = ifelse(is.na(pys), NA, "reported"))

  # 1. From the 95% interval of the incidence rate
  dat <- dat %>%
    mutate(se = (upper - ir) / 1.96,
           pys_source = ifelse(is.na(pys) & !is.na(upper), "interval", pys_source),
           pys = ifelse(is.na(pys) & !is.na(upper), (ir * (1 - ir)) / (se^2), pys))

  # 2. From the number HIV-negative at baseline (or total participants), using the
  #    log-linear relationship among observations reporting both
  fit <- lme4::lmer(log(pys) ~ log(N_neg) + (1 | doi),
                    data = filter(dat, !str_detect(doi, "PHIA"), StudyClass != "CSI",
                                  !is.na(N_neg), !is.na(pys), pys != 0))
  beta <- lme4::fixef(fit)

  dat <- dat %>%
    mutate(pys_source = ifelse(is.na(pys) & !is.na(N_neg), "participants_negative", pys_source),
           N_neg = ifelse(is.na(N_neg), N_total - N_positive_baseline, N_neg),
           pys_source = ifelse(is.na(pys) & is.na(pys_source) & !is.na(N_neg), "participants_total_minus_positive", pys_source),
           N_neg = ifelse(is.na(N_neg), N_total, N_neg),
           pys_source = ifelse(is.na(pys) & is.na(pys_source) & !is.na(N_neg), "participants_total", pys_source),
           pys = ifelse(is.na(pys), exp(beta[[1]] + beta[[2]] * log(N_neg)), pys))

  # Cross-sectional incidence reported as counts of recent infections: mean duration
  # of recent infection (MDRI) from the study, else 130 days
  csi <- dat %>%
    filter(is.na(ir), StudyClass == "CSI") %>%
    left_join(inclusion %>% filter(StudyClass == "CSI") %>% select(doi, mdri), by = "doi") %>%
    mutate(mdri = ifelse(is.na(mdri), 130, mdri),
           pys = N_neg * mdri / 365,
           pys_source = "recency",
           ir = N_RI / pys)

  dat %>%
    filter(!idx %in% csi$idx) %>%
    mutate(ir = ifelse(is.na(ir), N_Seroconvert / pys, ir)) %>%
    bind_rows(csi) %>%
    filter(!is.na(ir))
}

#' Lower age bounds rounded to five years (below 16 -> 15); missing upper bound -> 49
clean_age_groups <- function(dat) {
  ends <- seq(17, 97, 5)
  dat %>%
    mutate(sex = case_when(Sex == "F" ~ "female", Sex == "M" ~ "male", TRUE ~ "both")) %>%
    filter(AgeHigh > 16 | is.na(AgeHigh)) %>%
    mutate(age_group_start = ifelse(AgeLow < 16 | is.na(AgeLow), 15, plyr::round_any(AgeLow, 5)),
           AgeHigh = ifelse(is.na(AgeHigh), 49, AgeHigh),
           end_i = findInterval(AgeHigh, ends),
           age_group_end = ifelse(end_i == 0, NA, ends[pmax(end_i, 1)] + 2),
           age_group = paste0("Y0", age_group_start, "_0", age_group_end))
}

#' Midpoint calendar year of the study period
clean_year <- function(dat) {
  dat %>%
    separate(TimeLow, sep = "/", into = c(NA, NA, "year_start"), remove = FALSE, convert = TRUE) %>%
    separate(TimeHigh, sep = "/", into = c(NA, NA, "year_end"), remove = FALSE, convert = TRUE) %>%
    mutate(year = floor((year_start + year_end) / 2),
           year = case_when(is.na(year) & !is.na(year_end) ~ year_end,
                            is.na(year) & !is.na(year_start) ~ year_start,
                            TRUE ~ year))
}

################################################################################
#' ## Additional data from study teams, PHIA surveys and key population reviews
################################################################################

read_additional_data <- function(inclusion, search_areas, dir = file.path(path_raw, "additional"),
                                 phia_file = file.path(path_raw, "phia-incid-extract_INTERNAL_2024-05-07.csv")) {

  # HPTN071/PopART: Arm C (standard of care) communities
  popart_areas <- tribble(
    ~demcom, ~area_name, ~area_id,
    3, "Kabwe", "ZMB_2_5",
    4, "Ndola", "ZMB_2_21",
    7, "Lusaka", "ZMB_2_53",
    12, "Choma", "ZMB_2_78",
    15, "Cape Town MM", "ZAF_2_CPT",
    17, "Cape Town MM", "ZAF_2_CPT",
    21, "Cape Winelands DM", "ZAF_2_DC2"
  )
  popart <- read_csv(file.path(dir, "popart.csv"), show_col_types = FALSE) %>%
    filter(arm == "Arm C") %>%
    inner_join(popart_areas, by = "demcom") %>%
    mutate(sex = tolower(sexc), year = 2016, age_group = age_label_to_group(ageJeffcat)) %>%
    filter(!is.na(age_group)) %>%  # ages 18-19 are not a five-year group and are dropped
    select(area_name, area_id, year, sex, age_group, ir = Incidence, pys = Person_Years) %>%
    mutate(doi = "PopART", study_type = "CT")

  rakai <- read_csv(file.path(dir, "rakai.csv"), show_col_types = FALSE) %>%
    mutate(sex = tolower(Sex),
           low = as.integer(str_sub(Agecat, 2, 3)),
           high = as.integer(str_sub(Agecat, 5, 6)),
           age_group = age_group_code(low, high - 1),
           area_name = "Rakai", area_id = "UGA_3_121") %>%
    select(area_name, area_id, year = Year, sex, age_group, ir = crude_inc, pys = py) %>%
    mutate(doi = "RCCS", study_type = "PC")

  bcpp <- read_csv(file.path(dir, "bcpp.csv"), show_col_types = FALSE) %>%
    filter(random_arm == "Standard of Care") %>%
    group_by(sex, age_group) %>%
    mutate(idx = cur_group_id()) %>%
    ungroup() %>%
    mutate(doi = "BCCP/YaTsie", sex = tolower(sex), year = 2016,
           age_group = age_label_to_group(age_group),
           value = paste0(community, ", Botswana"),
           ir = infections / person_years,
           study_type = "CT") %>%
    rename(pys = person_years)

  phia <- read_csv(phia_file, show_col_types = FALSE) %>%
    filter(age_group %in% five_year_age_groups, sex != "both") %>%
    select(doi = survey_id, area_id = iso3, year = survey_year, sex, age_group, ir = incid, pys = n_pys_approx) %>%
    mutate(study_type = "CSI")

  manicaland <- read_csv(file.path(dir, "manicaland.csv"), show_col_types = FALSE) %>%
    mutate(ir = seroconversions / pys, doi = "Manicaland", sex = tolower(sex), study_type = "PC")

  # 2023 data were not to be used (study team)
  ahri <- bind_rows(readxl::read_excel(file.path(dir, "ahri.xlsx"), sheet = "Male"),
                    readxl::read_excel(file.path(dir, "ahri.xlsx"), sheet = "Female")) %>%
    mutate(sex = tolower(Female),
           age_group = age_group_code(as.integer(str_sub(AgeCat, 2, 3)), as.integer(str_sub(AgeCat, 5, 6)) - 1),
           area_id = "ZAF_2_DC27",
           ir = sero / pyrs) %>%
    select(area_id, year = Year, sex, age_group, ir, pys = pyrs) %>%
    filter(year != 2023) %>%
    mutate(doi = "AHRI", study_type = "PC")

  # SEARCH incidence is not sex-disaggregated (the "male" column is both sexes).
  # Cluster person-years unavailable: 8461, the average across clusters (study team).
  search <- read_csv(file.path(dir, "search.csv"), show_col_types = FALSE) %>%
    filter(intervention == 0) %>%
    left_join(search_areas, by = "district") %>%
    transmute(area_id, ir = incidence_male, sex = "both", age_group = "Y015_079", year = 2015, pys = 8461,
              doi = "SEARCH", study_type = "CT")

  kisumu <- read_kisumu(file.path(dir, "kisumu.csv")) %>%
    mutate(study_type = "PC")

  tasp_raw <- readRDS(file.path(dir, "tasp.rds"))
  tasp <- tasp_raw %>%
    filter(arm == "Control") %>%
    mutate(age_group = age_label_to_group(ifelse(age_group == "16-19", "15-19", as.character(age_group))),
           year = lubridate::year(date_midpoint),
           area_id = "ZAF_2_DC27",
           doi = "TasP") %>%
    select(area_id, sex, age_group, year, sero = n_seroconversions, pys = person_years_at_risk, doi) %>%
    mutate(ir = sero / pys, study_type = "CT")

  read_cohort <- function(file, doi, sero) {
    read_csv(file.path(dir, file), show_col_types = FALSE) %>%
      mutate(age_group = age_label_to_group(age_group_label), doi = doi, ir = .data[[sero]] / pys, study_type = "PC") %>%
      select(-age_group_label)
  }
  gpc <- read_cohort("gpc.csv", "Kyamulibwa GPC", "seroconversions") %>% filter(year < 2023)
  karonga <- read_cohort("karonga.csv", "Karonga", "seroconversion")
  # Kisesa incidence as supplied by the study team
  kisesa <- read_csv(file.path(dir, "kisesa.csv"), show_col_types = FALSE) %>%
    mutate(age_group = age_label_to_group(age_group_label), doi = "Kisesa", study_type = "PC") %>%
    select(-age_group_label)

  kp <- read_csv(file.path(dir, "kp_incidence.csv"), show_col_types = FALSE) %>%
    filter(year > 1995) %>%
    # Overlapping MSM observations (peer review, November 2026). Kamali 2015 reports
    # the IAVI Kilifi and Nairobi cohorts overall, by time on study and by calendar
    # period; Price 2012 is earlier follow-up of the same cohorts and Wahome 2018 the
    # same Kilifi cohort to 2016. Keep Kamali's calendar periods for Kilifi and its
    # overall estimate for Nairobi.
    filter(!ref %in% c("Price 2012", "Wahome 2018"),
           !(ref == "Kamali 2015" & area_id == "KEN_2_03" & !pys %in% c(304, 692)),
           !(ref == "Kamali 2015" & area_id == "KEN_1_8" & pys != 338),
           # Kimani 2019: the 15-34 row is the 15-24 and 25-34 rows combined
           !(ref == "Kimani 2019" & age_group == "Y015_034"),
           # Wirtz 2015 (doi:10.1097/QAI.0000000000000693): two sets of 2012-13 rows with the
           # same 7 seroconversions (all at the first follow-up) but different person-time.
           # Keep the set totalling ~147 py, consistent with 103 men followed Jan 2012 - May 2013
           !(ref == "Wirtz 2015" & round(pys, 3) %in% c(79.58, 23.125))) %>%
    mutate(kp = ifelse(kp == "MSM2", "MSM", kp),
           sex = recode(kp, FSW = "female", MSM = "male", PWID = "male", TGW = "female", CFSW = "male"),
           study_type = recode(method,
                               "prospective cohort" = "PC", "RCT" = "CT", "cohort" = "PC",
                               "retrospective cohort" = "RC", "prev backcalculation" = "RC",
                               "non-randomised plausability evaluation" = "CT", "cross-sectional" = "CSI - Other"),
           area_id = recode(area_id,
                            CMR_1_12 = "CMR_1_2",        # Yaounde merged into Centre
                            ZMB_1_14 = "ZMB_1_14fz",
                            KEN_3_3042le = "KEN_2_03"),
           doi = ifelse(kp == "FSW", paste0("Jones & Anderson | ", ref), paste0("Stannah | ", ref)),
           ir = infections / pys)

  kp_doi <- unique(kp$doi)
  add_data <- bind_rows(popart, rakai, manicaland, ahri, search, tasp, phia, kisumu, bcpp, gpc, karonga, kisesa, kp) %>%
    mutate(sexual_activity = ifelse(doi %in% kp_doi, 1, 0)) %>%
    select(doi, sexual_activity, study_type, value, area_id, age_group, year, sex, ir, pys) %>%
    filter(!age_group %in% c("Y080_084", "Y080_999", "Y085_089")) %>%
    mutate(pys_source = "reported")

  # Observations outside the original search strategy, extracted by hand
  new_data <- read_csv(file.path(dir, "../review/new_data.csv"), show_col_types = FALSE) %>%
    mutate(se = (upper - ir) / 1.96,
           pys_source = ifelse(is.na(pys) & !is.na(upper), "interval", "reported"),
           pys = ifelse(is.na(pys) & !is.na(upper), (ir * (1 - ir)) / (se^2), pys),
           year = round((year_start + year_end) / 2)) %>%
    left_join(inclusion %>% select(doi, sexual_activity), by = "doi") %>%
    mutate(sexual_activity = ifelse(is.na(sexual_activity), 0, 1)) %>%
    select(-source)

  bind_rows(add_data, new_data)
}

#' KEMRI-Siaya (Kisumu HDSS). No person-years for men in 2020-2022: imputed from
#' 2016-2019 male person-years and the female 2020-2022:2016-2019 ratio.
read_kisumu <- function(file) {
  kisumu <- read_csv(file, show_col_types = FALSE) %>%
    separate(year, into = c("year_start", "year_end"), convert = TRUE) %>%
    mutate(year = round((year_start + year_end) / 2),
           age_group = case_when(age_group_label == "35-44" ~ "Y035_044",
                                 age_group_label == "45-54" ~ "Y045_054",
                                 age_group_label == "55+" ~ "Y055_064", # cohort recruited ages 15-64
                                 TRUE ~ age_label_to_group(age_group_label))) %>%
    rename(sero = seroconversions)

  ratio <- kisumu %>%
    filter(year %in% c(2018, 2021), sex == "female") %>%
    select(year, age_group, pys) %>%
    pivot_wider(names_from = year, values_from = pys) %>%
    transmute(age_group, ratio = `2021` / `2018`)

  kisumu %>%
    filter(!(sex == "male" & year == 2021)) %>%
    bind_rows(kisumu %>%
                filter(sex == "male", year == 2021) %>%
                select(-pys) %>%
                left_join(kisumu %>%
                            filter(sex == "male", year == 2018) %>%
                            left_join(ratio, by = "age_group") %>%
                            transmute(sex, age_group, pys = pys * ratio),
                          by = c("sex", "age_group"))) %>%
    mutate(area_id = "KEN_2_41") %>%
    select(area_id, year, sex, age_group, sero, pys, ir, doi = source)
}

################################################################################
#' ## Matching observations to UNAIDS subnational areas (Supplementary Text S3)
################################################################################

#' One row per geocoding search string (study sites below province level)
make_geocode_searches <- function(dat, add_data) {

  havlir <- tibble(
    country = c(rep("Uganda", 20), rep("Kenya", 12)),
    site = c("Nsiika", "Bugamba", "Rugazi", "Mitooma", "Kitwe", "Rubaare", "Rwashamaire", "Ruhoko", "Kazo", "Nyamuyanja",
             "Nsiinze", "Nankoma", "Kiyunga", "Kamuge", "Bugono", "Muyembe", "Merikit", "Kiyeyi", "Kameke", "Kadama",
             "Nyatoto", "Nyamrisra", "Ogongo", "Kitare", "Magunga", "Kisegi", "Tom Mboya", "Sena", "Ongo", "Orthoro",
             "Sibouche", "Bware")
  ) %>%
    transmute(value = site, new_value = paste0(site, ", ", country), author = "Havlir", pubyear = 2019)

  searches <- dat %>%
    filter(level == "g") %>%
    mutate(area_name = case_when(is.na(area_name) & is.na(ISO2) & is.na(ISO1) ~ iso3,
                                 is.na(area_name) & is.na(ISO2) & !is.na(ISO1) ~ ISO1,
                                 is.na(area_name) & !is.na(ISO2) ~ ISO2,
                                 TRUE ~ area_name),
           area_name = str_replace(area_name, "\\,", "\\|")) %>%
    separate(area_name, into = paste0("split", 1:20), sep = "\\|", fill = "right") %>%
    pivot_longer(cols = split1:split20) %>%
    filter(!is.na(value)) %>%
    mutate(value = str_trim(value),
           value = case_when(str_detect(value, "Klerksdorp") ~ "Klerksdorp",
                             value == "Itended" ~ "Itende",
                             value == "Nsiinze" ~ "Nsinze",
                             value == "Bugono" ~ "Bugono Church of Uganda",
                             value == "Kiyeyi" ~ "Kiyeyi High School",
                             value == "Sena" ~ "Sena Centre",
                             value == "Sibouche" ~ "Sibouche Centre, Kenya",
                             value == "Kisegi" ~ "Kisigi Point",
                             value == "Ongo" ~ "Ongo Health Centre, Kenya",
                             value == "Bware" ~ "Bware Secondary School",
                             value == "Orthoro" ~ "Othoro Sub County Hospital",
                             value == "Tom Mboya" ~ "Tom Mboya Memorial Sub-County Hospital",
                             value == "Bugando Medical Centre (caters to Mwanza" ~ "Mwanza",
                             TRUE ~ value),
           iso3 = case_when(value == "Mbeya" & iso3 == "Uganda" & author == "Robb" ~ "Tanzania",
                            value == "Kericho" & iso3 == "Tanzania" & author == "Robb" ~ "Kenya",
                            TRUE ~ iso3),
           value = ifelse(!str_detect(iso3, "\\|"), paste0(value, ", ", iso3), value),
           value = case_when(value == "Epworth" ~ "Epworth, Zimbabwe",
                             value == "Jinja" ~ "Jinja, Uganda",
                             TRUE ~ value)) %>%
    left_join(havlir, by = c("value", "author", "pubyear")) %>%
    mutate(value = ifelse(!is.na(new_value), new_value, value))

  bcpp <- add_data %>%
    filter(doi == "BCCP/YaTsie") %>%
    mutate(idx = row_number() + max(searches$idx), study = "BCCP")

  bind_rows(searches, bcpp) %>%
    select(idx, study, value) %>%
    nest(data = idx) %>%
    mutate(google_idx = row_number())
}

#' Assign each review observation to one or more area_ids: national (level 0),
#' province (level 1, by name) or geocoded site (level "g"). Each area is counted
#' once per observation, however many map hits fall in it.
match_areas <- function(dat, areas, geocoded) {

  admin0 <- dat %>%
    filter(level == 0) %>%
    separate_longer_delim(iso3, "|") %>%
    mutate(area_id = countrycode::countrycode(iso3, "country.name", "iso3c")) %>%
    distinct(idx, area_id)

  admin1 <- dat %>%
    filter(level == 1) %>%
    mutate(ISO1 = ifelse(is.na(ISO1) & !is.na(ISO2), ISO2, ISO1),
           iso3_code = ifelse(!str_detect(iso3, "\\|"), countrycode::countrycode(iso3, "country.name", "iso3c"), NA),
           ISO1 = ifelse(iso3_code == "BWA" & ISO1 == "Northern", "North West|North East", ISO1)) %>%
    separate_longer_delim(ISO1, "|") %>%
    filter(!is.na(ISO1)) %>%  # province-level observations with no province recorded are dropped
    mutate(value = case_when(iso3_code == "RWA" & ISO1 == "Kigali Province" ~ "Kigali City",
                             iso3_code == "RWA" & ISO1 == "South" ~ "Southern Province",
                             iso3_code == "RWA" & ISO1 == "West" ~ "Western Province",
                             iso3_code == "RWA" & ISO1 == "North" ~ "Northern Province",
                             iso3_code == "RWA" & ISO1 == "East" ~ "Eastern Province",
                             iso3_code == "ZMB" & ISO1 == "North" ~ "Northern",
                             iso3_code == "BWA" & ISO1 == "North" ~ "Northern",
                             ISO1 == "KwaZulu Natal" ~ "KwaZulu-Natal",
                             ISO1 == "Harare Province" ~ "Harare",
                             ISO1 == "Nyanza Region" ~ "Nyanza",
                             TRUE ~ ISO1)) %>%
    left_join(areas %>% select(iso3, area_name, area_id, area_level),
              by = c("iso3_code" = "iso3", "value" = "area_name"), relationship = "many-to-many") %>%
    group_by(idx, value) %>%
    filter(area_level == min(area_level) | is.na(area_level)) %>%
    ungroup() %>%
    mutate(area_id = case_when(is.na(area_id) & value == "Kigali" ~ "RWA_1_1",
                               is.na(area_id) & value == "Nairobi" ~ "KEN_1_8",
                               TRUE ~ area_id)) %>%
    distinct(idx, value, area_id)

  if (any(is.na(admin1$area_id))) stop("Unassigned provinces")

  sites <- dat %>%
    filter(level == "g") %>%
    select(idx) %>%
    inner_join(geocoded, by = "idx", relationship = "one-to-many") %>%
    distinct(idx, area_id)

  bind_rows(admin0, select(admin1, idx, area_id), sites)
}

################################################################################
#' ## Matching to modelled incidence (Supplementary Text S4)
################################################################################

#' Join each observation-area to the modelled incidence for that area, year, sex
#' and age group, then collapse observations spanning several areas to the median
#' modelled incidence. `offset_col` selects the matched estimate (main analysis or
#' sensitivity variants).
match_modelled <- function(aligned, modelled, offset_col = "naomi") {

  matched <- aligned %>%
    left_join(modelled %>% select(area_id, year, sex, age_group, offset = all_of(offset_col)),
              by = c("area_id", "year", "sex", "age_group"))

  matched %>%
    filter(!is.na(offset),
           pys != 0) %>%  # pys = 0 arises from imputing pys from intervals when ir = 0
    group_by(doi, idx, year, sex, age_group) %>%
    mutate(area_id = toString(area_id)) %>%
    group_by(author, doi, study_type, idx, area_id, year, sex, age_group) %>%
    summarise(extrapolated_value = DescTools::Quantile(offset, pys, 0.5, names = FALSE, na.rm = TRUE),
              ir = median(ir),
              pys = median(pys),
              pys_source = first(pys_source),
              .groups = "drop") %>%
    mutate(extrapolated_value = ifelse(is.nan(extrapolated_value), NA, extrapolated_value))
}

################################################################################
#' ## Study classification
################################################################################

#' NB the submitted analysis classed Karonga as Other HHS (a typo in its
#' UTT/Cohort list); corrected here.
classify_studies <- function(matched, inclusion) {

  national_hhs_doi <- c("https://hsrc.ac.za/uploads/pageContent/10779/SABSSM%20V.pdf",
                        "https://hsrc.ac.za/uploads/pageContent/4565/SABSSM%20IV%20LEO%20final.pdf",
                        "https://hsrc.ac.za/wp-content/uploads/2024/07/SABSSM_VI_EXEC_REPORT_2PP.pdf")
  utt_cohort_doi <- c("PopART", "SEARCH", "BCCP/YaTsie", "TasP",
                      "RCCS", "Manicaland", "Karonga", "Kisesa", "Kisumu HDSS", "Kyamulibwa GPC", "AHRI")
  pregnant_doi <- filter(inclusion, sexual_activity == "Pregnant")$doi

  out <- matched %>%
    left_join(distinct(inclusion, doi, study_type, representative), by = c("doi", "study_type")) %>%
    mutate(infections = ir * pys,
           study_type = case_when(str_detect(doi, "Anderson \\|") ~ "FSW",
                                  str_detect(doi, "Stannah \\|") ~ "MSM",
                                  study_type == "CSI" ~ "Cross-sectional incidence",
                                  study_type %in% c("HIVPT", "HIVVT", "CT") ~ "Intervention Trials",
                                  study_type == "PC" ~ "Prospective Cohort",
                                  study_type == "RC" ~ "Retrospective Cohort",
                                  study_type == "STA" ~ "Secondary Trial Analysis",
                                  study_type == "SCS" ~ "Serial Cross Sectional",
                                  study_type == "NCC" ~ "Nested Case Control",
                                  TRUE ~ study_type),
           representative = ifelse(study_type %in% c("FSW", "MSM"), "Representative", representative),
           recruitment_type = case_when(
             doi %in% utt_cohort_doi ~ "UTT/Cohort",
             str_detect(doi, "PHIA") | doi %in% national_hhs_doi ~ "National HHS",
             study_type %in% c("Intervention Trials", "FSW", "MSM") ~ study_type,
             representative == "Representative" ~ "Other HHS",
             representative == "Specified inclusion criteria" ~ "Non-representative - Other",
             doi %in% pregnant_doi ~ "ANC",
             representative == "No inclusion criteria" ~ "No inclusion criteria"),
           recruitment_type = factor(recruitment_type, levels = recruitment_levels)) %>%
    group_by(doi) %>%
    mutate(study_idx = cur_group_id()) %>%
    ungroup()

  # SABSSM V reported no person-years: values from the survey report
  sabssm_v <- "https://hsrc.ac.za/uploads/pageContent/10779/SABSSM%20V.pdf"
  sabssm_pys <- tribble(
    ~sex, ~age_group, ~pys,
    "female", "Y015_024", 1009.777778,
    "male", "Y015_024", 468.3442937,
    "female", "Y025_079", 3180.727023,
    "male", "Y025_079", 1346.482577
  )
  out <- out %>%
    filter(doi != sabssm_v) %>%
    bind_rows(out %>%
                filter(doi == sabssm_v) %>%
                select(-pys) %>%
                left_join(sabssm_pys, by = c("sex", "age_group")) %>%
                # (the submitted analysis did not recompute infections here)
                mutate(infections = ir * pys, pys_source = "report")) %>%
    # Duplicates the age-stratified observations from the same study
    filter(!(doi == "10.1093/ofid/ofy060" & age_group == "Y015_064"))

  out
}

################################################################################
#' ## Whole pipeline up to observation-area pairs
################################################################################

#' Read, clean and area-match all observations. Returns one row per observation x
#' matched area (`aligned`), the exclusion tracker and the additional data counts
#' for Figure 1.
#'
#' Map hits are joined to searches by search text. (The submitted analysis joined
#' them by position in the geocode cache, which is offset by one for 20 searches,
#' so those sites took the neighbouring search's areas.)
align_observations <- function(geography = path_geography) {

  inclusion <- read_inclusion()
  match_level <- read_csv(file.path(path_raw, "review", "match_level.csv"), show_col_types = FALSE)
  areas <- read_csv(file.path(geography, "areas.csv"), show_col_types = FALSE)
  search_areas <- read_csv(file.path(geography, "search_areas.csv"), show_col_types = FALSE)
  geocoded_locations <- read_csv(file.path(geography, "geocoded_locations.csv"), show_col_types = FALSE)

  review <- clean_review_data(read_review_data(), inclusion, match_level)
  dat <- review$dat

  add_data <- read_additional_data(inclusion, search_areas)

  searches <- make_geocode_searches(dat, add_data)
  geocoded <- searches %>%
    inner_join(geocoded_locations %>% select(value, area_id), by = "value", relationship = "many-to-many") %>%
    unnest(data) %>%
    distinct(idx, study, value, area_id)

  review_areas <- match_areas(dat, areas, geocoded %>% filter(is.na(study)))

  aligned <- dat %>%
    select(author, doi, study_type = StudyClass, idx, pubyear, country = iso3, ISO1, ISO2, area_name,
           year, sex, age_group, ir, pys, pys_source) %>%
    inner_join(review_areas, by = "idx", relationship = "one-to-many")

  # Ya Tsie/BCPP communities are geocoded too
  bcpp <- add_data %>%
    filter(doi == "BCCP/YaTsie") %>%
    select(-area_id) %>%
    left_join(geocoded %>%
                filter(study == "BCCP") %>%
                distinct(value, area_id),
              by = "value", relationship = "many-to-many")

  add_data <- add_data %>%
    filter(doi != "BCCP/YaTsie") %>%
    bind_rows(bcpp) %>%
    mutate(idx = row_number() + max(aligned$idx))

  aligned <- bind_rows(aligned, add_data)
  if (any(is.na(aligned$area_id))) stop("Unmatched areas")

  # Figure 1: additional data sources (observations before matching). Rows of
  # new_data.csv for a study already in the review data (SABSSM V) are not a new study.
  # new_data.csv studies (primary trial reports for Joshi et al. secondary analyses, grey
  # literature survey reports and other published reports) form one box, "Additional
  # studies identified"; additional_sources gives the breakdown by report type
  new_data_source <- read_csv(file.path(path_raw, "review", "new_data.csv"), show_col_types = FALSE) %>%
    distinct(doi = str_trim(doi), report_source = source)
  primary_dois <- primary_report_dois()
  stopifnot(setequal(primary_dois, new_data_source$doi[new_data_source$report_source == "primary_report"]))
  additional_counts <- add_data %>%
    filter(!doi %in% dat$doi) %>%
    left_join(new_data_source, by = "doi") %>%
    mutate(source = case_when(str_detect(doi, "PHIA") ~ "PHIA surveys",
                              !is.na(report_source) ~ "Additional studies identified",
                              str_detect(doi, "Jones & Anderson") ~ "FSW: Jones and Anderson et al.",
                              str_detect(doi, "Stannah") ~ "MSM: Stannah et al.",
                              TRUE ~ doi)) %>%
    group_by(source) %>%
    summarise(studies = n_distinct(doi), observations = n_distinct(idx), .groups = "drop")
  additional_sources <- add_data %>%
    filter(!doi %in% dat$doi) %>%
    inner_join(new_data_source, by = "doi") %>%
    group_by(source = report_source) %>%
    summarise(studies = n_distinct(doi), observations = n_distinct(idx), .groups = "drop")

  list(aligned = aligned, tracker = review$tracker, additional_counts = additional_counts,
       additional_sources = additional_sources, primary_report_dois = primary_dois,
       review_max_idx = max(dat$idx),
       inclusion = inclusion, searches = searches)
}
