################################################################################
#' Copy the Box/OneDrive inputs the paper pipeline needs into private/
#'
#' private/ is gitignored. Run once from any working directory:
#'   source("~/Documents/GitHub/ssa-incidence-paper/prep/copy_private_inputs.R")
#' Each step is independent; a failed step is reported and the rest still run.
#'
#' The paths below are the authors' copies and are not public. Other users can
#' obtain the same inputs:
#'   - Spectrum model results (PJNZ files) can be requested from UNAIDS for
#'     research purposes: https://hivtools.unaids.org/spectrum-file-request/
#'   - Naomi subnational results for this estimates cycle can be downloaded from
#'     https://naomi-spectrum.unaids.org/
#'
#' Extracts derived from the Spectrum files are written to data/ and committed.
#' The Namibia 2026 file comes from "Non public Spectrum files", but its estimates
#' are on AIDSinfo and we have permission to release derived values. The Côte
#' d'Ivoire 2025 file ("2025 final shared") is public.
################################################################################

library(tidyverse)
library(sf)
library(naomi)

paper <- "~/Documents/GitHub/ssa-incidence-paper"
out <- file.path(paper, "private")
dir.create(out, showWarnings = FALSE, recursive = TRUE)

box <- "~/Library/CloudStorage/Box-Box/HIV Inference Group/Data/Spectrum files"
onedrive <- "~/OneDrive - Imperial College London/SSA incidence"

step <- function(name, expr) {
  message("== ", name)
  tryCatch(force(expr), error = function(e) message("FAILED ", name, ": ", conditionMessage(e)))
}

#' Naomi district outputs (areas + indicators), same vintages as ssa_incidence R/clean.R get_naomi()
get_naomi <- function(iso3) {
  message(iso3)
  if (iso3 %in% c("CAF", "SEN")) {
    path <- grep("zip", list.files(file.path(box, "2023 naomi final"), pattern = iso3, full.names = TRUE), value = TRUE)
    out <- read_output_package(path)
    areas <- out$meta_area
    ind <- filter(out$indicators, calendar_quarter == "CY2022Q4")
  } else if (iso3 == "BEN") {
    areas <- read_sf(file.path(box, "2024 naomi final/data-inputs/raw/BEN/6bef0c177daac4c1770e71e4814affc2/files/99770883A55DE9214123149A20A0F267.geojson"))
    ind <- read_output_package(file.path(box, "2024 naomi final/output-zip/BEN 2024 naomi_outputs.zip"))$indicators %>%
      filter(calendar_quarter == "CY2023Q4")
  } else if (iso3 %in% c("CIV", "ZWE")) {
    out <- read_output_package(list.files(file.path(box, "2025 naomi final/output-zips"), full.names = TRUE, pattern = iso3))
    areas <- out$meta_area
    ind <- filter(out$indicators, calendar_quarter == "CY2024Q4")
  } else {
    out <- read_output_package(list.files(file.path(box, "2026 naomi final/output-zips"), full.names = TRUE, pattern = iso3))
    areas <- out$meta_area
    ind <- filter(out$indicators, calendar_quarter == "CY2025Q4")
  }
  list(areas = areas,
       indicators = ind %>%
         filter(indicator %in% c("population", "plhiv", "prevalence", "incidence", "infections")) %>%
         mutate(year = calendar_quarter_to_year(calendar_quarter)) %>%
         select(area_level:age_group_label, year, indicator, mean))
}

naomi_iso3 <- c("RWA", "GNB", "KEN", "ZWE", "ZAF", "UGA", "NGA", "TZA", "SWZ", "MOZ", "ZMB", "MWI", "LSO",
                "CIV", "NAM", "BWA", "CMR", "LBR", "BFA", "ETH", "BEN", "SEN", "CAF", "MLI", "TGO")

step("Naomi areas and indicators", {
  naomi_res <- setNames(lapply(naomi_iso3, get_naomi), naomi_iso3)
  areas <- bind_rows(lapply(naomi_res, "[[", "areas"), .id = "iso3") %>% st_make_valid()
  saveRDS(areas, file.path(out, "naomi_areas.rds"))
  saveRDS(bind_rows(lapply(naomi_res, "[[", "indicators"), .id = "iso3"), file.path(out, "naomi_indicators.rds"))
})

step("Kenya 2023 areas (SEARCH community matching)", {
  path <- file.path(box, "2023 naomi final", "Subnational data inputs", "KEN")
  ids <- readRDS(file.path(path, "project_ids.rds"))
  inputs <- readRDS(file.path(path, ids$fit_id, "data.rds"))
  saveRDS(read_sf(file.path(path, ids$fit_id, "files", inputs$objects$data$shape$path)), file.path(out, "ken_areas_2023.rds"))
})

#' Subnational EPP fits from the Spectrum files, as R/funs.R calculate_spectrum_ratios()
spec_paths <- c(
  list.files(file.path(box, "2026 preliminary restricted/Public Spectrum files/mainland SSA"), recursive = TRUE, full.names = TRUE),
  file.path(box, "2026 final restricted/Non public Spectrum files/Namibia_050626_Final.pjnz"),
  file.path(box, "2025 final shared/WCA/cotedivoire_2025_Abidjan_13032025.pjnz")
)
no_phia_paths <- list.files(file.path(onedrive, "No PHIA incidence Spectrum"), full.names = TRUE)

read_subnational_epp <- function(paths) {
  lapply(paths, function(x) {
    message(basename(x))
    tryCatch(
      epp::read_spt(x) %>%
        lapply(function(d) mutate(d, year = attr(d, "row.names"))) %>%
        bind_rows(.id = "area_name") %>%
        mutate(iso3 = countrycode::countrycode(eppasm::read_country(x), "country.name", "iso3c"),
               file = basename(x)),
      error = function(e) NULL)
  }) %>%
    bind_rows()
}

step("Subnational EPP (published Spectrum files)", {
  saveRDS(read_subnational_epp(spec_paths), file.path(out, "subnational_epp.rds"))
  write_csv(tibble(file = basename(spec_paths)), file.path(out, "spectrum_files_used.csv"))
})

step("Subnational EPP (Spectrum refit without PHIA incidence)", {
  saveRDS(read_subnational_epp(no_phia_paths), file.path(out, "subnational_epp_no_phia.rds"))
  write_csv(tibble(file = basename(no_phia_paths)), file.path(out, "no_phia_files.csv"))
})

step("National Spectrum, refit without PHIA incidence", {
  saveRDS(bind_rows(lapply(no_phia_paths, naomi::extract_pjnz_naomi)), file.path(out, "spectrum_no_phia_national.rds"))
})

step("OneDrive supplementary figure/table files (reference copies)", {
  dir.create(file.path(out, "reference"), showWarnings = FALSE)
  figs <- file.path(onedrive, "Figs")
  file.copy(list.files(figs, full.names = TRUE, recursive = TRUE), file.path(out, "reference"), overwrite = TRUE)
})

step("match_level.csv from Downloads (read by clean_joshi)", {
  file.copy("~/Downloads/match_level.csv", file.path(out, "match_level_downloads.csv"), overwrite = TRUE)
})

message("Done. Contents of ", out, ":")
print(list.files(out, recursive = TRUE))
