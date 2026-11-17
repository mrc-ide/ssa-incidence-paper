################################################################################
#' Build the area lookups in data/geography/ from the Naomi area boundaries
#'
#' Needs private/ (prep/copy_private_inputs.R). Run from the repository root.
#' Writes:
#'   areas.csv               area hierarchy (no boundaries) for the 25 countries
#'                           with data, and the Spectrum area each area's
#'                           modelled incidence trend comes from
#'   geocoded_locations.csv  study-site search strings -> area_id (Google Maps
#'                           API results, Supplementary Text S3, joined to the
#'                           lowest-level Naomi areas)
#'   search_areas.csv        SEARCH trial districts -> area_id
################################################################################

library(tidyverse)
library(sf)

source("R/functions.R")
source("R/data_functions.R")

dir.create(path_geography, showWarnings = FALSE, recursive = TRUE)

areas_sf <- readRDS(file.path(path_private, "naomi_areas.rds")) %>% st_make_valid()

#' Level of the Spectrum (or Thembisa) estimates whose trends are applied to each
#' country's areas: 0 national, otherwise EPP subpopulation / province level
spectrum_level <- tribble(
  ~iso3, ~spectrum_level,
  "BEN", 0, "BFA", 0, "BWA", 0, "CIV", 0, "CMR", 0, "GNB", 0, "KEN", 1, "LSO", 0, "MOZ", 1, "MWI", 0,
  "NAM", 1, "NGA", 0, "RWA", 0, "SEN", 0, "SWZ", 1, "TZA", 2, "UGA", 0, "ZAF", 1, "ZMB", 0, "ZWE", 1,
  "CAF", 0, "MLI", 0, "LBR", 0, "TGO", 0, "ETH", 0
)

areas <- areas_sf %>%
  st_drop_geometry() %>%
  select(iso3, area_id, area_name, area_level, parent_area_id, spectrum_region_code) %>%
  left_join(spectrum_level, by = "iso3")

#' Ancestor at the Spectrum level; areas at or above it map to themselves.
#' (Eswatini: every area maps to itself.)
ancestor_at <- function(area_id, level, areas) {
  parent <- setNames(areas$parent_area_id, areas$area_id)
  lvl <- setNames(areas$area_level, areas$area_id)
  map2_chr(area_id, level, function(a, l) {
    while (!is.na(a) && lvl[[a]] > l) a <- parent[[a]]
    a
  })
}

areas <- areas %>%
  mutate(spectrum_area_id = ifelse(iso3 == "SWZ", area_id, ancestor_at(area_id, spectrum_level, areas)))

write_csv(areas, file.path(path_geography, "areas.csv"), na = "")

################################################################################
#' ## Geocoded study sites
################################################################################

inclusion <- read_inclusion()
match_level <- read_csv(file.path(path_raw, "review", "match_level.csv"), show_col_types = FALSE)
review <- clean_review_data(read_review_data(), inclusion, match_level)

#' SEARCH districts: Ugandan districts by name (any level), Kenyan sub-counties
#' from the 2023 boundaries assigned to their county; Suba is now part of Homa Bay
ken_2023 <- readRDS(file.path(path_private, "ken_areas_2023.rds")) %>% st_drop_geometry()
ken_low <- ken_2023 %>%
  filter(area_level == 3) %>%
  select(area_name, area_id = parent_area_id)

search_areas <- read_csv(file.path(path_raw, "additional", "search.csv"), show_col_types = FALSE) %>%
  filter(intervention == 0) %>%
  distinct(district) %>%
  left_join(bind_rows(areas %>% filter(iso3 == "UGA") %>% select(area_name, area_id), ken_low),
            by = c("district" = "area_name"), relationship = "many-to-many") %>%
  mutate(area_id = ifelse(district == "Suba", "KEN_2_43", area_id))

write_csv(search_areas, file.path(path_geography, "search_areas.csv"))

add_data <- read_additional_data(inclusion, search_areas)
searches <- make_geocode_searches(review$dat, add_data)

#' Google Maps API results, one list element per search (cached by ssa_incidence)
gps <- readRDS(file.path(path_private, "matched_gps.rds")) %>%
  bind_rows(.id = "google_idx") %>%
  mutate(google_idx = as.integer(google_idx))

stopifnot(setequal(gps$google_search, searches$value))

#' Keep hits that resolved to an address; drop hits checked against the source
#' as being in the wrong district
hits <- gps %>%
  filter(!is.na(value), str_detect(value, ",")) %>%
  filter(!value %in% c("Itende, Manyoni, Singida Region, Tanzania",
                       "Bugando, Bukombe, Geita Region, Tanzania",
                       "Bugando, Kwimba, Mwanza REgion, Tanzania"),
         !str_detect(value, "Lugazi")) %>%
  st_as_sf()

lowest_areas <- areas_sf %>%
  group_by(iso3) %>%
  filter(area_level == max(area_level)) %>%
  ungroup() %>%
  select(area_id)

geocoded <- hits %>%
  st_join(lowest_areas) %>%
  st_drop_geometry() %>%
  filter(!is.na(area_id)) %>%
  distinct(value = google_search, google_result = value, area_id)

#' Searches whose hits fall in more than one district (kept, as in the analysis)
geocoded %>%
  distinct(value, area_id) %>%
  count(value) %>%
  filter(n > 1) %>%
  print()

write_csv(geocoded, file.path(path_geography, "geocoded_locations.csv"))
