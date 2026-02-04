library(tidyverse)
library(INLA)

## Age and time trends from cohort data

u25_dat <- final_data %>%
  filter(doi %in% c("RCCS", "AHRI", "Manicaland", "Kyamulibwa GPC", "Kisesa", "Karonga")) %>%
  bind_rows(final_data %>% filter(doi == "Kisumu HDSS") %>% mutate(doi = "KEMRI")) %>%
  filter(age_group %in% c("Y015_019", "Y020_024", "Y025_029", "Y030_034", "Y035_039", "Y035_044", "Y040_044", "Y045_049", "Y045_054")) %>%
  mutate(infections = ir * pys,
         age_group = ifelse(age_group %in% c("Y015_019", "Y020_024"), "Y015_024", "Y025_049")) %>%
  group_by(age_group, area_id, year, doi, sex) %>%
  summarise(infections = sum(infections),
            pys = sum(pys)) %>%
  ungroup() %>%
  mutate(ir = infections/pys)

prec.prior <- list(prec= list(prior = "normal", param = c(2,1)))

cohort_u25_inla <-
  crossing(year = unique(u25_dat$year),
           distinct(u25_dat, doi, area_id),
           sex = c("male", "female"),
           age_group = c("Y015_024", "Y025_049")) %>%
  mutate(pys = 1) %>%
  bind_rows(u25_dat) %>%
  mutate(year_min = year - 2015,
         year_min2 = year_min,
         year_min3 = year_min,
         year_min_f = ifelse(sex == "female", year_min, NA),
         year_min_m = ifelse(sex == "male", year_min, NA),
         year_min_f2 = year_min_f,
         year_min_m2 = year_min_m,
         doi2 = factor(doi),
         doi3 = doi,
         doi4 = doi,
         doi_int = as.numeric(doi2),
         sex_int = as.numeric(factor(sex)),
         age_int = as.numeric(factor(age_group)),
         age_group2 = fct_rev(factor(age_group)),
         age_int = as.numeric(age_group2)
  ) %>%
  mutate(interaction = group_indices(., doi, age_group)) %>%
  ungroup()


u25_t_mod_test <- inla(infections ~
                         year_min*sex
                       + age_group2*sex
                       + age_group2*year_min
                       + f(year_min2, model = "ar", order = 2, replicate = age_int, constr = T, hyper = prec.prior)
                       + f(year_min_f, model = "ar", order = 2, group = doi_int, control.group = list(model = "iid")) +
                         + f(year_min_m, model = "ar", order = 2, group = doi_int, control.group = list(model = "iid")) +
                         + f(doi3, replicate = age_int)
                       ,
                       data = cohort_u25_inla %>% mutate(age_group2 = factor(age_group2, levels = c("Y015_024", "Y025_049"))),
                       family = "xpoisson",
                       control.compute=list(config = TRUE),
                       E = pys)


cohort_u25_time_res <- function(u25_t_mod, cohort_u25_inla) {


  df <- cohort_u25_inla %>%
    filter(across(all_of("infections"), ~!is.na(.x)))

  samples <- inla.posterior.sample(1000, u25_t_mod)
  contents = u25_t_mod$misc$configs$contents
  effect = "Predictor"
  id.effect = which(contents$tag==effect)
  ind.effect = contents$start[id.effect]-1 + (1:contents$length[id.effect])

  ind.effect <- 1:(nrow(cohort_u25_inla) - nrow(df))

  samples.effect = lapply(samples, function(x) x$latent[ind.effect])

  u25_samples <- matrix(sapply(samples.effect, cbind), ncol=1000)
  ident <- cohort_u25_inla[ind.effect, ]

  u25_samples <- ident %>% cbind(exp(u25_samples))

  # u25_samples <- u25_samples %>%
  #   ungroup() %>%
  #   group_by(year, sex, doi) %>%
  #   summarise(across(as.character(1:1000), ~sum(.x[age_group == "Y015_024"])/(sum(.x)/sum(pys)))) %>%
  #   ungroup()

  mean <- rowMeans(u25_samples[,as.character(1:1000)])
  qtls <- apply(u25_samples[,as.character(1:1000)], 1, quantile, c(0.025, 0.5, 0.975))

  u25_time_res <-
    # u25_samples %>%
    # select(year, sex, doi) %>%
    ident %>%
    # distinct(year, doi) %>%
    mutate(
      lower = qtls[1,],
      mean = mean,
      upper = qtls[3,]
    )

}
pool_u25_time_res <- function(u25_t_mod, cohort_u25_inla_pool) {


  df <- cohort_u25_inla %>%
    filter(across(all_of("infections"), ~!is.na(.x)))

  samples <- inla.posterior.sample(1000, u25_t_mod)
  contents = u25_t_mod$misc$configs$contents
  effect = "Predictor"
  id.effect = which(contents$tag==effect)
  ind.effect = contents$start[id.effect]-1 + (1:contents$length[id.effect])

  ind.effect <- 1:(nrow(cohort_u25_inla_pool) - nrow(df))

  samples.effect = lapply(samples, function(x) x$latent[ind.effect])

  u25_samples <- matrix(sapply(samples.effect, cbind), ncol=1000)
  ident <- cohort_u25_inla_pool[ind.effect, ]

  u25_samples <- ident %>% cbind(exp(u25_samples))

  # u25_samples <- u25_samples %>%
  #   ungroup() %>%
  #   group_by(year, sex) %>%
  #   summarise(across(as.character(1:1000), ~sum(.x[age_group == "Y015_024"])/(sum(.x)/sum(pys)))) %>%
  #   ungroup()

  mean <- rowMeans(u25_samples[,as.character(1:1000)])
  qtls <- apply(u25_samples[,as.character(1:1000)], 1, quantile, c(0.025, 0.5, 0.975))

  u25_time_res <-
    # u25_samples %>%
    # select(year, age_group, sex) %>%
    ident %>%
    # distinct(year, doi) %>%
    mutate(
      lower = qtls[1,],
      mean = mean,
      upper = qtls[3,]
    )

}

u25_time_res <- cohort_u25_time_res(u25_t_mod_test, cohort_u25_inla) %>%
  mutate(
    sex = str_to_sentence(sex))


### PHIA age-time trends

u25_phia <- read_csv("inputs/phia-incid-extract_2025-11-13b.csv", show_col_types = F) %>%
  filter(age_group %in% c("Y015_024", "Y025_049"),
         sex != "both",
         # iso3 %in% c("ZWE", "ZMB", "SWZ", "MWI", "LSO", "UGA", "TZA")
  ) %>%
  mutate(infections = incid * n_pys_approx) %>%
  select(iso3, survey_id, year = survey_year, age_group, sex, incid, lower = ci_l.incid, upper = ci_u.incid, pys = n_pys_approx)

u25_d <- crossing(year = 2016:2022,
                  sex = c("male", "female"),
                  age_group = c("Y015_024", "Y025_049")) %>%
  mutate(pys = 1) %>%
  bind_rows(
    u25_phia
    # filter(age_group == "Y015_024")
  ) %>%
  mutate(infections = incid*pys,
         year_min = year - min(year),
         year_min2 = year_min,
         iso3 = factor(iso3),
         iso32 = iso3,
         age_group = fct_rev(age_group))


u25_phia_mod <- inla(infections ~
                       year_min*sex +
                       age_group*sex +
                       age_group*year_min +
                       f(iso32) +
                       f(iso3, year_min, model = "iid", hyper = prec.prior),
                     data = u25_d ,
                     family = "xpoisson",
                     control.compute=list(config = TRUE),
                     E = pys)

summary(u25_phia_mod)
