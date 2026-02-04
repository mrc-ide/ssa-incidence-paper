library(tidyverse)
library(INLA)
library(patchwork)

sample_model <-  function (inla_mod, inla_df, col)
  {
    df <- inla_df %>% filter(across(all_of(col), ~!is.na(.x)))
    samples <- inla.posterior.sample(1000, inla_mod)
    contents = inla_mod$misc$configs$contents
    effect = "Predictor"
    id.effect = which(contents$tag == effect)
    ind.effect = contents$start[id.effect] - 1 + (1:contents$length[id.effect])
    ind.effect <- 1:(nrow(inla_df) - nrow(df))
    samples.effect = lapply(samples, function(x) x$latent[ind.effect])
    samples <- matrix(sapply(samples.effect, cbind), ncol = 1000)
    mean <- rowMeans(samples)
    qtls <- apply(samples, 1, quantile, c(0.025, 0.5, 0.975))
    ident <- inla_df[ind.effect, ]
    samples <- ident %>% cbind(samples)
    res <- ident %>% mutate(mean = mean, lower = qtls[1, ],
                            upper = qtls[3, ])
  }

final_data <- read_csv("data/submission_data.csv", show_col_types = F)

inla_dat <- crossing(
  year = 1990:2025,
  sex = c("female", "male"),
  recruitment_type = unique(final_data$recruitment_type),
  pys = 1
) %>%
  filter(!(sex == "male" & recruitment_type == "ANC"),
         !(sex == "male" & recruitment_type == "FSW"),
         !(sex == "female" & recruitment_type == "MSM")) %>%
  bind_rows(final_data) %>%
  filter(year %in% 1990:2025) %>%
  mutate(year2 = year,
         recruitment_type = factor(recruitment_type, levels = c("UTT/Cohort", "National HHS", "Other HHS", "ANC", "Intervention Trials",  "Non-representative - Other", "No inclusion criteria", "FSW", "MSM")),
         sex_int = case_when(
           sex == "female" ~ 0,
           sex == "both" ~ 0.5,
           sex == "male" ~ 1
         ),
         csi_dummy = as.numeric(study_type == "Cross-sectional incidence"),
         binary_sex = ifelse(sex == "female", 1, 2),
         std_year = year - 2015,
         id.year = group_indices(., year)) %>%
  ungroup() %>%
  mutate(id.year.recruitment = group_indices(., year, recruitment_type))

A1 <- matrix(as.numeric(unique(inla_dat$year) == 2015), ncol = length(1990:2025)) ## Forces 2015 to be 0
A2 <- matrix(scale(1990:2025 - 2015), nrow = 1) ## Forces the slope in 2015 to be 0 [identify the linear time trend too]
A <- rbind(A1, A2)

R_period <- dfertility::make_rw_structure_matrix(36, 2) ## Makes 34x34 RW2 structure
R_recruitment <- as(diag(1, 9), "sparseMatrix")

Q <- kronecker(R_recruitment, R_period)

e2 <- matrix(0, nrow = 2)


########
n_years <- 36
n_studies <- 9
n_interactions <- n_years * n_studies

# Slope and coefficient constraint for year 25
A_slope_coef <- matrix(0, nrow = 2, ncol = n_interactions)

# Slope constraint between year 24 and 25 (zero slope)
for (study in 1:n_studies) {
  A_slope_coef[1, (24 - 1) * n_studies + study] <- -1  # Year 24 for each study
  A_slope_coef[1, (25 - 1) * n_studies + study] <- 1   # Year 25 for each study
}

# Coefficient value for year 25 is zero for all studies
for (study in 1:n_studies) {
  A_slope_coef[2, (25 - 1) * n_studies + study] <- 1  # Year 25 for each study
}

# Sum-to-zero constraint within each year for study effects
A_sum_study <- matrix(0, nrow = n_years, ncol = n_interactions)

# Loop over each year to impose sum-to-zero constraint across the 5 study groups
for (year in 1:n_years) {
  A_sum_study[year, ((year - 1) * n_studies + 1):(year * n_studies)] <- 1
}

# Combine slope/coef constraints with sum-to-zero constraints
A_combined <- rbind(A_slope_coef, A_sum_study)

# Constraints vector (zeros for both slope, coefficient, and sum-to-zero)
e3 <- matrix(0, nrow(A_combined), nrow = 1)


#########
mod <- INLA::inla(infections ~ 0
                  + sex_int
                  + std_year
                  + recruitment_type
                  # + csi_dummy

                  + f(year2, model = "rw2", replicate = binary_sex, extraconstr = list(A = A, e = e2))
                  + f(id.year.recruitment, model = "generic0", Cmatrix = Q, extraconstr = list(A = A_combined, e = e3), rankdef = n_years + n_studies -1)


                  + f(study_idx, model = "iid", constr = T)
                  + offset(log(extrapolated_value)),

                  family = "xpoisson",
                  data = inla_dat,
                  E = pys,
                  control.compute=list(config = TRUE))


summary(mod)

res <- sample_model(mod, inla_dat, "infections") %>%
  left_join(
    final_data %>%
      filter(sex != "both",
             !(sex == "male" & recruitment_type == "ANC")) %>%
      group_by(sex, recruitment_type) %>%
      filter(year %in% range(year)) %>%
      distinct(recruitment_type, sex, year) %>%
      arrange(recruitment_type, sex) %>%
      mutate(nam = c('low', 'high')) %>%
      pivot_wider(names_from = nam, values_from = year) %>%
      mutate(year = list(low:high)) %>%
      select(-c(low, high)) %>%
      unnest(year) %>%
      mutate(has = T)
  ) %>%
  filter(sex != "both",
         has == T) %>%
  mutate(
    sex = str_to_sentence(sex),
    sex = factor(sex, levels = c("Female", "Male", "Both"))) %>%
  filter(!(sex == "Female" & recruitment_type == "MSM"),
         !(sex == "Male" & recruitment_type == "FSW"),
         recruitment_type != "No inclusion criteria",
         sex != "Both") %>%
  left_join(count(final_data, recruitment_type))

