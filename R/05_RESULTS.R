### 01. MONTHLY BLEACHING RESULTS ####

#### A. Primary warming-season comparisons ####

## Overall treatment probabilities 

genus_weights <- c(ACRO = 632 / 800, POR = 168 / 800)

bleach_warming_overall <- emmeans(
  bleach_warming_model,
  ~ treatment * genus | warming_season
) %>%
  regrid(transform = "response") %>%
  emmeans(~ treatment | warming_season, weights = genus_weights)

summary(bleach_warming_overall, infer = c(TRUE, FALSE))

##  Overall shallow versus deep 

# Differences are on the probability scale.
summary(
  pairs(bleach_warming_overall, reverse = TRUE),
  infer = TRUE
)

## Genus-specific treatment contrasts

bleach_warming_emm <- emmeans(
  bleach_warming_model,
  ~ treatment | warming_season * genus
)

summary(bleach_warming_emm, type = "response", infer = c(TRUE, FALSE))

summary(
  pairs(bleach_warming_emm, reverse = TRUE),
  type = "response", infer = TRUE
)

## Change in treatment effect between years

# Change in shallow/deep odds ratio, identical across genera
# under the selected fixed-effect structure.
emmeans(bleach_warming_model, ~ treatment * warming_season) %>%
  contrast(interaction = c("revpairwise", "revpairwise")) %>%
  summary(infer = TRUE)


#### B. Population-marginal probabilities: parametric bootstrap ####

# Integrate over coral, frame and survey-date random intercepts.
# Refit simulated responses to incorporate fixed-effect and variance-component uncertainty.

marginal_bleaching <- function(model) {
  re_variance <- sum(vapply(
    VarCorr(model)$cond,
    function(x) as.numeric(x[1, 1]), numeric(1)
  ))
  
  marginal_prob <- function(eta) {
    integrate(
      function(z) plogis(eta + sqrt(re_variance) * z) * dnorm(z),
      lower = -Inf, upper = Inf
    )$value
  }
  
  emmeans(model, ~ treatment | warming_season * genus) %>%
    summary(type = "link") %>%
    as_tibble() %>%
    transmute(
      treatment, warming_season, genus,
      probability = vapply(emmean, marginal_prob, numeric(1))
    )
}

## Original model estimates 

bleach_marginal_est <- marginal_bleaching(bleach_warming_model)

bleach_marginal_diff <- bleach_marginal_est %>%
  pivot_wider(names_from = treatment, values_from = probability) %>%
  mutate(difference = shallow - deep) %>%
  select(warming_season, genus, difference)

##### Bootstrap refitting ####

# Simulate the original model's complete random-effects structure.
# Failed/non-positive-definite fits are excluded and counted.

n_boot <- 500
bleach_boot <- vector("list", n_boot)

for (b in seq_len(n_boot)) {
  boot_data <- bleach_warming
  boot_data$bleached <- simulate(bleach_warming_model, nsim = 1)[[1]]
  
  boot_fit <- tryCatch(
    suppressWarnings(update(bleach_warming_model, data = boot_data)),
    error = function(e) NULL
  )
  
  if (is.null(boot_fit)) next
  if (boot_fit$fit$convergence != 0 || !isTRUE(boot_fit$sdr$pdHess)) next
  
  bleach_boot[[b]] <- tryCatch(
    marginal_bleaching(boot_fit) %>% mutate(replicate = b),
    error = function(e) NULL
  )
}

bleach_boot_draws <- bind_rows(bleach_boot)

cat("Successful bootstrap replicates:",
    n_distinct(bleach_boot_draws$replicate), "/", n_boot, "\n")

stopifnot(n_distinct(bleach_boot_draws$replicate) >= 475)

## Bootstrap confidence intervals

bleach_marginal_ci <- bleach_boot_draws %>%
  group_by(treatment, warming_season, genus) %>%
  summarise(
    lower = quantile(probability, 0.025),
    upper = quantile(probability, 0.975),
    .groups = "drop"
  )

bleach_marginal_emms <- bleach_marginal_est %>%
  left_join(bleach_marginal_ci,
            by = c("treatment", "warming_season", "genus"))

bleach_marginal_diff_ci <- bleach_boot_draws %>%
  select(replicate, treatment, warming_season, genus, probability) %>%
  pivot_wider(names_from = treatment, values_from = probability) %>%
  mutate(difference = shallow - deep) %>%
  group_by(warming_season, genus) %>%
  summarise(
    lower = quantile(difference, 0.025),
    upper = quantile(difference, 0.975),
    .groups = "drop"
  )

bleach_marginal_contrasts <- bleach_marginal_diff %>%
  left_join(bleach_marginal_diff_ci,
            by = c("warming_season", "genus")) %>%
  mutate(across(c(difference, lower, upper), ~ .x * 100))

print(bleach_marginal_emms, n = Inf)
print(bleach_marginal_contrasts, n = Inf)

## Save marginal estimates 

monthly_stats_dir <- file.path(stats_dir, "01_monthly_bleaching")
dir.create(monthly_stats_dir, recursive = TRUE, showWarnings = FALSE)

write_csv(bleach_marginal_emms,
          file.path(monthly_stats_dir, "primary_marginal_probabilities.csv"))

write_csv(bleach_marginal_contrasts,
          file.path(monthly_stats_dir, "primary_marginal_differences.csv"))

# Compact bootstrap record for reproducing confidence intervals.
saveRDS(bleach_boot_draws,
        file.path(monthly_stats_dir, "primary_marginal_bootstrap.rds"))

#### C. Table 1: primary population-marginal estimates ####

## Genus-specific marginal estimates and contrasts

bleach_table <- bleach_marginal_emms %>%
  mutate(
    year = stringr::str_extract(as.character(warming_season), "\\d{4}"),
    estimate = sprintf("%.1f (%.1f–%.1f)",
                       100 * probability, 100 * lower, 100 * upper)
  ) %>%
  select(genus, year, treatment, estimate) %>%
  pivot_wider(names_from = treatment, values_from = estimate) %>%
  left_join(
    bleach_marginal_contrasts %>%
      mutate(
        year = stringr::str_extract(as.character(warming_season), "\\d{4}"),
        contrast = sprintf("%.1f (%.1f–%.1f)", difference, lower, upper)
      ) %>%
      select(genus, year, contrast),
    by = c("genus", "year")
  ) %>%
  arrange(genus, year) %>%
  select(genus, year, deep, shallow, contrast)

## Overall population-marginal estimates 

# Standardise across genera using the original planting proportions.
overall_est <- bleach_marginal_emms %>%
  mutate(weight = genus_weights[as.character(genus)]) %>%
  group_by(warming_season, treatment) %>%
  summarise(probability = weighted.mean(probability, weight), .groups = "drop")

overall_boot <- bleach_boot_draws %>%
  mutate(weight = genus_weights[as.character(genus)]) %>%
  group_by(replicate, warming_season, treatment) %>%
  summarise(probability = weighted.mean(probability, weight), .groups = "drop")

overall_ci <- overall_boot %>%
  group_by(warming_season, treatment) %>%
  summarise(
    lower = quantile(probability, 0.025),
    upper = quantile(probability, 0.975),
    .groups = "drop"
  )

overall_contrasts <- overall_boot %>%
  select(replicate, warming_season, treatment, probability) %>%
  pivot_wider(names_from = treatment, values_from = probability) %>%
  mutate(difference = shallow - deep) %>%
  group_by(warming_season) %>%
  summarise(
    lower = quantile(difference, 0.025),
    upper = quantile(difference, 0.975),
    .groups = "drop"
  )

overall_table <- overall_est %>%
  left_join(overall_ci, by = c("warming_season", "treatment")) %>%
  mutate(
    year = stringr::str_extract(as.character(warming_season), "\\d{4}"),
    estimate = sprintf("%.1f (%.1f–%.1f)",
                       100 * probability, 100 * lower, 100 * upper)
  ) %>%
  select(year, treatment, estimate) %>%
  pivot_wider(names_from = treatment, values_from = estimate) %>%
  left_join(
    overall_est %>%
      select(warming_season, treatment, probability) %>%
      pivot_wider(names_from = treatment, values_from = probability) %>%
      mutate(
        year = stringr::str_extract(as.character(warming_season), "\\d{4}"),
        difference = shallow - deep
      ) %>%
      left_join(overall_contrasts, by = "warming_season") %>%
      transmute(
        year,
        contrast = sprintf("%.1f (%.1f–%.1f)",
                           100 * difference, 100 * lower, 100 * upper)
      ),
    by = "year"
  ) %>%
  mutate(genus = "Overall", .before = 1)

## Final manuscript table 

bleach_table <- bleach_table %>%
  mutate(genus = as.character(genus)) %>%
  bind_rows(overall_table) %>%
  mutate(genus = factor(genus, levels = c("ACRO", "POR", "Overall"))) %>%
  arrange(genus, year)

print(bleach_table, n = Inf)


write_csv(
  bleach_table,
  file.path(tables_dir, "T1_monthly_bleaching_primary_table.csv")
)

rm(boot_data, boot_fit, bleach_boot, bleach_boot_draws, bleach_marginal_ci,
   bleach_marginal_diff, bleach_marginal_diff_ci, marginal_bleaching,
   overall_est, overall_boot, overall_ci, overall_contrasts, overall_table,
   monthly_stats_dir, n_boot, b)

#### D. Sensitivity: full monitoring period ####

## Overall treatment probabilities 

bleach_monthly_overall <- emmeans(
  bleach_monthly_model,
  ~ treatment * genus | period
) %>%
  regrid(transform = "response") %>%
  emmeans(~ treatment | period, weights = genus_weights)

summary(bleach_monthly_overall, infer = c(TRUE, FALSE))

summary(pairs(bleach_monthly_overall, reverse = TRUE),infer = TRUE)

## Genus-specific treatment contrasts

bleach_monthly_emm <- emmeans(
  bleach_monthly_model,
  ~ treatment | period * genus)

summary(pairs(bleach_monthly_emm, reverse = TRUE),type = "response", infer = TRUE)

## Change in treatment effect between periods

emmeans(bleach_monthly_model, ~ treatment * period) %>%
  contrast(interaction = c("revpairwise", "revpairwise")) %>%
  summary(infer = TRUE)


#### E. Sensitivity: any visible bleaching ###

## Overall treatment comparison 

bleach_any_overall <- emmeans(
  bleach_any_model,
  ~ treatment * genus | warming_season
) %>%
  regrid(transform = "response") %>%
  emmeans(~ treatment | warming_season, weights = genus_weights)

summary(bleach_any_overall, infer = c(TRUE, FALSE))

summary(
  pairs(bleach_any_overall, reverse = TRUE),
  infer = TRUE
)

## Genus-specific treatment contrasts

bleach_any_emm <- emmeans(
  bleach_any_model,
  ~ treatment | warming_season * genus
)

summary(
  pairs(bleach_any_emm, reverse = TRUE),
  type = "response", infer = TRUE
)


### 02. WEEKLY BLEACHING RESULTS ####

#### A. Overall treatment effects ####

# Standardise across genera using original planting proportions.
# Predictions are evaluated at centred time (time_c = 0).
bleach_weekly_overall <- emmeans(
  bleach_weekly_linear,
  ~ treatment * genus,
  at = list(time_c = 0)
) %>%
  regrid(transform = "response") %>%
  emmeans(~ treatment, weights = genus_weights)

summary(bleach_weekly_overall, infer = c(TRUE, FALSE))

summary(
  pairs(bleach_weekly_overall, reverse = TRUE),
  infer = TRUE
)


#### B. Genus-specific treatment effects ####

# Shallow/deep odds ratios at centred time.
bleach_weekly_emm <- emmeans(
  bleach_weekly_linear,
  ~ treatment | genus,
  at = list(time_c = 0)
)

summary(bleach_weekly_emm, type = "response", infer = c(TRUE, FALSE))

summary(
  pairs(bleach_weekly_emm, reverse = TRUE),
  type = "response", infer = TRUE
)


#### C. Weekly temporal trajectories ####

# Predicted probabilities across the observed weekly survey dates.
# Confidence intervals are population-level, with random effects set to zero.
weekly_predictions <- emmeans(
  bleach_weekly_linear,
  ~ treatment * genus | time_c,
  at = list(time_c = sort(unique(weekly_model_data$time_c)))
) %>%
  summary(type = "response", infer = c(TRUE, FALSE)) %>%
  as_tibble() %>%
  mutate(
    survey_date = min(weekly_model_data$survey_date) +
      round((time_c + mean(unique(weekly_model_data$time_weeks))) * 7)
  )

print(weekly_predictions, n = Inf)


#### D. Change in treatment difference over time ####

# The treatment × time interaction tests whether log-odds slopes differ.
# Genus-specific treatment contrasts vary with time, but the interaction
# is common across genera under the selected model specification.
summary(bleach_weekly_linear)$coefficients$cond[
  "treatmentshallow:time_c", , drop = FALSE
] %>%
  print()

#### E. Weekly bleaching transitions ####

# Conditional proportions of onset and recovery between adjacent surveys.
# Descriptive only; repeated transitions within coral are not independent.

weekly_transition_summary <- weekly_transitions %>%
  group_by(genus, treatment, previous_bleached) %>%
  summarise(
    transitions = n(),
    changed = sum(bleached != previous_bleached),
    proportion = changed / transitions,
    .groups = "drop"
  ) %>%
  mutate(
    transition = if_else(
      previous_bleached == 0,
      "New substantial bleaching",
      "Recovery below threshold"
    ),
    transition = factor(
      transition,
      levels = c("New substantial bleaching", "Recovery below threshold")
    )
  )

print(weekly_transition_summary, n = Inf) 
### 03. SURVIVAL AND MORTALITY RESULTS ####

#### Overall treatment effects ####

survival_overall <- emmeans(
  survival_model, ~ treatment * genus
) %>%
  regrid(transform = "response") %>%
  emmeans(~ treatment, weights = genus_weights)

mortality_overall <- emmeans(
  mortality_other_model, ~ treatment * genus
) %>%
  regrid(transform = "response") %>%
  emmeans(~ treatment, weights = genus_weights)

summary(pairs(survival_overall, reverse = TRUE), infer = TRUE)
summary(pairs(mortality_overall, reverse = TRUE), infer = TRUE)

#### Genus-specific treatment effects ####

summary(
  pairs(emmeans(survival_model, ~ treatment | genus), reverse = TRUE),
  type = "response", infer = TRUE
)

summary(
  pairs(emmeans(mortality_other_model, ~ treatment | genus), reverse = TRUE),
  type = "response", infer = TRUE
)




### 04. POSITIVE SIZE EXTENSION RESULTS ####

#### A. Overall treatment comparison ####

# Standardise to initial genus composition (632 ACRO, 168 POR).
# Regridding gives response-scale means and absolute differences (cm/month).
growth_overall <- emmeans(growth_coral_gamma, ~ treatment * genus) %>%
  regrid(transform = "response") %>%
  emmeans(~ treatment, weights = genus_weights)

summary(growth_overall, infer = c(TRUE, FALSE))

growth_overall_contrast <- pairs(growth_overall, reverse = TRUE) %>%
  summary(infer = TRUE) %>%
  as_tibble()

print(growth_overall_contrast)

#### B. Treatment comparisons within genus ####

# Ratios are shallow/deep on the response scale.
growth_by_genus <- emmeans(growth_coral_gamma, ~ treatment | genus)

summary(growth_by_genus, type = "response", infer = c(TRUE, FALSE))

growth_genus_contrasts <- pairs(growth_by_genus, reverse = TRUE) %>%
  summary(type = "response", infer = TRUE) %>%
  as_tibble()

print(growth_genus_contrasts)

#### C. Treatment × genus interaction ####

summary(growth_coral_gamma)$coefficients$cond[
  "treatmentshallow:genusPOR", , drop = FALSE
]


### 05. SAVE MANUSCRIPT STATISTICS ####

#### Output directories ####

stats_sections <- c(
  monthly = "01_monthly_bleaching",
  weekly = "02_weekly_bleaching",
  survival = "03_survival_mortality",
  growth = "04_positive_extension"
)

purrr::walk(file.path(stats_dir, stats_sections),
            ~ dir.create(.x, recursive = TRUE, showWarnings = FALSE))

save_stats <- function(x, section, filename) {
  write_csv(as_tibble(x), file.path(stats_dir, stats_sections[[section]],
                                    paste0(filename, ".csv")))
}

#### A. Monthly bleaching: primary warming-season analysis ####

save_stats(summary(bleach_warming_overall, infer = c(TRUE, FALSE)),
           "monthly", "primary_overall_emms")

save_stats(summary(pairs(bleach_warming_overall, reverse = TRUE), infer = TRUE),
           "monthly", "primary_overall_contrasts")

save_stats(summary(bleach_warming_emm, type = "response", infer = c(TRUE, FALSE)),
           "monthly", "primary_genus_emms")

save_stats(summary(pairs(bleach_warming_emm, reverse = TRUE),
                   type = "response", infer = TRUE),
           "monthly", "primary_genus_contrasts")

save_stats(
  emmeans(bleach_warming_model, ~ treatment * warming_season) %>%
    contrast(interaction = c("revpairwise", "revpairwise")) %>%
    summary(infer = TRUE),
  "monthly", "primary_treatment_year_interaction"
)

#### B. Monthly bleaching: full-period sensitivity ####

save_stats(summary(bleach_monthly_overall, infer = c(TRUE, FALSE)),
           "monthly", "sensitivity_full_period_emms")

save_stats(summary(pairs(bleach_monthly_overall, reverse = TRUE), infer = TRUE),
           "monthly", "sensitivity_full_period_contrasts")

save_stats(summary(pairs(bleach_monthly_emm, reverse = TRUE),
                   type = "response", infer = TRUE),
           "monthly", "sensitivity_full_period_genus_contrasts")

save_stats(
  emmeans(bleach_monthly_model, ~ treatment * period) %>%
    contrast(interaction = c("revpairwise", "revpairwise")) %>%
    summary(infer = TRUE),
  "monthly", "sensitivity_full_period_interaction"
)

#### C. Monthly bleaching: any-bleaching sensitivity ####

save_stats(summary(bleach_any_overall, infer = c(TRUE, FALSE)),
           "monthly", "sensitivity_any_bleaching_emms")

save_stats(summary(pairs(bleach_any_overall, reverse = TRUE), infer = TRUE),
           "monthly", "sensitivity_any_bleaching_contrasts")

save_stats(summary(pairs(bleach_any_emm, reverse = TRUE),
                   type = "response", infer = TRUE),
           "monthly", "sensitivity_any_bleaching_genus_contrasts")


#### D. Weekly bleaching ####

save_stats(summary(bleach_weekly_overall, infer = c(TRUE, FALSE)),
           "weekly", "overall_emms")

save_stats(summary(pairs(bleach_weekly_overall, reverse = TRUE), infer = TRUE),
           "weekly", "overall_contrasts")

save_stats(summary(bleach_weekly_emm, type = "response", infer = c(TRUE, FALSE)),
           "weekly", "genus_emms")

save_stats(summary(pairs(bleach_weekly_emm, reverse = TRUE),
                   type = "response", infer = TRUE),
           "weekly", "genus_contrasts")
save_stats(weekly_transition_summary, "weekly", "bleaching_transitions")

# Temporal slopes on the log-odds scale.
save_stats(
  emtrends(bleach_weekly_linear, ~ treatment, var = "time_c") %>%
    summary(infer = TRUE),
  "weekly", "temporal_slopes"
)

save_stats(
  summary(bleach_weekly_linear)$coefficients$cond %>%
    as.data.frame() %>%
    tibble::rownames_to_column("term") %>%
    filter(term %in% c("treatmentshallow:time_c",
                       "treatmentshallow:genusPOR")),
  "weekly", "interaction_terms"
)

save_stats(weekly_predictions, "weekly", "temporal_predictions")


#### E. Survival and mortality ####

save_stats(summary(survival_overall, infer = c(TRUE, FALSE)),
           "survival", "survival_overall_emms")

save_stats(summary(pairs(survival_overall, reverse = TRUE), infer = TRUE),
           "survival", "survival_overall_contrasts")

save_stats(
  emmeans(survival_model, ~ treatment | genus) %>%
    summary(type = "response", infer = c(TRUE, FALSE)),
  "survival", "survival_genus_emms"
)

save_stats(
  pairs(emmeans(survival_model, ~ treatment | genus), reverse = TRUE) %>%
    summary(type = "response", infer = TRUE),
  "survival", "survival_genus_contrasts"
)

save_stats(summary(mortality_overall, infer = c(TRUE, FALSE)),
           "survival", "non_attachment_overall_emms")

save_stats(summary(pairs(mortality_overall, reverse = TRUE), infer = TRUE),
           "survival", "non_attachment_overall_contrasts")

save_stats(
  emmeans(mortality_other_model, ~ treatment | genus) %>%
    summary(type = "response", infer = c(TRUE, FALSE)),
  "survival", "non_attachment_genus_emms"
)

save_stats(
  pairs(emmeans(mortality_other_model, ~ treatment | genus), reverse = TRUE) %>%
    summary(type = "response", infer = TRUE),
  "survival", "non_attachment_genus_contrasts"
)

save_stats(rmst_group, "survival", "restricted_mean_survival")
save_stats(mortality_causes, "survival", "mortality_cause_distribution")

save_stats(
  survival_data %>%
    group_by(treatment, genus) %>%
    summarise(
      corals = n(),
      alive = sum(event == 0),
      deaths = sum(event == 1),
      survival_pct = 100 * mean(event == 0),
      .groups = "drop"
    ),
  "survival", "observed_endpoint_survival"
)


#### F. Positive size extension ####

save_stats(summary(growth_overall, infer = c(TRUE, FALSE)),
           "growth", "overall_emms")

save_stats(growth_overall_contrast, "growth", "overall_contrasts")

save_stats(summary(growth_by_genus, type = "response",
                   infer = c(TRUE, FALSE)),
           "growth", "genus_emms")

save_stats(growth_genus_contrasts, "growth", "genus_contrasts")

save_stats(
  summary(growth_coral_gamma)$coefficients$cond %>%
    as.data.frame() %>%
    tibble::rownames_to_column("term") %>%
    filter(term == "treatmentshallow:genusPOR"),
  "growth", "treatment_genus_interaction"
)

save_stats(
  growth_coral %>%
    group_by(treatment, genus) %>%
    summarise(
      corals = n(),
      median_measurements = median(measurements),
      mean_extension = mean(mean_positive_extension),
      median_extension = median(mean_positive_extension),
      .groups = "drop"
    ),
  "growth", "observed_extension_summary"
)

cat("\n=== MANUSCRIPT STATISTICS SAVED ===\n")
print(file.path(stats_dir, stats_sections))

rm(save_stats, stats_sections)
