
### 01. MONTHLY BLEACHING ####
#### A. Primary model: warming-season GLMM ####

# Tests treatment differences in substantial bleaching and whether
# the treatment effect changes between the 2025 and 2026 seasons.
# Random intercepts account for repeated corals, frames and survey dates.
bleach_warming_model <- glmmTMB(
  bleached ~ treatment * warming_season + treatment * genus +
    (1 | coral_id) + (1 | frame) + (1 | survey_date),
  family = binomial,
  data = bleach_warming
)

summary(bleach_warming_model)


#### B. Sensitivity: full monitoring-period GLMM ####

# Tests whether treatment effects remain consistent when all 16 monthly
# surveys are included, rather than restricting to warming seasons.
bleach_monthly_model <- glmmTMB(
  bleached ~ treatment * period + treatment * genus +
    (1 | coral_id) + (1 | frame) + (1 | survey_date),
  family = binomial,
  data = bleach_monthly
)

summary(bleach_monthly_model)


#### C. Sensitivity: alternative bleaching threshold ####

# Tests whether conclusions depend on defining bleaching as substantial
# (scores 3–5) rather than any visible bleaching (scores 1–5).
# Uses the same March–July observations as the primary model.
bleach_any_model <- glmmTMB(
  bleached_any ~ treatment * warming_season + treatment * genus +
    (1 | coral_id) + (1 | frame) + (1 | survey_date),
  family = binomial,
  data = bleach_warming
)

summary(bleach_any_model)


#### D. Model diagnostics ####

# All models use the same random-effects structure.
# Seed = 42 is defined in 00_SETUP.R.
monthly_models <- list(
  "Primary warming-season GLMM" = bleach_warming_model,
  "Full-period sensitivity" = bleach_monthly_model,
  "Any-bleaching sensitivity" = bleach_any_model
)

for (model_name in names(monthly_models)) {
  model <- monthly_models[[model_name]]
  simulated_res <- DHARMa::simulateResiduals(model, n = 1000)
  
  cat("\n===", model_name, "===\n")
  print(tibble(
    observations = nobs(model),
    AIC = AIC(model),
    convergence = model$fit$convergence,
    positive_definite_hessian = model$sdr$pdHess
  ))
  print(DHARMa::testUniformity(simulated_res))
  print(DHARMa::testDispersion(simulated_res))
  print(DHARMa::testOutliers(simulated_res))
}

rm(monthly_models, model_name, model, simulated_res)


#### E. Save model outputs ####

##### Model fits ####

saveRDS(bleach_warming_model, file.path(fits_dir, "bleach_warming_model.rds"))
saveRDS(bleach_monthly_model, file.path(fits_dir, "bleach_monthly_model.rds"))
saveRDS(bleach_any_model, file.path(fits_dir, "bleach_any_model.rds"))

##### Fixed-effect summaries ####

bind_rows(
  purrr::imap(
    list(
      warming_primary = bleach_warming_model,
      full_period_sensitivity = bleach_monthly_model,
      any_bleaching_sensitivity = bleach_any_model
    ),
    ~ as.data.frame(summary(.x)$coefficients$cond) %>%
      tibble::rownames_to_column("term") %>%
      mutate(model = .y, .before = 1)
  )
) %>%
  write_csv(file.path(summaries_dir, "bleaching_model_summaries.csv"))

### 02. WEEKLY BLEACHING 2026 ####
#### A. Primary model: linear GLMM ####

# Tests treatment differences in substantial bleaching and whether
# those differences change over the 11 weekly surveys.
# Random intercepts account for repeated corals, frames and survey dates.
bleach_weekly_linear <- glmmTMB(
  bleached ~ treatment * time_c + treatment * genus +
    (1 | coral_id) + (1 | frame) + (1 | survey_date),
  family = binomial,
  data = weekly_model_data
)

summary(bleach_weekly_linear)


#### B. Model selection: quadratic GLMM ####

# A quadratic time structure was evaluated because observed bleaching
# initially peaked before declining. However, additional curvature did
# not improve fit (LRT p = 0.143; AIC 2015.51 vs 2015.40).
# The linear model was retained for parsimony.
bleach_weekly_quadratic <- glmmTMB(
  bleached ~ treatment * (time_c + I(time_c^2)) + treatment * genus +
    (1 | coral_id) + (1 | frame) + (1 | survey_date),
  family = binomial,
  data = weekly_model_data
)

anova(bleach_weekly_linear, bleach_weekly_quadratic)
AIC(bleach_weekly_linear, bleach_weekly_quadratic)


#### C. Exploratory model: GAMM ####

# Nonlinear treatment-specific trajectories were also explored.
# Both smooths had EDF ~1, indicating approximately linear logit-scale
# trends. GAMM residuals showed underdispersion (~0.70), requiring
# caution; the simpler GLMM was preferred.
bleach_weekly_gamm <- mgcv::bam(
  bleached ~ treatment * genus +
    s(time_weeks, by = treatment, k = 5) +
    s(coral_id, bs = "re") +
    s(frame, bs = "re") +
    s(survey_factor, bs = "re"),
  family = binomial(link = "logit"),
  method = "fREML",
  discrete = TRUE,
  data = weekly_model_data
)

summary(bleach_weekly_gamm)
mgcv::k.check(bleach_weekly_gamm)


#### D. Primary model diagnostics ####

# Linear GLMM converged with a positive-definite Hessian.
# Previous diagnostics: uniformity p = 0.459; dispersion = 1.030
# (p = 0.226); outliers p = 0.500.
tibble(
  observations = nobs(bleach_weekly_linear),
  AIC = AIC(bleach_weekly_linear),
  convergence = bleach_weekly_linear$fit$convergence,
  positive_definite_hessian = bleach_weekly_linear$sdr$pdHess
) %>%
  print()

weekly_linear_res <- DHARMa::simulateResiduals(bleach_weekly_linear, n = 1000)

DHARMa::testUniformity(weekly_linear_res)
DHARMa::testDispersion(weekly_linear_res)
DHARMa::testOutliers(weekly_linear_res)


#### E. Save selected model outputs ####

# Alternative models are retained in code for reproducibility,
saveRDS(bleach_weekly_linear, file.path(fits_dir, "bleach_weekly_linear.rds"))

as.data.frame(summary(bleach_weekly_linear)$coefficients$cond) %>%
  tibble::rownames_to_column("term") %>%
  write_csv(file.path(summaries_dir, "weekly_bleaching_model_summary.csv"))

rm(bleach_weekly_quadratic, bleach_weekly_gamm, weekly_linear_res)




### 03. SURVIVAL AND MORTALITY ####
#### A. Primary model: overall survival GLMM ####

# Tests survival probability at day 504 by treatment and genus.
# Frame random intercept accounts for clustering of transplanted corals.
# Unlike Cox regression, this model does not require proportional hazards
# Because Cox Assumption of proportional hazard was violated)
survival_model <- glmmTMB(
  survived ~ treatment * genus + (1 | frame),
  family = binomial,
  data = survival_data
)

summary(survival_model)


#### B. Secondary model: non-attachment mortality GLMM ####

# Tests cumulative mortality from all causes except attachment failure.
# Attachment deaths remain in the denominator as competing outcomes.
# This estimates observed mortality risk, not cause-specific hazards.
mortality_other_model <- glmmTMB(
  other_mortality ~ treatment * genus + (1 | frame),
  family = binomial,
  data = survival_data
)

summary(mortality_other_model)


#### C. Descriptive survival and mortality analyses ####

##### Restricted mean survival time ####

# All surviving corals were observed at day 504, so RMST equals
# mean recorded survival duration capped at 504 days.
# Death dates reflect first detection rather than exact mortality timing.
rmst_group <- survival_data %>%
  group_by(treatment, genus) %>%
  summarise(
    corals = n(),
    deaths = sum(event),
    rmst_days = mean(survival_days),
    .groups = "drop"
  )

print(rmst_group)

##### Kaplan–Meier survival estimates ####

survival_km <- survival::survfit(
  survival::Surv(survival_days, event) ~ treatment + genus,
  data = survival_data
)

summary(survival_km, times = c(180, 365, 504))

##### Mortality cause distribution ####

# Each percentage uses all originally transplanted corals in its group.
mortality_causes <- survival_data %>%
  count(treatment, genus, death_cause, name = "corals") %>%
  group_by(treatment, genus) %>%
  mutate(percent_of_group = round(100 * corals / sum(corals), 1)) %>%
  ungroup()

print(mortality_causes, n = Inf)

##### Cause-specific cumulative incidence ####

# Attachment and non-attachment mortality are competing events.
# Aalen–Johansen estimates their cumulative incidence over time.
survival_cif <- survival::survfit(
  survival::Surv(survival_days, cause_group) ~ treatment + genus,
  data = survival_data
)

summary(survival_cif, times = c(180, 365, 504))


#### D. Model diagnostics ####

##### Cox proportional hazards assessment ####

# Conventional Cox regression was evaluated but not selected:
# proportional hazards were violated (global p = 0.00023).
survival_cox <- survival::coxph(
  survival::Surv(survival_days, event) ~ treatment * genus + cluster(frame),
  data = survival_data,
  ties = "efron",
  x = TRUE
)

survival_ph <- survival::cox.zph(survival_cox)
print(survival_ph)

##### Binomial GLMM diagnostics ####

survival_models <- list(
  "Overall survival" = survival_model,
  "Non-attachment mortality" = mortality_other_model
)

for (model_name in names(survival_models)) {
  model <- survival_models[[model_name]]
  simulated_res <- DHARMa::simulateResiduals(model, n = 1000)
  
  cat("\n===", model_name, "===\n")
  print(tibble(
    observations = nobs(model),
    AIC = AIC(model),
    convergence = model$fit$convergence,
    positive_definite_hessian = model$sdr$pdHess
  ))
  print(DHARMa::testUniformity(simulated_res))
  print(DHARMa::testDispersion(simulated_res))
  print(DHARMa::testOutliers(simulated_res))
}

rm(survival_models, model_name, model, simulated_res,
   survival_cox, survival_ph)


#### E. Save model outputs ####

##### Selected model fits ####

saveRDS(survival_model, file.path(fits_dir, "survival_model.rds"))
saveRDS(mortality_other_model, file.path(fits_dir, "mortality_other_model.rds"))
saveRDS(survival_km, file.path(fits_dir, "survival_km.rds"))
saveRDS(survival_cif, file.path(fits_dir, "survival_cif.rds"))

##### Fixed-effect summaries ####

bind_rows(
  purrr::imap(
    list(
      overall_survival = survival_model,
      non_attachment_mortality = mortality_other_model
    ),
    ~ as.data.frame(summary(.x)$coefficients$cond) %>%
      tibble::rownames_to_column("term") %>%
      mutate(model = .y, .before = 1)
  )
) %>%
  write_csv(file.path(summaries_dir, "survival_mortality_model_summaries.csv"))

##### Descriptive summaries ####

write_csv(rmst_group, file.path(summaries_dir, "survival_rmst_groups.csv"))
write_csv(mortality_causes, file.path(summaries_dir, "mortality_causes.csv"))



### 04. POSITIVE SIZE EXTENSION ####
#### A. Primary model ####

growth_coral_gamma <- glmmTMB(
  mean_positive_extension ~ treatment * genus + (1 | frame),
  family = Gamma(link = "log"),
  data = growth_coral
)

summary(growth_coral_gamma)

stopifnot(
  growth_coral_gamma$fit$convergence == 0,
  growth_coral_gamma$sdr$pdHess
)


#### B. Model selection: Gaussian versus Gamma ####

# Same coral-level response; AIC values directly comparable.
# Gaussian residuals showed nonuniformity and excess outliers.
# Gamma provided better distributional fit.

growth_coral_gaussian <- glmmTMB(
  mean_positive_extension ~ treatment * genus + (1 | frame),
  family = gaussian,
  data = growth_coral
)

AIC(growth_coral_gaussian, growth_coral_gamma)

growth_gaussian_res <- DHARMa::simulateResiduals(growth_coral_gaussian, n = 1000)

DHARMa::testUniformity(growth_gaussian_res)
DHARMa::testDispersion(growth_gaussian_res)
DHARMa::testOutliers(growth_gaussian_res)

rm(growth_coral_gaussian, growth_gaussian_res)

#### C. Model diagnostics ####

growth_coral_res <- DHARMa::simulateResiduals(growth_coral_gamma, n = 1000)

plot(growth_coral_res)
DHARMa::testUniformity(growth_coral_res)
DHARMa::testDispersion(growth_coral_res)
DHARMa::testOutliers(growth_coral_res)

#### E. Model outputs ####

saveRDS(growth_coral_gamma,
        file.path(fits_dir, paste0("growth_coral_gamma_", analysis_date, ".rds")))

summary(growth_coral_gamma)$coefficients$cond %>%
  as.data.frame() %>%
  tibble::rownames_to_column("term") %>%
  write_csv(file.path(summaries_dir,
                      paste0("growth_fixed_effects_", analysis_date, ".csv")))



