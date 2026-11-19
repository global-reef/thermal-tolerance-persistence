
### 01. EXPERIMENTAL DESIGN AND SAMPLING ####

# Formal exploratory analysis following Zuur et al. (2010):
# response distributions, outliers, heterogeneity, independence,
# predictor relationships and sampling structure.
# Model comparisons and diagnostics remain in 04_MODEL.R.

#### Analytical sample sizes ####

tibble(
  analysis = c("Monthly bleaching", "Weekly bleaching",
               "Survival", "Positive extension"),
  observations = c(nrow(bleach_monthly), nrow(bleach_weekly),
                   nrow(survival_analysis), nrow(growth_analysis)),
  corals = c(n_distinct(bleach_monthly$coral_id),
             n_distinct(bleach_weekly$coral_id),
             n_distinct(survival_analysis$coral_id),
             n_distinct(growth_analysis$coral_id)),
  frames = c(n_distinct(bleach_monthly$frame),
             n_distinct(bleach_weekly$frame),
             n_distinct(survival_analysis$frame),
             n_distinct(growth_analysis$frame)),
  surveys = c(n_distinct(bleach_monthly$survey_date),
              n_distinct(bleach_weekly$survey_date),
              NA_integer_, n_distinct(growth_analysis$survey_date))
) %>%
  print()

#### Treatment and genus balance ####

survival_analysis %>%
  count(treatment, genus, name = "corals") %>%
  print()

survival_analysis %>%
  count(frame, treatment, genus, name = "corals") %>%
  print(n = Inf)

#### Repeated-measurement structure ####

bind_rows(
  bleach_monthly %>%
    count(coral_id, name = "observations") %>%
    mutate(analysis = "Monthly bleaching"),
  bleach_weekly %>%
    count(coral_id, name = "observations") %>%
    mutate(analysis = "Weekly bleaching"),
  growth_analysis %>%
    count(coral_id, name = "observations") %>%
    mutate(analysis = "Positive extension")
) %>%
  group_by(analysis) %>%
  summarise(
    corals = n(),
    median_observations = median(observations),
    min_observations = min(observations),
    max_observations = max(observations),
    .groups = "drop"
  ) %>%
  print()

# Treatment assigned at frame level; frames nested within treatment.
# Repeated assessments within corals and shared survey dates
# require appropriate dependence structures in longitudinal models.


### 02. MONTHLY BLEACHING ####

#### Bleaching score distribution ####

# Upper severity classes are sparse; evaluate binary response suitability.
# Scores 3–5 correspond to substantial bleaching (50–100% classes).

bleach_monthly %>%
  filter(!is.na(warming_season)) %>%
  count(bleaching) %>%
  mutate(percent = 100 * n / sum(n)) %>%
  print(n = Inf)

bleach_monthly %>%
  filter(!is.na(warming_season)) %>%
  count(warming_season, treatment, genus, bleaching) %>%
  group_by(warming_season, treatment, genus) %>%
  mutate(percent = 100 * n / sum(n)) %>%
  ungroup() %>%
  print(n = Inf)

bleach_monthly %>%
  filter(!is.na(warming_season)) %>%
  ggplot(aes(factor(bleaching), fill = treatment)) +
  geom_bar(position = "fill") +
  facet_grid(genus ~ warming_season) +
  scale_fill_manual(values = treatment_palette) +
  scale_y_continuous(labels = scales::percent) +
  labs(x = "Bleaching score (0–5)", y = "Proportion") +
  theme_clean

#### Binary response and event representation ####

bleach_monthly %>%
  filter(!is.na(warming_season)) %>%
  group_by(warming_season, treatment, genus) %>%
  summarise(
    assessments = n(),
    any_bleaching = sum(bleached_any),
    substantial = sum(bleached),
    substantial_pct = 100 * mean(bleached),
    affected_corals = n_distinct(coral_id[bleached == 1]),
    affected_frames = n_distinct(frame[bleached == 1]),
    .groups = "drop"
  ) %>%
  print(n = Inf)

#### Monitoring coverage and seasonal balance ####

bleach_monthly %>%
  group_by(survey_date, warming_season, treatment, genus) %>%
  summarise(
    assessments = n(),
    corals = n_distinct(coral_id),
    substantial_pct = 100 * mean(bleached),
    .groups = "drop"
  ) %>%
  print(n = Inf)

#### Temporal bleaching patterns ####

bleach_monthly %>%
  group_by(survey_date, treatment, genus) %>%
  summarise(
    prevalence = mean(bleached),
    assessments = n(),
    .groups = "drop"
  ) %>%
  ggplot(aes(survey_date, prevalence, colour = treatment)) +
  geom_line() +
  geom_point(aes(size = assessments)) +
  facet_wrap(~genus) +
  scale_colour_manual(values = treatment_palette) +
  scale_y_continuous(labels = scales::percent) +
  labs(x = "Survey date", y = "Substantial bleaching prevalence",
       size = "Assessments") +
  theme_clean

#### Frame-level heterogeneity ####

bleach_monthly %>%
  filter(!is.na(warming_season)) %>%
  group_by(frame, treatment, warming_season) %>%
  summarise(
    assessments = n(),
    prevalence = mean(bleached),
    .groups = "drop"
  ) %>%
  ggplot(aes(treatment, prevalence, colour = treatment)) +
  geom_jitter(width = 0.12, height = 0) +
  facet_wrap(~warming_season) +
  scale_colour_manual(values = treatment_palette) +
  scale_y_continuous(labels = scales::percent) +
  labs(x = "Treatment", y = "Frame-level bleaching prevalence") +
  theme_clean

# Model rationale:
# Substantial bleaching retained as primary binomial response.
# Any bleaching retained as threshold sensitivity.
# Matched March–July windows support comparison between warming years.
# Coral, frame and survey-date intercepts account for dependence.


### 03. WEEKLY BLEACHING ####

#### Response distribution and event sparsity ####

bleach_weekly %>%
  count(treatment, genus, bleaching) %>%
  group_by(treatment, genus) %>%
  mutate(percent = 100 * n / sum(n)) %>%
  ungroup() %>%
  print(n = Inf)

bleach_weekly %>%
  group_by(treatment, genus) %>%
  summarise(
    assessments = n(),
    substantial = sum(bleached),
    prevalence_pct = 100 * mean(bleached),
    affected_corals = n_distinct(coral_id[bleached == 1]),
    affected_frames = n_distinct(frame[bleached == 1]),
    .groups = "drop"
  ) %>%
  print()

#### Temporal prevalence and trajectory shape ####

weekly_prevalence <- bleach_weekly %>%
  group_by(survey_date, treatment, genus) %>%
  summarise(
    assessments = n(),
    substantial = sum(bleached),
    prevalence = mean(bleached),
    .groups = "drop"
  )

print(weekly_prevalence, n = Inf)

weekly_prevalence %>%
  ggplot(aes(survey_date, prevalence, colour = treatment)) +
  geom_line() +
  geom_point(aes(size = assessments)) +
  facet_wrap(~genus) +
  scale_colour_manual(values = treatment_palette) +
  scale_y_continuous(labels = scales::percent) +
  labs(x = "Survey date", y = "Substantial bleaching prevalence",
       size = "Assessments") +
  theme_clean

#### Within-coral dependence ####

bleach_weekly %>%
  arrange(coral_id, survey_date) %>%
  group_by(coral_id) %>%
  mutate(
    previous_bleaching = lag(bleached),
    previous_date = lag(survey_date)
  ) %>%
  ungroup() %>%
  filter(!is.na(previous_bleaching)) %>%
  group_by(treatment, genus) %>%
  summarise(
    paired_observations = n(),
    lag1_correlation = cor(previous_bleaching, bleached),
    .groups = "drop"
  ) %>%
  print()

#### Frame-level heterogeneity ####

bleach_weekly %>%
  group_by(frame, treatment) %>%
  summarise(
    assessments = n(),
    prevalence = mean(bleached),
    .groups = "drop"
  ) %>%
  print(n = Inf)

#### Thermal predictor coverage ####

weekly_environment <- bleach_weekly %>%
  distinct(survey_date, temp_mean, hotspot, dhw_local, dhw_complete) %>%
  arrange(survey_date) %>%
  mutate(days_since_start = as.numeric(survey_date - min(survey_date)))

print(weekly_environment, n = Inf)

weekly_environment %>%
  summarise(
    surveys = n(),
    temperature_available = sum(!is.na(temp_mean)),
    dhw_available = sum(!is.na(dhw_local)),
    complete_dhw = sum(dhw_complete, na.rm = TRUE)
  ) %>%
  print()

#### Temporal and thermal predictor associations ####

# Date-level correlations avoid pseudoreplication from shared exposure.
weekly_environment %>%
  select(days_since_start, temp_mean, hotspot, dhw_local) %>%
  cor(use = "pairwise.complete.obs", method = "spearman") %>%
  round(3) %>%
  print()

# Model rationale:
# Temporal trajectories evaluated using linear, quadratic and GAMM fits.
# Linear GLMM retained following model comparisons and diagnostics.
# Incomplete thermal exposure limits temperature-based inference.

rm(weekly_prevalence, weekly_environment)


### 04. SURVIVAL AND MORTALITY ####

#### Endpoint survival ####

survival_analysis %>%
  group_by(treatment, genus) %>%
  summarise(
    corals = n(),
    alive = sum(event == 0),
    deaths = sum(event == 1),
    survival_pct = 100 * mean(event == 0),
    frames = n_distinct(frame),
    .groups = "drop"
  ) %>%
  print()

#### Mortality cause distribution ####

survival_analysis %>%
  filter(event == 1) %>%
  count(treatment, genus, death_cause) %>%
  group_by(treatment, genus) %>%
  mutate(percent_of_deaths = 100 * n / sum(n)) %>%
  ungroup() %>%
  print(n = Inf)

#### Cause-specific incidence at study endpoint ####

survival_analysis %>%
  mutate(
    outcome = case_when(
      event == 0 ~ "Alive",
      death_cause == "METH" ~ "Attachment",
      TRUE ~ "Other mortality"
    )
  ) %>%
  count(treatment, genus, outcome) %>%
  group_by(treatment, genus) %>%
  mutate(percent_of_group = 100 * n / sum(n)) %>%
  ungroup() %>%
  print(n = Inf)

#### Event timing and censoring ####

survival_analysis %>%
  group_by(treatment, genus, event) %>%
  summarise(
    corals = n(),
    median_days = median(survival_days),
    min_days = min(survival_days),
    max_days = max(survival_days),
    .groups = "drop"
  ) %>%
  print(n = Inf)

survival_analysis %>%
  filter(event == 1) %>%
  count(endpoint_date, name = "deaths") %>%
  arrange(endpoint_date) %>%
  mutate(cumulative_deaths = cumsum(deaths)) %>%
  print(n = Inf)

#### Frame-level mortality variation ####

survival_analysis %>%
  group_by(frame, treatment) %>%
  summarise(
    corals = n(),
    mortality_pct = 100 * mean(event == 1),
    .groups = "drop"
  ) %>%
  print(n = Inf)

# Model rationale:
# One survival endpoint per coral; no repeated-response coral intercept.
# Attachment failure (METH) separated from other mortality.
# Frame clustering retained.
# Cox PH assumption assessed in 04_MODEL.R; endpoint GLMMs selected.
# KM, CIF and RMST retained for descriptive time-to-event summaries.


### 05. POSITIVE SIZE EXTENSION ####

#### Size measurement coverage ####

coral_monthly %>%
  group_by(survey_date) %>%
  summarise(measured = sum(!is.na(size_cm)), .groups = "drop") %>%
  left_join(
    growth_analysis %>%
      group_by(survey_date) %>%
      summarise(
        positive_extensions = n(),
        median_extension = median(extension_cm_month),
        .groups = "drop"
      ),
    by = "survey_date"
  ) %>%
  mutate(month = format(survey_date, "%B %Y")) %>%
  select(month, measured, positive_extensions, median_extension) %>%
  print(n = Inf)

#### Measurement intervals ####

growth_analysis %>%
  count(days_elapsed, sort = TRUE) %>%
  print(n = Inf)

# Extension already standardised to cm/month in 02_PREP.R.
# Missing measurements not interpolated; negative changes excluded.

#### Interval-level distribution ####

growth_analysis %>%
  group_by(treatment, genus) %>%
  summarise(
    observations = n(),
    corals = n_distinct(coral_id),
    mean_extension = mean(extension_cm_month),
    median_extension = median(extension_cm_month),
    sd_extension = sd(extension_cm_month),
    p99 = quantile(extension_cm_month, 0.99),
    maximum = max(extension_cm_month),
    .groups = "drop"
  ) %>%
  print()

growth_analysis %>%
  ggplot(aes(extension_cm_month, fill = treatment)) +
  geom_histogram(bins = 40, alpha = 0.5, position = "identity") +
  facet_wrap(~genus, scales = "free_y") +
  scale_fill_manual(values = treatment_palette) +
  labs(x = "Positive extension (cm/month)", y = "Observations") +
  theme_clean

growth_analysis %>%
  ggplot(aes(log(extension_cm_month), fill = treatment)) +
  geom_histogram(bins = 40, alpha = 0.5, position = "identity") +
  facet_wrap(~genus, scales = "free_y") +
  scale_fill_manual(values = treatment_palette) +
  labs(x = "Log positive extension", y = "Observations") +
  theme_clean

#### Temporal variation ####

growth_analysis %>%
  group_by(survey_date, treatment, genus) %>%
  summarise(
    mean_extension = mean(extension_cm_month),
    median_extension = median(extension_cm_month),
    observations = n(),
    .groups = "drop"
  ) %>%
  ggplot(aes(survey_date, mean_extension, colour = treatment)) +
  geom_line() +
  geom_point(aes(size = observations)) +
  facet_wrap(~genus) +
  scale_colour_manual(values = treatment_palette) +
  labs(x = "Survey date", y = "Mean positive extension (cm/month)",
       size = "Observations") +
  theme_clean

#### Coral-level response structure ####

# Aggregation considered to reduce sensitivity to interval-level noise
# and provide one positive-extension outcome per coral.

growth_coral_explore <- growth_coral

growth_coral_explore %>%
  group_by(treatment, genus) %>%
  summarise(
    corals = n(),
    mean_extension = mean(mean_positive_extension),
    median_extension = median(mean_positive_extension),
    median_measurements = median(measurements),
    min_measurements = min(measurements),
    max_measurements = max(measurements),
    single_measurement_corals = sum(measurements == 1),
    .groups = "drop"
  ) %>%
  print()

#### Coral-level distribution ####

growth_coral_explore %>%
  ggplot(aes(mean_positive_extension, fill = treatment)) +
  geom_histogram(bins = 30, alpha = 0.5, position = "identity") +
  facet_wrap(~genus, scales = "free_y") +
  scale_fill_manual(values = treatment_palette) +
  labs(x = "Mean positive extension (cm/month)", y = "Corals") +
  theme_clean

growth_coral_explore %>%
  ggplot(aes(treatment, mean_positive_extension, fill = treatment)) +
  geom_boxplot(outlier.alpha = 0.3) +
  facet_wrap(~genus) +
  scale_fill_manual(values = treatment_palette) +
  labs(x = "Treatment", y = "Mean positive extension (cm/month)") +
  theme_clean

# Model rationale:
# Raw and log extension distributions examined.
# Strong temporal heterogeneity; interval-level GLMMs and GAMMs explored.
# Coral-level Gaussian and Gamma alternatives compared in 04_MODEL.R.
# Gamma retained for distributional fit and interpretable treatment effects.
# Unequal positive-extension histories and survivorship remain limitations.

rm(growth_coral_explore)


### 06. LOCAL TEMPERATURE AND THERMAL STRESS ####

#### Annual environmental summary ####

temp_daily %>%
  mutate(year = lubridate::year(date)) %>%
  group_by(year) %>%
  summarise(
    days_observed = sum(complete_day, na.rm = TRUE),
    mean_temp = mean(temp_mean, na.rm = TRUE),
    max_daily_mean = max(temp_mean, na.rm = TRUE),
    peak_dhw = if (any(!is.na(dhw_local)))
      max(dhw_local, na.rm = TRUE) else NA_real_,
    .groups = "drop"
  ) %>%
  print()

#### Monthly coverage and thermal stress ####

temp_daily %>%
  mutate(month = lubridate::floor_date(date, "month")) %>%
  group_by(month) %>%
  summarise(
    days_observed = sum(complete_day, na.rm = TRUE),
    missing_days = sum(!complete_day),
    mean_temp = mean(temp_mean, na.rm = TRUE),
    peak_dhw = if (any(!is.na(dhw_local)))
      max(dhw_local, na.rm = TRUE) else NA_real_,
    .groups = "drop"
  ) %>%
  print(n = Inf)

#### Observed temperature trajectory ####

temp_daily %>%
  ggplot(aes(date, temp_mean)) +
  geom_line(na.rm = TRUE) +
  geom_hline(yintercept = c(mmm, mmm + 1), linetype = "dotted") +
  labs(x = "Date", y = "Daily mean temperature (°C)") +
  theme_clean

#### Investigating peaks ####
temp_daily %>%
  filter(complete_day, temp_mean >= 31) %>%
  select(date, temp_mean) %>%
  arrange(date) %>%
  print(n = Inf)

#### max readout #### 
temp_daily %>%
  filter(complete_day) %>%
  group_by(year = lubridate::year(date)) %>%
  summarise(
    max_daily_mean = max(temp_mean),
    max_recorded_temp = max(temp_max),
    .groups = "drop"
  )
#### Growth interval temperature ####

# Descriptive only; no temperature-growth model selected.
# Exposure shared across corals measured within each interval.

growth_analysis %>%
  distinct(previous_date, survey_date) %>%
  rowwise() %>%
  mutate(
    temp_readings = sum(temp_clean$date > previous_date &
                          temp_clean$date <= survey_date &
                          !is.na(temp_clean$temperature)),
    mean_temp = mean(temp_clean$temperature[
      temp_clean$date > previous_date &
        temp_clean$date <= survey_date
    ], na.rm = TRUE)
  ) %>%
  ungroup() %>%
  arrange(survey_date) %>%
  print(n = Inf)

# Environmental context retained for interpretation of warming periods.
# No interpolated temperatures or DHW values used for inference.
