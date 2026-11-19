### 01. MONTHLY BLEACHING ####

#### Monitoring periods and warming seasons ####

# March–July 2025 and 2026 provide matched seasonal windows.
# Full monitoring periods are retained for sensitivity analysis.
bleach_monthly <- coral_monthly %>%
  filter(!is.na(bleached)) %>%
  mutate(
    period = factor(
      if_else(survey_date < as.Date("2025-12-01"),
              "First monitoring period", "Second monitoring period"),
      levels = c("First monitoring period", "Second monitoring period")
    ),
    warming_season = factor(
      case_when(
        between(survey_date, as.Date("2025-03-01"), as.Date("2025-07-31")) ~ "2025 warming",
        between(survey_date, as.Date("2026-03-01"), as.Date("2026-07-31")) ~ "2026 warming"
      ),
      levels = c("2025 warming", "2026 warming")
    )
  ) %>%
  rename(treatment = origin) %>%
  select(coral_id, frame, treatment, genus, survey_date,
         period, warming_season, bleaching, bleached, bleached_any)


### 02. WEEKLY BLEACHING ####

#### Binary bleaching and local thermal conditions ####

bleach_weekly <- coral_weekly %>%
  filter(!is.na(bleached)) %>%
  left_join(
    temp_daily %>%
      select(date, temp_mean, hotspot, dhw_local, dhw_complete),
    by = c("survey_date" = "date")
  ) %>%
  rename(treatment = origin) %>%
  select(coral_id, frame, treatment, genus, survey_date,
         bleaching, bleached, bleached_any,
         temp_mean, hotspot, dhw_local, dhw_complete)


#### Consecutive weekly bleaching transitions ####

# Compare adjacent survey rounds only; retain complete assessment pairs.
weekly_transitions <- coral_weekly %>%
  arrange(coral_id, survey_index) %>%
  group_by(coral_id) %>%
  mutate(
    previous_bleached = lag(bleached),
    previous_index = lag(survey_index)
  ) %>%
  ungroup() %>%
  filter(
    survey_index - previous_index == 1,
    !is.na(previous_bleached),
    !is.na(bleached)
  ) %>%
  rename(treatment = origin) %>%
  select(coral_id, frame, treatment, genus, survey_date,
         previous_bleached, bleached)

stopifnot(
  all(weekly_transitions$previous_bleached %in% 0:1),
  all(weekly_transitions$bleached %in% 0:1)
)

### 03. SURVIVAL ####

#### Time to first observed complete mortality ####

# Event dates represent first detection; survivors are censored at study end.
survival_analysis <- coral_mortality %>%
  left_join(
    coral_clean %>% distinct(coral_id, baseline_date),
    by = "coral_id"
  ) %>%
  mutate(
    endpoint_date = coalesce(death_date, study_end),
    survival_days = as.numeric(endpoint_date - baseline_date)
  ) %>%
  rename(treatment = origin) %>%
  select(coral_id, frame, treatment, genus, baseline_date,
         endpoint_date, survival_days, event, death_cause)


### 04. POSITIVE SIZE EXTENSION ####

#### Consecutive monthly measurements ####

# Negative changes are excluded, not converted to zero.
# Missing measurements are not interpolated or bridged.
growth_analysis <- coral_monthly %>%
  filter(survey_date > as.Date("2025-03-28"),
         survey_date <= as.Date("2026-03-27")) %>%
  select(coral_id, frame, origin, genus, survey_date, size_cm) %>%
  bind_rows(
    coral_clean %>%
      distinct(coral_id, frame, origin, genus, baseline_date, baseline_size_cm) %>%
      transmute(coral_id, frame, origin, genus,
                survey_date = baseline_date, size_cm = baseline_size_cm)
  ) %>%
  arrange(coral_id, survey_date) %>%
  group_by(coral_id) %>%
  mutate(
    previous_date = lag(survey_date),
    previous_size = lag(size_cm),
    days_elapsed = as.numeric(survey_date - previous_date),
    size_change = size_cm - previous_size
  ) %>%
  ungroup() %>%
  filter(!is.na(size_change), days_elapsed > 0, size_change > 0) %>%
  mutate(extension_cm_month = size_change / days_elapsed * 30.44) %>%
  rename(treatment = origin) %>%
  select(coral_id, frame, treatment, genus, previous_date, survey_date,
         days_elapsed, size_change, extension_cm_month)



### 05. MODEL-READY DATASETS ####

#### Matched warming-season bleaching ####

# Five monthly surveys per warming season, 2025 and 2026.
bleach_warming <- bleach_monthly %>%
  filter(!is.na(warming_season)) %>%
  mutate(warming_season = droplevels(warming_season))

#### Weekly bleaching with elapsed time ####

# Centre elapsed weeks for interpretable treatment comparisons.
weekly_model_data <- bleach_weekly %>%
  mutate(
    coral_id = factor(coral_id),
    time_weeks = as.numeric(survey_date - min(survey_date)) / 7,
    time_c = time_weeks - mean(unique(time_weeks)),
    survey_factor = factor(survey_date)
  )

#### Survival and non-attachment mortality ####

# Attachment failure (METH) separated from other mortality causes.
# Attachment deaths remain competing outcomes for non-attachment mortality.
survival_data <- survival_analysis %>%
  mutate(
    survived = 1L - event,
    other_mortality = as.integer(event == 1 & death_cause != "METH"),
    cause_group = factor(
      case_when(
        event == 0 ~ "Censored",
        death_cause == "METH" ~ "Attachment",
        TRUE ~ "Other mortality"
      ),
      levels = c("Censored", "Attachment", "Other mortality")
    )
  )

#### Coral-level positive extension ####

# One mean positive, monthly-standardised extension rate per coral.
# Retain measurement counts to document unequal observation histories.
growth_coral <- growth_analysis %>%
  group_by(coral_id, frame, treatment, genus) %>%
  summarise(
    mean_positive_extension = mean(extension_cm_month),
    measurements = n(),
    .groups = "drop"
  )

### 06. ANALYSIS DATA VALIDATION ####

#### Structural and response validation ####

stopifnot(
  nrow(bleach_monthly) == 8986,
  nrow(bleach_weekly) == 4493,
  all(bleach_monthly$bleached == as.integer(bleach_monthly$bleaching >= 3)),
  all(bleach_weekly$bleached == as.integer(bleach_weekly$bleaching >= 3)),
  !anyNA(bleach_monthly$bleached),
  !anyNA(bleach_weekly$bleached),
  nrow(survival_analysis) == 800,
  !anyNA(survival_analysis$survival_days),
  all(survival_analysis$survival_days > 0),
  nrow(growth_analysis) == 4352,
  all(growth_analysis$extension_cm_month > 0),
  all(growth_analysis$days_elapsed > 0),
  !anyNA(growth_analysis$extension_cm_month)
)

#### Model-ready dataset validation ####

stopifnot(
  nrow(bleach_warming) == 5618,
  n_distinct(bleach_warming$survey_date) == 10,
  nrow(weekly_model_data) == 4493,
  all(is.finite(weekly_model_data$time_c)),
  nrow(survival_data) == 800,
  sum(survival_data$event) == 439,
  all(survival_data$survival_days[survival_data$event == 0] == 504),
  !anyNA(survival_data$other_mortality),
  !anyNA(survival_data$cause_group),
  nrow(growth_coral) == 726,
  n_distinct(growth_coral$coral_id) == 726,
  all(growth_coral$mean_positive_extension > 0),
  all(is.finite(growth_coral$mean_positive_extension))
)

#### Warming-season survey coverage ####

bleach_monthly %>%
  filter(!is.na(warming_season)) %>%
  distinct(survey_date, warming_season) %>%
  count(warming_season, name = "survey_dates") %>%
  print()

stopifnot(
  n_distinct(bleach_monthly$survey_date[!is.na(bleach_monthly$warming_season)]) == 10,
  nrow(filter(bleach_monthly, !is.na(warming_season))) == 5618
)



#### Weekly survey temperatures ####

bleach_weekly %>%
  distinct(survey_date, temp_mean, hotspot, dhw_local, dhw_complete) %>%
  arrange(survey_date) %>%
  print(n = Inf)



### 06. SAVE ANALYSIS DATA ####

write_csv(bleach_monthly, file.path(data_clean_dir, "bleach_monthly_analysis.csv"))
write_csv(bleach_weekly, file.path(data_clean_dir, "bleach_weekly_analysis.csv"))
write_csv(survival_analysis, file.path(data_clean_dir, "survival_analysis.csv"))
write_csv(growth_analysis, file.path(data_clean_dir, "growth_analysis.csv"))
write_csv(growth_coral, file.path(data_clean_dir, "growth_coral_analysis.csv"))

cat("\n=== ANALYSIS PREPARATION COMPLETE: 5 DATASETS VALIDATED AND SAVED ===\n")
