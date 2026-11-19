
### 01. IMPORTS AND CHECKS ####


#### Import data ####

monthly <- readr::read_csv(file.path(data_raw_dir, "20260727 data monthly.csv"), show_col_types = FALSE) %>%
  mutate(survey_month = lubridate::make_date(year, month, 1)) %>%
  arrange(fragment_ID, survey_month)

weekly <- readr::read_csv(file.path(data_raw_dir, "20260730 weekly data.csv"), show_col_types = FALSE) %>%
  mutate(
    survey_date = lubridate::dmy(date),
    survey_month = lubridate::floor_date(survey_date, "month")
  ) %>%
  arrange(fragment_ID, survey_date)

#### Check structure ####

glimpse(monthly)
glimpse(weekly)

#### Check sampling and identifiers ####

monthly %>% summarise(
  fragments = n_distinct(fragment_ID),
  frames = n_distinct(frame_number),
  surveys = n_distinct(survey_month),
  observations = n()
)

weekly %>% summarise(
  fragments = n_distinct(fragment_ID),
  frames = n_distinct(frame_number),
  surveys = n_distinct(survey_date),
  observations = n()
)

#### Check duplicates ####

monthly %>% count(fragment_ID, survey_month) %>% filter(n > 1)
weekly %>% count(fragment_ID, survey_date) %>% filter(n > 1)

#### Check missing values ####

monthly %>% summarise(across(everything(), ~ sum(is.na(.x))))
weekly %>% summarise(across(everything(), ~ sum(is.na(.x))))

#### Check treatment and genus consistency ####

bind_rows(
  monthly %>% select(fragment_ID, frame_number, condition, genus),
  weekly %>% select(fragment_ID, frame_number, condition, genus)
) %>%
  distinct() %>%
  count(fragment_ID) %>%
  filter(n > 1)

#### Check survey coverage ####

monthly %>% count(survey_month)
weekly %>% count(survey_date)


#### Monthly and weekly agreement ####

weekly_last <- weekly %>%
  group_by(fragment_ID, survey_month) %>%
  slice_max(survey_date, n = 1, with_ties = FALSE) %>%
  ungroup()

monthly %>%
  filter(survey_month %in% weekly_last$survey_month) %>%
  select(fragment_ID, survey_month, bleaching, survival, mortality) %>%
  inner_join(
    weekly_last %>%
      select(fragment_ID, survey_month, survey_date, bleaching, survival, mortality),
    by = c("fragment_ID", "survey_month"),
    suffix = c("_monthly", "_weekly")
  ) %>%
  group_by(survey_month, survey_date) %>%
  summarise(
    bleaching_diff = sum(bleaching_monthly != bleaching_weekly),
    survival_diff = sum(survival_monthly != survival_weekly),
    mortality_diff = sum(coalesce(mortality_monthly, "NA") !=
                           coalesce(mortality_weekly, "NA")),
    .groups = "drop"
  )


#### April survival discrepancy ####

monthly %>%
  filter(survey_month == as.Date("2026-04-01")) %>%
  select(fragment_ID, survival_monthly = survival) %>%
  inner_join(
    weekly %>%
      filter(survey_month == as.Date("2026-04-01")) %>%
      select(fragment_ID, survey_date, survival_weekly = survival),
    by = "fragment_ID"
  ) %>%
  group_by(survey_date) %>%
  summarise(
    monthly_dead = sum(survival_monthly),
    weekly_dead = sum(survival_weekly),
    disagreements = sum(survival_monthly != survival_weekly),
    .groups = "drop"
  )


#### Weekly mortality coding ####

weekly %>%
  group_by(survey_date) %>%
  summarise(
    dead = sum(survival == 1),
    alive = sum(survival == 0),
    .groups = "drop"
  )

monthly %>%
  group_by(survey_month) %>%
  summarise(dead = sum(survival == 1), .groups = "drop")

## AUDIT: 
# 105 WEEKLY SURVIVAL
# - raw weekly records clearly contain mortality = 5 throughout Apr
# - cleaned weekly survival incorrectly codes all 800 alive until May 22
# - binary recoding/processing error, exact cause unknown
# - cannot reconstruct accurate weekly death dates from monthly data
# - need audited source + cleaning scripts before correcting


#### Bleaching coding ####

monthly %>%
  count(survey_month, condition, bleaching) %>%
  tidyr::pivot_wider(names_from = bleaching, values_from = n, names_prefix = "bleaching_")

weekly %>%
  count(survey_date, condition, bleaching) %>%
  tidyr::pivot_wider(names_from = bleaching, values_from = n, names_prefix = "bleaching_")

#### Bleaching consistency ####

monthly %>%
  count(bleaching)

weekly %>%
  count(bleaching)

#### Bleaching recorded after death ####

monthly %>%
  group_by(fragment_ID) %>%
  mutate(dead_from_first_event = cummax(survival == 1)) %>%
  ungroup() %>%
  summarise(
    observations_after_first_death = sum(dead_from_first_event),
    bleaching_1_after_death = sum(bleaching == 1 & dead_from_first_event),
    bleaching_0_after_death = sum(bleaching == 0 & dead_from_first_event)
  )

# 3,733 observations (29.2%) occur at or after the coral's first recorded death.

#### Bleaching among living fragments ####

monthly %>%
  group_by(fragment_ID) %>%
  mutate(alive = !cummax(survival == 1)) %>%
  ungroup() %>%
  group_by(survey_month, condition) %>%
  summarise(
    n_alive = sum(alive),
    bleached = sum(bleaching == 1 & alive),
    prevalence = bleached / n_alive,
    .groups = "drop"
  ) %>%
  print(n = Inf)


### 01. A1_5 DATA AUDIT ####

#### Monthly records ####

monthly %>%
  filter(fragment_ID == "A1_5") %>%
  select(survey_month, bleaching, survival, mortality, growth_rate) %>%
  print(n = Inf)

#### Weekly records ####

weekly %>%
  filter(fragment_ID == "A1_5") %>%
  select(survey_date, bleaching, survival, mortality) %>%
  print(n = Inf)

#### Check against known death date ####

weekly %>%
  filter(fragment_ID == "A1_5") %>%
  summarise(
    false_alive = sum(survival == 0),
    first_cleaned_death = min(survey_date[survival == 1]),
    mortality_codes = paste(unique(mortality), collapse = ", ")
  )



### 02. CROSS-VARIABLE AUDIT ####

#### Monthly data ####

monthly %>%
  arrange(fragment_ID, survey_month) %>%
  group_by(fragment_ID) %>%
  mutate(
    previously_dead = lag(cummax(survival == 1), default = FALSE),
    ever_dead = cummax(survival == 1),
    first_death = survival == 1 & !previously_dead
  ) %>%
  summarise(
    reversals = sum(previously_dead & survival == 0),
    bleach_after_death = sum(ever_dead & bleaching == 1),
    unbleached_after_death = sum(ever_dead & bleaching == 0),
    growth_after_death = sum(ever_dead & !is.na(growth_rate)),
    cause_changes_after_death = n_distinct(mortality[ever_dead & !is.na(mortality)]) > 1,
    .groups = "drop"
  ) %>%
  summarise(
    fragments_with_reversals = sum(reversals > 0),
    reversal_records = sum(reversals),
    bleach_after_death = sum(bleach_after_death),
    unbleached_after_death = sum(unbleached_after_death),
    growth_after_death = sum(growth_after_death),
    fragments_with_cause_changes = sum(cause_changes_after_death)
  )

