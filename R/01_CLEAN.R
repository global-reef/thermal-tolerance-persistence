
### 01. IMPORT RAW DATA ####

#### Workbook structure ####

raw_file <- file.path(data_raw_dir, "2026.11.16_ThermalTolerance_MASTER.xlsx")

sheets <- c(
  A = "A (shallow)", B = "B (deep)",
  C = "C (shallow)", D = "D (deep)"
)

# Eight identification/baseline columns followed by 27 survey blocks
column_names <- c(
  "coral_id", "frame", "quadrat", "position", "genus", "physical_damage",
  "baseline_date", "baseline_size_cm",
  paste0("s01_", c("date", "bleaching", "mortality", "cause")),
  purrr::map(2:27, ~ paste0(sprintf("s%02d_", .x),
                            c("date", "size_cm", "bleaching", "mortality", "cause"))) %>%
    unlist()
)

#### Import and reshape ####

# Fixed range excludes formatting rows; text import preserves original entries.
raw_wide <- purrr::imap_dfr(sheets, ~ {
  readxl::read_excel(
    raw_file, sheet = .x, range = "A2:EL201",
    col_names = column_names, col_types = "text"
  ) %>%
    mutate(source_sheet = .y, .before = 1)
})

raw_long <- raw_wide %>%
  pivot_longer(
    cols = matches("^s\\d{2}_"),
    names_to = c("survey_index", ".value"),
    names_pattern = "s(\\d{2})_(.*)",
    values_drop_na = FALSE
  )


### 02. STANDARDISE AND CORRECT VARIABLES ####

#### Mortality cause definitions ####

# Confirmed transcription corrections and primary-cause allocation rules
cause_dictionary <- c(
  PRED = "PRD", PR = "PRD",
  ALG = "OG", AL = "OG",
  BLER = "BLE", UKNN = "UKN",
  MRYH = "METH",
  "BLE + UKN" = "BLE",
  "BLE + PHY" = "BLE",
  "DIS + BLE" = "DIS",
  "UKN + PHY" = "UKN",
  N = NA_character_,
  "NA" = NA_character_
)

#### Variable types and coding ####

coral_clean <- raw_long %>%
  mutate(
    source_sheet = factor(source_sheet, levels = names(sheets)),
    coral_id = str_trim(coral_id),
    frame = factor(str_trim(frame)),
    origin = factor(if_else(source_sheet %in% c("A", "C"), "shallow", "deep"),
                    levels = c("deep", "shallow")),
    genus = factor(if_else(coral_id == "D3_15", "ACRO", genus),
                   levels = c("ACRO", "POR")),
    position = as.integer(as.numeric(position)),
    survey_index = as.integer(survey_index),
    baseline_date = as.Date(as.numeric(baseline_date), origin = "1899-12-30"),
    survey_date = as.Date(as.numeric(date), origin = "1899-12-30"),
    bleaching = suppressWarnings(as.integer(as.numeric(bleaching))),
    mortality = suppressWarnings(as.integer(as.numeric(mortality))),
    baseline_size_cm = na_if(suppressWarnings(as.numeric(baseline_size_cm)), 0),
    size_cm = size_cm %>%
      na_if("NA") %>%
      na_if("-") %>%
      na_if(".") %>%
      na_if("O") %>%
      str_remove_all("\\*") %>%
      str_replace_all(fixed("18/.22"), "18.22") %>%
      str_remove("^`") %>%
      as.numeric() %>%
      na_if(0),
    cause = factor(
      recode(str_squish(cause), !!!cause_dictionary),
      levels = c("ACC", "BLE", "COMP", "DIS", "METH",
                 "OG", "PHY", "PRD", "UKN")
    )
  ) %>%
  select(-date, -physical_damage, -quadrat)

#### Verified photographic corrections ####

# Seven mortality scores and associated causes verified against field photographs
mortality_corrections <- tibble(
  coral_id = c("A2_18", "C3_24", "C3_25", "C5_28",
               "C5_30", "C5_40", "D5_1"),
  survey_date = as.Date(c("2025-06-27", rep("2025-10-23", 5), "2026-04-24")),
  mortality_fix = c(0L, 1L, 0L, 0L, 0L, 0L, 1L),
  cause_fix = c(NA, "METH", NA, NA, NA, NA, "OG")
)

stopifnot(
  nrow(inner_join(coral_clean, mortality_corrections,
                  by = c("coral_id", "survey_date"))) == 7
)

coral_clean <- coral_clean %>%
  left_join(mortality_corrections, by = c("coral_id", "survey_date")) %>%
  mutate(
    cause = factor(if_else(!is.na(mortality_fix), cause_fix, as.character(cause)),
                   levels = levels(cause)),
    mortality = coalesce(mortality_fix, mortality)
  ) %>%
  select(-mortality_fix, -cause_fix)


### 03. RECONSTRUCT MORTALITY ####

#### Irreversible complete mortality ####

# Score 5 is irreversible; partial tissue mortality (1–4) may recover.
# Retain the first complete-death cause, assigning UKN if missing.
# Bleaching and size are ineligible at and after complete mortality.
coral_clean <- coral_clean %>%
  arrange(coral_id, survey_index) %>%
  group_by(coral_id) %>%
  mutate(
    dead = cumany(coalesce(mortality == 5L, FALSE)),
    cause = if_else(
      dead,
      coalesce(cause[match(5L, mortality)], factor("UKN", levels = levels(cause))),
      cause
    ),
    mortality = replace(mortality, dead, 5L),
    bleaching = replace(bleaching, dead, NA_integer_),
    size_cm = replace(size_cm, dead, NA_real_)
  ) %>%
  ungroup()


### 04. VALIDATE SIZE MEASUREMENTS ####

#### Flag implausible size reversals ####

# Flag a >10 cm change followed by an opposite >10 cm change
# returning within 5 cm of the original measurement.
# Remove only the middle measurement; genuine fragmentation is retained.
# Apply once, without interpolation or iterative removal.
size_errors <- coral_clean %>%
  filter(survey_index %in% 1:12) %>%
  select(coral_id, survey_date, size_cm) %>%
  bind_rows(
    coral_clean %>%
      distinct(coral_id, baseline_date, baseline_size_cm) %>%
      transmute(coral_id, survey_date = baseline_date,
                size_cm = baseline_size_cm)
  ) %>%
  filter(!is.na(size_cm)) %>%
  arrange(coral_id, survey_date) %>%
  group_by(coral_id) %>%
  mutate(
    size_change = size_cm - lag(size_cm),
    previous_change = lag(size_change),
    middle_date = lag(survey_date)
  ) %>%
  filter(
    abs(size_change) > 10,
    abs(previous_change) > 10,
    sign(size_change) != sign(previous_change),
    abs(previous_change + size_change) < 5
  ) %>%
  ungroup() %>%
  transmute(coral_id, survey_date = middle_date) %>%
  distinct()

coral_clean <- coral_clean %>%
  left_join(size_errors %>% mutate(size_error = TRUE),
            by = c("coral_id", "survey_date")) %>%
  mutate(size_cm = if_else(coalesce(size_error, FALSE), NA_real_, size_cm)) %>%
  select(-size_error)


### 05. DERIVE ANALYTICAL VARIABLES ####

#### Binary bleaching and mortality ####

# Primary response: substantial bleaching (scores 3–5).
# Retain any bleaching (scores 1–5) for sensitivity analysis.
# Complete mortality is irreversible; missing assessments remain NA.
coral_clean <- coral_clean %>%
  mutate(
    bleached = case_when(
      is.na(bleaching) ~ NA_integer_,
      bleaching %in% 0:2 ~ 0L,
      bleaching %in% 3:5 ~ 1L
    ),
    bleached_any = case_when(
      is.na(bleaching) ~ NA_integer_,
      bleaching == 0L ~ 0L,
      bleaching %in% 1:5 ~ 1L
    ),
    mortality_binary = case_when(
      dead ~ 1L,
      is.na(mortality) ~ NA_integer_,
      TRUE ~ 0L
    )
  ) %>%
  select(-dead)


### 06. CONSTRUCT ANALYTICAL DATASETS ####

#### Study period ####

study_end <- as.Date("2026-07-31")

#### Monthly observations ####

# April–June 2026 use the final weekly survey of each month.
coral_monthly <- coral_clean %>%
  filter(survey_index %in% c(1:12, 16, 20, 23, 24),
         survey_date <= study_end) %>%
  arrange(coral_id, survey_date)

#### Weekly observations ####

coral_weekly <- coral_clean %>%
  filter(survey_index %in% 13:23,
         survey_date <= study_end) %>%
  arrange(coral_id, survey_date)

#### Fragment-level mortality ####

# One record per fragment. Death dates indicate first observed complete mortality.
coral_mortality <- coral_clean %>%
  filter(survey_date <= study_end) %>%
  arrange(coral_id, survey_index) %>%
  group_by(coral_id, frame, origin, genus) %>%
  summarise(
    first_damage_cause = cause[match(TRUE, mortality %in% 1:4)],
    death_date = survey_date[match(5L, mortality)],
    death_cause = cause[match(5L, mortality)],
    last_alive_date = dplyr::last(
      survey_date[mortality %in% 0:4], default = as.Date(NA)
    ),
    last_observed_date = dplyr::last(
      survey_date[!is.na(mortality)], default = as.Date(NA)
    ),
    .groups = "drop"
  ) %>%
  mutate(event = as.integer(!is.na(death_date)))

### 07. NOAA CLIMATOLOGY AND LOCAL TEMPERATURE ####

#### Extract Chalok Bay climatological MMM ####
# Rounded from NOAA Coral Reef Watch Chalok Bay pixel MMM (30.15 °C).
mmm <- 30

# NOAA Coral Reef Watch 5 km v3.1 climatology; 30 °C is the rounded
# analysis baseline defined in 00_SETUP.R, not the exact NOAA pixel value.
noaa_file <- list.files(data_raw_dir, pattern = "\\.nc$", full.names = TRUE)
stopifnot(length(noaa_file) == 1)
noaa_climatology <- terra::rast(noaa_file)
stopifnot("sst_clim_mmm" %in% names(noaa_climatology))
chalok_coords <- terra::vect(data.frame(lon = 99.82727, lat = 10.066069),
                             geom = c("lon", "lat"), crs = "EPSG:4326")
chalok_mmm <- terra::extract(noaa_climatology[["sst_clim_mmm"]], chalok_coords)[[2]]
stopifnot(length(chalok_mmm) == 1, is.finite(chalok_mmm))
cat("\n=== CHALOK NOAA MMM (EXACT) AND ANALYSIS BASELINE ===\n")
print(tibble(noaa_mmm = chalok_mmm, analysis_mmm = mmm))
rm(noaa_climatology, chalok_coords)



#### Import and validate local logger ####

temp_clean <- read_csv(file.path(data_raw_dir, "temp.csv"), show_col_types = FALSE) %>%
  mutate(
    datetime = as.POSIXct(paste(date, time), format = "%d/%m/%Y %H:%M:%S",
                          tz = "Asia/Bangkok"),
    date = as.Date(datetime, tz = "Asia/Bangkok"),
    temperature = as.numeric(temperature)
  ) %>%
  select(datetime, date, temperature) %>%
  arrange(datetime)

# Exclude logger observations after retrieval on 14 July 2026 at 12:00.
temp_clean <- temp_clean %>%
  filter(datetime <= as.POSIXct("2026-07-14 12:00:00", tz = "Asia/Bangkok"))

stopifnot(!anyNA(temp_clean$datetime), !anyNA(temp_clean$temperature),
          !anyDuplicated(temp_clean$datetime))

#### Daily temperature statistics and missing days ####


# Six four-hourly readings constitute a complete day. Incomplete or
# missing days do not contribute to the observed-only DHW calculation.
temp_daily <- temp_clean %>%
  group_by(date) %>%
  summarise(
    n_readings = n(),
    temp_mean = mean(temperature),
    temp_min = min(temperature),
    temp_max = max(temperature),
    .groups = "drop"
  ) %>%
  tidyr::complete(date = seq.Date(min(date), max(date), by = "day"),
                  fill = list(n_readings = 0L)) %>%
  arrange(date) %>%
  mutate(
    complete_day = n_readings == 6,
    hotspot = pmax(temp_mean - mmm, 0),
    heat_daily = if_else(complete_day, if_else(hotspot >= 1, hotspot, 0), NA_real_)
  )

stopifnot(all(temp_daily$n_readings %in% 0:6),
          !anyDuplicated(temp_daily$date))

#### Rolling 84-day observed-only local DHW ####

# Only HotSpots >= 1 °C accumulate; 84 fully observed days are required.
# Incomplete windows remain NA instead of silently treating gaps as zero.
temp_daily <- temp_daily %>%
  mutate(
    n_complete_84d = slider::slide_int(complete_day, sum, .before = 83, .complete = TRUE),
    dhw_complete = n_complete_84d == 84,
    dhw_local = slider::slide_dbl(heat_daily, ~ sum(.x, na.rm = TRUE) / 7,
                                 .before = 83, .complete = TRUE),
    dhw_local = if_else(dhw_complete, dhw_local, NA_real_)
  ) %>%
  select(-heat_daily)

#### Temperature integrity summary ####

cat("\n=== LOCAL LOGGER AND DHW COVERAGE ===\n")
temp_daily %>%
  summarise(
    first_date = min(date), last_date = max(date),
    days = n(), missing_or_incomplete_days = sum(!complete_day),
    complete_dhw_windows = sum(dhw_complete, na.rm = TRUE),
    peak_complete_dhw = if (any(dhw_complete, na.rm = TRUE))
      max(dhw_local, na.rm = TRUE) else NA_real_
  ) %>%
  print()

### 08. VALIDATION AND SUMMARIES ####

#### Structural integrity ####

stopifnot(
  nrow(raw_wide) == 800,
  n_distinct(raw_wide$coral_id) == 800,
  nrow(raw_long) == 21600,
  nrow(coral_clean) == 21600,
  nrow(coral_monthly) == 12800,
  nrow(coral_weekly) == 8800,
  nrow(coral_mortality) == 800,
  nrow(distinct(coral_clean, coral_id, survey_index)) == nrow(coral_clean),
  n_distinct(coral_clean$frame) == 20,
  all(as.integer(table(coral_clean$frame, coral_clean$survey_index)) == 40),
  !anyNA(coral_clean$genus),
  all(!is.na(coral_clean$survey_date[coral_clean$survey_index <= 24])),
  all(is.na(coral_clean$bleaching) | coral_clean$bleaching %in% 0:5),
  all(is.na(coral_clean$mortality) | coral_clean$mortality %in% 0:5),
  all(is.na(coral_clean$bleached) | coral_clean$bleached %in% 0:1),
  all(is.na(coral_clean$mortality_binary) | coral_clean$mortality_binary %in% 0:1),
  all(is.na(coral_clean$size_cm) | coral_clean$size_cm > 0),
  all(coral_clean$mortality[coral_clean$mortality_binary == 1L] == 5L),
  all(is.na(coral_clean$bleaching[coral_clean$mortality_binary == 1L])),
  all(is.na(coral_clean$size_cm[coral_clean$mortality_binary == 1L])),
  all(is.na(coral_clean$bleached[coral_clean$mortality_binary == 1L])),
  all(is.na(coral_clean$bleached_any) | coral_clean$bleached_any %in% 0:1),
  all(is.na(coral_clean$bleached_any[coral_clean$mortality_binary == 1L])),
  all(is.na(coral_clean$bleached) |
        coral_clean$bleached == as.integer(coral_clean$bleaching >= 3))
)

cat("\n=== STRUCTURAL INTEGRITY: PASSED ===\n")

#### Experimental balance ####

cat("\n=== EXPERIMENTAL DESIGN BY FRAME AND GENUS ===\n")
coral_clean %>%
  distinct(coral_id, frame, origin, genus) %>%
  count(frame, origin, genus) %>%
  arrange(frame, genus) %>%
  print(n = Inf)

#### Survival and mortality causes ####

cat("\n=== SURVIVAL BY EXPERIMENTAL GROUP ===\n")
coral_mortality %>%
  group_by(origin, genus) %>%
  summarise(
    fragments = n(),
    deaths = sum(event),
    survival_pct = round(100 * mean(event == 0), 1),
    .groups = "drop"
  ) %>%
  print(n = Inf)

cat("\n=== CAUSES OF COMPLETE MORTALITY ===\n")
coral_mortality %>%
  filter(event == 1) %>%
  count(death_cause, sort = TRUE) %>%
  mutate(percent = round(100 * n / sum(n), 1)) %>%
  print(n = Inf)

#### Mortality chronology ####

cat("\n=== MORTALITY CHRONOLOGY ===\n")
coral_mortality %>%
  filter(event == 1) %>%
  count(death_date, name = "new_deaths") %>%
  arrange(death_date) %>%
  mutate(cumulative_deaths = cumsum(new_deaths)) %>%
  print(n = Inf)

#### Observation completeness ####

cat("\n=== MONTHLY AND WEEKLY OBSERVATION COMPLETENESS ===\n")
bind_rows(
  coral_monthly %>% mutate(dataset = "Monthly"),
  coral_weekly %>% mutate(dataset = "Weekly")
) %>%
  group_by(dataset, survey_date) %>%
  summarise(
    alive = sum(mortality_binary == 0, na.rm = TRUE),
    dead = sum(mortality_binary == 1, na.rm = TRUE),
    unknown = sum(is.na(mortality_binary)),
    bleaching_observed = sum(!is.na(bleached)),
    size_observed = sum(!is.na(size_cm)),
    .groups = "drop"
  ) %>%
  print(n = Inf)

#### Bleaching distributions ####

cat("\n=== BLEACHING SEVERITY DISTRIBUTIONS ===\n")
bind_rows(
  coral_monthly %>% mutate(dataset = "Monthly"),
  coral_weekly %>% mutate(dataset = "Weekly")
) %>%
  count(dataset, bleaching) %>%
  group_by(dataset) %>%
  mutate(percent = round(100 * n / sum(n), 1)) %>%
  ungroup() %>%
  print(n = Inf)

#### Bleaching prevalence by experimental group ####

cat("\n=== MONTHLY BLEACHING PREVALENCE AND SEVERITY BY GROUP ===\n")
coral_monthly %>%
  group_by(survey_date, origin, genus) %>%
  summarise(
    assessed = sum(!is.na(bleached)),
    bleached = sum(bleached == 1, na.rm = TRUE),
    prevalence_pct = round(100 * bleached / na_if(assessed, 0), 1),
    mean_severity = round(mean(bleaching, na.rm = TRUE), 2),
    .groups = "drop"
  ) %>%
  print(n = Inf)

#### Missingness and mortality consistency ####

cat("\n=== MISSINGNESS AND MORTALITY CONSISTENCY ===\n")
coral_clean %>%
  filter(survey_date <= study_end) %>%
  summarise(
    missing_mortality = sum(is.na(mortality_binary)),
    missing_bleaching_alive = sum(mortality_binary == 0 & is.na(bleached)),
    bleaching_after_death = sum(mortality_binary == 1 & !is.na(bleached)),
    size_after_death = sum(mortality_binary == 1 & !is.na(size_cm))
  ) %>%
  print()

#### Size screening summary ####

cat("\n=== SIZE MEASUREMENT SCREENING ===\n")
tibble(
  flagged_measurements = nrow(size_errors),
  affected_fragments = n_distinct(size_errors$coral_id),
  remaining_size_observations = sum(!is.na(coral_monthly$size_cm))
) %>%
  print()


### 09. SAVE CLEAN DATA ####

write_csv(coral_clean, file.path(data_clean_dir, "coral_clean.csv"), na = "")
write_csv(coral_monthly, file.path(data_clean_dir, "coral_monthly.csv"), na = "")
write_csv(coral_weekly, file.path(data_clean_dir, "coral_weekly.csv"), na = "")
write_csv(coral_mortality, file.path(data_clean_dir, "coral_mortality.csv"), na = "")
write_csv(temp_clean, file.path(data_clean_dir, "temp_clean.csv"))
write_csv(temp_daily, file.path(data_clean_dir, "temp_daily.csv"))

cat("\n=== CLEANING COMPLETE: 4 CORAL AND 2 TEMPERATURE DATASETS SAVED ===\n")


### 10. SAMPLING EFFORT FOR PUBLICATION ####

#### Experimental design and monitoring effort ####

cat("\n=== PUBLICATION: EXPERIMENTAL DESIGN AND MONITORING EFFORT ===\n")
coral_clean %>%
  filter(survey_date <= study_end) %>%
  summarise(
    fragments = n_distinct(coral_id),
    frames = n_distinct(frame),
    shallow_fragments = n_distinct(coral_id[origin == "shallow"]),
    deep_fragments = n_distinct(coral_id[origin == "deep"]),
    acropora = n_distinct(coral_id[genus == "ACRO"]),
    porites = n_distinct(coral_id[genus == "POR"]),
    first_survey = min(survey_date),
    last_survey = max(survey_date),
    survey_dates = n_distinct(survey_date),
    fragment_observations = n()
  ) %>%
  print(width = Inf)

#### Analytical sample sizes ####

cat("\n=== PUBLICATION: ANALYTICAL SAMPLE SIZES ===\n")
bind_rows(
  coral_monthly %>% mutate(dataset = "Monthly"),
  coral_weekly %>% mutate(dataset = "Weekly")
) %>%
  group_by(dataset) %>%
  summarise(
    surveys = n_distinct(survey_date),
    fragments = n_distinct(coral_id),
    observations = n(),
    bleaching_assessments = sum(!is.na(bleaching)),
    .groups = "drop"
  ) %>%
  print(width = Inf)

#### Mortality and size measurement effort ####

cat("\n=== PUBLICATION: MORTALITY AND SIZE MEASUREMENT EFFORT ===\n")
tibble(
  confirmed_deaths = sum(coral_mortality$event),
  surviving_fragments = sum(coral_mortality$event == 0),
  baseline_sizes = sum(!is.na(
    coral_clean %>% distinct(coral_id, baseline_size_cm) %>% pull(baseline_size_cm)
  )),
  followup_sizes = sum(!is.na(coral_monthly$size_cm)),
  fragments_with_followup_sizes = n_distinct(
    coral_monthly$coral_id[!is.na(coral_monthly$size_cm)]
  )
) %>%
  print(width = Inf)


