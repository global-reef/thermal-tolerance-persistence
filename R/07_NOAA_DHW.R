### 01. SETUP ####

# NOAA CRW v3.1 daily 5 km products, nearest grid cell to Chalok Bay.
# Satellite data represent regional surface conditions, not coral-level exposure.
noaa_products <- list(
  dhw_noaa = c("CRW_dhw_v1_0", "degree_heating_week"),
  sst_noaa = c("CRW_sst_v3_1", "analysed_sst"),
  hotspot_noaa = c("CRW_hs_v1_0", "hotspot")
)

expected_dates <- seq.Date(as.Date("2025-01-01"), as.Date("2026-07-31"), by = "day")


### 02. DOWNLOAD AND VALIDATE NOAA DATA ####

#### Retrieve original satellite records ####

get_noaa <- function(product, variable, name) {
  file <- file.path(data_raw_dir, paste0(product, "_Chalok_2025_2026.csv"))
  url <- paste0(
    "https://oceanwatch.pifsc.noaa.gov/erddap/griddap/",
    product, ".csv0?", variable,
    "[(2025-01-01T12:00:00Z):1:(2026-07-31T12:00:00Z)]",
    "[(10.075)][(99.825)]"
  )
  
  if (!file.exists(file)) {
    download.file(url, file, mode = "wb", method = "libcurl")
  }
  
  raw <- read_csv(
    file, col_names = c("time", "latitude", "longitude", "value"),
    show_col_types = FALSE
  ) %>%
    mutate(date = as.Date(substr(time, 1, 10)))
  
  stopifnot(
    nrow(raw) > 0,
    all(!is.na(raw$date) & raw$date %in% expected_dates),
    !anyDuplicated(raw$date),
    all(abs(raw$latitude - 10.075) < 1e-4),
    all(abs(raw$longitude - 99.825) < 1e-4),
    all(is.finite(raw$value))
  )
  
  raw %>% transmute(date, !!name := as.numeric(value))
}

noaa_daily <- purrr::imap(noaa_products, ~ get_noaa(.x[1], .x[2], .y)) %>%
  purrr::reduce(full_join, by = "date") %>%
  right_join(tibble(date = expected_dates), by = "date") %>%
  arrange(date)

# Missing satellite dates remain NA; no interpolation.
stopifnot(
  identical(noaa_daily$date, expected_dates),
  all(noaa_daily$dhw_noaa >= 0, na.rm = TRUE),
  all(noaa_daily$sst_noaa > 0, na.rm = TRUE)
)

cat("\n=== NOAA DATA COVERAGE ===\n")
noaa_daily %>%
  summarise(
    days = n(),
    missing_dhw = sum(is.na(dhw_noaa)),
    missing_sst = sum(is.na(sst_noaa)),
    missing_hotspot = sum(is.na(hotspot_noaa))
  ) %>%
  print()

write_csv(noaa_daily, file.path(data_clean_dir, "noaa_daily.csv"))


### 03. SEASONAL NOAA AND IN SITU COMPARISON ####

#### Import local logger summaries ####

temp_daily <- read_csv(
  file.path(data_clean_dir, "temp_daily.csv"),
  show_col_types = FALSE
)

#### Summarise matched warming periods ####

thermal_summary <- noaa_daily %>%
  left_join(temp_daily %>% select(date, temp_mean, dhw_local), by = "date") %>%
  mutate(year = lubridate::year(date)) %>%
  filter(lubridate::month(date) %in% 3:7) %>%
  group_by(year) %>%
  summarise(
    noaa_dhw_peak = max(dhw_noaa, na.rm = TRUE),
    noaa_dhw_peak_date = date[which.max(dhw_noaa)],
    noaa_sst_peak = max(sst_noaa, na.rm = TRUE),
    noaa_hotspot_days_ge_1 = sum(hotspot_noaa >= 1, na.rm = TRUE),
    noaa_dhw_days_ge_4 = sum(dhw_noaa >= 4, na.rm = TRUE),
    noaa_dhw_missing_days = sum(is.na(dhw_noaa)),
    local_temp_peak = max(temp_mean, na.rm = TRUE),
    local_dhw_peak = if (any(!is.na(dhw_local)))
      max(dhw_local, na.rm = TRUE) else NA_real_,
    paired_sst_days = sum(!is.na(temp_mean) & !is.na(sst_noaa)),
    mean_local_minus_satellite = mean(temp_mean - sst_noaa, na.rm = TRUE),
    paired_dhw_days = sum(!is.na(dhw_local) & !is.na(dhw_noaa)),
    .groups = "drop"
  )

print(thermal_summary, width = Inf)

# daily compmarison 
temp_daily %>%
  filter(between(date, as.Date("2026-04-16"), as.Date("2026-05-14"))) %>%
  select(date, temp_mean, complete_day) %>%
  left_join(noaa_daily %>% select(date, sst_noaa, hotspot_noaa), by = "date") %>%
  summarise(
    paired_days = sum(complete_day & !is.na(sst_noaa)),
    local_mean = mean(temp_mean[complete_day], na.rm = TRUE),
    satellite_mean = mean(sst_noaa, na.rm = TRUE),
    local_days_ge_31 = sum(complete_day & temp_mean >= 31, na.rm = TRUE),
    satellite_days_ge_31.15 = sum(sst_noaa >= 31.15, na.rm = TRUE),
    local_max = max(temp_mean, na.rm = TRUE),
    satellite_max = max(sst_noaa, na.rm = TRUE),
    mean_difference = mean(temp_mean - sst_noaa, na.rm = TRUE)
  )


write_csv(
  thermal_summary,
  file.path(summaries_dir, "noaa_local_thermal_comparison.csv")
)

rm(get_noaa, noaa_products, expected_dates)
