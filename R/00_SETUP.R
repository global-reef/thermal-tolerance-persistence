
### 00. SETUP ####
set.seed(42)
#### Packages ####

suppressPackageStartupMessages({
  library(readxl)
  library(tidyr)
  library(janitor)
  library(lubridate)
  library(glmmTMB)
  library(survival)
  library(survminer)
  library(DHARMa)
  library(performance)
  library(emmeans)
  library(tidyverse)
  library(patchwork)
  library(terra)
  library(dplyr)
  library(beepr)
})

#### Analysis date ####

analysis_date <- "2026.11.17" # Update for each new analysis version

#### Directories ####

data_raw_dir   <- "data/raw"
data_clean_dir <- "data/clean"
r_dir          <- "R"
docs_dir       <- "docs"

purrr::walk(c(data_raw_dir, data_clean_dir, docs_dir),
            ~ dir.create(.x, showWarnings = FALSE, recursive = TRUE))

#### Output structure ####

output_dir    <- paste0("Analysis_", analysis_date)

fits_dir      <- file.path(output_dir, "fits")
plots_dir     <- file.path(output_dir, "plots")
tables_dir    <- file.path(output_dir, "tables")
stats_dir     <- file.path(output_dir, "stats")
summaries_dir <- file.path(output_dir, "summaries")
eda_dir       <- file.path(output_dir, "eda")

purrr::walk(c(output_dir, fits_dir, plots_dir, tables_dir,
              stats_dir, summaries_dir, eda_dir),
            ~ dir.create(.x, showWarnings = FALSE, recursive = TRUE))

#### Plot theme ####

theme_clean <- theme_minimal(base_family = "sans") +
  theme(
    legend.position = "bottom",
    panel.grid = element_blank(),
    plot.title = element_blank(),
    panel.background = element_rect(fill = "transparent", colour = NA),
    plot.background = element_rect(fill = "transparent", colour = NA)
  )


#### Colour palettes ####

treatment_palette <- c(
  shallow = "#E2A76F",
  deep    = "#1F4E6B"
)

condition_palette <- c(
  "0" = "#1F4E6B",
  "1" = "#477F91",
  "2" = "#7FA9A8",
  "3" = "#B9C7B5",
  "4" = "#D1C8B5",
  "5" = "#E8D8B5",
  "Dead" = "#90766C"
)

mortality_palette <- c(
  METH = "#1F4E6B",
  OG   = "#477F91",
  DIS  = "#B9C7B5",
  BLE  = "#E8D8B5",
  PHY  = "#9B8579",
  ACC  = "gray70",
  UKN  = "gray30"
)

mortality_labels <- c(
  METH = "Attachment failure",
  OG   = "Overgrowth",
  DIS  = "Disease",
  BLE  = "Bleaching",
  PHY  = "Physical damage",
  ACC  = "Acclimation",
  UKN  = "Unknown"
)


#### Shared helper functions ####


## citations 
packages <- c("mclogit", "emmeans", "vegan", "ggalluvial", "tidyverse", "DHARMa", "vegan", "glmmTMB", "terra")

for (pkg in packages) {
  cat("\n\n###", pkg, "###\n")
  print(citation(pkg))
}

sessionInfo()
