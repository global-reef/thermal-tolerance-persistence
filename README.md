# Coral thermal tolerance following transplantation

Analysis repository for a Global Reef study investigating whether differences in coral thermal tolerance associated with donor depth persist after transplantation to a common reef environment at Chalok Bay, Koh Tao, Thailand.

## Overview

In March 2025, 800 fragments of *Acropora* and *Porites* were transplanted from shallow (1–2 m) and deep (approximately 10 m) donor environments onto 20 experimental frames. Coral condition, survival and size were monitored through July 2026, encompassing the 2025 and 2026 warming periods. Weekly bleaching assessments were also conducted during the 2026 warming period.

## Research questions

1. Do shallow- and deep-origin corals differ in substantial bleaching prevalence, and does the treatment difference change between the 2025 and 2026 warming periods?
2. How do bleaching severity and prevalence change over the 2026 warming period?
3. Does donor treatment affect survival and complete mortality, including mortality attributed to different causes?
4. Does positive size extension differ by donor treatment and coral genus?

## Analysis framework

The analyses use R and distinguish donor **treatment** (shallow or deep origin) from the shared transplantation environment.

- **Monthly bleaching:** Binomial mixed-effects models compare substantial bleaching (scores 3–5, representing at least 50% bleached tissue) across treatments, warming seasons and genera. Population-marginal probabilities and treatment differences are estimated with parametric-bootstrap confidence intervals.
- **Weekly bleaching:** Observed severity distributions and prevalence are presented alongside a binomial mixed-effects model of temporal treatment differences.
- **Survival and mortality:** Binomial mixed-effects models assess endpoint survival and non-attachment mortality. Kaplan–Meier survival curves and Aalen–Johansen cumulative incidence estimates describe event timing and competing causes of mortality.
- **Positive size extension:** A Gamma mixed-effects model compares coral-level mean positive extension between treatments and genera. Zero and negative size changes are excluded from this specific response.

## Repository structure

```text
R/
├── 00_SETUP.R
├── 01_CLEAN.R
├── 02_PREP.R
├── 03_EXPLORE.R
├── 04_MODEL.R
├── 05_RESULTS.R
├── 06_PLOT.R
└── 99_scratch.R

data/
├── raw/
└── clean/
docs/

Analysis_YYYY.MM.DD/
├── eda/
├── fits/
├── plots/
├── stats/
├── summaries/
└── tables/
```

### Scripts

| Script | Purpose |
|---|---|
| `00_SETUP.R` | Loads packages; defines paths, analysis date, plotting settings and random seed. |
| `01_CLEAN.R` | Imports and validates coral monitoring and local temperature records. |
| `02_PREP.R` | Constructs analysis-ready bleaching, survival and positive-extension datasets. |
| `03_EXPLORE.R` | Examines sampling structure, response distributions and analytical assumptions. |
| `04_MODEL.R` | Fits primary models, sensitivity analyses and model diagnostics. |
| `05_RESULTS.R` | Generates estimated effects, contrasts, confidence intervals and manuscript tables. |
| `06_PLOT.R` | Produces manuscript figures and supplementary visualisations. |
| `99_scratch.R` | Informal exploratory work outside the reproducible analysis pipeline. |

## Reproducibility

Run scripts `00_SETUP.R` through `06_PLOT.R` sequentially from the repository root. The source data are included under `data/raw/`, and the required R packages are listed in `00_SETUP.R`. The setup script uses `set.seed(42)`.

Results are written to `Analysis_YYYY.MM.DD/`, allowing dated runs to be retained independently. Manuscript figures are exported as PNG and PDF; detailed estimates, model objects and tables are saved in the corresponding output subdirectories.

## Interpretation notes

Bleaching analyses concern **substantial bleaching**, not all visible paling. Observations after complete mortality are excluded from bleaching and size responses. Mortality dates indicate the first survey at which complete mortality was recorded, rather than the exact date of death. Temperature measurements provide environmental context; incomplete daily records are not interpolated for thermal-stress calculations.

## Project

Global Reef, Koh Tao, Thailand. For research enquiries, contact [scarlett@global-reef.com](mailto:scarlett@global-reef.com).
