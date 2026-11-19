### 01. FIGURE SETUP ####

#### Publication dimensions ####

figure_width <- 7.2
half_page_height <- 4.5
full_page_height <- 8
figure_dpi <- 600

#### Figure export helper ####
save_figure <- function(plot,
                        filename,
                        size = c("half", "full"),
                        width = figure_width,
                        dir = plots_dir,
                        dpi = figure_dpi,
                        bg = "white") {
  
  size <- match.arg(size)
  height <- if (size == "half") half_page_height else full_page_height
  
  png_path <- file.path(dir, paste0(filename, ".png"))
  pdf_path <- file.path(dir, paste0(filename, ".pdf"))
  
  ggsave(
    filename = png_path,
    plot = plot,
    width = width,
    height = height,
    units = "in",
    dpi = dpi,
    bg = bg)
  
  ggsave(
    filename = pdf_path,
    plot = plot,
    width = width,
    height = height,
    units = "in",
    bg = bg)
}


### 02. FIGURE 2: MONTHLY BLEACHING AND TEMPERATURE ####

#### A. Observed substantial bleaching prevalence ####

# Scores 3–5 represent >=50% tissue bleached. Colour scale is capped at 40%.
f2_heatmap <- bleach_warming %>%
  mutate(
    year = factor(lubridate::year(survey_date)),
    month = factor(format(survey_date, "%b"), levels = month.abb[3:7]),
    treatment = factor(treatment, levels = c("deep", "shallow"),
                       labels = c("Deep", "Shallow"))
  ) %>%
  group_by(year, month, genus, treatment) %>%
  summarise(prevalence = mean(bleached), .groups = "drop") %>%
  ggplot(aes(month, treatment, fill = prevalence)) +
  geom_tile(colour = "white", linewidth = 0.8) +
  facet_grid(genus ~ year,
             labeller = labeller(genus = c(ACRO = "Acropora", POR = "Porites"))) +
  scale_fill_gradientn(
    colours = unname(condition_palette[as.character(0:5)]),
    limits = c(0, 0.4), breaks = seq(0, 0.4, 0.1),
    labels = scales::percent, oob = scales::squish,
    guide = guide_colourbar(
      direction = "vertical", title.position = "top", title.hjust = 0.5,
      barheight = unit(3.2, "cm"), barwidth = unit(0.35, "cm"))
  ) +
  labs(x = "Month", y = NULL, fill = "Bleaching\nprevalence") +
  theme_clean +
  theme(
    strip.text.x = element_text(size = 9),
    strip.text.y = element_text(size = 9, face = "italic"),
    axis.text = element_text(size = 8),
    legend.position = "left", legend.direction = "vertical",
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 8),
    panel.spacing = unit(0.7, "lines"),
    legend.margin = margin(0, 2, 0, 0),
    legend.box.margin = margin(0, 0, 0, 0)
  )

#### B. Local temperature during warming periods ####

# Restrict to complete days and preserve gaps in logger coverage.
f2_temperature <- temp_daily %>%
  filter(lubridate::month(date) %in% 3:7) %>%
  mutate(
    year = factor(lubridate::year(date)),
    day_of_season = as.numeric(date - as.Date(paste0(year, "-03-01"))),
    segment = cumsum(!complete_day)
  ) %>%
  filter(complete_day) %>%
  ggplot(aes(day_of_season, group = interaction(year, segment))) +
  geom_ribbon(aes(ymin = temp_min, ymax = temp_max),
              fill = condition_palette[["2"]], alpha = 0.5) +
  geom_line(aes(y = temp_mean), colour = condition_palette[["0"]],
            linewidth = 0.7) +
  geom_hline(yintercept = c(mmm, mmm + 1),
             linetype = "dotted", colour = "gray50") +
  facet_wrap(~ year, nrow = 1) +
  scale_x_continuous(
    breaks = as.numeric(as.Date(paste0("2025-", sprintf("%02d", 3:7), "-01")) -
                          as.Date("2025-03-01")),
    labels = month.abb[3:7]
  ) +
  scale_y_continuous(breaks = 27:33) +
  labs(x = "Month", y = "Temperature (°C)") +
  theme_clean +
  theme(strip.text = element_text(size = 9), axis.text = element_text(size = 8),
        plot.margin = margin(5, 5, 5, 75))

#### C. Assemble and export ####

fig2 <- (f2_heatmap / f2_temperature) +
  patchwork::plot_layout(heights = c(1.1, 1)) +
  patchwork::plot_annotation(
    tag_levels = "a", tag_prefix = "(", tag_suffix = ")",
    theme = theme(plot.tag = element_text(face = "bold", size = 11))
  )

print(fig2)
save_figure(fig2, "F2_bleaching_temperature", size = "half")
rm(f2_heatmap, f2_temperature)


### 03. FIGURE 3: WEEKLY BLEACHING ####

#### A. Prepare weekly bleaching-score distributions ####

# Follow corals assessed at the first weekly survey.
# Each bar represents the initial treatment × genus cohort.
# Missing bleaching scores indicate complete mortality.

f3_stacked_data <- coral_weekly %>%
  semi_join(
    coral_weekly %>%
      filter(survey_index == 13, !is.na(bleaching)) %>%
      distinct(coral_id),
    by = "coral_id"
  ) %>%
  mutate(
    treatment = factor(origin, levels = c("deep", "shallow"),
                       labels = c("Deep", "Shallow")),
    genus = factor(genus, levels = c("ACRO", "POR"),
                   labels = c("Acropora", "Porites")),
    state = if_else(is.na(bleaching), "Dead", as.character(bleaching)),
    state = factor(state, levels = c(as.character(0:5), "Dead")),
    survey = factor(
      survey_index,
      levels = 13:23,
      labels = format(
        coral_weekly %>%
          filter(survey_index %in% 13:23) %>%
          distinct(survey_index, survey_date) %>%
          arrange(survey_index) %>%
          pull(survey_date),
        "%d %b"
      )
    )
  ) %>%
  count(survey, genus, treatment, state, name = "corals") %>%
  group_by(survey, genus, treatment) %>%
  mutate(proportion = corals / sum(corals)) %>%
  ungroup()

#### B. Observed weekly bleaching-score distributions ####

f3_stacked <- ggplot(f3_stacked_data,
                     aes(survey, proportion, fill = state)) +
  geom_col(width = 0.88) +
  facet_grid(genus ~ treatment) +
  scale_fill_manual(
    values = condition_palette,
    breaks = c(as.character(0:5), "Dead"),
    drop = FALSE
  ) +
  guides(fill = guide_legend(nrow = 1, byrow = TRUE)) +
  scale_x_discrete(
    breaks = levels(f3_stacked_data$survey)[c(2, 6, 10)],
    labels = c("Apr", "May", "Jun")
  ) +
  scale_y_continuous(
    labels = scales::percent,
    breaks = seq(0, 1, 0.25),
    expand = expansion(mult = c(0, 0.01))
  ) +
  labs(x = "Survey date", y = "Proportion of corals",
       fill = "Bleaching score") +
  theme_clean +
  theme(
    strip.text.x = element_text(size = 9),
    strip.text.y = element_text(size = 9, face = "italic"),
    axis.text.x = element_text(size = 8),
    legend.position = "bottom",
    legend.title = element_text(size = 9),
    legend.text = element_text(size = 8),
    legend.key.size = unit(0.35, "cm"),
    panel.spacing = unit(0.6, "lines")
  )


#### C. Observed weekly bleaching trajectories ####

# Wilson 95% CIs describe observed proportions, without accounting
# for repeated assessments or frame-level clustering.

f3_trajectories <- bleach_weekly %>%
  group_by(survey_date, genus, treatment) %>%
  summarise(
    assessed = n(),
    events = sum(bleached),
    prevalence = events / assessed,
    .groups = "drop"
  ) %>%
  mutate(
    z = qnorm(0.975),
    centre = (prevalence + z^2 / (2 * assessed)) / (1 + z^2 / assessed),
    margin = z * sqrt(prevalence * (1 - prevalence) / assessed +
                        z^2 / (4 * assessed^2)) / (1 + z^2 / assessed),
    lower = pmax(0, centre - margin),
    upper = pmin(1, centre + margin)
  ) %>%
  ggplot(aes(survey_date, prevalence, colour = treatment,
             fill = treatment, group = treatment)) +
  geom_ribbon(aes(ymin = lower, ymax = upper),
              alpha = 0.15, colour = NA) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 2.1) +
  facet_grid(
    genus ~ .,
    labeller = labeller(genus = c(ACRO = "Acropora", POR = "Porites"))
  ) +
  scale_colour_manual(values = treatment_palette,
                      labels = c(deep = "Deep", shallow = "Shallow")) +
  scale_fill_manual(values = treatment_palette, guide = "none") +
  scale_x_date(
    breaks = as.Date(c("2026-04-17", "2026-05-15", "2026-06-12")),
    labels = c("Apr", "May", "Jun")
  ) +
  scale_y_continuous(
    breaks = seq(0, 0.8, 0.1),
    labels = scales::percent,
    expand = expansion(mult = c(0, 0.02))
  ) +
  coord_cartesian(ylim = c(0, 0.6)) +
  labs(x = "Survey date", y = "Bleaching prevalence", colour = NULL) +
  theme_clean +
  theme(
    strip.text.y = element_text(size = 9, face = "italic"),
    axis.text.x = element_text(size = 8),
    legend.position = "bottom",
    legend.key.size = unit(0.35, "cm"),
    panel.spacing.y = unit(0.8, "lines")
  )


#### D. Assemble and export ####

fig3 <- (f3_stacked + f3_trajectories) +
  patchwork::plot_layout(ncol = 2, widths = c(2, 1)) +
  patchwork::plot_annotation(
    tag_levels = "a",
    tag_prefix = "(",
    tag_suffix = ")",
    theme = theme(plot.tag = element_text(face = "bold", size = 11))
  ) &
  theme(legend.position = "bottom")

print(fig3)
save_figure(fig3, "F3_weekly_bleaching", size = "half")
rm(f3_stacked_data, f3_stacked, f3_trajectories)




### 04. FIGURE 4: SURVIVAL AND MORTALITY ####

#### A. Setup ####

transplant_date <- as.Date("2025-03-14")
f4_date_breaks <- seq(transplant_date, as.Date("2026-06-14"), by = "3 months")
f4_day_breaks <- as.numeric(f4_date_breaks - transplant_date)
f4_date_labels <- format(f4_date_breaks, "%b\n%Y")

#### B. Cause-specific cumulative incidence ####

# First complete-mortality event per coral.
# Aalen–Johansen estimates account for competing mortality causes.

f4_events <- survival_data %>%
  mutate(
    cause = factor(
      if_else(event == 0, "Censored", as.character(death_cause)),
      levels = c("Censored", names(mortality_palette), "COMP", "PRD")
    )
  )

stopifnot(!anyNA(f4_events$cause))

f4_cif_fit <- survival::survfit(
  survival::Surv(survival_days, cause) ~ 1,
  data = f4_events
)

f4_cif_summary <- summary(f4_cif_fit)

f4_cif_data <- as_tibble(f4_cif_summary$pstate) %>%
  mutate(day = f4_cif_summary$time) %>%
  pivot_longer(-day, names_to = "cause", values_to = "incidence") %>%
  filter(cause %in% names(mortality_palette)) %>%
  group_by(cause) %>%
  filter(max(incidence, na.rm = TRUE) > 0) %>%
  ungroup() %>%
  mutate(cause = factor(cause, levels = names(mortality_palette))) %>%
  bind_rows(
    tibble(
      day = 0,
      cause = factor(names(mortality_palette),
                     levels = names(mortality_palette)),
      incidence = 0
    )
  ) %>%
  arrange(cause, day)

#### C. Plot cumulative mortality ####

f4_cif_stacked <- ggplot(
  f4_cif_data,
  aes(day, incidence, fill = cause)
) +
  geom_area(alpha = 0.95, colour = "white", linewidth = 0.2) +
  scale_fill_manual(
    values = mortality_palette,
    breaks = names(mortality_palette),
    labels = mortality_labels,
    drop = TRUE
  ) +
  scale_x_continuous(
    breaks = f4_day_breaks,
    labels = f4_date_labels,
    limits = c(0, 504),
    expand = expansion(mult = c(0, 0.01))
  ) +
  scale_y_continuous(
    breaks = seq(0, 0.6, 0.1),
    labels = scales::percent,
    expand = expansion(mult = c(0, 0))
  ) +
  coord_cartesian(ylim = c(0, 0.6)) +
  guides(fill = guide_legend(ncol = 1)) +
  labs(x = "Date", y = "Cumulative mortality", fill = "Mortality cause") +
  theme_clean +
  theme(
    legend.position = c(0.03, 0.97),
    legend.justification = c(0, 1),
    legend.direction = "vertical",
    legend.background = element_rect(fill = "white", colour = NA),
    legend.title = element_text(size = 9),
    legend.text = element_text(size = 7),
    legend.key.size = unit(0.5, "cm"),
    axis.text = element_text(size = 8)
  )

#### D. Kaplan–Meier survival curves ####

# Descriptive all-cause survival with 95% confidence intervals.
# Censoring marks and risk sets are handled by ggsurvplot.

f4_km_acro <- survminer::ggsurvplot(
  survival::survfit(
    survival::Surv(survival_days, event) ~ treatment,
    data = filter(survival_data, genus == "ACRO")
  ),
  data = filter(survival_data, genus == "ACRO"),
  conf.int = TRUE,
  conf.int.alpha = 0.20,
  censor = TRUE,
  risk.table = FALSE,
  palette = unname(treatment_palette[c("deep", "shallow")]),
  legend.title = " ",
  legend.labs = c("Deep", "Shallow"),
  xlim = c(0, 504),
  ylim = c(0, 1),
  ggtheme = theme_clean
)

f4_km_por <- survminer::ggsurvplot(
  survival::survfit(
    survival::Surv(survival_days, event) ~ treatment,
    data = filter(survival_data, genus == "POR")
  ),
  data = filter(survival_data, genus == "POR"),
  conf.int = TRUE,
  conf.int.alpha = 0.20,
  censor = TRUE,
  risk.table = FALSE,
  palette = unname(treatment_palette[c("deep", "shallow")]),
  legend.title = " ",
  legend.labs = c("Deep", "Shallow"),
  xlim = c(0, 504),
  ylim = c(0, 1),
  ggtheme = theme_clean
)

#### E. Format Kaplan–Meier panels ####

f4_km_acro$plot <- f4_km_acro$plot +
  scale_x_continuous(
    breaks = f4_day_breaks,
    labels = f4_date_labels,
    limits = c(0, 504)
  ) +
  scale_y_continuous(labels = scales::percent) +
  labs(title = "Acropora", x = "Date", y = "Survival probability") +
  theme(plot.title = element_text(face = "italic", size = 10))

f4_km_por$plot <- f4_km_por$plot +
  scale_x_continuous(
    breaks = f4_day_breaks,
    labels = f4_date_labels,
    limits = c(0, 504)
  ) +
  scale_y_continuous(labels = scales::percent) +
  labs(title = "Porites", x = "Date", y = "Survival probability") +
  theme(plot.title = element_text(face = "italic", size = 10))


# Standardise figure typography

f4_km_theme <- theme(
  plot.title = element_text(size = 10, face = "italic"),
  axis.title = element_text(size = 10),
  axis.text = element_text(size = 7),
  legend.title = element_text(size = 9),
  legend.text = element_text(size = 8),
  legend.key.size = unit(0.4, "cm")
)

f4_km_acro$plot <- f4_km_acro$plot + f4_km_theme
f4_km_por$plot <- f4_km_por$plot + f4_km_theme

#### F. Assemble and export ####

# Shared treatment legend for KM panels only.
f4_km_right <- (f4_km_acro$plot / f4_km_por$plot) +
  patchwork::plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

fig4 <- (f4_cif_stacked + f4_km_right) +
  patchwork::plot_layout(widths = c(2, 1)) +
  patchwork::plot_annotation(
    tag_levels = "a",
    tag_prefix = "(",
    tag_suffix = ")",
    theme = theme(plot.tag = element_text(face = "bold", size = 11))
  )

print(fig4)
save_figure(fig4, "F4_survival_mortality", size = "half")
rm(transplant_date, f4_date_breaks, f4_day_breaks, f4_date_labels,
   f4_events, f4_cif_fit, f4_cif_summary, f4_cif_data, f4_cif_stacked,
   f4_km_acro, f4_km_por, f4_km_theme, f4_km_right)

### 05. FIGURE 5: POSITIVE SIZE EXTENSION ####
#### A. Observed temporal variation in positive extension ####
f5_temporal_data <- growth_analysis %>%
  mutate(survey_month = factor(
    format(survey_date, "%b %Y"),
    levels = format(sort(unique(survey_date)), "%b %Y")
  ))

f5_temporal_means <- f5_temporal_data %>%
  group_by(survey_month, genus, treatment) %>%
  summarise(mean_extension = mean(extension_cm_month), .groups = "drop")

f5_temporal <- ggplot(
  f5_temporal_data,
  aes(survey_month, extension_cm_month, colour = treatment)
) +
  geom_point(
    position = position_jitter(width = .2, height = 0),
    size = 0.5, alpha = 0.12
  ) +
  geom_line(
    data = f5_temporal_means,
    aes(y = mean_extension, group = treatment),
    linewidth = 1.5
  ) +
  facet_grid(
    genus ~ .,
    labeller = labeller(genus = c(ACRO = "Acropora", POR = "Porites"))
  ) +
  scale_colour_manual(
    values = treatment_palette,
    breaks = c("deep", "shallow"),
    labels = c("Deep", "Shallow")
  ) +
  scale_x_discrete(labels = \(x) substr(x, 1, 3)) +
  scale_y_continuous(breaks = 0:6) +
  coord_cartesian(ylim = c(0, 6)) +
  labs(x = "Survey date (2025–2026)",
       y = "Positive extension (cm/month)", colour = NULL) +
  theme_clean +
  theme(
    strip.text.y = element_text(size = 9, face = "italic"),
    axis.text.x = element_text(size = 8, hjust = 1),
    legend.position = "none"
  )


#### B. Coral-level distributions and model estimates ####

# One mean positive-extension rate per coral.
# Model estimates are response-scale means from the primary Gamma GLMM.

f5_model_estimates <- emmeans(
  growth_coral_gamma, ~ treatment | genus
) %>%
  summary(type = "response", infer = c(TRUE, FALSE)) %>%
  as_tibble()

f5_distributions <- ggplot(
  growth_coral,
  aes(treatment, mean_positive_extension, fill = treatment)
) +
  geom_violin(
    alpha = 0.25, width = 0.75,
    colour = NA, trim = TRUE
  ) +
  geom_jitter(
    aes(colour = treatment),
    width = 0.12, height = 0,
    size = 0.65, alpha = 0.25
  ) +
  geom_pointrange(
    data = f5_model_estimates,
    aes(x = treatment, y = response,
        ymin = asymp.LCL, ymax = asymp.UCL),
    inherit.aes = FALSE,
    colour = "black",
    linewidth = 0.7,
    fatten = 2.5
  ) +
  facet_grid(
    genus ~ .,
    labeller = labeller(genus = c(ACRO = "Acropora", POR = "Porites"))
  ) +
  scale_colour_manual(values = treatment_palette, guide = "none") +
  scale_fill_manual(values = treatment_palette, guide = "none") +
  scale_y_log10(
    breaks = c(0.25, 0.5, 1, 2, 5, 10),
    labels = c("0.25", "0.5", "1", "2", "5", "10")
  ) +
  scale_x_discrete(labels = c(deep = "Deep", shallow = "Shallow")) +
  labs(x = NULL, y = "Mean positive extension (cm/month, log scale)") +
  theme_clean +
  theme(
    strip.text = element_text(size = 9, face = "italic"),
    axis.text = element_text(size = 8),
    panel.spacing = unit(0.7, "lines")
  )

#### C. Assemble and export ####

fig5 <- (f5_temporal + f5_distributions) +
  patchwork::plot_layout(widths = c(2, 1)) +
  patchwork::plot_annotation(
    tag_levels = "a",
    tag_prefix = "(",
    tag_suffix = ")",
    theme = theme(plot.tag = element_text(face = "bold", size = 11))
  )

print(fig5)
save_figure(fig5, "F5_positive_extension", size = "half")
rm(f5_temporal_data, f5_temporal_means, f5_temporal, f5_model_estimates,
   f5_distributions)

### 06. SUPPLEMENTARY FIGURE S1: CAUSE-SPECIFIC MORTALITY ####

#### A. Cumulative incidence by donor treatment ####

# First observed cause of complete mortality; competing events retained.
s1_events <- survival_data %>%
  mutate(cause = factor(
    if_else(event == 0, "Censored", as.character(death_cause)),
    levels = c("Censored", names(mortality_palette), "COMP", "PRD")
  ))

stopifnot(!anyNA(s1_events$cause))

s1_cif <- survival::survfit(
  survival::Surv(survival_days, cause) ~ treatment,
  data = s1_events
)

s1_cif_data <- summary(s1_cif) %>%
  (\(s) {
    as_tibble(s$pstate) %>%
      mutate(day = s$time, treatment = sub("treatment=", "", s$strata)) %>%
      pivot_longer(all_of(names(mortality_palette)),
                   names_to = "cause", values_to = "incidence")
  })() %>%
  bind_rows(tidyr::expand_grid(
    treatment = c("deep", "shallow"),
    cause = names(mortality_palette), day = 0, incidence = 0
  )) %>%
  mutate(cause = factor(cause, levels = names(mortality_palette)))

#### B. Assemble and export ####

figS1 <- ggplot(s1_cif_data, aes(day, incidence, fill = cause)) +
  geom_area(colour = "white", linewidth = 0.2) +
  facet_wrap(~ treatment, ncol = 2) +
  scale_fill_manual(values = mortality_palette, labels = mortality_labels) +
  scale_y_continuous(labels = scales::percent, limits = c(0, 0.6)) +
  labs(x = "Days since transplantation", y = "Cumulative mortality",
       fill = "Mortality cause") +
  theme_clean

print(figS1)
save_figure(figS1, "S1_mortality_by_treatment", size = "half")
rm(s1_events, s1_cif, s1_cif_data)

### print figures #### 
print(fig2)
print(fig3)
print(fig4)
print(fig5)
print(figS1)
 