
    
    ### 01. CUMULATIVE MORTALITY BY TREATMENT ####
    
    scratch_cif <- survival_data %>%
      mutate(
        cause = factor(
          if_else(event == 0, "Censored", as.character(death_cause)),
          levels = c("Censored", names(mortality_palette), "COMP", "PRD")
        )
      )
    
    scratch_fit <- survival::survfit(
      survival::Surv(survival_days, cause) ~ treatment,
      data = scratch_cif
    )
    
    scratch_cif_plot <- summary(scratch_fit) %>%
      (\(s) {
        as_tibble(s$pstate) %>%
          mutate(day = s$time, treatment = sub("treatment=", "", s$strata)) %>%
          pivot_longer(all_of(names(mortality_palette)),
                       names_to = "cause", values_to = "incidence")
      })() %>%
      bind_rows(
        tidyr::expand_grid(
          treatment = c("deep", "shallow"),
          cause = names(mortality_palette),
          day = 0, incidence = 0
        )
      ) %>%
      mutate(cause = factor(cause, levels = names(mortality_palette))) %>%
      ggplot(aes(day, incidence, fill = cause)) +
      geom_area(colour = "white", linewidth = 0.2) +
      facet_wrap(~ treatment, ncol = 2) +
      scale_fill_manual(values = mortality_palette, labels = mortality_labels) +
      scale_y_continuous(labels = scales::percent, limits = c(0, 0.6)) +
      labs(x = "Days since transplantation", y = "Cumulative mortality",
           fill = "Mortality cause") +
      theme_clean
    
    print(scratch_cif_plot)
    
    
    survival_data %>%
      group_by(treatment) %>%
      summarise(
        corals = n(),
        bleaching_deaths = sum(event == 1 & death_cause == "BLE"),
        pct_all_corals = 100 * bleaching_deaths / corals,
        pct_of_deaths = 100 * bleaching_deaths / sum(event == 1),
        .groups = "drop"
      ) %>%
      print()
    
    
    
    #### A. Observed warming-season bleaching heatmap ####
    
    # Substantial bleaching (scores 3–5); darker colours indicate lower prevalence.
    
    f3_heat_data <- bleach_warming %>%
      mutate(
        year = factor(lubridate::year(survey_date)),
        month = factor(format(survey_date, "%b"), levels = month.abb[3:7]),
        treatment = factor(treatment, levels = c("deep", "shallow"),
                           labels = c("Deep", "Shallow"))
      ) %>%
      group_by(year, month, genus, treatment) %>%
      summarise(
        assessed = n(),
        prevalence = mean(bleached),
        .groups = "drop"
      )
    
    f3_heatmap <- ggplot(
      f3_heat_data,
      aes(month, treatment, fill = prevalence)
    ) +
      geom_tile(colour = "white", linewidth = 0.8) +
      facet_grid(
        genus ~ year,
        labeller = labeller(genus = c(ACRO = "Acropora", POR = "Porites"))
      ) +
      scale_fill_gradientn(
        colours = unname(condition_palette[as.character(0:5)]),
        limits = c(0, 0.4),
        breaks = seq(0, 0.4, 0.1),
        labels = scales::percent,
        oob = scales::squish
      ) +
      labs(x = "Month", y = NULL, fill = "Bleaching prevalence") +
      theme_clean +
      theme(
        strip.text.x = element_text(size = 9),
        strip.text.y = element_text(size = 9, face = "italic"),
        axis.text = element_text(size = 8),
        legend.position = "bottom",
        legend.title = element_text(size = 9),
        legend.text = element_text(size = 8),
        panel.spacing = unit(0.7, "lines")
      )
    
    print(f3_heatmap)
    
    ### 02. FIGURE 3: MONTHLY BLEACHING ####
    #### A. Observed warming-season bleaching heatmap ####
    
    # Observed prevalence of substantial bleaching (scores 3–5).
    # Dark = no bleaching; light = higher bleaching prevalence.
    
    f3_heat_data <- bleach_warming %>%
      mutate(
        year = factor(lubridate::year(survey_date)),
        month = factor(format(survey_date, "%b"), levels = month.abb[3:7]),
        treatment = factor(treatment, levels = c("deep", "shallow"),
                           labels = c("Deep", "Shallow"))
      ) %>%
      group_by(year, month, genus, treatment) %>%
      summarise(
        assessed = n(),
        prevalence = mean(bleached),
        .groups = "drop"
      )
    
    f3_heatmap <- ggplot(
      f3_heat_data,
      aes(month, treatment, fill = prevalence)
    ) +
      geom_tile(colour = "white", linewidth = 0.8) +
      facet_grid(
        genus ~ year,
        labeller = labeller(genus = c(ACRO = "Acropora", POR = "Porites"))
      ) +
      scale_fill_gradientn(
        colours = unname(condition_palette[as.character(0:5)]),
        limits = c(0, 0.4),
        breaks = seq(0, 0.4, 0.1),
        labels = scales::percent,
        oob = scales::squish,
        guide = guide_colourbar(
          direction = "vertical",
          title.position = "top",
          title.hjust = 0.5,
          barheight = unit(3.2, "cm"),
          barwidth  = unit(0.35, "cm")
        )
      ) +
      labs(x = "Month", y = NULL, fill = "Bleaching\nprevalence") +
      theme_clean +
      theme(
        strip.text.x = element_text(size = 9),
        strip.text.y = element_text(size = 9, face = "italic"),
        axis.text = element_text(size = 8),
        legend.position = "left",
        legend.direction = "vertical",
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 8),
        panel.spacing = unit(0.7, "lines"),
        legend.margin = margin(0, 2, 0, 0),
        legend.box.margin = margin(0, 0, 0, 0)
      )
    
    print(f3_heatmap)
    
    
    
    
    
    
    
    #### B. Temperature during warming periods ####
    
    f3_temperature <- temp_daily %>%
      filter(month(date) %in% 3:7) %>%
      mutate(
        year = factor(year(date)),
        day_of_season = as.numeric(date - as.Date(paste0(year, "-03-01"))),
        segment = cumsum(!complete_day)
      ) %>%
      filter(complete_day) %>%
      ggplot(aes(day_of_season, group = interaction(year, segment))) +
      geom_ribbon(
        aes(ymin = temp_min, ymax = temp_max),
        fill = condition_palette[["2"]], alpha = 0.5
      ) +
      geom_line(aes(y = temp_mean),
                colour = condition_palette[["0"]], linewidth = 0.7) +
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
      theme(
        strip.text = element_text(size = 9),
        axis.text = element_text(size = 8)
      )
    
    print(f3_temperature)
    
    
    
    
    #### C. Population-marginal bleaching probabilities ####
    
    f3_probabilities <- bleach_marginal_emms %>%
      mutate(year = stringr::str_extract(as.character(warming_season), "\\d{4}")) %>%
      ggplot(aes(year, probability, colour = treatment)) +
      geom_pointrange(
        aes(ymin = lower, ymax = upper),
        position = position_dodge(width = 0.45),
        linewidth = 0.7
      ) +
      facet_wrap(
        ~ genus,
        labeller = labeller(genus = c(ACRO = "Acropora", POR = "Porites"))
      ) +
      scale_colour_manual(values = treatment_palette) +
      scale_y_continuous(
        labels = scales::percent,
        limits = c(0, NA),
        expand = expansion(mult = c(0, 0.08))
      ) +
      labs(x = NULL, y = "Bleaching probability") +
      theme_clean +
      theme(
        strip.text = element_text(size = 9, face = "italic"),
        legend.position = "none",
        axis.text.x = element_text(size = 7, angle = 45, hjust = 1),
        panel.spacing = unit(0.5, "lines")
      )
    
    
    #### D. Marginal treatment differences ####
    
    # Shallow minus deep, expressed in percentage points.
    # Zero indicates no difference in marginal bleaching probability.
    
    f3_contrasts <- bleach_marginal_contrasts %>%
      mutate(
        label = paste(
          recode(as.character(genus), ACRO = "Acropora\n", POR = "Porites\n"),
          stringr::str_extract(as.character(warming_season), "\\d{4}")
        ),
        label = factor(label, levels = rev(unique(label)))
      ) %>%
      ggplot(aes(difference, label)) +
      geom_vline(xintercept = 0, linetype = "dashed",
                 colour = "grey60", linewidth = 0.5) +
      geom_segment(
        aes(x = lower, xend = upper, yend = label),
        linewidth = 0.8, colour = treatment_palette[["deep"]]
      ) +
      geom_point(size = 2.7, colour = treatment_palette[["deep"]]) +
      scale_x_continuous(
        labels = scales::label_number(accuracy = 1)
      ) +      
      labs(x = "Shallow − deep (pp)", y = NULL) +
      scale_y_discrete(labels = function(x) {
        genus <- stringr::str_extract(x, "Acropora|Porites")
        year <- stringr::str_extract(x, "\\d{4}")
        paste0("*", genus, "*<br>", year)
      }) +
      theme_clean +
      theme(
        axis.text.y.left = ggtext::element_markdown(
          size = 8, lineheight = 0.85, hjust = 1
        )
      )
    
    f3_contrasts 
    
    
    
    #### E. Assemble and export ####
    
    f3_left <- f3_heatmap / f3_temperature +
      patchwork::plot_layout(heights = c(1, 1))
    
    f3_right <- f3_probabilities / f3_contrasts +
      patchwork::plot_layout(heights = c(1, 1))
    
    fig3 <- (f3_left | f3_right) +
      patchwork::plot_layout(widths = c(3, 1)) +
      patchwork::plot_annotation(
        tag_levels = "a",
        tag_prefix = "(",
        tag_suffix = ")",
        theme = theme(plot.tag = element_text(face = "bold", size = 11))
      )
    
    print(fig3)
    
    
    ##  #or #### 
    #### E. Assemble and export ####
    
    fig2 <- (f3_heatmap / f3_temperature) +
      patchwork::plot_layout(heights = c(1.1, 1)) +
      patchwork::plot_annotation(
        tag_levels = "a",
        tag_prefix = "(",
        tag_suffix = ")",
        theme = theme(plot.tag = element_text(face = "bold", size = 11))
      )
    
    print(fig2)
    save_figure(fig2, "F2_bleaching_temperature", size = "half")
    