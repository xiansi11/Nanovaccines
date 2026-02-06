##------------------------------------------------------###
##---------------------Fig.S4-----------------------------####
##------------------------------------------------------###

library(readr)
library(here)
library(dplyr)
library(ggplot2)
library(lubridate)
library(xtable)
library(tidyr)
library(lme4)
OCR_data <- read_csv("SHSplenicDC_cleanedOCR data.csv")



library(purrr)
library(patchwork)  # For combining plots
# Split data into treatment and control groups
treatment_data <- OCR_data %>% 
  filter(Treatment != "Control")

control_data <- OCR_data %>% 
  filter(Treatment == "Control")

# Find matching controls for each treatment
matched_pairs <- treatment_data %>%
  inner_join(control_data,
             by = c("Age", "ExpID", "Date", "PlateDay", "Plate"),
             suffix = c("_treatment", "_control")) %>%
  select(
    # Treatment information
    starts_with("Treatment"),
    starts_with("Age"),
    starts_with("ExpID"),
    starts_with("Date"),
    starts_with("Plate"),
    # Keep other columns if needed
    matches("_treatment|_control"))


# Split into list by treatment
split_list <- split(matched_pairs, 
                    matched_pairs$Treatment_treatment)


split_list_log <- lapply(split_list, function(treatment_df) {
  min_spare <- min(c(treatment_df$Spare_treatment[treatment_df$Spare_treatment > 0],
                     treatment_df$Spare_control[treatment_df$Spare_control > 0])) / 2
  
  min_leak <- min(c(treatment_df$Leak_treatment[treatment_df$Leak_treatment > 0],
                    treatment_df$Leak_control[treatment_df$Leak_control > 0])) / 2
  treatment_df <- treatment_df %>%
    mutate(
      Spare_treatment = ifelse(Spare_treatment == 0, min_spare, Spare_treatment+min_spare),
      Spare_control = ifelse(Spare_control == 0, min_spare, Spare_control+min_spare),
      Leak_treatment = ifelse(Leak_treatment == 0, min_leak, Leak_treatment+min_leak),
      Leak_control = ifelse(Leak_control == 0, min_leak, Leak_control+min_leak)
    )
  treatment_df %>%
    mutate(
      log_Basal_ratio = log2(Basal_treatment / Basal_control),
      log_ATP_ratio = log2(ATP_treatment / ATP_control),
      log_Max_ratio = log2(Max_treatment / Max_control),
      log_Spare_ratio = log2(Spare_treatment / Spare_control),
      log_Leak_ratio = log2(Leak_treatment / Leak_control),
      log_NonMito_ratio = log2(NonMito_treatment / NonMito_control)
      # log_ECARBasal_ratio = log2(ECARBasal_treatment / ECARBasal_control)
    )
})


treatment_plots <- imap(split_list_log, function(treatment_df, treatment_name) {
  long_data <- treatment_df %>%
    pivot_longer(
      cols = starts_with("log_"),
      names_to = "metric",
      values_to = "log_ratio"
    ) %>%
    mutate(
      metric = factor(metric,
                      levels = c("log_Basal_ratio", "log_ATP_ratio",
                                 "log_Max_ratio", "log_Spare_ratio",
                                 "log_Leak_ratio", "log_NonMito_ratio"),
                      labels = c("Basal", "ATP", "Max", 
                                 "Spare", "Leak", "NonMito")),
      Age = factor(Age)
    )
  
  # Create individual plots for each metric
  plot_list <- map(levels(long_data$metric), \(m) {
    data_subset <- filter(long_data, metric == m)
    
    ggplot(data_subset, aes(x = Age, y = log_ratio, color = Age, shape = Age)) +
      geom_jitter(width = 0, size = 2, alpha = 0.8) +
      labs(
        title = m,  # Add metric name as the subplot title
        y = "",
        x = NULL
      ) +
      theme_bw(base_size = 8) +
      theme(
        legend.position = "none",
        panel.grid.minor = element_blank(),
        plot.title = element_text(size = 8)  # Adjust subplot title size
      ) +
      scale_color_manual(values = c("Aged" = "blue", "Young" = "red")) +
      scale_shape_manual(values = c("Aged" = 1, "Young" = 2))
  })
  
  # Combine all 6 plots into one, add overall treatment title
  combined_plot <- wrap_plots(plot_list, nrow = 1) +  # Arrange in 1 row
    plot_annotation(
      title = paste("Treatment:", treatment_name),
      theme = theme(plot.title = element_text(hjust = 0.5, size = 12))
    ) 
  return(combined_plot)
})

# View the combined plots
treatment_plots[[1]]