#### Bootstrapping Data on Percentage of Cells Per Zone 
#### Created by Aurelia Leona, 2025-01-28
#### Adapted from [original source if applicable]

# Load required libraries
required_packages <- c("readr","dplyr","ggplot2","purrr")

# Install missing packages
new_packages <- required_packages[!(required_packages %in% installed.packages()[,"Package"])]
if(length(new_packages)) install.packages(new_packages)

# Load all packages
lapply(required_packages, library, character.only = TRUE)

# Load necessary libraries
library(readr)
library(dplyr)
library(ggplot2)
library(purrr)
library(tidyr)
library(rstatix)  # for tidy statistical tests


# Initialize directories
getwd()
rawfileDir <- "/projects/b1042/GoyalLab/aleona/Zonation_Analysis/rawData/"
plotDir <- "/projects/b1042/GoyalLab/aleona/Zonation_Analysis/plots/"
extrdataDir <- "/projects/b1042/GoyalLab/aleona/Zonation_Analysis/extractedData/"
plotscriptsdir <- "/projects/b1042/GoyalLab/aleona/Zonsation_Analysis/plotScripts/"
extScriptsdir <- "/projects/b1042/GoyalLab/aleona/Zonation_Analysis/extractScripts/"

source(paste0(extScriptsdir, "etazonefunction.R"))

# Import the combined eta zones data
etawithzones_directory <- paste0(extrdataDir, "etawithzones_data/")
etazones <- read_csv(paste0(etawithzones_directory, "combined_etawithzones.csv"))

#Define new Directories 
svg_directory <- paste0(plotDir, "svg/")
png_directory <- paste0(plotDir, "png/")
# Perform tests for each zone
bootstrapdir <- paste0(extrdataDir, "bootstrap_results/")

# Create directories
create_dir_if_not_exists(svg_directory)
create_dir_if_not_exists(png_directory)
create_dir_if_not_exists(bootstrapdir)

# Define conditions and zones
list_name <- c("Normal", "AC", "AH")
list_of_zones <- c("Zone 1", "Zone 2", "Zone 3")
subsampling <- c(20,40,60)


######Bootstrapping functions

## Wilcoxon test for Zone distribution between conditions
# Function to perform bootstrap analysis with different subsample size 
bootstrap_zone_comparison <- function(data, n_iterations = 500, subsample_prop = 0.7) {
  bootstrap_results <- list()
  
  for(i in 1:n_iterations) {
    # Sample with replacement using specified proportion
    boot_data <- data %>%
      group_by(Condition) %>%
      slice_sample(prop = subsample_prop, replace = TRUE) %>%
      ungroup()
    
    # Calculate zone percentages for each condition
    zone_props <- boot_data %>%
      group_by(Condition, Zone) %>%
      summarise(count = n(), .groups = "drop") %>%
      group_by(Condition) %>%
      mutate(percentage = (count / sum(count)) * 100) %>%
      ungroup() %>%
      select(Condition, Zone, percentage)
    
    # Store results with iteration number
    bootstrap_results[[i]] <- zone_props %>%
      mutate(iteration = i)
  }
  
  # Combine all results
  all_results <- bind_rows(bootstrap_results)
  
  # Reshape data for comparison (wide format by Zone)
  comparison_data <- all_results %>%
    pivot_wider(
      id_cols = c(iteration, Zone),
      names_from = Condition,
      values_from = percentage
    )
  
  # Perform Wilcoxon tests for each zone with NA handling
  wilcoxon_results <- comparison_data %>%
    group_by(Zone) %>%
    summarise(
      normal_vs_ac_p = wilcox.test(Normal, AC, paired = TRUE, na.rm = TRUE)$p.value,
      normal_vs_ah_p = wilcox.test(Normal, AH, paired = TRUE, na.rm = TRUE)$p.value,
      normal_mean = mean(Normal, na.rm = TRUE),
      ac_mean = mean(AC, na.rm = TRUE),
      ah_mean = mean(AH, na.rm = TRUE),
      normal_ci_lower = quantile(Normal, 0.025, na.rm = TRUE),
      normal_ci_upper = quantile(Normal, 0.975, na.rm = TRUE),
      ac_ci_lower = quantile(AC, 0.025, na.rm = TRUE),
      ac_ci_upper = quantile(AC, 0.975, na.rm = TRUE),
      ah_ci_lower = quantile(AH, 0.025, na.rm = TRUE),
      ah_ci_upper = quantile(AH, 0.975, na.rm = TRUE),
      .groups = "drop"
    )
  
  return(list(
    all_iterations = all_results,
    comparison_data = comparison_data,
    summary = wilcoxon_results
  ))
}

# Create empty dataframe to store final results
combined_bootstrap <- data.frame(
  Condition = character(),      # Will contain "Normal", "AC", "AH"
  Zone = character(),          # Will contain "Zone 1", "Zone 2", "Zone 3"
  mean_prop = numeric(),       # Mean percentage for each condition/zone
  lower_ci = numeric(),        # Lower bound of confidence interval
  upper_ci = numeric(),        # Upper bound of confidence interval
  p_value = numeric(),         # P-value comparing to Normal condition
  Subsample = character(),     # Will contain "20%", "40%", "60%"
  stringsAsFactors = FALSE     # Prevent automatic factor conversion
)


# Loop through each subsampling percentage (20%, 40%, 60%)
for (prop in subsampling) {
  # Get bootstrap results for this subsample size
  results <- bootstrap_zone_comparison(etazones, subsample_prop = prop/100)
  
  # Create filenames with prop value
  csvfile_iterations <- paste0(bootstrapdir, "bootstrap_iterations_", prop, ".csv")
  csvfile_comparison <- paste0(bootstrapdir, "bootstrap_comparisondata_", prop, ".csv")
  csvfile_summary <- paste0(bootstrapdir, "bootstrap_summary_", prop, ".csv")
  
  # Write results
  write.csv(results$all_iterations, csvfile_iterations, row.names = FALSE)
  write.csv(results$comparison_data, csvfile_comparison, row.names = FALSE)
  write.csv(results$summary, csvfile_summary, row.names = FALSE)
  
  # Format results for combining:
  summary_formatted <- results$summary %>%
    
    # Step 1: Convert wide to long format for means
    pivot_longer(
      cols = c(normal_mean, ac_mean, ah_mean), # Column contain means
      names_to = "temp",                       # Temporary Column Name 
      values_to = "mean_prop"                  # Column for mean values
    ) %>%
    
    # Step 2: Create proper condition labels and match values
    mutate(
      Condition = case_when(
        # Convert temporary column names to condition labels
        temp == "normal_mean" ~ "Normal",     
        temp == "ac_mean" ~ "AC",
        temp == "ah_mean" ~ "AH"
      ),
      # Match lower CI bounds to correct condition
      lower_ci = case_when(
        Condition == "Normal" ~ normal_ci_lower,
        Condition == "AC" ~ ac_ci_lower,
        Condition == "AH" ~ ah_ci_lower
      ),
      # Match upper CI bounds to correct condition
      upper_ci = case_when(
        Condition == "Normal" ~ normal_ci_upper,
        Condition == "AC" ~ ac_ci_upper,
        Condition == "AH" ~ ah_ci_upper
      ),
      # Add p-values (comparing each treatment to Normal)
      p_value = case_when(
        Condition == "AC" ~ normal_vs_ac_p,
        Condition == "AH" ~ normal_vs_ah_p,
        TRUE ~ NA_real_  # Normal gets NA (reference)
      ),
      # Add subsample percentage label
      Subsample = paste0(prop, "%")
      
    ) %>%
    # Step 3: Select only the columns we want in final output
    select(Condition, Zone, mean_prop, lower_ci, upper_ci, p_value, Subsample)
  
  # Add to combined_bootstrap
  combined_bootstrap <- bind_rows(combined_bootstrap, summary_formatted)
}


# results70 <- bootstrap_zone_comparison(etazones)
# write.csv(results70$all_iterations, "bootstrap_iterations70.csv", row.names = FALSE)
# write.csv(results70$comparison_data, "bootstrap_comparisondata70.csv", row.names = FALSE)
# write.csv(results70$summary, "bootstrap_summary70.csv", row.names = FALSE)
# 
# results50 <- bootstrap_zone_comparison(etazones, subsample_prop = 0.5)
# write.csv(results50$all_iterations, "bootstrap_iterations50.csv", row.names = FALSE)
# write.csv(results50$comparison_data, "bootstrap_comparisondata50.csv", row.names = FALSE)
# write.csv(results50$summary, "bootstrap_summary50.csv", row.names = FALSE)
# 
# results90 <- bootstrap_zone_comparison(etazones, subsample_prop = 0.9)
# write.csv(results50$all_iterations, "bootstrap_iterations90.csv", row.names = FALSE)
# write.csv(results90$comparison_data, "bootstrap_comparisondata90.csv", row.names = FALSE)
# write.csv(results50$summary, "bootstrap_summary90.csv", row.names = FALSE)



##Plot all the data above
# Combine all bootstrap results into one dataframe
# Order by Normal, AC and AH 
combined_bootstrap$Condition <- factor(combined_bootstrap$Condition, 
                                       levels = c("Normal", "AC", "AH"))

#Plot the data 
ggplot(combined_bootstrap, 
       aes(x = Zone, y = mean_prop, 
           color = Condition, 
           shape = Subsample)) +
  geom_point(position = position_dodge(width = 0.5),
             size = 3) +
  geom_errorbar(aes(ymin = lower_ci, ymax = upper_ci),
                position = position_dodge(width = 0.5),
                width = 0.2) +
  scale_color_manual(values = c("Normal" = "#224b5e", 
                                "AC" = "#94b594", 
                                "AH" = "#edc775")) +
  labs(title = "Zone Distribution Across Conditions",
       subtitle = "Bootstrap Results with Different Sampling Proportions",
       x = "Hepatic Zone",
       y = "Percentage of Cells (%)") +
  theme_classic() +
  theme(legend.position = "top",
        plot.title = element_text(size = 14),
        plot.subtitle = element_text(size = 12),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 12))

# Save the plot
ggsave(filename = file.path(png_directory, "bootstrapped_zone_distribution_dots.png"), 
       width = 10, height = 6, dpi = 300)
ggsave(filename = file.path(svg_directory, "bootstrapped_zone_distribution_dots.svg"), 
       width = 8, height = 6, dpi = 300)


