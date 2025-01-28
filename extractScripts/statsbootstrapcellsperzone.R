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


######Bootstrapping functions 
# Function to perform bootstrap with different subsample sizes

## Wilcoxon test for Zone distribution between conditions
# Function to perform bootstrap analysis with zone percentages
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
  
  
  # Perform Wilcoxon tests for each zone
  wilcoxon_results <- comparison_data %>%
    group_by(Zone) %>%
    summarise(
      normal_vs_ac_p = wilcox.test(Normal, AC, paired = TRUE)$p.value,
      normal_vs_ah_p = wilcox.test(Normal, AH, paired = TRUE)$p.value,
      normal_mean = mean(Normal),
      ac_mean = mean(AC),
      ah_mean = mean(AH),
      normal_ci_lower = quantile(Normal, 0.025),
      normal_ci_upper = quantile(Normal, 0.975),
      ac_ci_lower = quantile(AC, 0.025),
      ac_ci_upper = quantile(AC, 0.975),
      ah_ci_lower = quantile(AH, 0.025),
      ah_ci_upper = quantile(AH, 0.975),
      .groups = "drop"
    )
  
  return(list(
    all_iterations = all_results,
    comparison_data = comparison_data,
    summary = wilcoxon_results
  ))
}


# plot_zone_distributions <- function(bootstrap_results) {
#   ggplot(bootstrap_results$all_iterations, 
#          aes(x = percentage, fill = Condition)) +
#     geom_density(alpha = 0.5) +
#     facet_wrap(~Zone, scales = "free") +
#     theme_minimal() +
#     labs(title = "Distribution of Zone Percentages Across Bootstrap Iterations",
#          x = "Percentage",
#          y = "Density")
# }
# 


#Indicate subsampling data 
subsampling <- c(20, 40, 60) 
combined_bootstrap <- 


for (prop in subsampling){
  decimal_ver <- prop/100
  results <- bootstrap_zone_comparison(etazones, subsample_prop = decimal_ver )
  
  csvfile_name <- paste0(bootstrapdir, "bootstrap_iterations", subsampling, ".csv")
  write.csv(results$all_iterations, csvfile_name, row.names = FALSE)
  csvfile_name <- paste0(bootstrapdir, "bootstrap_comparisondata", subsampling, ".csv")
  write.csv(results$comparison_data, csvfile_name, row.names = FALSE)
  csvfile_name <- paste0(bootstrapdir, "bootstrap_summary", subsampling, ".csv")
  write.csv(results$summary, csvfile_name, row.names = FALSE)
  
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
combined_bootstrap <- bind_rows(
  # 50% results
  mutate(all_results$`0.5`$ci_results, Subsample = "50%"),
  # 70% results
  mutate(all_results$`0.7`$ci_results, Subsample = "70%"),
  # 90% results
  mutate(all_results$`0.9`$ci_results, Subsample = "90%")
) %>%
  # Convert proportions to percentages
  mutate(
    mean_prop = mean_prop * 100,
    lower_ci = lower_ci * 100,
    upper_ci = upper_ci * 100
  )


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
ggsave(filename = file.path(png_directory, "zone_distribution_dots.png"), 
       width = 10, height = 6, dpi = 300)
ggsave(filename = file.path(svg_directory, "zone_distribution_dots.svg"), 
       width = 8, height = 6, dpi = 300)



##EXTRA -----------------
# bootstrap_with_subsampling <- function(data, subsample_prop = 0.7, n_iterations = 1000) {
#   bootstrap_results <- list()
#   all_data <- list()  # New list to store all iteration data
#   
#   # Get original sample sizes for reference
#   original_sizes <- data %>%
#     group_by(Condition) %>%
#     summarise(n = n())
#   
#   print("Original sample sizes:")
#   print(original_sizes)
#   
#   for(i in 1:n_iterations) {
#     # Sample with replacement using specified proportion
#     boot_data <- data %>%
#       group_by(Condition) %>%
#       slice_sample(prop = subsample_prop, replace = TRUE) %>%
#       ungroup()
#     
#     # Store the complete bootstrap sample
#     all_data[[i]] <- boot_data
#     
#     # Calculate proportions for this bootstrap sample
#     props <- boot_data %>%
#       group_by(Condition, Zone) %>%
#       summarise(Count = n(), .groups = 'drop') %>%
#       group_by(Condition) %>%
#       mutate(CountNorm = Count / sum(Count)) %>%
#       ungroup() %>%
#       mutate(
#         Condition = factor(Condition, levels = list_name),
#         Zone = factor(Zone, levels = list_of_zones)
#       )
#     
#     bootstrap_results[[i]] <- props
#   }
#   
#   # Combine all bootstrap results
#   bootstrap_df <- bind_rows(bootstrap_results, .id = "iteration")
#   
#   # Combine all raw data with iteration numbers
#   all_data_df <- bind_rows(all_data, .id = "iteration")
#   
#   # Calculate confidence intervals
#   ci_results <- bootstrap_df %>%
#     group_by(Condition, Zone) %>%
#     summarise(
#       mean_prop = mean(CountNorm),
#       lower_ci = quantile(CountNorm, 0.025),
#       upper_ci = quantile(CountNorm, 0.975),
#       .groups = 'drop'
#     )
#   
#   return(list(
#     ci_results = ci_results,
#     subsample_sizes = original_sizes %>% 
#       mutate(subsample_size = round(n * subsample_prop)),
#     all_bootstrap_data = all_data_df,     # Return all bootstrapped data
#     bootstrap_summaries = bootstrap_df     # Return all summary statistics
#   ))
# }
# 


# Create a directory for text output if it doesn't exist
txt_directory <- paste0(extrdataDir, "bootstrap_results/")
create_dir_if_not_exists(txt_directory)

# Function to write results to file
write_bootstrap_results <- function(results, subsample_prop, file_path) {
  # Open connection to file
  con <- file(file_path, "w")
  
  # Write header
  writeLines(sprintf("Bootstrap Results for %d%% Subsampling\n", subsample_prop * 100), con)
  
  # Write original sample sizes
  writeLines("Original Sample Sizes:", con)
  capture.output(results$subsample_sizes, file = con, append = TRUE)
  writeLines("\n", con)
  
  # Write confidence intervals
  writeLines("Confidence Intervals:", con)
  ci_formatted <- results$ci_results %>%
    arrange(Condition, Zone) %>%
    mutate(across(where(is.numeric), ~round(.*100, 1))) %>%
    mutate(
      CI = sprintf("%.1f%% (%.1f%% - %.1f%%)", 
                   mean_prop, lower_ci, upper_ci)
    )
  capture.output(ci_formatted, file = con, append = TRUE)
  
  # Close connection
  close(con)
}

# Try different subsampling proportions
subsample_proportions <- c(0.5, 0.7, 0.9)
all_results <- list()

for(prop in subsample_proportions) {
  set.seed(123) # for reproducibility
  results <- bootstrap_with_subsampling(etazones, subsample_prop = prop)
  all_results[[as.character(prop)]] <- results
  
  # Write results to file
  file_path <- paste0(txt_directory, sprintf("bootstrap_results_%dpercent.txt", prop * 100))
  write_bootstrap_results(results, prop, file_path)
  
  # Also save as CSV for potential further analysis
  csv_path <- paste0(txt_directory, sprintf("bootstrap_ci_%dpercent.csv", prop * 100))
  write_csv(results$ci_results, csv_path)
}

# Create summary file with comparison across all proportions
summary_file <- paste0(txt_directory, "bootstrap_summary_comparison.txt")
con <- file(summary_file, "w")

writeLines("Comparison of Bootstrap Results Across Different Subsampling Proportions\n", con)

for(prop in subsample_proportions) {
  results <- all_results[[as.character(prop)]]
  
  writeLines(sprintf("\n%d%% Subsampling Results:", prop * 100), con)
  writeLines("\nSample Sizes:", con)
  capture.output(results$subsample_sizes, file = con, append = TRUE)
  
  writeLines("\nConfidence Intervals:", con)
  ci_formatted <- results$ci_results %>%
    arrange(Condition, Zone) %>%
    mutate(across(where(is.numeric), ~round(.*100, 1))) %>%
    mutate(
      CI = sprintf("%.1f%% (%.1f%% - %.1f%%)", 
                   mean_prop, lower_ci, upper_ci)
    )
  capture.output(ci_formatted, file = con, append = TRUE)
  writeLines("\n-----------------------------------\n", con)
}

close(con)

# Compare results visually
plot_subsample_comparison <- function(all_results, subsample_prop) {
  results <- all_results[[as.character(subsample_prop)]]
  
  ggplot(results$ci_results, aes(x = Condition, y = mean_prop, fill = Zone)) +
    geom_col(position = "stack") +
    geom_errorbar(aes(ymin = lower_ci, ymax = upper_ci),
                  position = position_stack(vjust = 0.5),
                  width = 0.2) +
    geom_text(aes(label = sprintf("%.1f%%", mean_prop*100)),
              position = position_stack(vjust = 0.5),
              color = "black", size = 4) +
    labs(x = "Conditions", 
         y = "Percentage of Cells", 
         title = sprintf("Distribution with %d%% Subsampling", subsample_prop * 100)) +
    scale_fill_manual(values = palette) +
    scale_y_continuous(labels = scales::percent) +
    theme_classic() +
    theme(legend.position = "right")
}
# Create plots as before
for(prop in subsample_proportions) {
  p <- plot_subsample_comparison(all_results, prop)
  ggsave(filename = paste0(png_directory, 
                           sprintf("StackedPlot_bootstrap_%dpercent.png", prop * 100)), 
         plot = p, width = 8, height = 6, dpi = 300)
}



# First, reshape the bootstrap results into a longer format for plotting
bootstrap_long <- all_results[["0.5"]]$ci_results %>%
  mutate(
    mean_percent = mean_prop * 100,
    lower_percent = lower_ci * 100,
    upper_percent = upper_ci * 100,
    ci_width = upper_percent - lower_percent
  ) %>%
  arrange(Zone, Condition)

# Add percentage labels to the plot
p <- ggplot(bootstrap_long, aes(x = Zone, y = mean_percent, fill = Condition)) +
  geom_bar(stat = "identity", 
           position = position_dodge(width = 0.9),
           width = 0.8) +
  geom_errorbar(aes(ymin = lower_percent, ymax = upper_percent),
                position = position_dodge(width = 0.9),
                
                width = 0.25) +
  # Add percentage labels on bars
  geom_text(aes(label = sprintf("%.1f%%", mean_percent)),
            position = position_dodge(width = 0.9),
            vjust = -0.5,  # Adjust this value to position labels above bars
            size = 3.5) +
  scale_fill_manual(values = c("Normal" = "#224b5e", 
                               "AC" = "#94b594", 
                               "AH" = "#edc775")) +
  labs(title = "Zone Distribution Across Conditions",
       subtitle = "Means with 95% Confidence Intervals",
       x = "Hepatic Zone",
       y = "Percentage of Cells (%)") +
  theme_classic() +
  theme(legend.position = "top",
        plot.title = element_text(size = 14),
        plot.subtitle = element_text(size = 12),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 12))


# Save the plot
ggsave(filename = file.path(png_directory, "zone_distribution_bars.png"), 
       width = 10, height = 6, dpi = 300)
ggsave(filename = file.path(svg_directory, "zone_distribution_bars.svg"), 
       width = 8, height = 6, dpi = 300)


# Create tile plot showing confidence interval widths
p2 <- ggplot(combined_bootstrap, aes(x = Zone, y = Condition)) +
  geom_tile(aes(fill = ci_width)) +
  scale_fill_gradient(low = "white", high = "steelblue",
                      name = "CI Width (%)") +
  geom_text(aes(label = sprintf("%.1f", ci_width)), 
            color = "black", size = 3) +
  labs(title = "Precision of Estimates",
       subtitle = "Width of 95% Confidence Intervals") +
  theme_classic() +
  theme(legend.position = "right")

ggsave(filename = file.path(png_directory, "bootstrap_precision_heatmap.png"), 
       plot = p2, width = 8, height = 5, dpi = 300)

# Create summary table with formatted text
summary_table <- bootstrap_long %>%
  mutate(
    Estimate = sprintf("%.1f%% (%.1f%% - %.1f%%)", 
                       mean_percent, lower_percent, upper_percent),
    `Sample Size` = n,
    `CI Width` = sprintf("%.2f%%", ci_width)
  ) %>%
  select(Zone, Condition, Estimate, `Sample Size`, `CI Width`) %>%
  arrange(Zone, Condition)

# Write summary to file
write.csv(summary_table, 
          file.path(txt_directory, "bootstrap_summary_formatted.csv"),
          row.names = FALSE)

# Create detailed comparison text file
sink(file.path(txt_directory, "bootstrap_analysis_report.txt"))
cat("Bootstrap Analysis of Hepatic Zone Distribution\n")
cat("=============================================\n\n")

# Overall summary
cat("Sample Sizes and Zone Distribution\n")
cat("--------------------------------\n")
for(cond in unique(combined_bootstrap$Condition)) {
  zone_stats <- combined_bootstrap %>% 
    filter(Condition == cond)
  cat(sprintf("\n%s Condition:\n", cond))
  cat(sprintf("Total cells: %d\n", sum(zone_stats$n)))
  cat("Zone distribution:\n")
  for(i in 1:nrow(zone_stats)) {
    cat(sprintf("%s: %.1f%% (%.1f%% - %.1f%%)\n", 
                zone_stats$Zone[i],
                zone_stats$mean_percent[i],
                zone_stats$lower_percent[i],
                zone_stats$upper_percent[i]))
  }
}

# Key findings
cat("\nKey Findings\n")
cat("------------\n")

# Compare conditions for each zone
for(zone in unique(bootstrap_long$Zone)) {
  cat(sprintf("\n%s:\n", zone))
  zone_data <- bootstrap_long %>% 
    filter(Zone == zone) %>%
    arrange(desc(mean_percent))
  
  # Calculate max difference
  max_diff <- max(zone_data$mean_percent) - min(zone_data$mean_percent)
  cat(sprintf("- Range: %.1f percentage points\n", max_diff))
  
  # Report ordering
  cat("- Ordering: ")
  cat(paste(sprintf("%s (%.1f%%)", 
                    zone_data$Condition, 
                    zone_data$mean_percent), 
            collapse = " > "))
  cat("\n")
}

sink()


# Calculate gaps/overlaps between confidence intervals
check_ci_overlap <- function(data) {
  overlap_results <- data.frame()
  
  # For each zone
  for(zone in unique(data$Zone)) {
    zone_data <- data %>% filter(Zone == zone)
    
    # Normal vs AC
    normal <- zone_data %>% filter(Condition == "Normal")
    ac <- zone_data %>% filter(Condition == "AC")
    ah <- zone_data %>% filter(Condition == "AH")
    
    overlap_results <- rbind(overlap_results, data.frame(
      Zone = zone,
      Comparison = c("Normal-AC", "Normal-AH"),
      Gap = c(
        ac$lower_percent - normal$upper_percent,
        ah$lower_percent - normal$upper_percent,
        ah$lower_percent - ac$upper_percent
      )
    ))
  }
  return(overlap_results)
}


# Calculate and print the overlap analysis
overlap_analysis <- check_ci_overlap(bootstrap_long)
print("Confidence Interval Gap Analysis:")
print(overlap_analysis %>%
        mutate(
          Status = ifelse(Gap < 0,
                          sprintf("Overlaps by %.2f%%", abs(Gap)),
                          sprintf("Gap of %.2f%%", Gap))
        ))
