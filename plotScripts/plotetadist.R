#### Plotting the Distribution of Eta Across all Cells in Each Condition
#### Created by Aurelia Leona, 2024-07-16
#### Edited by Aurelia Leona, 2025-01-25
# Load required libraries
required_packages <- c("readr","dplyr","ggplot2")

# Install missing packages
new_packages <- required_packages[!(required_packages %in% installed.packages()[,"Package"])]
if(length(new_packages)) install.packages(new_packages)

# Load all packages
lapply(required_packages, library, character.only = TRUE)

# Load necessary libraries
library(ggplot2)
library(dplyr)
library(readr)

# Set random seed for reproducibility
set.seed(23)

# Initialize directories
getwd()
rawfileDir <- "/projects/b1042/GoyalLab/aleona/Zonation_Analysis/rawData/"
plotDir <- "/projects/b1042/GoyalLab/aleona/Zonation_Analysis/plots/"
extrdataDir <- "/projects/b1042/GoyalLab/aleona/Zonation_Analysis/extractedData/"
plotscriptsdir <- "/projects/b1042/GoyalLab/aleona/Zonation_Analysis/plotScripts/"
extScriptsdir <- "/projects/b1042/GoyalLab/aleona/Zonation_Analysis/extractScripts/"

source(paste0(extScriptsdir, "etazonefunction.R"))

# Import the combined eta zones data
etazones_directory <- paste0(extrdataDir, "etawithzones_data/")
etazones <- read_csv(paste0(etazones_directory, "combined_etawithzones.csv"))

#Define new Directories 
svg_directory <- paste0(plotDir, "svg/")
png_directory <- paste0(plotDir, "png/")

# Create directories
create_dir_if_not_exists(svg_directory)
create_dir_if_not_exists(png_directory)

# Define color palette for conditions
palette <- c("Normal" = "#224b5e", "AC" = "#94b594", "AH" = "#edc775")

# Factor the Condition with a specific order
etazones$Condition <- factor(etazones$Condition, levels = c("Normal", "AC", "AH"))

# Create histogram plot of eta distribution
p <- ggplot(etazones, aes(x = eta, fill = Condition)) +
  geom_histogram(binwidth = 0.015, color = "black") +
  labs(x = "Eta Values", 
       y = "Frequency", 
       title = "Distribution of Eta Values by Condition") +
  scale_fill_manual(values = palette) +  # Apply custom color palette
  theme_classic() +
  theme(legend.position = "right")

print(p)

# Save the plot
ggsave(filename = paste0(png_directory, "Distribution_of_Eta.png"), plot = p, width = 10, height = 6, dpi = 300)
ggsave(filename =  paste0(svg_directory, "Distribution_of_Eta.svg"), plot = p, width = 8, height = 6, dpi = 300)

# Create histogram plot of eta distribution
p2 <- ggplot(etazones, aes(x = eta, fill = Condition)) +
  geom_histogram(binwidth = 0.015, color = "black") +
  labs(x = "Eta Values", 
       y = "Frequency (log scale)", 
       title = "Distribution of Eta Values by Condition") +
  scale_fill_manual(values = palette) +  # Apply custom color palette
  scale_y_log10() +  # Add log scale transformation to y-axis
  theme_classic() +
  theme(legend.position = "right")


# Save the plot
ggsave(filename = paste0(png_directory, "Distribution_of_Eta_log.png"), plot = p2, width = 10, height = 6, dpi = 300)
ggsave(filename =  paste0(svg_directory, "Distribution_of_Eta_log.svg"), plot = p2, width = 8, height = 6, dpi = 300)

# Print confirmation message
print(paste("Plot saved as", paste0(png_directory, "Distribution_of_Eta.png")))


#### Statistical Significance of the distribution of eta across each condition-----------------
# Perform pairwise Kolmogorov-Smirnov tests
# First, separate data by condition
normal_eta <- etazones %>% 
  filter(Condition == "Normal") %>% 
  pull(eta)

ac_eta <- etazones %>% 
  filter(Condition == "AC") %>% 
  pull(eta)

ah_eta <- etazones %>% 
  filter(Condition == "AH") %>% 
  pull(eta)

# Perform K-S tests
ks_normal_ac <- ks.test(normal_eta, ac_eta)
ks_normal_ah <- ks.test(normal_eta, ah_eta)
ks_ac_ah <- ks.test(ac_eta, ah_eta)

# Create a data frame with results
ks_results <- data.frame(
  Comparison = c("Normal vs AC", "Normal vs AH", "AC vs AH"),
  D_statistic = c(ks_normal_ac$statistic, 
                  ks_normal_ah$statistic, 
                  ks_ac_ah$statistic),
  P_value = c(ks_normal_ac$p.value,
              ks_normal_ah$p.value,
              ks_ac_ah$p.value)
)

# Print results
print("Kolmogorov-Smirnov Test Results:")
print("--------------------------------")
print(ks_results)

#Save the resulst 
write.csv(ks_results, paste0(etazones_directory, "KSanalysis.csv"), row.names = FALSE)

# 
# # Add statistical annotation to the plot
# p_with_stats <- p +
#   annotate("text", x = max(etazones$eta), y = Inf,
#            label = paste("K-S test results:\n",
#                          "Normal vs AC: ", format_pvalue(ks_normal_ac$p.value), "\n",
#                          "Normal vs AH: ", format_pvalue(ks_normal_ah$p.value), "\n",
#                          "AC vs AH: ", format_pvalue(ks_ac_ah$p.value)),
#            hjust = 1, vjust = 1,
#            size = 3)
# 
# # Save the updated plot
# ggsave(filename = paste0(png_directory, "Distribution_of_Eta_with_KS.png"), 
#        plot = p_with_stats, width = 10, height = 6, dpi = 300)
# ggsave(filename = paste0(svg_directory, "Distribution_of_Eta_with_KS.svg"), 
#        plot = p_with_stats, width = 8, height = 6, dpi = 300)

## Compare different sample sizes on a quantile quantile plot ---------------------------------
create_qq_plot <- function(data1, data2, label1, label2, output_dir = bootstrapdir) {
  # Get the quantiles for both datasets
  n_points <- min(length(data1), length(data2))  # Use the smaller sample size
  probs <- seq(0, 1, length.out = n_points)
  
  # Calculate quantiles
  quantiles1 <- quantile(data1, probs, type = 1)
  quantiles2 <- quantile(data2, probs, type = 1)
  
  # Create data frame for plotting
  df <- data.frame(
    theoretical = quantiles1,
    sample = quantiles2,
    probability = probs
  )
  
  # Create descriptive names for the files
  comparison_name <- paste0(label1, "_vs_", label2)
  comparison_name <- gsub(" ", "_", comparison_name)  # Replace spaces with underscores
  
  # Save the quantile data
  write.csv(df, 
            file = file.path(etazones_directory, paste0("qq_data_", comparison_name, ".csv")),
            row.names = FALSE)
  
  # Save summary statistics
  summary_stats <- data.frame(
    Comparison = comparison_name,
    Label1 = label1,
    Label2 = label2,
    n1 = length(data1),
    n2 = length(data2),
    mean1 = mean(data1),
    mean2 = mean(data2),
    sd1 = sd(data1),
    sd2 = sd(data2)
  )
  
  write.csv(summary_stats,
            file = file.path(etazones_directory, paste0("qq_stats_", comparison_name, ".csv")),
            row.names = FALSE)
  
  # Create the plot
  ggplot(df, aes(x = theoretical, y = sample)) +
    geom_point(alpha = 0.6) +
    geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed") +
    labs(
      title = paste("Q-Q Plot:", label1, "vs", label2),
      x = paste("Quantiles of", label1),
      y = paste("Quantiles of", label2)
    ) +
    theme_classic() +
    # Add sample sizes to the plot
    annotate("text", 
             x = min(quantiles1), 
             y = max(quantiles2),
             hjust = 0, vjust = 1,
             label = sprintf("%s n=%d\n%s n=%d", 
                             label1, length(data1),
                             label2, length(data2)))
}

# Create and save Q-Q plots
# Normal vs AC
normal_ac_qq <- create_qq_plot(normal_eta, ac_eta, "Normal", "AC")
ggsave(filename = file.path(png_directory, "QQnormalac.png"), 
       plot = normal_ac_qq, 
       width = 10, height = 6, dpi = 300)
ggsave(filename = file.path(svg_directory, "QQnormalac.svg"), 
       plot = normal_ac_qq, 
       width = 8, height = 6, dpi = 300)

# Normal vs AH
normal_ah_qq <- create_qq_plot(normal_eta, ah_eta, "Normal", "AH")
ggsave(filename = file.path(png_directory, "QQnormalah.png"), 
       plot = normal_ah_qq, 
       width = 10, height = 6, dpi = 300)
ggsave(filename = file.path(svg_directory, "QQnormalah.svg"), 
       plot = normal_ah_qq, 
       width = 8, height = 6, dpi = 300)

# AC vs AH
ac_ah_qq <- create_qq_plot(ac_eta, ah_eta, "AC", "AH")
ggsave(filename = file.path(png_directory, "QQacah.png"), 
       plot = ac_ah_qq, 
       width = 10, height = 6, dpi = 300)
ggsave(filename = file.path(svg_directory, "QQacah.svg"), 
       plot = ac_ah_qq, 
       width = 8, height = 6, dpi = 300)

#Subsample Normal v Normal 

# Randomly split Normal data into two halves
n_samples <- floor(length(normal_eta)/2)  # This will be 4425 (half of 8851)
indices <- sample(1:length(normal_eta), n_samples)

normal_half1 <- normal_eta[indices]
normal_half2 <- normal_eta[-indices]

# Create QQ plot using existing function
normal_comparison <- create_qq_plot(normal_half1, normal_half2, "Normal Half 1", "Normal Half 2")

# Save the plot
ggsave(filename = paste0(png_directory, "QQ_normal_halves.png"), 
       plot = normal_comparison, width = 10, height = 6, dpi = 300)
ggsave(filename = paste0(svg_directory, "QQ_normal_halves.svg"), 
       plot = normal_comparison, width = 8, height = 6, dpi = 300)



##QQ metrics 
# Function to calculate Q-Q metrics
calculate_qq_metrics <- function(data1, data2, label1, label2) {
  # Get quantiles
  n_points <- min(length(data1), length(data2))
  probs <- seq(0, 1, length.out = n_points)
  q1 <- quantile(data1, probs = probs)
  q2 <- quantile(data2, probs = probs)
  
  # Calculate correlation
  cor_coef <- cor(q1, q2)
  
  # Calculate MSE
  mse <- mean((q1 - q2)^2)
  
  # Return results
  data.frame(
    Comparison = paste(label1, "vs", label2),
    Correlation = cor_coef,
    MSE = mse
  )
}

# Calculate metrics for all your comparisons
qq_metrics <- rbind(
  calculate_qq_metrics(normal_eta, ac_eta, "Normal", "AC"),
  calculate_qq_metrics(normal_eta, ah_eta, "Normal", "AH"),
  calculate_qq_metrics(normal_half1, normal_half2, "Normal1", "Normal2")
)

# Format results nicely
qq_metrics_formatted <- qq_metrics %>%
  mutate(
    Correlation = sprintf("%.4f", Correlation),
    MSE = sprintf("%.6f", MSE)
  )

# Write to CSV
write.csv(qq_metrics_formatted, 
          file = file.path(etazones_directory, "qq_plot_metrics.csv"), 
          row.names = FALSE)

# Print results
print(qq_metrics_formatted)
