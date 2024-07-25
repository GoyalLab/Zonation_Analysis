#### Script for Visualizing Euclidean Distances Between Liver Zones Across Conditions
#### Created by Aurelia Leona, 2024-07-16
#### Adapted from [original source if applicable]
# Load required libraries
required_packages <- c("tidyverse","dplyr","tidyr", "ggplot2")

# Install missing packages
new_packages <- required_packages[!(required_packages %in% installed.packages()[,"Package"])]
if(length(new_packages)) install.packages(new_packages)

# Load all packages
lapply(required_packages, library, character.only = TRUE)

library(tidyverse)
library(ggplot2)
library(dplyr)
library(tidyr)

# Initialize directories
getwd()
rawfileDir <- "/projects/b1042/GoyalLab/aleona/github_uploads/Zonation_Analysis/rawData/"
plotDir <- "/projects/b1042/GoyalLab/aleona/github_uploads/Zonation_Analysis/plots/"
extrdataDir <- "/projects/b1042/GoyalLab/aleona/github_uploads/Zonation_Analysis/extractedData/"
plotscriptsdir <- "/projects/b1042/GoyalLab/aleona/github_uploads/Zonation_Analysis/plotScripts/"
extScriptsdir <- "/projects/b1042/GoyalLab/aleona/github_uploads/Zonation_Analysis/extractScripts/"

source(paste0(extScriptsdir, "etazonefunction.R"))


#Define new Directories 
svg_directory <- paste0(plotDir, "svg/")
png_directory <- paste0(plotDir, "png/")

# Read the plot_zones.csv file
sorteddata_directory <- paste0(extrdataDir, "sorted_data/")
plot_zones <- read_csv(paste0(sorteddata_directory, "plot_zones.csv"))

# Calculate summary statistics
summary_df <- plot_zones %>%
  group_by(Condition) %>%
  summarize(
    Mean_eucdist = mean(eucdist),
    SD_eucdist = sd(eucdist)
  )

# Join summary statistics with the main dataframe
dfplot <- left_join(plot_zones, summary_df, by = "Condition")

# Ensure Condition is a factor with the correct order
dfplot$Condition <- factor(dfplot$Condition, levels = c("Normal", "AC", "AH"))

# Define color palette
palette <- c("#F4A582", "#67001F", "#4393C3","black")

# Create the plot
plot <- ggplot(dfplot, aes(x = Condition)) +
  geom_line(aes(y = eucdist, color = Zones, group = Zones), linewidth = 1) +
  geom_point(aes(y = eucdist, color = Zones), size = 5, shape = 15) +
  geom_point(aes(y = Mean_eucdist, color = "Mean"), size = 5, shape = 18) +
  geom_line(aes(y = Mean_eucdist, group = 1, color = "Mean"), linewidth = 1, linetype = "dashed") +
  scale_color_manual(values = palette, name = "Compared Zones", 
                     breaks = c("1_2", "1_3", "2_3", "Mean"),
                     labels = c("Zone 1 vs 2", "Zone 1 vs 3", "Zone 2 vs 3", "Mean")) + 
  theme_classic() +
  labs(
    y = "Euclidean Distance",
    x = "Condition",
    title = "Euclidean distances between Zones of the same Condition"
  ) +
  theme(
    legend.position = "right",
    plot.title = element_text(size = 11),
    legend.title = element_text(size = 12),
    strip.text = element_text(size = 10, face = "bold"),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    legend.text = element_text(size = 10)
  )

# Save the plot as PNG
ggsave(filename = paste0(png_directory, "PlotEucDistances_Zones.png"), 
       plot = plot, 
       width = 8, 
       height = 6, 
       dpi = 300)

# Save the plot as SVG
ggsave(filename = paste0(svg_directory, "PlotEucDistances_Zones.svg"), 
       plot = plot, 
       width = 8, 
       height = 6, 
       dpi = 300)


# Optionally, you can save the dfplot dataframe which includes the summary statistics
write_csv(dfplot, paste0(extrdataDir, "dfplot_with_summary.csv"))
