#### Script for Plot Pairwise Euclidean Distances between the Cells in Each Zones in Each Condition
#### Created by Aurelia Leona, 2024-07-16
#### Adapted from [original source if applicable]

# Load required libraries
required_packages <- c("readr","dplyr","ggplot2","purrr")

# Install missing packages
new_packages <- required_packages[!(required_packages %in% installed.packages()[,"Package"])]
if(length(new_packages)) install.packages(new_packages)

# Load all packages
lapply(required_packages, library, character.only = TRUE)

library(tidyverse)
library(readr)
library(ggplot2)

# Set random seed for reproducibility
set.seed(23)

# Initialize directories
getwd()
rawfileDir <- "/projects/b1042/GoyalLab/aleona/github_uploads/Zonation_Analysis/rawData/"
plotDir <- "/projects/b1042/GoyalLab/aleona/github_uploads/Zonation_Analysis/plots/"
extrdataDir <- "/projects/b1042/GoyalLab/aleona/github_uploads/Zonation_Analysis/extractedData/"
plotscriptsdir <- "/projects/b1042/GoyalLab/aleona/github_uploads/Zonation_Analysis/plotScripts/"
extScriptsdir <- "/projects/b1042/GoyalLab/aleona/github_uploads/Zonation_Analysis/extractScripts/"

compare_zones_dir <- paste0(extrdataDir, "Compare_Zones/")
svg_directory <- paste0(plotDir, "svg/")
png_directory <- paste0(plotDir, "png/")

# Read the combined distances CSV file
dplot <- read_csv(paste0(compare_zones_dir,"combined_distances.csv"))

# Calculate summary statistics and add them to the dataframe
dfplot <- dplot %>%
  group_by(Zone, Condition) %>%
  mutate(
    Mean_eucdist = mean(eucdist),
    SD_eucdist = sd(eucdist)
  ) %>%
  ungroup()

dfplot$Condition <- factor(dplot$Condition, levels = c("Normal", "AC", "AH")) 
dfplot$Zone <- factor(dplot$Zone , levels = c("Zone 1", "Zone 2", "Zone 3")) 


# Define color palette
palette <- c("Normal" = "#224b5e", "AC" = "#94b594", "AH" = "#edc775", "Mean" = "black")

# Create the plot
plot <- ggplot(dfplot, aes(x = Condition, y = eucdist)) +
  geom_jitter(aes(color = Condition), alpha = 0.7, size = 1.5, shape = 16, height = 0.01) +
  geom_pointrange(aes(y = Mean_eucdist, 
                      ymin = Mean_eucdist - SD_eucdist * 2, 
                      ymax = Mean_eucdist + SD_eucdist * 2, 
                      color = "Mean"),
                  size = 1, fatten = 3) +
  scale_color_manual(values = palette, 
                     breaks = c("Normal", "AC", "AH", "Mean"),
                     labels = c("Normal", "AC", "AH", "Mean")) + 
  facet_wrap(~ Zone, scales = "free_y") +
  theme_classic() +
  labs(y = "Euclidean Distance",
       x = "Condition",
       title = "Euclidean distances between Cells within the same Zones") +
  theme(
    legend.position = "right",
    plot.title = element_text(size = 11),
    legend.title = element_blank(),
    strip.text = element_text(size = 10, face = "bold")
  )

# Save the plot
ggsave(filename = paste0(png_directory, "Euclidean_Distances_by_Zone_and_Condition.png"), 
       plot = plot, width = 12, height = 8, dpi = 300)
ggsave(filename = paste0(svg_directory, "Euclidean_Distances_by_Zone_and_Condition.svg"), 
       plot = plot, width = 12, height = 8, dpi = 300)

print("Processing and plotting complete. Plots saved in the plots directory.")

# Save dplot as CSV
write_csv(dfplot, paste0(extrdataDir, "dfplot_euclidean_distances.csv"))
print(paste("dfplot saved as:", paste0(extrdataDir, "dplot_euclidean_distances.csv")))
