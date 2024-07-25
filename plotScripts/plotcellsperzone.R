#### Plotting the Stack Plot of Normalized Number of Cells in Each Condition
#### Created by Aurelia Leona, 2024-07-16
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

# Initialize directories
getwd()
rawfileDir <- "/projects/b1042/GoyalLab/aleona/github_uploads/Zonation_Analysis/rawData/"
plotDir <- "/projects/b1042/GoyalLab/aleona/github_uploads/Zonation_Analysis/plots/"
extrdataDir <- "/projects/b1042/GoyalLab/aleona/github_uploads/Zonation_Analysis/extractedData/"
plotscriptsdir <- "/projects/b1042/GoyalLab/aleona/github_uploads/Zonation_Analysis/plotScripts/"
extScriptsdir <- "/projects/b1042/GoyalLab/aleona/github_uploads/Zonation_Analysis/extractScripts/"

source(paste0(extScriptsdir, "etazonefunction.R"))

# Import the combined eta zones data
etawithzones_directory <- paste0(extrdataDir, "etawithzones_data/")
etazones <- read_csv(paste0(etawithzones_directory, "combined_etawithzones.csv"))

#Define new Directories 
svg_directory <- paste0(plotDir, "svg/")
png_directory <- paste0(plotDir, "png/")

# Create directories
create_dir_if_not_exists(svg_directory)
create_dir_if_not_exists(png_directory)

# Define conditions and zones
list_name <- c("Normal", "AC", "AH")
list_of_zones <- c("Zone 1", "Zone 2", "Zone 3")

# Calculate the number of cells per zone per condition
etazonesplot <- etazones %>%
  group_by(Condition, Zone) %>%
  summarise(Count = n(), .groups = 'drop') %>%
  group_by(Condition) %>%
  mutate(CountNorm = Count / sum(Count)) %>%
  ungroup() %>%
  mutate(
    Condition = factor(Condition, levels = list_name),
    Zone = factor(Zone, levels = list_of_zones)
  )

# Save etazonesplot as a CSV
write_csv(etazonesplot, paste0(extrdataDir, "cells_perZC.csv"))

# Define color palette
palette <- c("Zone 1" ="#f5c34d", "Zone 2" = "#e67424", "Zone 3" = "#8d1c06")

# Create stacked bar plot
p <- ggplot(etazonesplot, aes(x = Condition, y = CountNorm, fill = Zone)) +
  geom_col(position = "stack") +
  # geom_text(aes(label = sprintf("%.1f%%", CountNorm*100)), 
  #           position = position_stack(vjust = 0.5), 
  #           color = "black", size = 4) +
  labs(x = "Conditions", 
       y = "Percentage of Cells", 
       title = "Percentage of Cells per Zone Across Conditions") +
  scale_fill_manual(values = palette) +
  scale_y_continuous(labels = scales::percent) +
  theme_classic() +
  theme(legend.position = "right")

# Save the plot
ggsave(filename =  paste0(png_directory, "StackedPlotCells3Zones.png"), plot = p, width = 8, height = 6, dpi = 300)
ggsave(filename =  paste0(svg_directory, "StackedPlotCells3Zones.svg"), plot = p, width = 8, height = 6, dpi = 300)

## Plot with texts for visualization 
p2 <- p + geom_text(aes(label = sprintf("%.1f%%", CountNorm*100)),
            position = position_stack(vjust = 0.5),
            color = "black", size = 4)

ggsave(filename =  paste0(png_directory, "StackedPlotCells3Zoneslabels.png"), plot = p2, width = 8, height = 6, dpi = 300)

# Print confirmation messages
print(paste("Etazonesplot data saved as", paste0(extrdataDir, "cells_perZC.csv")))
print(paste("Plot saved as", paste0(png_directory, "StackedPlotCells3Zones.png"))
