#### Plotting the Distribution of Eta Across all Cells in Each Condition
#### Created by Aurelia Leona, 2024-07-16
#### Adapted from [original source if applicable]

# Load necessary libraries
library(ggplot2)
library(dplyr)
library(readr)

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

# Save the plot
ggsave(filename = paste0(png_directory, "Distribution_of_Eta.png"), plot = p, width = 10, height = 6, dpi = 300)
ggsave(filename =  paste0(svg_directory, "Distribution_of_Eta.svg"), plot = p, width = 8, height = 6, dpi = 300)

# Print confirmation message
print(paste("Plot saved as", paste0(png_directory, "Distribution_of_Eta.png")))
