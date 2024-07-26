#### Script for Plotting Zonation Genes in Liver Single-Cell Data (HeatMap)
#### Created by Aurelia Leona, 2024-07-16
#### Adapted from [original source if applicable]

# Load required libraries
required_packages <- c("readr","pheatmap","ggplotify","purrr","dplyr","ggplot2")

# Install missing packages
new_packages <- required_packages[!(required_packages %in% installed.packages()[,"Package"])]
if(length(new_packages)) install.packages(new_packages)

# Load all packages
lapply(required_packages, library, character.only = TRUE)

library(tidyverse)
library(pheatmap)
library(ggplotify)
library(purrr)
library(dplyr)
library(ggplot2)


# Initialize directories
getwd()
rawfileDir <- "/projects/b1042/GoyalLab/aleona/github_uploads/Zonation_Analysis/rawData/"
plotDir <- "/projects/b1042/GoyalLab/aleona/github_uploads/Zonation_Analysis/plots/"
extrdataDir <- "/projects/b1042/GoyalLab/aleona/github_uploads/Zonation_Analysis/extractedData/"
plotscriptsdir <- "/projects/b1042/GoyalLab/aleona/github_uploads/Zonation_Analysis/plotScripts/"
extScriptsdir <- "/projects/b1042/GoyalLab/aleona/github_uploads/Zonation_Analysis/extractScripts/"

sorteddata_directory <- paste0(extrdataDir, "sorted_data/")
svg_directory <- paste0(plotDir, "svg/")
png_directory <- paste0(plotDir, "png/")

# Import necessary functions
source(paste0(extScriptsdir, "etazonefunction.R"))  

# Import the combined gene availability file
genes_to_include <- read_csv(paste0(extrdataDir, "ZonationGenesavailable.csv"))


# Import sorted data from CSV files
conditions <- c("Normal", "AC", "AH")
sorteddata_directory <- paste0(extrdataDir, "sorted_data/")
list_sorted_data <- map(conditions, function(condition) {
  read_csv(paste0(sorteddata_directory, condition, "_sorted_data.csv"))
}) %>% set_names(conditions)

# Create a new list of filtered sorted data with all the new information added
filtered_sorted_data <- imap(list_sorted_data, function(data, condition_name) {
  data %>%
    filter(gene %in% new_dataset$gene) %>%
    left_join(genes_to_include, by = "gene") %>%
    mutate(
      Gene_type = factor(gene_type, levels = c("central", "portal")),
      Condition = condition_name
    )
})



# Define your custom color palette and breaks
diverging_colors <- colorRampPalette(c("#053061", "#2166AC", "#4393C3", "#92C5DE", "#D1E5F0", "#F7F7F7", "#FDDBC7", "#F4A582", "#D6604D", "#B2182B", "#67001F"))(100)
custom_breaks <- c(0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0)  # Example breaks, adjust as needed

# Create heatmaps for each condition
heatmap_plots <- imap(filtered_sorted_data, function(data, condition) {
  # Reshape the data from wide to long format
  data_long <- data %>%
    pivot_longer(cols = c("Zone 1", "Zone 2", "Zone 3"), names_to = "zone", values_to = "expression")
  
  # Create the heatmap
  ggplot(data_long, aes(x = zone, y = gene, fill = expression)) +
    geom_tile() +
    scale_fill_gradientn(colors = diverging_colors,
                         limits = c(min(custom_breaks), max(custom_breaks)),
                         breaks = custom_breaks,
                         labels = custom_breaks) +
    facet_grid(gene_type ~ ., scales = "free_y", space = "free_y") +
    theme_minimal() +
    theme(axis.text.y = element_text(size = 6),
          axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = paste("Gene Expression Heatmap -", condition),
         x = "Zones",
         y = "Genes",
         fill = "Expression")
})

# View the heatmaps
print(heatmap_plots$Normal)
print(heatmap_plots$AC)
print(heatmap_plots$AH)

# Save the heatmaps
for (condition in names(heatmap_plots)) {
  gg_heatmap <- heatmap_plots[[condition]]
  naming <- condition
  
  ggsave(gg_heatmap, file = paste0(svg_directory, 'HeatmapZonationGenes_', naming, '.svg'), width = 8, height = 6, dpi = 300)
  ggsave(gg_heatmap, file = paste0(png_directory, 'HeatmapZonationGenes_', naming, '.png'), width = 8, height = 6, dpi = 300)
}


# Combine all data into a single dataframe
all_data <- bind_rows(filtered_sorted_data, .id = "condition")

# Save the combined data to a single CSV file
write_csv(all_data, file = paste0(sorteddata_directory, "all_filtered_sorted_data.csv"))
cat("Saved: all_filtered_sorted_data.csv\n")
