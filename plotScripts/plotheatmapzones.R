#### Script for Calculating and Visualizing Zonation Distances in Liver Single-Cell Data (HeatMap)
#### Created by Aurelia Leona, 2024-07-16
#### Adapted from [original source if applicable]

# Load required libraries
required_packages <- c("readr","pheatmap","ggplotify")

# Install missing packages
new_packages <- required_packages[!(required_packages %in% installed.packages()[,"Package"])]
if(length(new_packages)) install.packages(new_packages)

# Load all packages
lapply(required_packages, library, character.only = TRUE)

library(tidyverse)
library(pheatmap)
library(ggplotify)

# Initialize directories
getwd()
rawfileDir <- "/projects/b1042/GoyalLab/aleona/github_uploads/Zonation_Analysis/rawData/"
plotDir <- "/projects/b1042/GoyalLab/aleona/github_uploads/Zonation_Analysis/plots/"
extrdataDir <- "/projects/b1042/GoyalLab/aleona/github_uploads/Zonation_Analysis/extractedData/"
plotscriptsdir <- "/projects/b1042/GoyalLab/aleona/github_uploads/Zonation_Analysis/plotScripts/"
extScriptsdir <- "/projects/b1042/GoyalLab/aleona/github_uploads/Zonation_Analysis/extractScripts/"


svg_directory <- paste0(plotDir, "svg/")
png_directory <- paste0(plotDir, "png/")
sorteddata_directory <- paste0(extrdataDir, "sorted_data/")

# Import necessary functions
source(paste0(extScriptsdir, "etazonefunction.R"))  

# Import the combined Zonation genes available
genes_to_exclude <- read_csv(paste0(extrdataDir, "ZonationGenesavailable.csv"))

print(paste("Number of genes to exclude:", length(genes_to_exclude$gene)))

# Import sorted data from CSV files
conditions <- c("Normal", "AC", "AH")

list_sorted_data <- map(conditions, function(condition) {
  read_csv(paste0(sorteddata_directory, condition, "_sorted_data.csv"))
}) %>% set_names(conditions)

plot_zones <- data.frame(
  eucdist = numeric(0), 
  Condition = character(0),
  Zones = character(0),
  stringsAsFactors = FALSE
)

diverging_colors <- colorRampPalette(c("#053061", "#2166AC", "#4393C3", "#92C5DE", "#D1E5F0", "#F7F7F7", "#FDDBC7", "#F4A582", "#D6604D", "#B2182B", "#67001F"))(100)
list_zonedist <- list()
for (naming in names(list_sorted_data)){
  n <- 3
  high_variance_data <- list_sorted_data[[naming]]
  randall <- high_variance_data %>%
    filter(!gene %in% genes_to_exclude) %>%
    filter(Variance >= 0.1) %>%
    arrange(desc(Variance))
  
  rownames(randall) <- randall$gene
  randall$gene <- NULL
  
  distances <- matrix(NA, n, n)  # Initialize empty matrix for distances
  
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      print(paste(i, "and", j))
      dist_ij <- sqrt(sum((randall[,i] - randall[,j])^2))  # Calculate Euclidean distance
      distances[i, j] <- dist_ij
      distances[j, i] <- dist_ij  # Since it's symmetric, fill both sides
      plot_zones= rbind (plot_zones, data.frame(
        eucdist = dist_ij,
        Condition = naming, 
        Zones = paste0(i, "_", j)))
    }
  }
  
  diag (distances) <- sqrt(sum((randall[,i] - randall[,i])^2))
  
  result_df <- as.data.frame(distances)
  colnames(result_df) <- c( "Zone 1", "Zone 2", "Zone 3")
  rownames(result_df) <- c( "Zone 1", "Zone 2", "Zone 3")
  
  title_p = paste0("Distances Between Zones for ", naming, " Condition")
  hplot <- pheatmap(
    result_df,
    main = title_p,
    color = diverging_colors,
    breaks = seq(0, 20, length.out = 101),
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    scale = "none",  # Important: we don't want to scale the rows now
    show_rownames = TRUE,  # Set to TRUE if you want to show gene names
    show_colnames = TRUE,
    legend_breaks = seq(0, 30, by = 5),  # This sets the tick marks on the legend
    legend_labels = as.character(seq(0, 30, by = 5))  # This sets the labels for the tick marks
    
  )
  gg_heatmap <- ggplotify::as.ggplot(hplot$gtable)
  ggsave(gg_heatmap, file = paste0(svg_directory, 'HeatmapZonation_',naming,'.svg'), width = 8, height = 6, dpi = 300)
  ggsave(gg_heatmap, file = paste0(png_directory, 'HeatmapZonation_',naming,'.png'), width = 8, height = 6, dpi = 300)
  
  
  list_zonedist[[naming]] <- result_df
}

names(list_zonedist) <- names(list_sorted_data)

# Save plot_zones
write_csv(plot_zones, paste0(sorteddata_directory, "plot_zones.csv"))

print("Processing complete. Results and heatmaps saved.")



