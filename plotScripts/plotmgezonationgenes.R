#### Scripts to Plot Mean Gene Expression of Zonation Genes for Liver Zonation Analysis in Single-Cell RNA Sequencing Data 
#### Created by Aurelia Leona, 2024-07-16
#### Adapted from [original source if applicable]

required_packages <- c("tidyverse","dplyr","tidyr", "ggplot2", "ggrepel")

# Install missing packages
new_packages <- required_packages[!(required_packages %in% installed.packages()[,"Package"])]
if(length(new_packages)) install.packages(new_packages)

# Load all packages
lapply(required_packages, library, character.only = TRUE)

library(tidyverse)
library(ggplot2)
library(ggrepel)
library(dplyr)
library(tidyr)

# Set random seed for reproducibility
set.seed(23)

# Initialize directories
getwd()
rawfileDir <- "/projects/b1042/GoyalLab/aleona/github_uploads/Zonation_Analysis/rawData/"
plotDir <- "/projects/b1042/GoyalLab/aleona/github_uploads/Zonation_Analysis/plots/"
extrdataDir <- "/projects/b1042/GoyalLab/aleona/github_uploads/Zonation_Analysis/extractedData/"
plotscriptsdir <- "/projects/b1042/GoyalLab/aleona/github_uploads/Zonation_Analysis/plotScripts/"
extScriptsdir <- "/projects/b1042/GoyalLab/aleona/github_uploads/Zonation_Analysis/extractScripts/"

# Import custom functions for eta calculation
source(paste0(extScriptsdir, "etazonefunction.R"))
source(paste0(extScriptsdir, "gvefunction.R"))

# Define directories
x_directory <- paste0(extrdataDir, "x_data/")
xwithzones_directory <- paste0(extrdataDir, "xwithzones_data/")
seurat_obj_directory <- paste0(extrdataDir, "seurat_objects/")
mge_directory <- paste0(extrdataDir, "mge/")
svg_directory <- paste0(plotDir, "svg/")
png_directory <- paste0(plotDir, "png/")

# Import Seurat objects (assuming they're stored as RDS files)
list_subset <- list(
  Normal = readRDS(paste0(seurat_obj_directory, "hepnormalwithpc.rds")),
  AH = readRDS(paste0(seurat_obj_directory, "hepahwithpc.rds"))
)

# Import the combined ratios file
combined_x <- read_csv(paste0(xwithzones_directory, "combined_xwithzones.csv"))

#Import the central and portal genes 
genes_to_include <- read_csv(paste0(extrdataDir, "ZonationGenesavailable.csv"))

mge_list <- list()

for (condition in names(list_subset)) {
  # Extract etazones for the current condition
  x <- combined_x %>% filter(Condition == condition)
  
  # Compute Pmat and get subsampled cells
  result <- compute_pmat_mge(x)
  Pmat <- result$Pmat
  subsampled_cells <- result$subsampled_cells
  
  # Compute MGE
  heps <- list_subset[[condition]]
  mat_norm <- as.matrix(heps@assays[["RNA"]]@data)
  # Convert the gene column to a vector
  genes_to_keep <- genes_to_include$gene
  
  # Subset matnorm to keep only the rows (genes) that are in genes_to_keep
  mat_norm_subsampled <- mat_norm[rownames(mat_norm) %in% genes_to_keep, ]
  mat_norm_subsampled <- t(mat_norm_subsampled)
  
  mat_normdf <- as.data.frame(mat_norm_subsampled) %>%
    rownames_to_column(var = "cell_id") %>%
    left_join(combined_x, by = "cell_id")
  
  # Calculate means of non-zero values for each column of filtered Z1 
  filterZ1 <- filter(mat_normdf, Zone %in% "Zone 1") 
  nmeandfZ1 <- filterZ1 %>%
    sapply( function(x) mean(x[x != 0], na.rm = TRUE))%>%
    t() %>%
    as.data.frame() %>%
    mutate(
      cell_id = "Mean",
      Zone = "Zone 1", 
      Condition = condition
    )
  filterZ1 <- rbind(filterZ1, nmeandfZ1)
  
  # Calculate means of non-zero values for each column of filtered Z2 
  filterZ2 <- filter(mat_normdf, Zone %in% "Zone 2") 
  nmeandfZ2 <- filterZ2 %>%
    sapply( function(x) mean(x[x != 0], na.rm = TRUE))%>%
    t() %>%
    as.data.frame() %>%
    mutate(
      cell_id = "Mean",
      Zone = "Zone 2", 
      Condition = condition
    )
  filterZ2 <- rbind(filterZ2, nmeandfZ2)
  
  # Calculate means of non-zero values for each column of filtered Z3 
  filterZ3 <- filter(mat_normdf, Zone %in% "Zone 3")
  nmeandfZ3 <- filterZ3 %>%
    sapply( function(x) mean(x[x != 0], na.rm = TRUE))%>%
    t() %>%
    as.data.frame() %>%
    mutate(
      cell_id = "Mean",
      Zone = "Zone 3", 
      Condition = condition
    )
  filterZ3 <- rbind(filterZ3, nmeandfZ3)
  
  nmeandfZ3 <- subset(nmeandfZ3, select = -c(Condition, cell_id))
  nmeandfZ2 <- subset(nmeandfZ2, select = -c(Condition, cell_id))
  nmeandfZ1 <- subset(nmeandfZ1, select = -c(Condition, cell_id))
  mean_data <- rbind(nmeandfZ1,nmeandfZ2, nmeandfZ3)
  mean_datadf <- mean_data%>%
    column_to_rownames("Zone")%>%
    subset(select = -c(x))%>%
    t()%>%
    as.data.frame()%>%
    rownames_to_column( var = "gene") %>%
    left_join(genes_to_include, by = "gene") %>%
    mutate(
      gene_type = factor(gene_type, levels = c("central", "portal")),
      Condition = condition
    )
  
  # Save results as CSV
  write_csv(mean_datadf, paste0(mge_directory, condition, "_MGEZ123.csv"))
  
  print(paste("Processed and saved data for", condition))

  mge_list[[condition]]  <- mean_datadf
}


# # Define your custom color palette and breaks
# diverging_colors <- colorRampPalette(c("#053061", "#2166AC", "#4393C3", "#92C5DE", "#D1E5F0", "#F7F7F7", "#FDDBC7", "#F4A582", "#D6604D", "#B2182B", "#67001F"))(100)
# custom_breaks <- c(0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0)  # Example breaks, adjust as needed
# 
# # Create heatmaps for each condition
# heatmap_plots <- imap(mge_list, function(data, condition) {
#   # Reshape the data from wide to long format
#   data_long <- data %>%
#     pivot_longer(cols = c("Zone 1","Zone 2", "Zone 3"), names_to = "zone", values_to = "expression")
#   
#   # Create the heatmap
#   plot <- ggplot(data_long, aes(x = zone, y = gene, fill = expression)) +
#     geom_tile() +
#     scale_fill_gradientn(colors = diverging_colors,
#                          limits = c(min(custom_breaks), max(custom_breaks)),
#                          breaks = custom_breaks,
#                          labels = custom_breaks) +
#     facet_grid(gene_type ~ ., scales = "free_y", space = "free_y") +
#     theme_minimal() +
#     theme(axis.text.y = element_text(size = 6),
#           axis.text.x = element_text(angle = 45, hjust = 1)) +
#     labs(title = paste("Mean Gene Expression of Zonation Genes -", condition),
#          x = "Zones",
#          y = "Genes",
#          fill = "Expression")
#   # Save as PNG
#   ggsave(filename = paste0(png_directory, "MGEZonationGenesHeatmapZ123_", condition, ".png"),
#          plot = plot,
#          width = 8,
#          height = 6,
#          dpi = 300)
#   
#   # Save as SVG
#   ggsave(filename = paste0(svg_directory, "MGEZonationGenesHeatmapZ123_", condition, ".svg"),
#          plot = plot,
#          width = 8,
#          height = 6,
#          dpi = 300)
#   
#   # Return the plot object
#   return(plot)
# })


# Define colors for central, portal, and mean
palette <- c("Normal" = "#224b5e", "AH" = "#edc775", "Mean" = "black")

# Combine the list of dataframes into one
combined_df <- bind_rows(mge_list)

# Reshape the data from wide to long format
long_df <- combined_df %>%
  pivot_longer(cols = c("Zone 1", "Zone 2", "Zone 3"), 
               names_to = "Zone", 
               values_to = "Expression") %>%
  mutate(Condition = factor(Condition, levels = c("Normal", "AH")))

# Determine the overall range of Expression values
y_min <- min(long_df$Expression, na.rm = TRUE)
y_max <- max(long_df$Expression, na.rm = TRUE)

# Calculate means and standard deviations
mean_data <- long_df %>%
  group_by(Zone, Condition, gene_type) %>%
  summarise(
    Mean_Expression = mean(Expression, na.rm = TRUE),
    SD_Expression = sd(Expression, na.rm = TRUE),
    .groups = 'drop'
  )

mean_data$Condition <- factor(mean_data$Condition, levels = c("Normal", "AH")) 
mean_data$gene_type <- factor(mean_data$gene_type, levels = c("central", "portal")) 
mean_data$Zone <- factor(mean_data$Zone, levels = c("Zone 1", "Zone 2","Zone 3")) 

line_data <- long_df %>%
  filter(!is.na(Expression)) %>%
  group_by(Zone, gene, gene_type) %>%
  filter(n() == 2 & n_distinct(Condition) == 2)  # Keep only genes that appear in both conditions


#Separate the dataset into only central genes
central_long <- filter(long_df, gene_type == "central")
central_mean <- filter(mean_data, gene_type == "central")
central_gene_name <- filter(line_data, gene_type == "central")

plot3 <- ggplot() +
  # Individual data points
  geom_point(data = central_long, aes(x = Condition, y = Expression, color = Condition), position = position_dodge(width = 1.0),
             size = 2.5, alpha = 0.7) +
  # Mean points with error bars (black)
  geom_pointrange(data = central_mean,
                  aes(x = Condition, y = Mean_Expression, 
                      ymin = Mean_Expression - SD_Expression, 
                      ymax = Mean_Expression + SD_Expression,
                      group = Condition),
                  position = position_dodge(width = 0.8),
                  color = "black",
                  size = 1, fatten = 2) +
  #Line between the same Gene
  geom_line(data = central_long, aes(x=Condition, y= Expression, group= interaction(gene)),
            color = "antiquewhite3", size = 0.3) +
  
  #Line between the same Mean
  geom_line(data = central_mean, aes(x=Condition, y= Mean_Expression, group = Zone),
            color = "black", size = 0.9, linetype = "dashed" ) +
  
  # Gene labels
  geom_text_repel(data = central_gene_name %>% filter(Condition == "AH"),
                  aes(x = Condition, y = Expression, label = gene),
                  position = position_dodge(width = 0.8),
                  size = 2.5, 
                  box.padding = unit(0.2, "lines"),
                  point.padding = unit(0.2, "lines"),
                  max.overlaps = 20) +
  
  scale_color_manual(values = palette,
                     name = "Condition") +
  facet_wrap(~Zone) +
  scale_fill_manual(values = palette,
                    name = "Condition") +
  theme_classic() +
  labs(title = "Gene Expression of Central Genes for Normal and AH",
       x = "Condition",
       y = "Gene Expression") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "right")

print(plot3)

# Save as PNG
ggsave(filename = paste0(png_directory, "CentralMGENormalAH.png"),
       plot = plot3,
       width = 8,
       height = 6,
       dpi = 300)

# Save as SVG
ggsave(filename = paste0(svg_directory, "CentralMGENormalAH.svg"),
       plot = plot3,
       width = 8,
       height = 6,
       dpi = 300)


#Separate the dataset into only portal genes
portal_long <- filter(long_df, gene_type == "portal")
portal_mean <- filter(mean_data, gene_type == "portal")
portal_gene_name <- filter(line_data, gene_type == "portal")

plot4 <- ggplot() +
  # Individual data points
  geom_point(data = portal_long, aes(x = Condition, y = Expression, color = Condition), position = position_dodge(width = 1.0),
             size = 2.5, alpha = 0.7) +
  # Mean points with error bars (black)
  geom_pointrange(data = portal_mean,
                  aes(x = Condition, y = Mean_Expression, 
                      ymin = Mean_Expression - SD_Expression, 
                      ymax = Mean_Expression + SD_Expression,
                      group = Condition),
                  position = position_dodge(width = 0.8),
                  color = "black",
                  size = 1, fatten = 2) +
  #Line between the same Gene
  geom_line(data = portal_long, aes(x=Condition, y= Expression, group= interaction(gene)),
            color = "antiquewhite3", size = 0.3) +
  
  #Line between Means
  geom_line(data = portal_mean, aes(x=Condition, y= Mean_Expression, group = Zone),
            color = "black", size = 0.9, linetype = "dashed" ) +
  # Gene labels
  geom_text_repel(data = portal_gene_name %>% filter(Condition == "AH"),
                  aes(x = Condition, y = Expression, label = gene),
                  position = position_dodge(width = 0.8),
                  size = 2.5, 
                  box.padding = unit(0.2, "lines"),
                  point.padding = unit(0.2, "lines"),
                  max.overlaps = 20) +
  
  scale_color_manual(values = palette,
                     name = "Condition") +
  facet_wrap(~Zone) +
  scale_fill_manual(values = palette,
                    name = "Condition") +
  theme_classic() +
  labs(title = "Gene Expression of Portal Genes for Normal and AH",
       x = "Condition",
       y = "Gene Expression") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "right")

print(plot4)




# Save as PNG
ggsave(filename = paste0(png_directory, "PortalMGENormalAH.png"),
       plot = plot4,
       width = 8,
       height = 6,
       dpi = 300)

# Save as SVG
ggsave(filename = paste0(svg_directory, "PortalMGENormalAH.svg"),
       plot = plot4,
       width = 8,
       height = 6,
       dpi = 300)




