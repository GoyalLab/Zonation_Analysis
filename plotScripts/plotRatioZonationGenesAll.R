#### Scripts to Plot Overall Ratio of Zonation Genes for Liver Zonation Analysis in Single-Cell RNA Sequencing Data 
#### Created by Aurelia Leona, 2024-07-16
#### Adapted from [original source if applicable]

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

# Import the combined ratios file
combined_x <- read_csv(paste0(xwithzones_directory, "combined_xwithzones.csv"))

# Import sorted data from CSV files
conditions <- c("Normal", "AC", "AH")
list_x <- map(conditions, function(condition) {
  read_csv(paste0(xwithzones_directory,"xwithzones_", condition,".csv"))
}) %>% set_names(conditions)


xzonesplot <- combined_x %>%
  group_by(Zone, Condition) %>%
  mutate(
    Mean_eucdist = mean(x),
    SD_eucdist = sd(x)
  ) %>%
  ungroup()

# Factor the Condition with a specific order
xzonesplot$Condition <- factor(xzonesplot$Condition, levels = c("Normal", "AC", "AH"))
# Factor the Condition with a specific order
xzonesplot$Zone <- factor(xzonesplot$Zone, levels = c("Zone 1", "Zone 2", "Zone 3"))

filtered_data <- filter(xzonesplot,Condition %in% c("Normal","AH"))

# Define color palette for conditions
palette <- c("Normal" = "#224b5e", "AH" = "#edc775", "Mean" = "black")

# Factor the Condition with a specific order
filtered_data$Condition <- factor(filtered_data$Condition, levels = c("Normal", "AH"))

# Create histogram plot of eta distribution
plot <- ggplot(filtered_data, aes(x = Condition, y = x)) +
  geom_jitter(aes(color = Condition), alpha = 0.7, size = 1.5, shape = 16, height = 0.01) +
  scale_color_manual(values = palette, 
                     breaks = c("Normal", "AH", "Mean"),
                     labels = c("Normal", "AH", "Mean")) + 
  theme_classic() +
  labs(y = "Ratio of Zonation Genes",
       x = "Condition",
       title = "Ratio of Zonation Genes of Cells") +
  theme(
    legend.position = "right",
    plot.title = element_text(size = 11),
    legend.title = element_blank(),
    strip.text = element_text(size = 10, face = "bold")
  )

# Save the plot
ggsave(filename = paste0(png_directory, "Ratio_Zonation_Genes_NormalAHAll.png"), 
       plot = plot, width = 8, height = 8, dpi = 300)
ggsave(filename = paste0(svg_directory, "Ratio_Zonation_Genes_NormalAHAll.svg"), 
       plot = plot, width = 8, height = 8, dpi = 300)

print("Processing and plotting complete. Plots saved in the plots directory.")

# Save data as CSV
write_csv(xzonesplot, paste0(xwithzones_directory, "xzoneswithsummary.csv"))
print(paste("xzonesplot saved as:", paste0(xwithzones_directory, "xzoneswithsummary.csv")))

write_csv(filtered_data, paste0(xwithzones_directory, "NormalAH_xzonesAll.csv"))
print(paste("filtered_data saved as:", paste0(xwithzones_directory, "NormalAH_xzones.csv")))