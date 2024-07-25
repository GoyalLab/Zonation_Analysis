#### Script for Statistical Tests on Ratio of Zonation Genes Between Zones
#### Created by Aurelia Leona, 2024-07-16
#### Adapted from [original source if applicable]

# Initialize directories
getwd()
rawfileDir <- "/projects/b1042/GoyalLab/aleona/github_uploads/Zonation_Analysis/rawData/"
plotDir <- "/projects/b1042/GoyalLab/aleona/github_uploads/Zonation_Analysis/plots/"
extrdataDir <- "/projects/b1042/GoyalLab/aleona/github_uploads/Zonation_Analysis/extractedData/"
plotscriptsdir <- "/projects/b1042/GoyalLab/aleona/github_uploads/Zonation_Analysis/plotScripts/"
extScriptsdir <- "/projects/b1042/GoyalLab/aleona/github_uploads/Zonation_Analysis/extractScripts/"

#Import custom functions
source(paste0(extScriptsdir, "statsfunc.R"))
stats_directory <- paste0(extrdataDir, "stats_tests/")
create_dir_if_not_exists(stats_directory)


# Import data
dfplot <- read_csv(paste0(xwithzones_directory, "xzoneswithsummary.csv"))
colnames(dfplot)[[2]] <- "value"
dfplot$Condition <- fct_relevel(dfplot$Condition, "Normal", "AH")
list_of_zones <- c("Zone 1", "Zone 2", "Zone 3")
list_of_conditions <- c("Normal", "AH")

# Perform test
wtest <- list()
ttest <- list()
for (z in list_of_zones) {
  # Perform Wilcoxon tests on conditions
  wtest[[paste0(z, "Normal-AH")]] <- wilcoxon_test_condition(dfplot, z, "Normal", "AH")
  # Perform T tests on conditions
  ttest[[paste0(z, "Normal-AH")]] <- t_test_condition(dfplot, z, "Normal", "AH")
}

for (c in list_of_conditions) {
  # Perform Wilcoxon tests on zones
  wtest[[paste0(c, "Zone2_Zone3")]] <- wilcoxon_test_zones(dfplot, c, "Zone 2", "Zone 3")
  # Perform T tests on zones
  ttest[[paste0(c, "Zone2_Zone3")]] <- t_test_zones(dfplot, c, "Zone 2", "Zone 3")
  # Perform T tests on zones
  ttest[[paste0(c, "Zone1_Zone3")]] <- t_test_zones(dfplot, c, "Zone 1", "Zone 3")
  # Perform T tests on zones
  ttest[[paste0(c, "Zone1_Zone2")]] <- t_test_zones(dfplot, c, "Zone 1", "Zone 2")
}


# Extract and save results
wilcoxon_results <- map_dfr(names(wtest), ~ tibble(
  Comparison = .x,
  P_value = wtest[[.x]]$p.value,
  Statistic = wtest[[.x]]$statistic,
))
t_test_results <- map_dfr(names(ttest), ~ tibble(
  Comparison = .x,
  P_value = ttest[[.x]]$p.value,
  Statistic = ttest[[.x]]$statistic,
  Mean_g1 = ttest[[.x]]$estimate[[1]],
  Mean_g2 = ttest[[.x]]$estimate[[2]]
))

write_csv(wilcoxon_results, paste0(stats_directory, "wtestonNormalAHZ123.csv"))
write_csv(t_test_results, paste0(stats_directory, "ttestonNormalAHZ123.csv"))

print("Statistical tests completed. Results saved in the extractedData directory.")
