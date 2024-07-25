#### Script for Statistical Tests on Euclidean Distances Between Zones
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
dfplot <- read_csv(paste0(extrdataDir, "dfplot_euclidean_distances.csv"))
colnames(dfplot)[[3]] <- "value"
dfplot$Condition <- fct_relevel(dfplot$Condition, "Normal", "AC", "AH")
list_of_zones <- c("Zone 1", "Zone 2", "Zone 3")

# Perform Wilcoxon tests
wtest <- list()
for (z in list_of_zones) {
  wtest[[paste0(z, "Normal-AC")]] <- wilcoxon_test_condition(dfplot, z, "Normal", "AC")
  wtest[[paste0(z, "Normal-AH")]] <- wilcoxon_test_condition(dfplot, z, "Normal", "AH")
  wtest[[paste0(z, "AC-AH")]] <- wilcoxon_test_condition(dfplot, z, "AC", "AH")
}

# Perform t-tests
ttest <- list()
for (z in list_of_zones) {
  ttest[[paste0(z, "Normal-AC")]] <- t_test_condition(dfplot, z, "Normal", "AC")
  ttest[[paste0(z, "Normal-AH")]] <- t_test_condition(dfplot, z, "Normal", "AH")
  ttest[[paste0(z, "AC-AH")]] <- t_test_condition(dfplot, z, "AC", "AH")
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

write_csv(wilcoxon_results, paste0(stats_directory, "wilcoxon_test_eucdist.csv"))
write_csv(t_test_results, paste0(stats_directory, "t_test_eucdist.csv"))

print("Statistical tests completed. Results saved in the extractedData directory.")
