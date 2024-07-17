#### Script for Statistical Tests on Euclidean Distances Between Zones
#### Created by Aurelia Leona, 2024-07-16
#### Adapted from [original source if applicable]

# Initialize directories
getwd()
rawfileDir <- "/projects/b1042/GoyalLab/aleona/Zonation_Analysis/rawData/"
plotDir <- "/projects/b1042/GoyalLab/aleona/Zonation_Analysis/plots/"
extrdataDir <- "/projects/b1042/GoyalLab/aleona/Zonation_Analysis/extractedData/"
plotscriptsdir <- "/projects/b1042/GoyalLab/aleona/Zonation_Analysis/plotScripts/"
extScriptsdir <- "/projects/b1042/GoyalLab/aleona/Zonation_Analysis/extractScripts/"

#Import custom functions
source(extScriptsdir, "statsfunc.R")

# Import data
dfplot <- read_csv(paste0(extrdataDir, "dfplot_euclidean_distances.csv"))
dfplot$Condition <- fct_relevel(dfplot$Condition, "Normal", "AC", "AH")
list_of_zones <- c("Zone 1", "Zone 2", "Zone 3")

# Perform Wilcoxon tests
wtest <- list()
for (z in list_of_zones) {
  wtest[[paste0(z, "Normal-AC")]] <- perform_wilcoxon_test(dfplot, z, "Normal", "AC")
  wtest[[paste0(z, "Normal-AH")]] <- perform_wilcoxon_test(dfplot, z, "Normal", "AH")
  wtest[[paste0(z, "AC-AH")]] <- perform_wilcoxon_test(dfplot, z, "AC", "AH")
}

# Perform t-tests
ttest <- list()
for (z in list_of_zones) {
  ttest[[paste0(z, "Normal-AC")]] <- perform_t_test(dfplot, z, "Normal", "AC")
  ttest[[paste0(z, "Normal-AH")]] <- perform_t_test(dfplot, z, "Normal", "AH")
  ttest[[paste0(z, "AC-AH")]] <- perform_t_test(dfplot, z, "AC", "AH")
}

# Function to extract test results
extract_test_results <- function(test_list) {
  map_dfr(names(test_list), ~ tibble(
    Comparison = .x,
    P_value = test_list[[.x]]$p.value,
    Statistic = test_list[[.x]]$statistic,
    Mean_g1 = test_list[[.x]]$estimate[[1]],
    Mean_g2 = test_list[[.x]]$estimate[[2]]
  ))
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

write_csv(wilcoxon_results, paste0(extrdataDir, "wilcoxon_test_results.csv"))
write_csv(t_test_results, paste0(extrdataDir, "t_test_results.csv"))

print("Statistical tests completed. Results saved in the extractedData directory.")
