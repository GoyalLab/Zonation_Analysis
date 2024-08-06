#### Script for Statistical Tests on Gene Expression Between Zones
#### Created by Aurelia Leona, 2024-07-16
#### Adapted from [original source if applicable]

set.seed(23)
# Initialize directories
getwd()
rawfileDir <- "/projects/b1042/GoyalLab/aleona/github_uploads/Zonation_Analysis/rawData/"
plotDir <- "/projects/b1042/GoyalLab/aleona/github_uploads/Zonation_Analysis/plots/"
extrdataDir <- "/projects/b1042/GoyalLab/aleona/github_uploads/Zonation_Analysis/extractedData/"
plotscriptsdir <- "/projects/b1042/GoyalLab/aleona/github_uploads/Zonation_Analysis/plotScripts/"
extScriptsdir <- "/projects/b1042/GoyalLab/aleona/github_uploads/Zonation_Analysis/extractScripts/"

#Import custom functions
source(paste0(extScriptsdir, "statsfunc.R"))
source(paste0(extScriptsdir, "etazonefunction.R"))

#Define and Create Directories
mge_directory <- paste0(extrdataDir, "mge/")
stats_directory <- paste0(extrdataDir, "stats_tests/")
create_dir_if_not_exists(stats_directory)


# Import data
dfnormal <- read_csv(paste0(mge_directory,"Normal_MGEZ123.csv"))
dfAH <- read_csv(paste0(mge_directory,"AH_MGEZ123.csv"))
df<- rbind(dfnormal,dfAH)
list_of_zones <- c("Zone 1", "Zone 2","Zone 3")
list_of_conditions <- c("Normal", "AH")
df_reshaped <- df %>%
  pivot_longer(
    cols = list_of_zones,
    names_to = "Zone",
    values_to = "value"
  )

#Including all the zeros 
df_reshaped$Condition <- fct_relevel(df_reshaped$Condition, "Normal", "AH")

# Perform test
ttestp <- list()
ttest <- list()
for (z in list_of_zones) {
  # Perform T tests on conditions
  ttest[[paste0(z, "Normal-AH")]] <- t_test_condition(df_reshaped, z, "Normal", "AH")
}

for (c in list_of_conditions) {
  # Perform two-sided T tests on zones
  ttest[[paste0(c, "Zone1_Zone3")]] <- t_test_zones(df_reshaped, c, "Zone 1", "Zone 3")
  ttest[[paste0(c, "Zone1_Zone2")]] <- t_test_zones(df_reshaped, c, "Zone 1", "Zone 2")
  ttest[[paste0(c, "Zone2_Zone3")]] <- t_test_zones(df_reshaped, c, "Zone 2", "Zone 3")
}


for (z in list_of_zones){
  subset_data <- subset(df_reshaped, gene_type == "central" & Zone == z & Condition %in% c("Normal", "AH"))
  x <- subset(subset_data, Condition =="Normal") %>% select(Expression)
  y <- subset(subset_data, Condition =="AH") %>% select(Expression)
  x = x$Expression
  y = y$Expression
  ttest[[paste0("centralNormal_AH",z)]] <- t.test(x,y, alternative = "two.sided", 
                                          mu = 0, paired = FALSE, var.equal = FALSE, conf.level = 0.95)
  ttestp[[paste0("centralNormal_AH",z)]] <- t.test(value ~ Condition, data = subset_data, alternative = "two.sided", 
                                                  mu = 0, paired = TRUE, var.equal = FALSE, conf.level = 0.95)
  
  subset_data <- subset(df_reshaped, gene_type == "portal" & Zone == z & Condition %in% c("Normal", "AH"))
  x <- subset(subset_data, Condition =="Normal") %>% select(Expression)
  y <- subset(subset_data, Condition =="AH") %>% select(Expression)
  x = x$Expression
  y = y$Expression
  ttest[[paste0("portalNormal_AH",z)]] <- t.test(value ~ Condition, data = subset_data, alternative = "two.sided", 
                                                   mu = 0, paired = FALSE, var.equal = FALSE, conf.level = 0.95)
  ttestp[[paste0("portalNormal_AH",z)]] <- t.test(value ~ Condition, data = subset_data, alternative = "two.sided", 
                                                 mu = 0, paired = TRUE, var.equal = FALSE, conf.level = 0.95)
  # for (c in list_of_conditions) {
  #   subset_data <- subset(df_reshaped, Condition == c & Zone == z & gene_type %in% c("central", "portal"))
  #   ttest[[paste0("genes",c,z)]] <- t.test(value ~ gene_type, data = subset_data, alternative = "two.sided", 
  #                         mu = 0, paired = FALSE, var.equal = FALSE, conf.level = 0.95)
  # }
}

# Import data for combined gene
dfallnormal <- read_csv(paste0(mge_directory,"Normal_MGEall.csv"))
dfallAH <- read_csv(paste0(mge_directory,"AH_MGEall.csv"))
dfall<- rbind(dfallnormal,dfallAH)
dfall$Condition <- fct_relevel(dfall$Condition, "Normal", "AH")

subset_data <- subset(dfall, gene_type == "central"& Condition %in% c("Normal", "AH"))
x <- subset(subset_data, Condition =="Normal") %>% select(Expression)
y <- subset(subset_data, Condition =="AH") %>% select(Expression)
x = x$Expression
y = y$Expression
# ttest[[paste0("central")]] <- t.test(Expression ~ Condition, data = subset_data, alternative = "less", 
#                                      mu = 0, paired = FALSE, var.equal = FALSE, conf.level = 0.95)

ttestp[[paste0("central")]] <- t.test(x,y, alternative = "greater",
                                                mu = 0, paired = TRUE, var.equal = FALSE, conf.level = 0.95)

subset_data <- subset(dfall, gene_type == "portal" & Condition %in% c("Normal", "AH"))
x <- subset(subset_data, Condition =="Normal") %>% select(Expression)
y <- subset(subset_data, Condition =="AH") %>% select(Expression)

x = x$Expression
y = y$Expression
# ttest[[paste0("portal")]] <- t.test(Expression ~ Condition, data = subset_data, alternative = "less", 
#                                      mu = 0, paired = FALSE, var.equal = FALSE, conf.level = 0.95)
ttestp[[paste0("portal")]] <- t.test(x,y, alternative = "greater", 
                                    mu = 0, paired = TRUE, var.equal = FALSE, conf.level = 0.95)
                                               
# Extract and save results
t_test_unpaired<- map_dfr(names(ttest), ~ tibble(
  Comparison = .x,
  P_value = ttest[[.x]]$p.value,
  Statistic = ttest[[.x]]$statistic,
  Mean_g1 = ttest[[.x]]$estimate[[1]],
  Mean_g2 = ttest[[.x]]$estimate[[2]]
))

t_test_paired<- map_dfr(names(ttestp), ~ tibble(
  Comparison = .x,
  P_value = ttestp[[.x]]$p.value,
  Statistic = ttestp[[.x]]$statistic,
  Mean_diff = ttestp[[.x]]$estimate[[1]],
))


write_csv(t_test_paired, paste0(stats_directory, "ttestoncentralportalpaired.csv"))
write_csv(t_test_unpaired, paste0(stats_directory, "ttestoncentralportalunpaired.csv"))
print("Statistical tests completed. Results saved in the extractedData directory.")
