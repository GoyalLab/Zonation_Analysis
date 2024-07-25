#### Functions for Statistical Tests on Euclidean Distances Between Zones
#### Created by Aurelia Leona, 2024-07-16
#### Adapted from [original source if applicable]

# Load required libraries
required_packages <- c("tidyverse","tibble")

# Install missing packages
new_packages <- required_packages[!(required_packages %in% installed.packages()[,"Package"])]
if(length(new_packages)) install.packages(new_packages)

# Load all packages
lapply(required_packages, library, character.only = TRUE)

# Load required libraries
library(tidyverse)
library(tibble)

#' Perform Wilcoxon Rank Sum Test
#'
#' This function performs a Wilcoxon rank sum test on the Euclidean distances
#' between two conditions for a given zone.
#'
#' @param data A dataframe containing the Euclidean distances.
#' @param zone The zone to be tested.
#' @param condition1 The first condition to compare.
#' @param condition2 The second condition to compare.
#'
#' @return A list containing the test results.
wilcoxon_test_condition <- function(data, zone, condition1, condition2) {
  subset_data <- subset(data, Zone == zone & Condition %in% c(condition1, condition2))
  test_result <- wilcox.test(value ~ Condition, data = subset_data)
  return(test_result)
}

#' Perform T-Test
#'
#' This function performs a t-test on the Euclidean distances
#' between two conditions for a given zone.
#'
#' @param data A dataframe containing the Euclidean distances.
#' @param zone The zone to be tested.
#' @param condition1 The first condition to compare.
#' @param condition2 The second condition to compare.
#'
#' @return A list containing the test results.
t_test_condition <- function(data, zone,  condition1, condition2) {
  subset_data <- subset(data, Zone == zone & Condition %in% c(condition1, condition2))
  test_result <- t.test(value ~ Condition, data =subset_data, alternative = "two.sided", 
                        mu = 0, paired = FALSE, var.equal = FALSE, conf.level = 0.95)
  return(test_result)
}

wilcoxon_test_zones <- function(data, condition, zone1, zone2){
  subset_data <- subset(data, Condition == condition & Zone %in% c(zone1, zone2))
  test_result <- wilcox.test(value ~ Zone, data = subset_data)
  return(test_result)
}

t_test_zones <- function(data, condition, zone1, zone2) {
  subset_data <- subset(data, Condition == condition & Zone %in% c(zone1, zone2))
  test_result <- t.test(value ~ Zone, data = subset_data, alternative = "two.sided", 
                        mu = 0, paired = FALSE, var.equal = FALSE, conf.level = 0.95)
  return(test_result)
}


wilcoxon_test_genes <- function(data, condition, zone1, zone2){
  subset_data <- subset(data, Condition == condition & Zone %in% c(zone1, zone2))
  test_result <- wilcox.test(value ~ Zone, data = subset_data)
  return(test_result)
}

t_test_genes <- function(data, zone, condition, gene1 = "central", gene2 = "portal") {
  subset_data <- subset(data, Condition == condition & Zone == zone & gene_type %in% c(gene1, gene2))
  test_result <- t.test(value ~ Zone, data = subset_data, alternative = "two.sided", 
                        mu = 0, paired = FALSE, var.equal = FALSE, conf.level = 0.95)
  return(test_result)
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

