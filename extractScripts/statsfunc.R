#### Functions for Statistical Tests on Euclidean Distances Between Zones
#### Created by Aurelia Leona, 2024-07-16
#### Adapted from [original source if applicable]

# Load required libraries
library(tidyverse)

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
perform_wilcoxon_test <- function(data, zone, condition1, condition2) {
  subset_data <- subset(data, Zone == zone & Condition %in% c(condition1, condition2))
  test_result <- wilcox.test(eucdist ~ Condition, data = subset_data)
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
perform_t_test <- function(data, zone, condition1, condition2) {
  subset_data <- subset(data, Zone == zone & Condition %in% c(condition1, condition2))
  test_result <- t.test(eucdist ~ Condition, data = subset_data, 
                        var.equal = FALSE, alternative = "less")
  return(test_result)
}