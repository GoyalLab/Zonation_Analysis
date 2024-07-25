#### Functions for Porto-Central Coordinate Calculation and Gene Variance Analysis in Liver Zonation
#### Created by Aurelia Leona, 2024-07-16
#### Adapted from [original source if applicable]

# Load required libraries
required_packages <- c("Seurat", "Matrix", "tidyverse", "tibble")

# Install missing packages
new_packages <- required_packages[!(required_packages %in% installed.packages()[,"Package"])]
if(length(new_packages)) install.packages(new_packages)

# Load all packages
lapply(required_packages, library, character.only = TRUE)

# Load required libraries
library(tidyverse)  # For data manipulation and tidying
library(Seurat)     # For working with single-cell data
library(Matrix)     # For sparse matrix operations

#' Calculate Porto-Central Matrix (Pmat) and Subsample Cells
#'
#' This function computes the porto-central matrix (Pmat) and subsamples cells
#' based on their zonation.
#'
#' @param etazones A dataframe containing cell IDs, zones, and other relevant information.
#' @param n_subsample Number of cells to subsample from each zone (default: 28).
#'
#' @return A list containing:
#'   - Pmat: The computed porto-central matrix
#'   - subsampled_cells: Vector of subsampled cell IDs
#'
#' @examples
#' result <- compute_pmat_mge(etazones_data)
#' Pmat <- result$Pmat
#' subsampled_cells <- result$subsampled_cells
#' 
compute_pmat_mge <- function(etazones, n_subsample = 28) {
  # Create a unique list of all cell_ids and zones
  all_cell_ids <- unique(etazones$cell_id)
  all_zones <- c("Zone 1", "Zone 2", "Zone 3")
  
  # Create an empty matrix with cell_ids as rows and zones as columns
  zone_matrix <- matrix(0, nrow = length(all_cell_ids), ncol = length(all_zones))
  rownames(zone_matrix) <- all_cell_ids
  colnames(zone_matrix) <- all_zones
  
  # Fill the matrix with 1s based on the zone assignments
  for (i in 1:nrow(etazones)) {
    cell <- etazones$cell_id[i]
    zone <- etazones$Zone[i]
    zone_matrix[cell, zone] <- 1
  }
  
  # Subsampling step
  subsampled_cells <- c()
  total_samples <- n_subsample
  
  for (zone in all_zones) {
    set.seed(123)  # For reproducibility
    zone_cells <- rownames(zone_matrix)[zone_matrix[, zone] == 1]
    sampled_cells <- sample(zone_cells, n_subsample, replace = FALSE)
    subsampled_cells <- c(subsampled_cells, sampled_cells)
  }
  
  # Subset the zone_matrix to include only the subsampled cells
  zone_matrix_subsampled <- zone_matrix[subsampled_cells, ]
  
  # Compute pmat
  Pmat <- zone_matrix_subsampled / colSums(zone_matrix_subsampled)

  return(list(Pmat = Pmat, subsampled_cells = subsampled_cells))
}


# This function calculates the variance for each row in a dataframe
RowVar <- function(df) {
  # Calculate row-wise variances
  variances <- apply(df, 1, function(row) {
    var <- var(row)
    if (is.na(var)) var <- NA  # Handle cases where variance is NA
    return(var)
  })
  
  # Create dataframe with rownames and variances
  result_df <- data.frame(Variance = variances, row.names = rownames(df))
  
  return(result_df)
}

#' Compute Variance and Sort Genes
#'
#' This function computes the variance for each gene in the MGE data
#' and sorts the genes based on their variance in descending order.
#'
#' @param MGE A matrix of gene expression data (genes as rows, cells as columns).
#' @param hep_seurat A Seurat object (not used in the current implementation, 
#'                   but kept for potential future use or compatibility).
#'
#' @return A dataframe of genes sorted by variance in descending order, 
#'         including gene names and their corresponding variances.
#'
#' @examples
#' sorted_data <- compute_var(MGE_data, seurat_object)
compute_var <- function(MGE, hep_seurat) {
  data <- as.data.frame(MGE)
  datarowvar <- RowVar(data)
  data <- rownames_to_column(data, var="gene")
  datarowvar <- rownames_to_column(as.data.frame(datarowvar), var="gene")
  data <- left_join(data, datarowvar, by ="gene")
  sorted_data <- data[order(data$Variance, decreasing=TRUE), , drop=FALSE]
  return(sorted_data)
}

