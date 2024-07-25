#### Functions for finding PortoCentral Coordinates and Zones for Liver Zonation Analysis in Single-Cell RNA Sequencing Data 
#### Created by Aurelia Leona, 2024-07-16
#### Adapted from [original source if applicable]


# Load required libraries
required_packages <- c("Seurat", "R.matlab", "dplyr", "Matrix", "tidyverse", "tibble")

# Install missing packages
new_packages <- required_packages[!(required_packages %in% installed.packages()[,"Package"])]
if(length(new_packages)) install.packages(new_packages)

# Load all packages
lapply(required_packages, library, character.only = TRUE)

library(Seurat)        # For handling single-cell RNA-seq data (Seurat objects)
library(R.matlab)      # For reading .mat files (readMat function)
library(dplyr)         # For data manipulation (often used with Seurat)
library(Matrix) 
library(tidyverse)
library(tibble)

###Prepare Seurat Object for Zonation Analysis---------------------------
#'
#' This function prepares a Seurat object for zonation analysis by performing
#' normalization, feature selection, scaling, and PCA.
#'
#' @param seurat_Obj A Seurat object containing single-cell RNA sequencing data.
#' @param assaytorun The name of the assay to use for analysis (e.g., "RNA").
#'
#' @return A processed Seurat object with normalized "RNA" assays.
#'
#' @details The function performs the following steps:
#'          1. Sets the default assay
#'          2. Normalizes data using LogNormalize method
#'          3. Finds variable features
#'          4. Scales data
#'          5. Runs PCA
prepzone<- function(seurat_Obj, assaytorun) {
  DefaultAssay(object = seurat_Obj) <- assaytorun
  seurat_obj<- NormalizeData(seurat_Obj, normalization.method = "LogNormalize") %>%
    FindVariableFeatures( selection.method = "vst", nfeatures =  2000, clip.max = "auto", binning.method = "equal_width")%>% 
    ScaleData(model.use = "linear", min.cells.to.block =  3000, block.size =  1000, scale.max= 10)%>%
    RunPCA(npcs = 50, weight.by.var = TRUE)
  print("normalized ")
  
  return(seurat_obj)
}

### Calculate Porto-Central Coordinates (Eta) for Liver Zonation--------------------

#' This function computes the porto-central coordinates (eta) for each cell in a liver single-cell RNA sequencing dataset.
#' Eta represents the relative position of a cell along the porto-central axis of the liver lobule,
#' based on the ratio of the expression of zone-specific marker genes.
#'
#' @param hepseurat A Seurat object containing hepatocyte cells of the liver single-cell RNA sequencing data.
#'                  The object should have an "RNA" assay with normalized expression data.
#'
#' @param ZonationParams A list containing two elements:
#'                       - genes.cv: A list of central vein (pericentral) marker genes
#'                       - genes.pn: A list of portal node (periportal) marker genes
#'                       Each element should be a nested list where gene names are stored as character vectors.
#'
#' @return A numeric vector of eta values, one for each cell in the dataset.
#'         Eta ranges from 0 (most pericentral) to 1 (most periportal).
#'
#' @details The function performs the following steps:
#'          1. Extracts and uppercases central vein and portal node marker genes from ZonationParams.
#'          2. Calculates the sum of expression for central vein and portal node genes for each cell.
#'          3. Computes a raw coordinate (x) as the ratio of portal node gene expression to total marker gene expression.
#'          4. Normalizes x to create eta, which ranges from 0 to 1.
#'
#' @seealso
#' Relevant papers or resources on liver zonation and spatial reconstruction in single-cell data.
findeta <- function(hepseurat, ZonationParams) {
  
  genes.cv <- list()
  for (i in 1: length(ZonationParams[["genes.cv"]]) ){
    genes.cv[[i]] <- toupper(ZonationParams[["genes.cv"]][[i]][[1]])
  }
  
  genes.pn <- list()
  for (i in 1: length(ZonationParams[["genes.pn"]]) ){
    genes.pn[[i]] <- toupper(ZonationParams[["genes.pn"]][[i]][[1]])
  }
  
  bool_cv <- hepseurat@assays[["RNA"]]@counts@Dimnames[[1]] %in% genes.cv
  cv_vec <- colSums(hepseurat@assays[["RNA"]]@data[bool_cv, ])
  
  bool_pn <- hepseurat@assays[["RNA"]]@counts@Dimnames[[1]] %in% genes.pn
  pn_vec <- colSums(hepseurat@assays[["RNA"]]@data[bool_pn, ])
  
  # Compute spatial coordinate (eta)
  x <- pn_vec / (cv_vec + pn_vec)
  X <- as.data.frame(x)
  X$x <- replace(X$x, which(is.nan(X$x)), 0 )
  x2 <- as.matrix(X)
  eta <- (x2 - min(x2)) / (max(x2) - min(x2))
  
  print("eta found")
  return(list(X = X, eta = eta))
}

### Check Availability of Zonation Marker Genes in Dataset-----------------------
#'
#' This function checks which of the specified central vein and portal node
#' marker genes are present in the single-cell RNA sequencing dataset.
#'
#' @param hepseurat A Seurat object containing hepatocyte cells of the liver single-cell RNA sequencing data.
#'                  The object should have an "RNA" assay with normalized expression data.
#' @param ZonationParams A list containing central vein and portal node marker genes.
#'
#' @return A list containing:
#'         - available_cv: Central vein marker genes present in the dataset
#'         - available_pn: Portal node marker genes present in the dataset
#'         - missing_cv: Central vein marker genes not found in the dataset
#'         - missing_pn: Portal node marker genes not found in the dataset
#'
check_available_genes <- function(hepseurat, ZonationParams) {
  # Extract genes from ZonationParams
  genes.cv <- unlist(lapply(ZonationParams[["genes.cv"]], function(x) toupper(x[[1]])))
  genes.pn <- unlist(lapply(ZonationParams[["genes.pn"]], function(x) toupper(x[[1]])))
  
  # Get available genes in the dataset
  available_genes <- hepseurat@assays[["RNA"]]@counts@Dimnames[[1]]
  
  # Find intersections and differences
  available_cv <- intersect(genes.cv, available_genes)
  available_pn <- intersect(genes.pn, available_genes)
  missing_cv <- setdiff(genes.cv, available_genes)
  missing_pn <- setdiff(genes.pn, available_genes)
  
  # Return results
  return(list(
    available_cv = available_cv,
    available_pn = available_pn,
    missing_cv = missing_cv,
    missing_pn = missing_pn
  ))
}

#' Create Directory if It Doesn't Exist
#'
#' This function creates a directory at the specified path if it doesn't already exist.
#'
#' @param dir_path A character string specifying the path of the directory to be created.
#'
#' @return No return value. The function creates the directory and prints a message if successful.
#'
#' @examples
#' create_dir_if_not_exists("path/to/new/directory")
create_dir_if_not_exists <- function(dir_path) {
  if (!dir.exists(dir_path)) {
    dir.create(dir_path, recursive = TRUE)
    print(paste("Created directory:", dir_path))
  }
}

#' Calculate Euclidean Distances
#'
#' This function calculates the Euclidean distances between all pairs of columns in a given data matrix.
#'
#' @param data A numeric matrix where columns represent different observations or entities.
#'
#' @return A symmetric matrix of Euclidean distances between all pairs of columns in the input data.
#'
#' @examples
#' data_matrix <- matrix(rnorm(20), nrow = 5, ncol = 4)
#' distances <- euc_distances(data_matrix)

euc_distances <- function(data) {
  n <- ncol(data)
  eucdist <- matrix(0, n, n)
  for (i in 1:(n-1)) {
    for (j in (i+1):n) {
      dist_ij <- sqrt(sum((data[,i] - data[,j])^2))
      eucdist[i,j] <- eucdist[j,i] <- dist_ij
    }
  }

  return(eucdist)
}

#' Categorize Gene
#'
#' This function categorizes a gene as "central", "portal", or "unknown" based on its presence in predefined lists of central vein and portal node genes.
#'
#' @param gene A character string representing the gene name to be categorized.
#' @param ZonationParams A list containing two elements: "genes.cv" and "genes.pn", each a list of central vein and portal node genes, respectively.
#'
#' @return A character string: "central" if the gene is in the central vein list, "portal" if it's in the portal node list, or "unknown" if it's in neither.
#'
#' @examples
#' ZonationParams <- list(
#'   genes.cv = list(c("GENE1", "GENE2")),
#'   genes.pn = list(c("GENE3", "GENE4"))
#' )
categorize_gene <- function(gene) {
  if (gene %in% unlist(lapply(ZonationParams[["genes.cv"]], function(x) toupper(x[[1]])))) {
    return("central")
  } else if (gene %in% unlist(lapply(ZonationParams[["genes.pn"]], function(x) toupper(x[[1]])))) {
    return("portal")
  } else {
    return("unknown")
  }
}
