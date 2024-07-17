#### Script for Subsampling and Calculating Pairwise Euclidean Distances between the Cells in each Zones in each Conditions
#### Created by Aurelia Leona, 2024-07-16
#### Adapted from [original source if applicable]

# Load necessary libraries
library(Seurat)
library(dplyr)
library(tidyr)

# Initialize directories
getwd()
rawfileDir <- "/projects/b1042/GoyalLab/aleona/Zonation_Analysis/rawData/"
plotDir <- "/projects/b1042/GoyalLab/aleona/Zonation_Analysis/plots/"
extrdataDir <- "/projects/b1042/GoyalLab/aleona/Zonation_Analysis/extractedData/"
plotscriptsdir <- "/projects/b1042/GoyalLab/aleona/Zonation_Analysis/plotScripts/"
extScriptsdir <- "/projects/b1042/GoyalLab/aleona/Zonation_Analysis/extractScripts/"


# Import custom functions for eta calculation
source(paste0(extScriptsdir, "etazonefunction.R"))

# Define new directories
compare_zones_dir <- paste0(extrdataDir, "Compare_Zones/")

# Create directories
create_dir_if_not_exists(compare_zones_dir)

# Import Seurat objects (assuming they're stored as RDS files)
seurat_obj_directory <- paste0(extrdataDir, "seurat_objects/")
list_subset <- list(
  Normal = readRDS(paste0(seurat_obj_directory, "hepnormalwithpc.rds")),
  AC = readRDS(paste0(seurat_obj_directory, "hepacwithpc.rds")),
  AH = readRDS(paste0(seurat_obj_directory, "hepahwithpc.rds"))
)

# Import the combined eta zones data
etazones_directory <- paste0(extrdataDir, "etawithzones_data/")
etazones <- read_csv(paste0(etazones_directory, "combined_etawithzones.csv"))

# Initialize lists to collect information
hds_list1 <- list()
hds_list2 <- list()
list_of_zones <- c("Zone 1", "Zone 2", "Zone 3")

# Based on the max number of dataset we want to subset for each condition 
# In this case, Normal has a max of 28 datapoints 
num_subsample <- c(14, 50, 50)


# Subsample random cells from each zone in each condition
for (cond in names(list_subset)) {
  counter <- 1
  for (whichzone in list_of_zones) {
    # Filter cells for the current zone and condition
    zone3 <- etazones %>%
      filter(Zone == whichzone & Condition == cond)
    naming <- paste0(whichzone, "_", cond)
    
    # First subsample
    set.seed(123)  # For reproducibility of random sampling
    zone3_1 <- zone3 %>%
      sample_n(num_subsample[[counter]])
    hds_list1[[naming]] <- zone3_1
    
    # Second subsample from remaining cells
    remaining_cells <- etazones %>%
      anti_join(zone3_1, by = "cell_id") %>%
      filter(Zone == whichzone & Condition == cond)
    
    set.seed(456)  # For reproducibility of random sampling
    zone3_2 <- remaining_cells %>%
      sample_n(num_subsample[[counter]])
    hds_list2[[naming]] <- zone3_2
    
    counter <- counter + 1
  }
}


# Initialize an empty list to store all distance dataframes
all_distances <- list()

# Calculate pairwise distances within each zone and condition, and save results
for (sample_name in names(list_subset)) {
  sample <- list_subset[[sample_name]]
  for (whichzone in list_of_zones) {
    naming <- paste0(whichzone, "_", sample_name)
    pca_data <- as.data.frame(sample@reductions$pca@cell.embeddings)
    
    # Combine subsampled cells
    random_1 <- pca_data[hds_list1[[naming]][["cell_id"]], ]
    random_2 <- pca_data[hds_list2[[naming]][["cell_id"]], ]
    randall <- rbind(random_1, random_2)
    
    # Calculate pairwise Euclidean distances using the new function
    eucdist <- euc_distances(t(randall))
    
    # Convert distance matrix to dataframe
    result_df <- as.data.frame(eucdist)
    colnames(result_df) <- rownames(randall)
    rownames(result_df) <- rownames(randall)
    
    # Save the distance matrix for this zone and condition
    savepath <- paste0(compare_zones_dir, naming, ".csv")
    write.csv(result_df, savepath)
    
    #Convert distance matrix to long format dataframe
    result_df2 <- as.data.frame(as.table(eucdist)) %>%
      rename(cell1 = Var1, cell2 = Var2, eucdist = Freq) %>%
      filter(as.character(cell1) < as.character(cell2))  # Keep only lower triangle
    
    # Add zone and condition information
    result_df2$Zone <- whichzone
    result_df2$Condition <- sample_name
    
    # Store in list
    all_distances[[naming]] <- result_df2
    
    print(paste("Saved distance matrix for", naming))
  }
}

# Combine all distance dataframes into a single dataframe
combined_distances <- bind_rows(all_distances)

# Print the dimensions of the combined dataframe
print(paste("Combined dataframe dimensions:", 
            nrow(combined_distances), "rows by", 
            ncol(combined_distances), "columns"))

# Save this combined dataframe 
write_csv(combined_distances, paste0(compare_zones_dir, "combined_distances.csv"))

print("Processing complete. Distance matrices saved in the Compare_Zones directory.")
