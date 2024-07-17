#### Scripts for Porto-Central Coordinate Calculation and Gene Variance Analysis in Liver Zonation
#### Created by Aurelia Leona, 2024-07-16
#### Adapted from [original source if applicable]

# Initialize directories
getwd()
rawfileDir <- "/projects/b1042/GoyalLab/aleona/Zonation_Analysis/rawData/"
plotDir <- "/projects/b1042/GoyalLab/aleona/Zonation_Analysis/plots/"
extrdataDir <- "/projects/b1042/GoyalLab/aleona/Zonation_Analysis/extractedData/"
ascriptsdir <- "/projects/b1042/GoyalLab/aleona/Zonation_Analysis/AnalysisScripts/"
plotscriptsdir <- "/projects/b1042/GoyalLab/aleona/Zonation_Analysis/plotScripts/"
extScriptsdir <- "/projects/b1042/GoyalLab/aleona/Zonation_Analysis/extractScripts/"

# Import necessary functions
source(paste0(extScriptsdir, "gvefunction.R"))  
source(paste0(extScriptsdir, "etazonefunction.R"))  

# Import the combined eta zones data
etawithzones_directory <- paste0(extrdataDir, "etawithzones_data/")
etazones <- read_csv(paste0(etawithzones_directory, "combined_etawithzones.csv"))

# Import Seurat objects (assuming they're stored as RDS files)
seurat_obj_directory <- paste0(extrdataDir, "seurat_objects/")

list_subset <- list(
  Normal = readRDS(paste0(seurat_obj_directory, "hepnormalwithpc.rds")),
  AC = readRDS(paste0(seurat_obj_directory, "hepacwithpc.rds")),
  AH = readRDS(paste0(seurat_obj_directory, "hepahwithpc.rds"))
)

#Define new directories 
mge_directory <- paste0(extrdataDir, "test/mge/")
pmat_directory <- paste0(extrdataDir, "test/pmat/")
sorteddata_directory <- paste0(extrdataDir, "test/sorted_data/")

create_dir_if_not_exists(mge_directory)
create_dir_if_not_exists(pmat_directory)
create_dir_if_not_exists(sorteddata_directory)
# 
# mge_directory <- paste0(extrdataDir, "mge/")
# pmat_directory <- paste0(extrdataDir, "pmat/")
# sorteddata_directory <- paste0(extrdataDir, "sorted_data/")

#Create directories 
create_dir_if_not_exists(mge_directory)
create_dir_if_not_exists(pmat_directory)
create_dir_if_not_exists(sorteddata_directory)

##To double check make a list 
sorted_list <- list()
mge_list <- list()
pmat_list <- list()

# Process each condition
for (condition in names(list_subset)) {
  # Extract etazones for the current condition
  etz <- etazones %>% filter(Condition == condition)
  
  # Compute Pmat and get subsampled cells
  result <- compute_pmat_mge(etz)
  Pmat <- result$Pmat
  subsampled_cells <- result$subsampled_cells
  
  # Compute MGE
  heps <- list_subset[[condition]]
  mat_norm <- as.matrix(heps@assays[["RNA"]]@data)
  mat_norm_subsampled <- mat_norm[, subsampled_cells]
  MGE <- mat_norm_subsampled %*% Pmat
  
  # Compute sorted data based on variance
  sorted_data <- compute_var(MGE, heps)
  
  # Save results as CSV
  # Save Pmat
  Pmat_with_rownames <- Pmat %>% 
    as.data.frame() %>% 
    rownames_to_column(var = "cell_id")
  write_csv(Pmat_with_rownames, paste0(pmat_directory, condition, "_Pmat.csv"))
  
  # Save MGE
  MGE_with_rownames <- MGE %>% 
    as.data.frame() %>% 
    rownames_to_column(var = "gene")
  write_csv(MGE_with_rownames, paste0(mge_directory, condition, "_MGE.csv"))
  
  # Save sorted_data (assuming it's already a dataframe with a 'gene' column)
  write_csv(sorted_data, paste0(sorteddata_directory, condition, "_sorted_data.csv"))
  
  print(paste("Processed and saved data for", condition))
  
  sorted_list[[condition]] <- sorted_data
  mge_list[[condition]]  <- MGE_with_rownames
  pmat_list[[condition]]  <- Pmat_with_rownames
}

print("Processing complete. Results saved as CSV files.")