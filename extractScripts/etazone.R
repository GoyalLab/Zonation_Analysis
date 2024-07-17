#### Scripts to find PortoCentral Coordinates and Zones for Liver Zonation Analysis in Single-Cell RNA Sequencing Data 
#### Created by Aurelia Leona, 2024-07-16
#### Adapted from [original source if applicable]

# Set random seed for reproducibility
set.seed(23)

# Initialize directories
getwd()
rawfileDir <- "/projects/b1042/GoyalLab/aleona/Zonation_Analysis/rawData/"
plotDir <- "/projects/b1042/GoyalLab/aleona/Zonation_Analysis/plots/"
extrdataDir <- "/projects/b1042/GoyalLab/aleona/Zonation_Analysis/extractedData/"
ascriptsdir <- "/projects/b1042/GoyalLab/aleona/Zonation_Analysis/AnalysisScripts/"
plotscriptsdir <- "/projects/b1042/GoyalLab/aleona/Zonation_Analysis/plotScripts/"
extScriptsdir <- "/projects/b1042/GoyalLab/aleona/Zonation_Analysis/extractScripts/"

# Import custom functions for eta calculation
source(paste0(extScriptsdir, "etazonefunction.R"))

# Define new directories
eta_directory <- paste0(extrdataDir, "eta_data/")
etawithzones_directory <- paste0(extrdataDir, "etawithzones_data/")
seurat_obj_directory <- paste0(extrdataDir, "seurat_objects/")

# Create directories
create_dir_if_not_exists(eta_directory)
create_dir_if_not_exists(etawithzones_directory)
create_dir_if_not_exists(seurat_obj_directory)


# Preprocessing
# Import initial dataset
LD <- readRDS(file = paste0(rawfileDir, '032722_final.rds'))

# Import zonation parameters from Shalev paper
ZonationParams <- readMat(paste0(rawfileDir,"Zonation_params.mat"))

# Subset data to include only hepatocytes
hepLD <- subset(x = LD, idents = "Hepatocyte")

# Divide dataset into three conditions: Normal, AC, and AH
hep_normal <- subset(hepLD, subset = Condition == "Normal")
hep_ac <- subset(hepLD, subset = Condition == "AC")
hep_ah <- subset(hepLD, subset = Condition == "AH")

# Initialize lists to collect information
list_subset <- list(hep_normal, hep_ac, hep_ah)
names(list_subset) <- c("Normal", "AC", "AH")


# Create a dataframe with all cell IDs from the original dataset
cell_id <- as.data.frame(LD@assays[["RNA"]]@counts@Dimnames[[2]])
colnames(cell_id) <- "cell_id"

# Process each condition
counter <- 1
for (cells in names(list_subset)) {
  hep <- list_subset[[cells]]
  
  # Normalize the dataset
  hep <- prepzone(hep, assaytorun = "RNA")
  
  # Calculate eta (porto-central coordinate) for each cell
  eta <- findeta(hep, ZonationParams)
  
  # Create dataframe with eta values
  etadf <- as.data.frame(eta)
  colnames(etadf) <- "eta"
  
  #Add the metadata to each of the subsetted cells 
  hep <- AddMetaData(hep, metadata = etadf, col.name = "porto-central_coord")
  
  #Change format 
  etadf <- rownames_to_column(etadf, var = "cell_id")
  # Save eta data
  eta_filesavepath <- paste0(eta_directory, "eta_", cells, ".csv")
  write.csv(etadf, eta_filesavepath, row.names = FALSE)
  print(paste("Saved eta data for", cells, "to", eta_filesavepath))
  
  # Prepare dataframes for further analysis
  etawithzones <- etadf
  condition <- data.frame(
    cell_id = hep@assays[["RNA"]]@data@Dimnames[[2]],
    Condition = rep(cells, length(hep@assays[["RNA"]]@data@Dimnames[[2]]))
  )
  
  # Assign zones based on eta values
  print(paste0("In eta there are ", colSums(!is.na(etadf)), " non-NA values"))
  etawithzones <- etawithzones %>%
    mutate(Zone = case_when(
      eta <= 1/3 ~ "Zone 1",
      eta > 1/3 & eta <= 2/3 ~ "Zone 2",
      eta > 2/3 ~ "Zone 3"
    ))
  etawithzones <- left_join(etawithzones, condition, by = "cell_id")
  
  # Save etawithzones data
  etawithzones_filesavepath <- paste0(etawithzones_directory, "etawithzones_", cells, ".csv")
  write.csv(etawithzones, etawithzones_filesavepath, row.names = FALSE)
  print(paste("Saved etawithzones data for", cells, "to", etawithzones_filesavepath))
  
  
  # Combine eta values across all conditions
  if (counter == 1) {
    etaall <- etadf
    etawithzonesall <- etawithzones
  } else {
    etaall <- rbind(etaall, etadf)
    etawithzonesall <- rbind(etawithzonesall, etawithzones)
  }
  print(paste0("In combined dataset there are ", colSums(!is.na(etaall)), " non-NA values"))
  
  counter <- counter + 1
}

# Prepare final dataframe with porto-central coordinates
cell_id <- left_join(cell_id, etaall, by = "cell_id")
cell_id <- column_to_rownames(cell_id, var = "cell_id")
print(paste0("In cell_id combined dataset there are ", colSums(!is.na(cell_id)), " non-NA values"))

# Add porto-central coordinates to the original Seurat object
LD <- AddMetaData(LD, metadata = cell_id, col.name = "porto-central_coord")

# Check availability of zonation marker genes in the dataset
list_genecheck <- check_available_genes(hep, ZonationParams)

# Save the combined results
genes_to_exclude <- do.call(rbind, lapply(names(list_genecheck), function(name) {
  data.frame(category = name, gene = list_genecheck[[name]])
}))
write.csv(genes_to_exclude, file = paste0(extrdataDir, "combined_gene_availability.csv"), row.names = FALSE)


# After the loop, save Seurat objects
saveRDS(LD, file = paste0(seurat_obj_directory, "LD_with_portocentral_coords.rds"))
saveRDS(hep_normal, file = paste0(seurat_obj_directory, "hepnormalwithpc.rds"))
saveRDS(hep_ac, file = paste0(seurat_obj_directory, "hepacwithpc.rds"))
saveRDS(hep_ah, file = paste0(seurat_obj_directory, "hepahwithpc.rds"))
print("Saved Seurat objects to seurat_objects directory")

# Save the combined eta with zones dataframe as CSV
write.csv(etawithzonesall, paste0(etawithzones_directory, "combined_etawithzones.csv"), row.names = FALSE)
print("Saved combined etawithzones data")
write.csv(etaall, paste0(eta_directory, "combined_eta.csv"), row.names = FALSE)
print("Saved combined eta data")

print("Processing complete. Results saved in the extractedData directory.")
