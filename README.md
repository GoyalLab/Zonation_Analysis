# Code Repository for Zonation Analysis of Liver Samples (Normal, Alcoholic Cirrhosis [AC], and Alcoholic Hepatitis [AH])
## Paper name .... 
Created by AL on 20240716

## Input files required
Seurat Object as an RDS file and Zonation_params.mat located in rawData file

## Steps Overview
1. Find the portocentral coordinates and zonation
2. Statistical Analysis of each Zone in each Condition
3. Generate the Normalized Mean Gene Expression Matrix and Sorted Gene Expression based on Variation across Zones 
4. Plotting the Datasets 

## Step 1: Find the portocentral coordinates and zonation using etazonefunc
### Scripts used:
"etazone.R" "etazonefunction.R" 
### Directories specified:
Input folder called rawData of .rds data and Zonation_params.mat

Output of extracted datasets including a folder of seurat objects with "porto-central_coord" as a new metadata, a folder of csv files of the porto-central coordinates, a csv file with both the coordinates and assigned zonation, and a csv files of all the zonation genes annotated in the dataset. 

### Instructions
1. Open etazone.R 
2. Modify the path to the directory that contains both the main and functions (extScripts).
3. Modify the directory to match where the .rds files are located and get the Zonation_params.mat from the rawData file (rawfiledir) and where you would like the output subfolders of .csv and .rds files to be located (extrDataDir). If you have multiple subfolders, specify those as well. For example:

```
      #Rename the .rds dataset
      LD <- readRDS(file = paste0(rawfileDir, '032722_final.rds'))
      
      #Define new directories 
      eta_directory <- paste0(extrdataDir, "eta_data/")

      ##Use the function in etafunctions.R to create the new directories or subfolders if they don't exist
      create_dir_if_not_exists(eta_directory)
```
3.Execute code cells that will source the function scripts as it will install and load all the packages needed. For example:
```
    source(paste0(extScriptsdir, "etazonefunction.R"))
```
4. Run the etazone.R script. This will subset Seurat object data based on the three conditions and run the seurat normalization pipeline. It will output multiple subfolders to be used later on. 

## Step 2: Statistical Analysis of each Zone in each Condition
### Scripts used:
"eucdist.R" "statstest.R" "etazonefunction.R" "statsfunc.R"
### Directories specified:
Input folder of .csv file with the marked zones and coordinates (combined_etawithzones.csv) and "seurat_objects" file. 

Output folder containing .csv files of the calculated dataset and reports of the analysis done. 

### Instructions
1. Open eucdist.R and statstest.R
2. Modify the directory to match where the combined_etawithzones.csv file  and seurat_objects folder are located (extrDatadir) and where you would like the output subfolders of .csv files to be located (extrDataDir). If you have multiple subfolders, specify those as well.
3. Execute code cells that will source the function scripts as it will install and load all the packages needed.
4. Run the eucdist.R before running the statstest.R as outputs from eucdist.R will be inputs for statstest.R.


## Step 3: Generate the Normalized Mean Gene Expression Matrix and Sorted Gene Expression based on Variation across Zones 
### Scripts used:
“gve.R" "gvefunction.R" "etazonefunction.R" 
### Directories specified
Input folder of .csv file with the marked zones and coordinates (combined_etawithzones.csv) and "seurat_objects" file. 

Output directory for holding .csv table of all matrices generated (Pmat, mge, sorted_data)

### Instructions:
1. Open “gve.R"
2. Modify the directory to match where the combined_etawithzones.csv file  and seurat_objects folder are located (extrDatadir) and where you would like the output subfolders of .csv files to be located (extrDataDir). If you have multiple subfolders, specify those as well.
3. There are code cells in the script that will source the function scripts as it will install and load all the packages needed.
4. Execute the script and it will output folders containing matrices of each condition of gene expression across all three zones and a sorted list of gene expression based on their variation across zones
   
## Step 4: Plotting Datasets 
### Scripts used:
all functions in the plotScripts file and "etafunctions.R"
### Directories specified:
Input are all the generated datasets from the previous steps, notably combined_etawithzones.csv, combined_distances.csv, combined_gene_availability.csv, and plot_zones.csv

Output are .svg and .png files of all the datasets generated above
### Instructions:
1. Modify the directory to match where the .csv files are located (extrDatadir) and where you would like the output subfolders of .csv files to be located (plotDir and extrDatadir). If you have multiple subfolders, specify those as well.
2. The script will also modified versions of the datasets for plotting. 


