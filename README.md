# Code Repository for Zonation Analysis of Liver Samples (Normal, Alcoholic Cirrhosis [AC], and Alcoholic Hepatitis [AH])
Created by AL on 20240716 and Last edited on 20240820

Inspired by Halpern, K., Shenhav, R., Massalha, H. et al. Paired-cell sequencing enables spatial gene expression mapping of liver endothelial cells. Nat Biotechnol 36, 962–970 (2018). (https://doi.org/10.1038/nbt.4231)

## Installation and Dependencies
1. Download R  using the link: https://cran.r-project.org/doc/manuals/r-patched/R-admin.html
2. Scripts in this repository has been implemented in R version 4.2.3 and can be installed using
```
    $> module load git
    $> git clone https://github.com/GoyalLab/Zonation_Analysis.git
```
3. Wishbone depends on a number of `R` packages available on CRAN and Bioconductor and these dependencies are listed in each of the function scripts as shown in the steps below.  All the dependencies will be automatically installed when running the steps below.
4. The runtime to download this github repository is less than 1 hour. The runtime to run all the scripts in the repository is around 1 hour.

## Disclaimer 
Please do start the code on a separate environment, in case there are conflicts occuring due to library differences. 

## Input files required
Seurat Object as an RDS file and Zonation_params.mat located in rawData file

## Steps Overview
1. Find the portocentral coordinates and zonation & euclidean distances between cells in each Zone in each Sample
3. Generate the Normalized Mean Gene Expression Matrix and Sorted Gene Expression based on Variation across Zones 
4. Plotting the Datasets
5. Statistical Analysis of each Zone in each Condition

## Step 1: Find the portocentral coordinates and zonation using etazonefunc & euclidean distances between cells in each Zone in each Sample
### Scripts used:
"etazone.R" "etazonefunction.R" "eucdist.R"
### Directories specified:
Input folder called rawData of .rds data and Zonation_params.mat

Output of extracted datasets including a folder of seurat objects with "porto-central_coord" as a new metadata, a folder of csv files of the porto-central coordinates and their assigned zones (eta_data and etawithzones_data), a csv file of all the zonation genes annotated in the dataset (ZonationGenesavailable.csv), and a folder of csv files of zonation ratios and their assigned zones (x_data and xwithzones_data). 

### Instructions
1. Open etazone.R and eucdist.R
2. Modify the path to the directory that contains both the main and functions (extScripts).
3. Modify the input directory path to match the location of your .rds files. Ensure that the 'Zonation_params.mat' file is located in the raw data directory (rawfiledir).
4. Specify the output directory (extrDataDir) where you want the generated .csv and .rds files to be saved. If you have multiple subfolders for output, make sure to specify those as well. For testing or as an example, you can use this [sample dataset](https://drive.google.com/drive/folders/1TEzlIAO7JGU2BAvdV0upSl1hMuyH05BR?usp=drive_link) as input for this step. If you have a more complex directory structure, adjust the paths in the script accordingly to match your specific setup. For instance, you might have separate folders for raw data, .rds files, and different types of output files.
```
rawData/
├── Zonation_params.mat
inputdata/
├── sample_data.rds
extData/
├── seurat_objects
│   └── seurat_with_coordinate.rds
├── eta_data
│   └── eta_data.csv
├── etawithzones_data
│   └── etawithzones_data.csv
├── x_data
│   └── x_data.csv
├── xwithzones_data
│   └── xwithzones_data.csv
├── ZonationGenesavailable.csv
```

```
      rawfileDir <- "/path/to/rawData"
      extrdataDir <- "/path/to/extData"
      input <- "/path/to/rds"
      
      #Rename the .rds dataset
      LD <- readRDS(file = paste0(input, '032722_final.rds'))
      
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
5. Run the eucdist.R script next.


## Step 2: Generate the Normalized Mean Gene Expression Matrix and Sorted Gene Expression based on Variation across Zones 
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
   
## Step 3: Plotting Datasets 
### Scripts used:
all functions in the plotScripts file and "etafunctions.R"
### Directories specified:
Input are all the generated datasets from the previous steps, notably combined_etawithzones.csv, combined_distances.csv, combined_gene_availability.csv, ZonationGenesavailable.csv, plot_zones.csv, and seurat_objects folder. 

Output are .svg and .png files of all the datasets generated above
### Instructions:
1. Modify the directory to match where the .csv files are located (extrDatadir) and where you would like the output subfolders of .csv files to be located (plotDir and extrDatadir). If you have multiple subfolders, specify those as well.
2. The script will also save the modified datasets used for plotting.

Here are a list of plots and their functions
- plotetadist.R -> Distribution of porto-central coordinates per cell in each condition,
    - p-value obtained from K-S test on the distribution comparing conditions to Normal.
    - Q-Q plot comparing conditions to Normal and Normal-Half 1 to Normal-Half 1 as control. 
- plotcellsperzone.R -> Percentage of Cells per Zone Across Conditions, use statsbootstrapcellsperzone.R below for statistical analysis on this data
- ploteucdist.R -> Euclidean Distances between Cells within Zones in PC Space
- plotheatmapzones.R -> Heatmap of Euclidean Distances between Zones of Each Condition
- plotZoneDistances.R -> Line Graph of Euclidean Distances between Zones of Each Condition
- plotmgezonationgenesall.R -> Scatter Plot of Overall Mean Expression of Zonation Genes  for Normal and SAH Condition

## Step 4: Statistical Analysis of each Zone in each Condition
### Scripts used:
"statstesteucdist.R" "statszonationgeneall.R" "statszonationgene_nonzero.R" "etazonefunction.R" "statsfunc.R" "statsbootstrapcellsperzone.R"
### Directories specified:
Input folder of .csv file with the marked zones and coordinates (combined_etawithzones.csv) and "seurat_objects" file. 

Output folder containing .csv files of the calculated dataset, plots and reports of the analysis done (stats_tests). 

### Instructions
1. Open statstesteucdist.R, statszonationgene.R, statsbootstrapcellsperzone.R and statsfunc.R
2. Modify the directory to match where the combined_etawithzones.csv file  and seurat_objects folder are located (extrDatadir) and where you would like the output subfolders of .csv files to be located (extrDataDir). If you have multiple subfolders, specify those as well.
3. Execute code cells that will source the function scripts as it will install and load all the packages needed.
4. Run the eucdist.R before running the statstest.R as outputs from eucdist.R will be inputs for statstest.R.
