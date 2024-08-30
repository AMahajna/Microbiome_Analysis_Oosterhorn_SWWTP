################################################################################
##Create folder for project organization

if(!dir.exists("input_data")){dir.create("input_data")}
if(!dir.exists("input_data/func_result")){dir.create("input_data/func_result")}
if(!dir.exists("input_data/org_result")){dir.create("input_data/org_result")}
if(!dir.exists("input_data/subsystems_result")){dir.create("input_data/subsystems_result")}

if(!dir.exists("output_data")){dir.create("output_data")}
if(!dir.exists("figures")){dir.create("figures")}
if(!dir.exists("scripts")){dir.create("scripts")}

################################################################################
##load packages 

#source(file = "scripts/install_load_packages.r")

install.packages("biomformat")
library(biomformat)

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

#BiocManager::install(version = "3.19")  # Ensure you have the latest Bioconductor
BiocManager::install()  # Update all installed packages to their latest versions

BiocManager::install("phyloseq")
BiocManager::install("biomformat")
BiocManager::install("TreeSummarizedExperiment")
BiocManager::install("microbiomeTree")
BiocManager::install("mia")
library(microbiomeTree)
library(TreeSummarizedExperiment)
library(phyloseq)
library(biom)
library(phyloseq)
library(treeio)
library(mia)
# Install the data.table package if it's not already installed
install.packages("data.table")

# Load the data.table package
library(data.table)
################################################################################

biom_data = biomformat::read_biom("input_data/Galaxy1760-[Bracken-biom_output_file_(including_metadata)].biom1")
tse <- makeTreeSEFromBiom(biom_data)
tse <- transformAssay(tse, method = "relabundance")

#assay(tse, "relabundance") |> head()
#check equal one all 
#colSums(assay(tse, "relabundance"))

################################################################################
#rowData is taxonom y-> change column names 
#rownames is NCBI Tax ID 
#colData is Sample information -> change rowname 
#colnames

# Change column names in rowData
row_data_new <- rowData(tse)
colnames(row_data_new) <- c('Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus','Species')

# Assign the modified rowData back to the SummarizedExperiment object
rowData(tse) <- row_data_new

# Change row names in colData
col_data_new <- colData(tse)
rownames(col_data_new) <- paste0("Sample_", 1:nrow(col_data_new))

# Assign the modified colData back to the SummarizedExperiment object
colData(tse) <- col_data_new


################################################################################
##Create bar plots

################################################################################
##Read samsa2 output 

# Define the path to your folder containing TSV files
folder_path_org <- "input_data/org_result"

# List all TSV files in the folder
file_list_org <- list.files(folder_path_org, pattern = "\\.tsv$", full.names = TRUE)

# Function to read and preprocess each file
read_and_preprocess <- function(file) {
  # Extract the file name without extension
  file_name <- tools::file_path_sans_ext(basename(file))
  
  # Read the TSV file
  df <- fread(file, header = FALSE, sep = "\t")
  
  # Rename columns
  colnames(df) <- c(paste0("rel_abundance_", file_name),
                    paste0("reads_", file_name),
                    "species")
  
  # Ensure 'species' column is character type for consistency
  df[, species := as.character(species)]
  
  # Return the data.table
  return(df)
}

# Read and preprocess all files
data_list_org <- lapply(file_list_org, read_and_preprocess)

# Merge all data frames by taxonomy column
merged_data_org <- Reduce(function(x, y) merge(x, y, by = "species", all = TRUE), data_list_org)

#generate a list of unique organisms 
unique_org = unique(merged_data_org$species)

# Replace NA with 0 in the merged data frame
merged_data_org[is.na(merged_data_org)] <- 0

#Check if transformation was done correctly
#sum of rel_abundance columns if 100
columns_sum_org = colSums(merged_data_org[, -1], na.rm = TRUE)

# Save the merged data to a new file (optional)
fwrite(merged_data_org, "output_data/merged_data_org.tsv", sep = "\t")

# Display the first few rows of the merged data
head(merged_data_org)




