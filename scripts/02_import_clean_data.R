##Read biom data  

#bioformat package 
biom_data = biomformat::read_biom("input_data/Galaxy1783-[Kraken-biom_output_file].biom1")

#mia package 
tse <- makeTreeSEFromBiom(biom_data)
#assays(tse)
#tse <- transformAssay(tse, method = "relabundance")
#assay(tse, "relabundance") |> head()
#check equal one all 
#colSums(assay(tse, "relabundance"))

#removal_efficiency respective to sample to be added to column data in TSE
removal_efficiency  = read_excel("input_data/Removal_Efficiency.xlsx")

#Relevant process data to be added to column data in TSE
process_data <- read_csv(file = "input_data/relevant_process_data.csv", show_col_types = FALSE)

################################################################################
##preprocess taxonomic data  

#rowData is taxonomy-> change column names 
#rownames is NCBI Tax ID 
#colData is Sample information -> change rowname 

# Change column names in rowData
row_data_new <- rowData(tse)
colnames(row_data_new) <- c('Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus','Species')
# Assign the modified rowData back to the SummarizedExperiment object
rowData(tse) <- row_data_new

col_data_new <- colData(tse)
# Convert the character column to Date
col_data_new$Date <- as.Date(col_data_new$Date, format = "%d-%m-%Y")

# Change row names in colData
rownames(col_data_new) <- col_data_new$Date
# Convert the year character column to numeric 
col_data_new$Year <- as.numeric(col_data_new$Year)
col_data_new = as.data.frame(col_data_new)
#rownames(col_data_new) <- paste0("Sample_", 1:nrow(col_data_new))

process_data <- subset(process_data, select = -rarity)
process_data[, 10] <- (process_data[, 10])/100
#Standardize environmental data - reflect on standardization method choice 
#process_data_normalized <- decostand(process_data, method ="standardize")
process_data_normalized <- process_data  # Create a copy to store the result
process_data_normalized[ , !names(process_data) %in% "Capacity_blowers_%"] <- decostand(process_data[ , !names(process_data) %in% "Capacity_blowers_%"], method = "standardize")
#summary(process_data_normalized)

col_data = DataFrame(cbind(col_data_new,process_data_normalized,(removal_efficiency[ ,3:6])/100))

# Assign the modified colData back to the SummarizedExperiment object
colData(tse) <- col_data
#tse

tse_bacteria <- tse[
  rowData(tse)$Kingdom %in% c("k__Bacteria"), ]
#Check
#unique(rowData(tse_bacteria)$Kingdom)

################################################################################
################################################################################
##Read samsa2 output 

################################################################################
##organisms: reading and cleaning 

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
#fwrite(merged_data_org, "output_data/merged_data_org.tsv", sep = "\t")

# Display the first few rows of the merged data
#head(merged_data_org)

# Remove columns that start with "rel"
merged_data_org <- merged_data_org %>%
  select(-starts_with("rel"))

# Remove the prefix "reads_" from column names
names(merged_data_org) <- gsub("^reads_", "", names(merged_data_org))

# Save the merged data to a new file (optional)
#fwrite(merged_data_org, "output_data/merged_data_org_reads.tsv", sep = "\t")

## Creating tse_active for metabolically active organisms 

#Same samples as in tse thus it's colData(tse)
#rowData is species column of merged_data_org
#assay in the numeric values of merged_data_org

#First, change column names in merged_data_org
name_changes <- as.data.frame(colData(tse))[ , c("SRA_accession", "Date")]
name_changes[] <- lapply(name_changes, as.character)


# Use mapvalues to rename columns
colnames(merged_data_org) <- mapvalues(colnames(merged_data_org), 
                                       from = name_changes$SRA_accession, 
                                       to = name_changes$Date, 
                                       warn_missing = FALSE)


# Rename the column
colnames(merged_data_org)[colnames(merged_data_org) == "species"] <- "Species"

# Specify the new order for the first four columns
new_order_first <- c("Species")

# Get the order of the remaining columns from df2
remaining_cols_order <- colnames(assay(tse))

# Reorder columns in df1
merged_data_org_reordered <- merged_data_org %>%
  select(all_of(new_order_first), all_of(remaining_cols_order))

rowData_active = as.data.frame(merged_data_org_reordered[ ,1])
#colnames(rowData_active) = "Species"
rownames(rowData_active) = rowData_active$Species

counts_active = as.matrix(merged_data_org_reordered[ ,2:33]) 
rownames(counts_active) = rowData_active$Species

#colData_active = as.data.frame(colData(tse))
colData_active = colData(tse)
#rownames(colData_active) = colnames(counts_active) 

#Condition:  
#colnames(counts_active) == rownames(colData_active)
#rownames(rowData_active) == rownames(counts_active) 

tse_active <- TreeSummarizedExperiment(assays = list(counts = counts_active ),
                                       colData = colData_active,
                                       rowData = rowData_active, 
                                       
)

################################################################################
##metabolism: reading and cleaning 

# Define the path to your folder containing TSV files
folder_path_metabolism <- "input_data/subsystems_result"

# List all TSV files in the folder
file_list_metabolism <- list.files(folder_path_metabolism, pattern = "\\.reduced$", full.names = TRUE)

# Function to read and preprocess each file
read_and_preprocess <- function(file) {
  # Extract the file name without extension
  file_name <- tools::file_path_sans_ext(basename(file))
  
  # Read the TSV file
  df <- fread(file, header = FALSE, sep = "\t", fill = TRUE)
  
  # Check if the data frame has 7 columns
  if (ncol(df) == 7) {
    # Drop the 7th column
    if (is.data.table(df)) {
      df <- df[, -7, with = FALSE]  # For data.table
    } else {
      df <- df[, -7]  # For data.frame
    }
  }
  
  # Rename columns
  names(df)[1:6] <- c(paste0("rel_abundance_", file_name),
                      paste0("reads_", file_name),
                      "enzyme",
                      "SEED_subsystem",
                      "SEED_subsystem_functional_category", 
                      "description")
  
  # Ensure 'protein' column is character type for consistency
  df[ , enzyme := as.character(enzyme)]
  
  # Return the data.table
  return(df)
}

# Read and preprocess all files
data_list_metabolism <- lapply(file_list_metabolism, read_and_preprocess)

# Merge all data frames by taxonomy column
merged_data_metabolism <- Reduce(function(x, y) merge(x, y, by = c("enzyme",
                                                                   "SEED_subsystem",
                                                                   "SEED_subsystem_functional_category", 
                                                                   "description"), all = TRUE), data_list_metabolism)

#generate a list of unique organisms 
#unique_metabolism = unique(merged_data_metabolism$enzyme)

#Replace NA with  in the merged data frame
merged_data_metabolism[is.na(merged_data_metabolism)] <- 0

#Check if transformation was done correctly
#sum of rel_abundance columns if 100
columns_sum_metabolism = colSums(merged_data_metabolism[, -c(1,2,3,4)], na.rm = TRUE)
columns_sum_metabolism

# Save the merged data to a new file (optional)
#fwrite(merged_data_metabolism, "output_data/merged_data_metabolism.tsv", sep = "\t")

# Display the first few rows of the merged data
#head(merged_data_metabolism)

# Remove columns that start with "rel"
merged_data_metabolism <- merged_data_metabolism %>%
  select(-starts_with("rel"))

# Remove the prefix "reads_" from column names
names(merged_data_metabolism) <- gsub("^reads_", "", names(merged_data_metabolism))

# Save the merged data to a new file (optional)
#fwrite(merged_data_metabolism, "output_data/merged_data_metabolism_reads.tsv", sep = "\t")

#enzyme lowest level 4
#SEED_subsystem level 3
#description level 2
#SEED_subsystem_functional_category highest level 1 

## Creating tse_pathway for metabolic pathways  

# Use mapvalues to rename columns
colnames(merged_data_metabolism) <- mapvalues(colnames(merged_data_metabolism), 
                                              from = name_changes$SRA_accession, 
                                              to = name_changes$Date, 
                                              warn_missing = FALSE)

# Specify the new order for the first four columns
new_order_first4 <- c("SEED_subsystem_functional_category", "description", "SEED_subsystem", "enzyme")

# Get the order of the remaining columns from df2
remaining_cols_order <- colnames(assay(tse))

# Reorder columns in df1
merged_data_metabolism_reordered <- merged_data_metabolism %>%
  select(all_of(new_order_first4), all_of(remaining_cols_order))

colnames(merged_data_metabolism_reordered)[1:4] <- c("Kingdom", "Phylum", "Class", "Order")
rowData_pathway = as.data.frame(merged_data_metabolism_reordered[ ,1:4])

#rownames(rowData_pathway) = rowData_pathway$enzyme

counts_pathway = as.matrix(merged_data_metabolism_reordered[ ,5:36]) 
rownames(counts_pathway) = rownames(rowData_pathway)

colData_pathway = colData(tse)
#rownames(colData_pathway) = colnames(counts_pathway) 

#Condition:  
#colnames(counts_active) == rownames(colData_active)
#rownames(rowData_active) == rownames(counts_active) 

tse_pathway <- TreeSummarizedExperiment(assays = list(counts = counts_pathway ),
                                        colData = colData_pathway,
                                        rowData = rowData_pathway, 
                                        
)

################################################################################
################################################################################
## MultiAssayExperiment (MAE) to combine 3 TSE 
#tse all organisms 
#tse_bacteria is bacteria from the total community 
#tse_active which is metabolically active community 
#tse_pathway which contains data regarding metabolic pathways 
################################################################################
##Transforming tse to relabundance 

tse <- transformAssay(tse, method = "relabundance")
tse_bacteria <- transformAssay(tse_bacteria, method = "relabundance")
tse_active <- transformAssay(tse_active, method = "relabundance")
tse_pathway <- transformAssay(tse_pathway, method = "relabundance")

################################################################################

##Create hierarchy tree 

#getHierarchyTree(tse)
#tse <- mia::addHierarchyTree(tse)
#tse_bacteria <- mia::addHierarchyTree(tse_bacteria)
#There is no hierarchy in tse_active as there is only species 
#tse_active <- addHierarchyTree(tse_active)
#tse_pathway <- mia::addHierarchyTree(tse_pathway)

