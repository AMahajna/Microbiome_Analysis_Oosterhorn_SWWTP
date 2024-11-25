##Reading input data

#read biom data using bioformat package 
biom_data = biomformat::read_biom("input_data/Galaxy1783-[Kraken-biom_output_file].biom1")

#convert biom to TSE usingmia package 
tse <- makeTreeSEFromBiom(biom_data)

#removal_efficiency respective to sample to be added to column data in TSE
removal_efficiency  = read_excel("input_data/Removal_Efficiency.xlsx")

#Relevant process data to be added to column data in TSE
process_data <- read_csv(file = "input_data/Weekly_Data.csv", show_col_types = FALSE)

################################################################################
################################################################################
##preprocess taxonomic data from Kraken2   

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

#Preprocess process and environmental data
process_data_normalized <- decostand(process_data[ ,1:13], method = "standardize")
#summary(process_data_normalized)

col_data = DataFrame(cbind(col_data_new,process_data_normalized,(removal_efficiency[ ,3:6])/100))

# Assign the modified colData back to the SummarizedExperiment object
colData(tse) <- col_data

tse <- transformAssay(tse, method = "relabundance")

#Generate a hierarchy tree on the fly
tse <- addHierarchyTree(tse)

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

# Replace NA with 0 in the merged data frame
merged_data_org[is.na(merged_data_org)] <- 0

#Check if transformation was done correctly
#sum of rel_abundance columns if 100
columns_sum_org = colSums(merged_data_org[, -1], na.rm = TRUE)

# Save the merged data to a new file (optional)
#fwrite(merged_data_org, "output_data/merged_data_org.tsv", sep = "\t")

# Remove columns that start with "rel"
merged_data_org_reads <- merged_data_org %>%
  select(-starts_with("rel"))

# Remove the prefix "reads_" from column names
names(merged_data_org_reads) <- gsub("^reads_", "", names(merged_data_org_reads))

# Save the merged data to a new file (optional)
#fwrite(merged_data_org, "output_data/merged_data_org_reads.tsv", sep = "\t")

#First, change column names in merged_data_org_reads
name_changes <- as.data.frame(colData(tse))[ , c("SRA_accession", "Date")]
name_changes[] <- lapply(name_changes, as.character)

# Use mapvalues to rename columns
colnames(merged_data_org_reads) <- mapvalues(colnames(merged_data_org_reads), 
                                       from = name_changes$SRA_accession, 
                                       to = name_changes$Date, 
                                       warn_missing = FALSE)

# Rename the column
colnames(merged_data_org_reads)[colnames(merged_data_org_reads) == "species"] <- "Species"

# Specify the new order for the first four columns
new_order_first <- c("Species")

# Get the order of the remaining columns from df2
remaining_cols_order <- colnames(assay(tse))

# Reorder columns in df1
merged_data_org_reordered <- merged_data_org_reads %>%
  select(all_of(new_order_first), all_of(remaining_cols_order))

#rowData structure
rowData_active = as.data.frame(merged_data_org_reordered[ ,1])
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
tse_active <- transformAssay(tse_active, method = "relabundance")


################################################################################
##gene: reading and cleaning 

# Define the path to your folder containing TSV files
folder_path_gene <- "input_data/func_result"

# List all TSV files in the folder
file_list_gene <- list.files(folder_path_gene, pattern = "\\.tsv$", full.names = TRUE)

# Function to read and preprocess each file
read_and_preprocess <- function(file) {
  # Extract the file name without extension
  file_name <- tools::file_path_sans_ext(basename(file))
  
  # Read the TSV file
  df <- fread(file, header = FALSE, sep = "\t")
  
  # Rename columns
  colnames(df) <- c(paste0("rel_abundance_", file_name),
                    paste0("reads_", file_name),
                    "Species")
  
  # Ensure 'Species' column is character type for consistency
  df[, Species := as.character(Species)]
  
  # Check for duplicates in the 'gene' column
  if (anyDuplicated(df$Species) > 0) {
    df <- df[, lapply(.SD, sum), by = Species, .SDcols = c(paste0("rel_abundance_", file_name), paste0("reads_", file_name))]
  }
  
  # Return the data.table
  return(df)
}

# Read and preprocess all files
data_list_gene <- lapply(file_list_gene, read_and_preprocess)

# Merge all data frames by taxonomy column
merged_data_gene <- Reduce(function(x, y) merge(x, y, by = "Species", all = TRUE), data_list_gene)

# Replace NA with 0 in the merged data frame
merged_data_gene[is.na(merged_data_gene)] <- 0

# unify 
merged_data_gene <- merged_data_gene %>%
  mutate(Species = case_when(
    Species %in% c("hypothetical protein", 
                   "hypothetical protein, partial", 
                   "hypothetical protein ") ~ "hypothetical protein", 
    TRUE ~ Species
  ))


merged_data_gene <- merged_data_gene %>%
  group_by(Species) %>%
  summarise(across(everything(), sum, na.rm = FALSE), .groups = 'drop')



#Check if transformation was done correctly
#sum of rel_abundance columns if 100
columns_sum_gene = colSums(merged_data_gene[, -1], na.rm = TRUE)

# Save the merged data to a new file (optional)
#fwrite(merged_data_gene, "output_data/merged_data_gene.tsv", sep = "\t")

# Remove columns that start with "rel"
merged_data_gene_reads <- merged_data_gene %>%
  select(-starts_with("rel"))

# Remove the prefix "reads_" from column names
names(merged_data_gene_reads) <- gsub("^reads_", "", names(merged_data_gene_reads))

# Save the merged data to a new file (optional)
#fwrite(merged_data_gene, "output_data/merged_data_gene_reads.tsv", sep = "\t")

#First, change column names in merged_data_gene_reads
name_changes <- as.data.frame(colData(tse))[ , c("SRA_accession", "Date")]
name_changes[] <- lapply(name_changes, as.character)

# Use mapvalues to rename columns
colnames(merged_data_gene_reads) <- mapvalues(colnames(merged_data_gene_reads), 
                                              from = name_changes$SRA_accession, 
                                              to = name_changes$Date, 
                                              warn_missing = FALSE)



# Specify the new order for the first four columns
new_order_first <- c("Species")

# Get the order of the remaining columns from df2
remaining_cols_order <- colnames(assay(tse))

# Reorder columns in df1
merged_data_gene_reordered <- merged_data_gene_reads %>%
  select(all_of(new_order_first), all_of(remaining_cols_order))

#rowData structure
rowData_gene = as.data.frame(merged_data_gene_reordered[ ,1])
rownames(rowData_gene) = rowData_gene$Species

counts_gene = as.matrix(merged_data_gene_reordered[ ,2:33]) 
rownames(counts_gene) = rowData_gene$Species

#colData_gene = as.data.frame(colData(tse))
colData_gene = colData(tse)
#rownames(colData_gene) = colnames(counts_gene) 

#Condition:  
#colnames(counts_gene) == rownames(colData_gene)
#rownames(rowData_gene) == rownames(counts_gene) 

tse_gene <- TreeSummarizedExperiment(assays = list(counts = counts_gene ),
                                     colData = colData_gene,
                                     rowData = rowData_gene, 
                                     
)
tse_gene <- transformAssay(tse_gene, method = "relabundance")

################################################################################
##pathway: reading and cleaning 

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


columns_sum_metabolism = colSums(merged_data_metabolism[, -c(1,2,3,4)], na.rm = TRUE)

merged_data_metabolism<- {
  # Directly access the columns by index without using 'within'
  
  # Check if column 4 (Order) contains "Retron-type reverse transcriptase"
  merged_data_metabolism[[1]] <- ifelse(merged_data_metabolism[[1]] == "Retron-type reverse transcriptase", 
                                        "NO HIERARCHY", 
                                        merged_data_metabolism[[1]])
  
  # Check column 1 (Kingdom) and update based on conditions
  merged_data_metabolism[[3]] <- ifelse(
    merged_data_metabolism[[3]] == "" & (merged_data_metabolism[[1]] == "" | merged_data_metabolism[[1]] == " " | merged_data_metabolism[[1]] == "NO HIERARCHY"), 
    "NO HIERARCHY", 
    ifelse(
      merged_data_metabolism[[3]] == "" & !(merged_data_metabolism[[1]] == "" | merged_data_metabolism[[1]] == " " | merged_data_metabolism[[1]] == "NO HIERARCHY"), 
      merged_data_metabolism[[4]], 
      merged_data_metabolism[[3]]
    )
  )
  
  # Ensure consistency for rows where Kingdom is set to "NO HIERARCHY"
  merged_data_metabolism[[4]] <- ifelse(merged_data_metabolism[[3]] == "NO HIERARCHY", "NO HIERARCHY", merged_data_metabolism[[4]])  # Phylum
  merged_data_metabolism[[2]] <- ifelse(merged_data_metabolism[[3]] == "NO HIERARCHY", "NO HIERARCHY", merged_data_metabolism[[2]])  # Class
  merged_data_metabolism[[1]] <- ifelse(merged_data_metabolism[[3]] == "NO HIERARCHY", "NO HIERARCHY", merged_data_metabolism[[1]])  # Order
  
  # Return the modified data frame
  merged_data_metabolism
}

#Replace NA with  in the merged data frame
merged_data_metabolism[is.na(merged_data_metabolism)] <- 0


merged_data_metabolism <- merged_data_metabolism %>%
  group_by(across(1:4)) %>%
  summarise(across(everything(), sum, na.rm = FALSE), .groups = 'drop')



#Check if transformation was done correctly
#sum of rel_abundance columns if 100
columns_sum_metabolism = colSums(merged_data_metabolism[, -c(1,2,3,4)], na.rm = TRUE)
#columns_sum_metabolism

# Save the merged data to a new file (optional)
#fwrite(merged_data_metabolism, "output_data/merged_data_metabolism.tsv", sep = "\t")

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

# Get the order of the remaining columns 
remaining_cols_order <- colnames(assay(tse))

# Reorder columns 
merged_data_metabolism_reordered <- merged_data_metabolism %>%
  select(all_of(new_order_first4), all_of(remaining_cols_order))

colnames(merged_data_metabolism_reordered)[1:4] <- c("Kingdom", "Phylum", "Class", "Order")
rowData_pathway = as.data.frame(merged_data_metabolism_reordered[ ,1:4])

counts_pathway = as.matrix(merged_data_metabolism_reordered[ ,5:36]) 
rownames(counts_pathway) = rownames(rowData_pathway)

colData_pathway = colData(tse)



tse_pathway <- TreeSummarizedExperiment(assays = list(counts = counts_pathway ),
                                        colData = colData_pathway,
                                        rowData = rowData_pathway, 
                                        
)
tse_pathway <- transformAssay(tse_pathway, method = "relabundance")

################################################################################

tse_enzyme


assay_enzyme <- merged_data_metabolism_reordered[ ,4:36] %>%
  group_by(Order) %>%
  summarise(across(everything(), sum, na.rm = TRUE))

colData_enzyme= colData(tse)
rowData_enzyme = as.data.frame(assay_enzyme[ ,1])
rownames(rowData_enzyme) <- rowData_enzyme[, 1]  
counts_enzyme = assay_enzyme[ ,2:33]
rownames(counts_enzyme)=rowData_enzyme[, 1]

tse_enzymes <- TreeSummarizedExperiment(assays = list(counts = counts_enzyme),
                                        colData = colData_enzyme,
                                        rowData = rowData_enzyme, 
                                        
)

tse_enzymes <- transformAssay(tse_enzymes, method = "relabundance")

################################################################################
################################################################################
##Create Multi-assay experiment mae  
# Create a list of TSE objects

tse_list <- list(microbiota = tse, bacteriota  = tse_bacteria, active = tse_active, 
                 functions = tse_pathway, enzymes = tse_enzyme, genes = tse_gene)

# Combine into a MultiAssayExperiment object
mae <- MultiAssayExperiment::MultiAssayExperiment(experiments = tse_list)

################################################################################
#Saving global varibales 

# Save global variable MAE (Multi Assay Experiment)
saveRDS(mae, file = "mae.rds")


