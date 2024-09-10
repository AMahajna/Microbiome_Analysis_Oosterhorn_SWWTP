source(file = "scripts/install_load_packages.R")
################################################################################
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
################################################################################
##Microbiome analysis tse_bacteria 

#Creating alternative experiment for esach taxonomic level
ranks = c('Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus','Species')
setTaxonomyRanks(ranks)
for (r in ranks) {
  altExp(tse_bacteria,r) <- agglomerateByRank(tse_bacteria, r, agglomerate.tree = TRUE)
}

#relative abundance for the top-10 phylum over a log-scaled axis
#png(filename="figures/density_plot_bacterial_phylum.png" ,units = 'in',width=9, height=6, res=1000)
plotAbundanceDensity(altExp(tse_bacteria, "Phylum"), layout = "jitter", 
                     assay.type = "relabundance",
                     n = 10, point_size=2, point_shape=19, colour_by="Season",
                     point_alpha=0.5) + 
  scale_x_log10(label=scales::percent)
#dev.off()

#relative abundance for the top-10 species over a log-scaled axis
#png(filename="figures/density_plot_bacteria.png" ,units = 'in',width=9, height=6, res=1000)
plotAbundanceDensity(altExp(tse_bacteria, "Species"), layout = "density", 
                     assay.type = "relabundance",
                     n = 10, colour_by="Season", 
                     point_alpha=1/10) +
  scale_x_log10(label=scales::percent)
#dev.off()

################################################################################
# Core community 

core_phylum = getPrevalence(
  tse, rank = "Phylum",detection = 1, sort = TRUE, assay.type = "counts",
  as.relative = FALSE) %>% head(8)

core_species =getPrevalence(
  tse, rank = "Species",detection = 1, sort = TRUE, assay.type = "counts",
  as.relative = FALSE) %>% head(14)

core_phylum_bacteria = getPrevalence(
  tse_bacteria, rank = "Phylum",detection = 1, sort = TRUE, assay.type = "counts",
  as.relative = FALSE) %>% head(5)

core_species_bacteria = getPrevalence(
  tse_bacteria, rank = "Species",detection = 1, sort = TRUE, assay.type = "counts",
  as.relative = FALSE) %>% head(13)

core_class_bacteria = getPrevalence(
  tse_bacteria, rank = "Class",detection = 1, sort = TRUE, assay.type = "counts",
  as.relative = FALSE) %>% head(12)

mapTaxonomy(tse_bacteria, taxa = "s__flagellatus")
#Homonyms in Taxonomy 

################################################################################
#Prevalence scater plot 
#Creating alternative experiment for esach taxonomic level
ranks = c('Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus','Species')
setTaxonomyRanks(ranks)
for (r in ranks) {
  altExp(tse,r) <- agglomerateByRank(tse, r, agglomerate.tree = TRUE)
}

#Prevalence of Phylum in total community labeled by kingdom 
rowData(altExp(tse,"Phylum"))$prevalence <- 
  getPrevalence(altExp(tse,"Phylum"), detection = 1/100, 
                sort = FALSE,
                assay.type = "counts", as.relative = TRUE)

#png(filename="figures/prevalence.png" ,units = 'in',width=9, height=6, res=1000)
plotRowData(altExp(tse,"Phylum"), "prevalence", point_size=5,
            colour_by = "Kingdom")
#dev.off()

#Core community is made of 8 Phyla (5 bacteria and 3 Eukaryota)
sum(getPrevalence(altExp(tse,"Phylum"), detection = 1/100, sort = FALSE,
                  assay.type = "counts", as.relative = TRUE, prevalence = 1) == 1)
################################################################################
#Prevalence taxonomic tree

tse_bacteria <- tse[
  rowData(tse)$Kingdom %in% c("k__Bacteria"), ]

altExps(tse_bacteria) <- lapply(
  altExps(tse_bacteria), function(y){
    rowData(y)$prevalence <- getPrevalence(
      y, detection = 1/100,
      sort = FALSE,
      assay.type = "counts",
      as.relative = TRUE)
    return(y)
  })

top_phyla <- getTopFeatures(
  altExp(tse_bacteria,"Phylum"),
  method="sum",
  top=5L,
  assay.type="counts")
top_phyla_mean <- getTopFeatures(
  altExp(tse_bacteria,"Phylum"),
  method="mean",
  top=5L,
  assay.type="counts")

#Number of ranks can also be decreased
x <- unsplitByRanks(tse_bacteria, ranks = taxonomyRanks(tse_bacteria)[1:7])
x <- addHierarchyTree(x)

#png(filename="figures/prevalence_tree_bacteria.png" ,units = 'in',width=9, height=6, res=1000)
plotRowTree(
  x[rowData(x)$Phylum %in% top_phyla,],
  edge_colour_by = "Phylum",
  tip_colour_by = "prevalence",
  node_colour_by = "prevalence")
#dev.off()

################################################################################
#Quality Control 
tse_bacteria <- addPerCellQC(tse_bacteria)

#png(filename="figures/library_size.png" ,units = 'in',width=9, height=6, res=1000)
plotColData(tse_bacteria,"sum","Season", colour_by = "Season", point_size = 5, 
            shape_by ="Season") + 
  theme(axis.text.x = element_text(angle = 45, hjust=1))
#dev.off()
################################################################################
#Composition barplot for core community based on class level

# Getting top taxa on a Phylum level
tse_class_bacteria <- agglomerateByRank(tse_bacteria, rank ="Class")
tse_class_bacteria <- transformAssay(tse_class_bacteria, assay.type = "counts", method = "relabundance")
top_taxa_bacteria <- getTopFeatures(tse_class_bacteria, top = 12, assay.type = "relabundance")

# Renaming the "Phylum" rank to keep only top taxa and the rest to "Other"
class_renamed <- lapply(rowData(tse_class_bacteria)$Class, function(x){
  if (x %in% top_taxa_bacteria) {x} else {"Other"}
})
rowData(tse_class_bacteria)$Class <- as.character(class_renamed)

tse_class_bacteria_sub <- agglomerateByRank(tse_class_bacteria, rank= "Class")


plots <- plotAbundance(tse_class_bacteria_sub, rank = "Class", 
                       assay.type = "relabundance", features = "Date")

# Modify the legend of the first plot to be smaller
plots[[1]] <- plots[[1]] +
  theme(
    legend.key.size = unit(0.3, 'cm'),
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 10))

# Modify the legend of the second plot to be smaller
plots[[2]] <- plots[[2]] +
  theme(
    legend.key.height = unit(0.3, 'cm'),
    legend.key.width = unit(0.3, 'cm'),
    legend.text = element_text(size = 6),
    legend.title = element_text(size = 8),
    legend.direction = "vertical")


# Combine legends
legend <- wrap_plots(
  as_ggplot(get_legend(plots[[1]])),
  as_ggplot(get_legend(plots[[2]])),
  ncol = 1)

# Removes legends from the plots
plots[[1]] <- plots[[1]] + theme(legend.position = "none")
plots[[2]] <- plots[[2]] +
  theme(legend.position = "none", axis.title.x=element_blank())

# Combine plots
plot <- wrap_plots(plots[[2]], plots[[1]], ncol = 1, heights = c(2, 10))
# Combine the plot with the legend

#png(filename="figures/barplot_class.png" ,units = 'in',width=9, height=6, res=1000)
wrap_plots(plot, legend, nrow = 1, widths = c(2, 1))
#dev.off()


# Getting top taxa on a Phylum level
tse_phylum_bacteria <- agglomerateByRank(tse_bacteria, rank ="Phylum")
tse_phylum_bacteria <- transformAssay(tse_phylum_bacteria, assay.type = "counts", method = "relabundance")
top_taxa_bacteria_phylum <- getTopFeatures(tse_phylum_bacteria, top = 10, assay.type = "relabundance")

# Renaming the "Phylum" rank to keep only top taxa and the rest to "Other"
phylum_renamed <- lapply(rowData(tse_phylum_bacteria)$Phylum, function(x){
  if (x %in% top_taxa_bacteria_phylum) {x} else {"Other"}
})
rowData(tse_phylum_bacteria)$Phylum <- as.character(phylum_renamed)

tse_phylum_bacteria_sub <- agglomerateByRank(tse_phylum_bacteria, rank= "Phylum")


plots <- plotAbundance(tse_phylum_bacteria_sub, rank = "Phylum", 
                       assay.type = "relabundance", features = "Date")

# Modify the legend of the first plot to be smaller
plots[[1]] <- plots[[1]] +
  theme(
    legend.key.size = unit(0.3, 'cm'),
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 10))

# Modify the legend of the second plot to be smaller
plots[[2]] <- plots[[2]] +
  theme(
    legend.key.height = unit(0.3, 'cm'),
    legend.key.width = unit(0.3, 'cm'),
    legend.text = element_text(size = 6),
    legend.title = element_text(size = 8),
    legend.direction = "vertical")


# Combine legends
legend <- wrap_plots(
  as_ggplot(get_legend(plots[[1]])),
  as_ggplot(get_legend(plots[[2]])),
  ncol = 1)

# Removes legends from the plots
plots[[1]] <- plots[[1]] + theme(legend.position = "none")
plots[[2]] <- plots[[2]] +
  theme(legend.position = "none", axis.title.x=element_blank())

# Combine plots
plot <- wrap_plots(plots[[2]], plots[[1]], ncol = 1, heights = c(2, 10))
# Combine the plot with the legend

#png(filename="figures/barplot_phylum.png" ,units = 'in',width=9, height=6, res=1000)
wrap_plots(plot, legend, nrow = 1, widths = c(2, 1))
#dev.off()





################################################################################
#Diversity  
#Visualizing significance in group-wise comparisons

index <- "shannon"
group_var <- "Season"

tse_bacteria <- estimateDiversity(tse_bacteria, assay.type = "counts", index = "shannon")


# Calculate p values
pvals <- pairwise.wilcox.test(
  tse_bacteria[[index]], tse_bacteria[[group_var]], p.adjust.method = "fdr")
# Put them to data.frame format
pvals <- pvals[["p.value"]] |>
  as.data.frame()
varname <- "group1"
pvals[[varname]] <- rownames(pvals)
# To long format
pvals <- reshape(
  pvals,
  direction = "long",
  varying = colnames(pvals)[ !colnames(pvals) %in% varname ],
  times = colnames(pvals)[ !colnames(pvals) %in% varname ],
  v.names = "p",
  timevar = "group2",
  idvar = "group1"
) |>
  na.omit()
# Add y-axis position
pvals[["y.position"]] <- apply(pvals, 1, function(x){
  temp1 <- tse[[index]][ tse[[group_var]] == x[["group1"]] ]
  temp2 <- tse[[index]][ tse[[group_var]] == x[["group2"]] ]
  temp <- max( c(temp1, temp2) )
  return(temp)
})
pvals[["y.position"]] <- max(pvals[["y.position"]]) +
  order(pvals[["y.position"]]) * 0.2
# Round values
pvals[["p"]] <- round(pvals[["p"]], 3)

# Create a boxplot
p <- plotColData(
  tse_bacteria, x = group_var, y = index, show_violin = TRUE, shape_by = "Season",
  point_size = 5, color_by = "Season") +
  theme(text = element_text(size = 10)) +
  stat_pvalue_manual(pvals) 

#png(filename="figures/diversity_plot.png" ,units = 'in',width=9, height=6, res=1000)
p
#dev.off()
################################################################################
## Plot species accumulation curve

ranks = c('Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus','Species')
for (r in ranks) {
  altExp(tse,r) <- agglomerateByRank(tse, r, agglomerate.tree = TRUE)
}

#Calculate species accumulation
species_accum <- specaccum(t(assay(tse)), method = 'random')
genera_accum <- specaccum(t(assay(altExp(tse,"Genus"))), method = 'random')
family_accum <- specaccum(t(assay(altExp(tse,"Family"))), method = 'random')

bacterial_species_accum <- specaccum(t(assay(tse_bacteria)), method = 'random')

active_species_accum <- specaccum(t(assay(tse_active)), method = 'random')

ranks = c('Kingdom', 'Phylum', 'Class', 'Order')
for (r in ranks) {
  altExp(tse_pathway,r) <- agglomerateByRank(tse_pathway, r, agglomerate.tree = TRUE)
}
enzyme_accum <- specaccum(t(assay(tse_pathway)), method = 'random')
pathway_accum <- specaccum(t(assay(altExp(tse_pathway,"Class"))), method = 'random')

#png(filename="figures/accum_curve.png" ,units = 'in',width=9, height=6, res=1000)
# Plot species accumulation curve
plot(species_accum, xlab = "Number of Samples", ylab = "Number of Features",
     col = "skyblue", lwd = 3,  xlim = c(0, 35), ci = 0, ylim = c(0,1500))

# Add the second curve using the lines() function
lines(
  genera_accum,
  col = "darkgreen",      # Different color for the second curve
  lwd = 3, # Line width
  ci = 0
)

# Add the second curve using the lines() function
lines(
  family_accum,
  col = "darkorange",      # Different color for the second curve
  lwd = 3, # Line width
  ci = 0
)

# Add the second curve using the lines() function
lines(
  bacterial_species_accum,
  col = "slategrey",      # Different color for the second curve
  lwd = 3, # Line width
  ci = 0
)

# Add the second curve using the lines() function
lines(
  pathway_accum,
  col = "firebrick",      # Different color for the second curve
  lwd = 3, # Line width
  ci = 0
)


# Add a legend to the plot
legend(
  "topleft",
  legend = c("Metabolic Pathways","Bacterial species","Species","Genera","Families"),  # Labels for each curve
  col = c("firebrick","slategrey","skyblue","darkgreen","darkorange"),         # Corresponding colors
  lwd = 3,
  cex = 1,# Line width in the legend,
  bty = "n"
  )
#dev.off()


################################################################################

#set taxonomy ranks 
colnames(rowData(tse_pathway))<- c("SEED_subsystem_functional_category", 
                                   "description", 
                                   "SEED_subsystem",
                                   "enzyme")
setTaxonomyRanks(colnames(rowData(tse_pathway)))
ranks =getTaxonomyRanks()

for (r in ranks) {
  altExp(tse_pathway,r) <- agglomerateByRank(tse_pathway, r, agglomerate.tree = TRUE)
}
#length(unique(rowData(tse_pathway)[ ,1]))

tse_functional_category <- agglomerateByRank(tse_pathway,rank = "SEED_subsystem_functional_category")

# Add clr-transformation on samples
tse_functional_category <- transformAssay(tse_functional_category, MARGIN = "samples", method = "clr", assay.type = "counts", pseudocount=1)

# Add standardize-transformation on features (taxa)
tse_functional_category <- transformAssay(tse_functional_category, assay.type = "clr",
                             MARGIN = "features", 
                             method = "standardize", name = "clr_z")

# Gets the assay table
mat <- assay(tse_functional_category, "clr_z")


#png(filename="figures/heatmap_functional_category.png" ,units = 'in',width=9, height=6, res=1000)
# Creates the heatmap
pheatmap(mat)
#dev.off()



core_description = getPrevalentFeatures(altExp(tse_pathway,"description"), detection = 0, prevalence = 99/100,
                                        rank = "description", sort = TRUE)

#check
sum((getPrevalence(altExp(tse_pathway,"description"), detection = 1/100, prevalence =1, 
                   sort = TRUE, assay.type = "counts",
                   as.relative = TRUE))== 1)


core_subsystem = getPrevalentFeatures(altExp(tse_pathway,"SEED_subsystem"), detection = 0, prevalence = 99/100,
                                        rank = "SEED_subsystem", sort = TRUE)


core_enzyme = getPrevalentFeatures(tse_pathway, detection = 0, prevalence = 0.999,
                                      rank = "enzyme", sort = TRUE)

sum((getPrevalence(tse_pathway, tank ="enzyme", detection = 1/100, prevalence =1,  
                   sort = TRUE, assay.type = "counts",
                   as.relative = TRUE))== 1)

tse_enzyme <- agglomerateByRank(tse_pathway,rank = "enzyme")

tse_enzyme  <- transformAssay(tse_enzyme , method = "relabundance")

# Add clr-transformation on samples
tse_enzyme <- transformAssay(tse_enzyme, MARGIN = "samples", method = "clr", assay.type = "counts", pseudocount=1)

# Add standardize-transformation on features (taxa)
tse_enzyme <- transformAssay(tse_enzyme, assay.type = "clr",
                                          MARGIN = "features", 
                                          method = "standardize", name = "clr_z")


  
tse_enzyme_subset <- tse_enzyme[core_enzyme, ]

# Add clr-transformation
tse_enzyme_subset <- transformAssay(tse_enzyme_subset, method = "clr",
                                    MARGIN="samples",
                                    assay.type = "counts")
# Does standardize-transformation
tse_enzyme_subset  <- transformAssay(tse_enzyme_subset , assay.type = "clr",
                                    MARGIN = "features", 
                                    method = "standardize", name = "clr_z")

# Gets the assay table

mat <- assay(tse_enzyme_subset, "clr_z")

png(filename="figures/heatmap_core_enzyme.png" ,units = 'in',width=9, height=6, res=1000)
# Creates the heatmap
pheatmap(mat)
dev.off()

# Hierarchical clustering
taxa_hclust <- hclust(dist(mat), method = "complete")

# Creates a phylogenetic tree
taxa_tree <- as.phylo(taxa_hclust)

# Plot taxa tree
taxa_tree <- ggtree(taxa_tree) + 
  theme(plot.margin=margin(0,0,0,0)) # removes margins

# Get order of taxa in plot
taxa_ordered <- get_taxa_name(taxa_tree)

# to view the tree, run
taxa_tree
# Creates clusters
taxa_clusters <- cutree(tree = taxa_hclust, k = 3)

# Converts into data frame
taxa_clusters <- data.frame(clusters = taxa_clusters)
taxa_clusters$clusters <- factor(taxa_clusters$clusters)

# Order data so that it's same as in phylo tree
taxa_clusters <- taxa_clusters[taxa_ordered, , drop = FALSE] 

# Prints taxa and their clusters
taxa_clusters

# Adds information to rowData
rowData(tse_enzyme_subset)$clusters <- taxa_clusters[order(match(rownames(taxa_clusters), rownames(tse_enzyme_subset))), ]

# Prints taxa and their clusters
rowData(tse_enzyme_subset)$clusters

# Hierarchical clustering
sample_hclust <- hclust(dist(t(mat)), method = "complete")

# Creates a phylogenetic tree
sample_tree <- as.phylo(sample_hclust)

# Plot sample tree
sample_tree <- ggtree(sample_tree) + layout_dendrogram() + 
  theme(plot.margin=margin(0,0,0,0)) # removes margins

# Get order of samples in plot
samples_ordered <- rev(get_taxa_name(sample_tree))

# to view the tree, run
sample_tree

# Creates clusters
sample_clusters <- factor(cutree(tree = sample_hclust, k = 3))

# Converts into data frame
sample_data <- data.frame(clusters = sample_clusters)

# Order data so that it's same as in phylo tree
sample_data <- sample_data[samples_ordered, , drop = FALSE] 



# Order data based on 
tse_enzyme_subset <- tse_enzyme_subset[ , rownames(sample_data)]

# Add sample type data
sample_data$sample_types <- unfactor(colData(tse_enzyme_subset)$Season)

############
# Use left_join to merge df1 and df2 by the shared column 'Date'
sample_data$Date = rownames(sample_data)
df2 = as.data.frame(colData(tse_enzyme_subset)[ , c("Date","Season")])
sample_data$Date = as.Date(sample_data$Date)
sample_data <- sample_data %>%
  left_join(df2, by = "Date")
rownames(sample_data) = sample_data$Date
sample_data$Date = NULL
############
png(filename="figures/heatmap_core_enzyme_clusters.png" ,units = 'in',width=9, height=6, res=1000)

breaks <- seq(-ceiling(max(abs(mat))), ceiling(max(abs(mat))), 
              length.out = ifelse( max(abs(mat))>5, 2*ceiling(max(abs(mat))), 10 ) )
colors <- colorRampPalette(c("darkblue", "blue", "white", "red", "darkred"))(length(breaks)-1)

pheatmap(mat, annotation_row = taxa_clusters, 
         annotation_col = sample_data,
         breaks = breaks,
         color = colors, border_color = grey)
dev.off()

################################################################################
#ls("package:mia")
#co-abundant groups as CAGs, which are clusters of taxa that co-vary across samples 
#getUniqueTaxa(tse_bacteria, rank = "Phylum")
