# Load the global variable 
mae <- readRDS("mae.rds")
tse = mae[[1]]
tse_bacteria = mae[[2]]
tse_active = mae[[3]]
tse_pathway = mae[[4]]
tse_enzyme = mae[[5]]
tse_gene = mae[[6]]

################################################################################
# Active Community 
tse_active_core <- agglomerateByRank(tse_active, rank ="Species")

##################################
##Core community 
core_active_species  = getPrevalence(
  tse_active_core, detection = 1, sort = TRUE, assay.type = "counts",
  as.relative = FALSE) %>% head(19)

#cat("There are 19 core active species in the microbiome")
#print(core_active_species)

core_active_species = names(core_active_species)
##################################

# Subset by feature
selected <- rowData(tse_active_core)$Species %in% core_active_species &
  !is.na(rowData(tse_active_core)$Species)
tse_active_core <- tse_active_core[selected, ]

################################################################################
tse_active_core <- transformAssay(tse_active_core,assay.type = "counts", method = "relabundance")
################################################################################
#Insert full taxonomy 
rowData_core_active = rowData(tse_active_core)

tax =read_excel("input_data//core_active_species_with_taxonomy.xlsx", sheet = "Sheet1")


rowData_core_active = merge(rowData_core_active,tax, by ="Species")

# Define the desired order of ranks
ranks <- c('Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species')

# Reorder the dataframe columns according to the 'ranks' vector
rowData_core_active <- rowData_core_active[, c(ranks, setdiff(names(rowData_core_active), ranks))]

rowData(tse_active_core) = rowData_core_active

################################################################################


#Generate a hierarchy tree on the fly
tse_active_core <- addHierarchyTree(tse_active_core)

tree <- rowTree(tse_active_core)

# Plot the tree with taxonomic labels and color

#png(filename="figures/tree.png" ,units = 'in',width=9, height=6, res=1000)
#ggtree(tree, layout = "rectangular", ladderize = TRUE, right = FALSE) +
#  geom_tiplab(hjust = 1,vjust = -0.5, size = 3, color = "black") 
#dev.off()
