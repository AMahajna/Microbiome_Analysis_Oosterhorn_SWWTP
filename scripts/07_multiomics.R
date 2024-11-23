#How does the microbial composition contribute to removal efficiency 
#How does the functional composition contribute to the removal efficiency 

#How does taxonomic composition contribute to functional composition 

# Load the global variable 
mae <- readRDS("mae.rds")
#Multi-Assay Experiment for Cross-Association Analysis 
mae_caa <- mae 
tse = mae[[1]]
tse_bacteria = mae[[2]]
tse_active = mae[[3]]
tse_pathway = mae[[4]]
tse_enzyme = mae[[5]]

unique_functions <- sapply(rowData(tse_function), function(x) length(unique(x)))

#core community from active species contributing to core function 

core_species_active = getPrevalence(
  tse_active, rank = "Species",detection = 1, sort = TRUE, assay.type = "counts",
  as.relative = FALSE) %>% head(19)

core_enzyme = getPrevalence(
  tse_pathway, rank = "Order",detection = 1, sort = TRUE, assay.type = "counts",
  as.relative = FALSE) %>% head(19)

core_pathway = as.data.frame(getPrevalence(
  tse_pathway, rank = "Class",detection = 1, sort = TRUE, assay.type = "counts",
  as.relative = FALSE) %>% head(34))

core_des = as.data.frame(getPrevalence(
  tse_pathway, rank = "Phylum",detection = 1, sort = TRUE, assay.type = "counts",
  as.relative = FALSE) %>% head(50))

core_cat = as.data.frame(getPrevalence(
  tse_pathway, rank = "Kingdom",detection = 1, sort = TRUE, assay.type = "counts",
  as.relative = FALSE) %>% head(50))

# Subset by feature
selected <- rowData(tse_pathway)$Kingdom %in% c("Stress Response", "Sulfur Metabolism") &
  !is.na(rowData(tse_pathway)$Kingdom)
tse_sub <- tse_pathway[selected, ]



# Getting top taxa on a species level - active from samsa2 
tse_pathway_core <- agglomerateByRank(tse_pathway, rank ="Class")

top_taxa_active <- getTopFeatures(tse_active_species, top = 19, assay.type = "relabundance")

# Renaming the "Phylum" rank to keep only top taxa and the rest to "Other"
species_renamed <- lapply(rowData(tse_active_species)$Species, function(x){
  if (x %in% top_taxa_active) {x} else {"Other"}
})
rowData(tse_active_species)$Species <- as.character(species_renamed)

################################################################################
#cross-Association Analysis 
mae_caa[[1]] <- agglomerateByPrevalence(mae[[1]], rank = "Family", na.rm = TRUE)

