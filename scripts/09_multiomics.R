#How does the microbial composition contribute to removal efficiency 
#How does the functional composition contribute to the removal efficiency 

#How does taxonomic composition contribute to functional composition 

# Load the global variable 
mae <- readRDS("mae.rds")
tse = mae[[1]]
tse_bacteria = mae[[2]]
tse_active = mae[[3]]
tse_pathway = mae[[4]]
tse_enzyme = mae[[5]]
tse_gene = mae[[6]]

################################################################################
# Stress Response 
selected_stress <- rowData(tse_pathway)$Kingdom %in% c("Stress Response") &
  !is.na(rowData(tse_pathway)$Kingdom)
tse_stress <- tse_pathway[selected_stress, ]

# # Active Community 
tse_active_barplot <- agglomerateByRank(tse_active, rank ="Species")

##########################
##Core community 
core_active_species = getPrevalence(
  tse_active, detection = 1, sort = TRUE, assay.type = "counts",
  as.relative = FALSE) %>% head(19)

#cat("There are 19 core active species in the microbiome")
#print(core_active_species)

core_active_species = names(core_active_species)
###########################

# Renaming the "Phylum" rank to keep only top taxa and the rest to "Other"
species_renamed <- lapply(rowData(tse_active_barplot)$Species, function(x){
  if (x %in% core_active_species) {x} else {"Other"}
})
rowData(tse_active_barplot)$Species <- as.character(species_renamed)
################################################################################
#Functional categories 
tse_func_cat <- agglomerateByRank(tse_pathway, rank ="Kingdom")
################################################################################

mae_active_list <- list(core_active = tse_active_barplot, func_cat = tse_func_cat, stress_response = tse_stress)

# Combine into a MultiAssayExperiment object
mae_active <- MultiAssayExperiment::MultiAssayExperiment(experiments = mae_active_list)


# Cross correlates data sets
res <- testExperimentCrossCorrelation(
  mae_active,
  experiment1 = 1,
  experiment2 = 2,
  assay.type1 = "relabundance",
  assay.type2 = "relabundance",
  method = "spearman",
  test.signif = TRUE,
  p_adj_threshold = NULL,
  cor_threshold = NULL,
  # Remove when mia is fixed
  mode = "matrix",
  sort = TRUE,
  show.warnings = FALSE)


# Function for marking significant correlations with "X"
add_signif <- function(j, i, x, y, width, height, fill) {
  # If the p-value is under threshold
  if( !is.na(res$p_adj[i, j]) & res$p_adj[i, j] < 0.05 ){
    # Print "X"
    grid.shadowtext(
      sprintf("%s", "X"), x, y, gp = gpar(fontsize = 8, col = "#f5f5f5"))
  }
}

# Create a heatmap
p <- Heatmap(res$cor,
             # Print values to cells
             cell_fun = add_signif,
             heatmap_legend_param = list(
               title = "correlation", legend_height = unit(5, "cm")),
             column_names_rot = 45
)
p





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
selected <- rowData(tse_pathway)$Phylum %in% c("Potassium metabolism") &
  !is.na(rowData(tse_pathway)$Phylum)
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

