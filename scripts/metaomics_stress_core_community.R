# Load the global variable 
mae <- readRDS("mae.rds")
tse = mae[[1]]
tse_bacteria = mae[[2]]
tse_active = mae[[3]]
tse_pathway = mae[[4]]
tse_enzyme = mae[[5]]
tse_gene = mae[[6]]

################################################################################
#All subsystems classes of stress 
# Stress Response 
selected_stress <- rowData(tse_pathway)$Kingdom %in% c("Stress Response") &
  !is.na(rowData(tse_pathway)$Kingdom)
tse_stress <- tse_pathway[selected_stress, ]

tse_stress_subsystems = agglomerateByRank(tse_stress, rank = "Class")

selected_stress_stress <- rowData(tse_pathway)$Phylum %in% c("Stress Response") &
  !is.na(rowData(tse_pathway)$Phylum)
tse_stress_stress <- tse_pathway[selected_stress_stress, ]

tse_stress_stress = agglomerateByRank(tse_stress_stress, rank = "Order")
#write_xlsx(as.data.frame(rowData(tse_stress_stress)), path = "output_data//stress_response_phylum.xlsx")

################################################################################
# Active Community 
tse_active_core <- agglomerateByRank(tse_active, rank ="Species")
##################################
##Core community 
core_active_species = getPrevalence(
  tse_active, detection = 1, sort = TRUE, assay.type = "counts",
  as.relative = FALSE) %>% head(19)

#cat("There are 19 core active species in the microbiome")
#print(core_active_species)

core_active_species = names(core_active_species)
##################################
# Renaming the "Phylum" rank to keep only top taxa and the rest to "Other"
species_renamed <- lapply(rowData(tse_active_core)$Species, function(x){
  if (x %in% core_active_species) {x} else {"Other"}
})
rowData(tse_active_core)$Species <- as.character(species_renamed)
tse_active_core = agglomerateByRank(tse_active_core, rank = "Species")
################################################################################
#Creat mae_stress 
tse_stress_list <- list(microbiota = tse_active_core, stress = tse_stress_subsystems)

# Combine into a MultiAssayExperiment object
mae_stress <- MultiAssayExperiment::MultiAssayExperiment(experiments = tse_stress_list)
  
################################################################################
# CLR transform mitigates compositional effects and is widely used in microbiome 
#studies. CLR-transformed data better represent underlying biological relationships.
mae_stress[[1]] <- transformAssay(mae_stress[[1]], method = "clr", pseudocount = TRUE)
mae_stress[[2]] <- transformAssay(mae_stress[[2]], method = "clr", pseudocount = TRUE)


# Cross correlates data sets

#Spearman correlation is a rank-based method that captures monotonic relationships, 
#making it more robust for microbiome data, which often do not meet normality assumptions.

# Cross correlates data sets
res <- testExperimentCrossCorrelation(
  mae_stress,
  experiment1 = 1,
  experiment2 = 2,
  assay.type1 = "clr",
  assay.type2 = "clr",
  method = "spearman",
  test.signif = TRUE,
  p_adj_threshold = NULL,    # Add significance threshold
  cor_threshold = NULL,      # Add correlation threshold
  # Remove when mia is fixed
  mode = "matrix",
  sort = TRUE,
  show.warnings = FALSE)


# Function for marking significant correlations with "X"
add_signif <- function(j, i, x, y, width, height, fill, p_adj_threshold = 0.05, cor_threshold = 0.3) {
  # Check if the p-value is under threshold and the correlation value exceeds the threshold
  if( !is.na(res$p_adj[i, j]) & res$p_adj[i, j] < p_adj_threshold & !is.na(res$cor[i, j]) & abs(res$cor[i, j]) > cor_threshold ){
    # Print "X" if both conditions are met
    grid.shadowtext(
      sprintf("%s", "X"), x, y, gp = gpar(fontsize = 8, col = "#f5f5f5"))
  }
}

#res = as.matrix(res)
# Create a heatmap
p <- Heatmap(res$cor,
             # Print values to cells
             cell_fun = add_signif,
             heatmap_legend_param = list(
               title = "correlation", legend_height = unit(5, "cm")),
             column_names_rot = -45,
             column_names_gp = gpar(fontsize = 8),  # Font size for column names
             row_names_gp = gpar(fontsize = 8)     # Font size for row names
)
             

#png(filename="figures/heatmap_stress_subsystems_core_active_genes_clr_spearman.png" ,units = 'in',width=9, height=6, res=1000)
p
#dev.off()
