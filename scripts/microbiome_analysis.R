suppressMessages({
  suppressWarnings({
    source(file = "scripts/01_install_load_packages.R")
  })
})

source(file = "scripts/02_import_clean_data.R")
#source(file = "scripts/03_feature_accumulation_curve.R")
#species function accumulation curve 
source(file = "scripts/04_exploration_quality_control_taxonomy.R")
source(file = "scripts/05_exploration_quality_control_function.R")
#source(file = "scripts/06_community_composition_taxonomy.R")
#source(file = "scripts/07_community_composition_function.R")
source(file = "scripts/08_diversity_taxonomy.R")

################################################################################
#ls("package:mia")

#Creating alternative experiment for esach taxonomic level
ranks = c('Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus','Species')
setTaxonomyRanks(ranks)

################################################################################
##Multi-omics: Cross-association
#This the only valid thing for our analysis out of meta-omics section 

setTaxonomyRanks( c('Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus','Species'))

# Agglomerate microbiome data at family level
mae[[1]] <- agglomerateByPrevalence(mae[[1]], rank = "Phylum")
# Does log10 transform for microbiome data
mae[[1]] <- transformAssay(mae[[1]], method = "log10", pseudocount = TRUE)

# Give unique names, so that we do not have problems when we are creating a plot
rownames(mae[[1]]) <- getTaxonomyLabels(mae[[1]])


mae[[2]] <- agglomerateByPrevalence(mae[[2]], rank = "Kingdom")
rownames(mae[[2]]) <- getTaxonomyLabels(mae[[2]])


# Cross correlates data sets
res <- testExperimentCrossCorrelation(
  mae,
  experiment1 = 1,
  experiment2 = 2,
  assay.type1 = "log10",
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

#res = as.matrix(res)
# Create a heatmap
p <- Heatmap(res$cor,
             # Print values to cells
             cell_fun = add_signif,
             heatmap_legend_param = list(
               title = "correlation", legend_height = unit(5, "cm")),
             column_names_rot = 45
)
p

