################################################################################
# Load the global variable 
mae <- readRDS("mae.rds")
tse = mae[[1]]
tse_bacteria = mae[[2]]
tse_active = mae[[3]]
tse_pathway = mae[[4]]
tse_enzyme = mae[[5]]
tse_gene = mae[[6]]


tse_phyla = agglomerateByRank(tse_bacteria, rank = "Phylum")
tse_func_cat = agglomerateByRank(tse_pathway, rank = "Kingdom")

################################################################################
#Creat mae_function
tse_function_list <- list(bacteriota = tse_phyla, functions = tse_func_cat)

# Combine into a MultiAssayExperiment object
mae_function <- MultiAssayExperiment::MultiAssayExperiment(experiments = tse_function_list)
################################################################################

mae_function[[1]] <- transformAssay(mae_function[[1]], method = "clr", pseudocount = TRUE)
mae_function[[2]] <- transformAssay(mae_function[[2]], method = "clr", pseudocount = TRUE)


# Cross correlates data sets
#Use Pearson between data is normally disributed 
# and we assume that an increase of a specific species will increase function linearly 
res <- testExperimentCrossCorrelation(
  mae_function,
  experiment1 = 1,
  experiment2 = 2,
  assay.type1 = "clr",
  assay.type2 = "relabundance",
  method = "pearson",
  test.signif = TRUE,
  p_adj_threshold = NULL,    # Add significance threshold
  cor_threshold = NULL,      # Add correlation threshold
  # Remove when mia is fixed
  mode = "matrix",
  sort = TRUE,
  show.warnings = FALSE)


# Function for marking significant correlations with "X"
add_signif <- function(j, i, x, y, width, height, fill, p_adj_threshold = 0.05, cor_threshold = 0.5) {
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
             row_names_gp = gpar(fontsize = 8),     # Font size for row names
             clustering_method_rows = "ward.D",
             clustering_method_columns = "ward.D"
)


png(filename="figures/heatmap_bacterial_phyla_clr_func_cat_relabundance_pearson_euclidean_wardd.png" ,units = 'in',width=9, height=6, res=1000)
print(p)
dev.off()
################################################################################


