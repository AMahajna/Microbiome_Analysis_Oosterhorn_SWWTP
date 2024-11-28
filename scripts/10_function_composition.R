# Load the global variable 
mae <- readRDS("mae.rds")
tse = mae[[1]]
tse_bacteria = mae[[2]]
tse_active = mae[[3]]
tse_pathway = mae[[4]]
tse_enzyme = mae[[5]]
tse_gene = mae[[6]]

#######################################################################
#Function 
# Agglomerate to phylum level
tse_pathway_cat <- agglomerateByRank(tse_pathway, rank = "Kingdom")

# Add clr-transformation on samples
tse_pathway_cat <- transformAssay(
  tse_pathway_cat, MARGIN = "samples", method = "clr", assay.type = "counts",
  pseudocount = 1)

# Add standardize-transformation on features (taxa)
tse_pathway_cat <- transformAssay(
  tse_pathway_cat, assay.type = "clr", MARGIN = "features",
  method = "standardize", name = "clr_z")

# Gets the assay table
mat <- assay(tse_pathway_cat, "clr_z")

# Create a heatmap with clusters

heatmap_func =pheatmap(mat ,fontsize_row = 7 ,fontsize_col = 8,
         scale = "none",
         clustering_distance_rows = "correlation",   
         clustering_distance_cols = "correlation",   
         clustering_method = "ward.D",
         annotation_legend = TRUE, heatmap_legend_param = list(title = "clr_z") )

png(filename="figures/heatmap_func_cat_clr_z.png" ,units = 'in',width=9, height=6, res=1000)
print(heatmap_func)
dev.off()

#png(filename="figures/heatmap_func_cat_clr_z_no_cluster.png" ,units = 'in',width=9, height=6, res=1000)
#pheatmap(mat, cluster_rows = FALSE,fontsize_row = 8 ,fontsize_col = 8,
#         cluster_cols = FALSE,annotation_legend = TRUE, heatmap_legend_param = list(title = "clr_z"))
#dev.off()