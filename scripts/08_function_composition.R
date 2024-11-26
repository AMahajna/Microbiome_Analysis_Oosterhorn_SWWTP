################################################################################
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

#png(filename="figures/heatmap_func_cat_clr_z.png" ,units = 'in',width=9, height=6, res=1000)
pheatmap(mat ,fontsize_row = 7 ,fontsize_col = 8,
         scale = "none",
         clustering_distance_rows = "correlation",   
         clustering_distance_cols = "correlation",   
         clustering_method = "ward.D",
         annotation_legend = TRUE, heatmap_legend_param = list(title = "clr_z") )
#dev.off()

#png(filename="figures/heatmap_func_cat_clr_z_no_cluster.png" ,units = 'in',width=9, height=6, res=1000)
#pheatmap(mat, cluster_rows = FALSE,fontsize_row = 8 ,fontsize_col = 8,
#         cluster_cols = FALSE,annotation_legend = TRUE, heatmap_legend_param = list(title = "clr_z"))
#dev.off()
#######################################################################
################################################################################

# Stress Response 
selected_stress <- rowData(tse_pathway)$Kingdom %in% c("Stress Response") &
  !is.na(rowData(tse_pathway)$Kingdom)
tse_stress <- tse_pathway[selected_stress, ]

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

#png(filename="figures/heatmap_func_cat_clr_z.png" ,units = 'in',width=9, height=6, res=1000)
pheatmap(mat ,fontsize_row = 7 ,fontsize_col = 8,
         scale = "none",
         clustering_distance_rows = "correlation",   
         clustering_distance_cols = "correlation",   
         clustering_method = "ward.D",
         annotation_legend = TRUE, heatmap_legend_param = list(title = "clr_z") )
#dev.off()






core_cat = getPrevalence(
  tse_pathway, rank = "Kingdom",detection = 1, sort = TRUE, assay.type = "counts",
  as.relative = FALSE) %>% head(30)

mapTaxonomy(tse_bacteria, taxa = "s__flagellatus")


# Identifying no hierarchy functional genes

unknown_genes <- c("NO HIERARCHY")
#"", " ", "hypothetical protein", 
#"hypothetical protein, partial","hypothetical protein",
#"Retron-type reverse transcriptase")

unique_genes = unique(rowData(tse_pathway)$Order)

known_genes <- setdiff(unique_genes, unknown_genes) 

tse_pathway_unknown = tse_pathway




# Subset by feature (core funationa categories)
selected <- rowData(tse_pathway)$Order %in% c("Retron-type reverse transcriptase") &
  !is.na(rowData(tse_pathway)$Order)


gene_renamed <- lapply(rowData(tse_pathway_unknown)$Order, function(x){
  if (x %in% unknown_genes) {x} else {"Other"}
})
rowData(tse_pathway_unknown)$Order <- as.character(gene_renamed)

plots <- plotAbundance(tse_pathway_unknown, rank = "Order", 
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
plot_class_relabundance <- wrap_plots(plots[[2]], plots[[1]], ncol = 1, heights = c(2, 10))
# Combine the plot with the legend

#png(filename="figures/barplot_class.png" ,units = 'in',width=9, height=6, res=1000)
wrap_plots(plot_class_relabundance, legend, nrow = 1, widths = c(2, 1))
#dev.off()



