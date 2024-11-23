#Composition barplot for core community based on class level

# Load the global variable 
mae <- readRDS("mae.rds")
tse = mae[[1]]
tse_bacteria = mae[[2]]
tse_active = mae[[3]]
tse_pathway = mae[[4]]
tse_enzyme = mae[[5]]

# Getting top taxa on a Phylum level
tse_class_bacteria <- agglomerateByRank(tse_bacteria, rank ="Class")

top_taxa_bacteria_class <- getTopFeatures(tse_class_bacteria, top = 12, assay.type = "relabundance")

# Renaming the "Phylum" rank to keep only top taxa and the rest to "Other"
class_renamed <- lapply(rowData(tse_class_bacteria)$Class, function(x){
  if (x %in% top_taxa_bacteria_class) {x} else {"Other"}
})
rowData(tse_class_bacteria)$Class <- as.character(class_renamed)

plots <- plotAbundance(tse_class_bacteria, rank = "Class", 
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

################################################################################
# Getting top taxa on a Phylum level
tse_phylum_bacteria <- agglomerateByRank(tse_bacteria, rank ="Phylum")

top_taxa_bacteria_phylum <- getTopFeatures(tse_phylum_bacteria, top = 10, assay.type = "relabundance")

# Renaming the "Phylum" rank to keep only top taxa and the rest to "Other"
phylum_renamed <- lapply(rowData(tse_phylum_bacteria)$Phylum, function(x){
  if (x %in% top_taxa_bacteria_phylum) {x} else {"Other"}
})
rowData(tse_phylum_bacteria)$Phylum <- as.character(phylum_renamed)

#tse_phylum_bacteria_sub <- agglomerateByRank(tse_phylum_bacteria, rank= "Phylum")


plots_phylum <- plotAbundance(tse_phylum_bacteria, rank = "Phylum", 
                       assay.type = "relabundance", features = "Date")

# Modify the legend of the first plot to be smaller
plots_phylum[[1]] <- plots_phylum[[1]] +
  theme(
    legend.key.size = unit(0.3, 'cm'),
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 10))

# Modify the legend of the second plot to be smaller
plots_phylum[[2]] <- plots_phylum[[2]] +
  theme(
    legend.key.height = unit(0.3, 'cm'),
    legend.key.width = unit(0.3, 'cm'),
    legend.text = element_text(size = 6),
    legend.title = element_text(size = 8),
    legend.direction = "vertical")

# Combine legends
legend_phylum <- wrap_plots(
  as_ggplot(get_legend(plots_phylum[[1]])),
  as_ggplot(get_legend(plots_phylum[[2]])),
  ncol = 1)

# Removes legends from the plots
plots_phylum[[1]] <- plots_phylum[[1]] + theme(legend.position = "none")
plots_phylum[[2]] <- plots_phylum[[2]] +
  theme(legend.position = "none", axis.title.x=element_blank())

# Combine plots
plot_phylum_relabundance <- wrap_plots(plots_phylum[[2]], plots_phylum[[1]], ncol = 1, heights = c(2, 10))
# Combine the plot with the legend

#png(filename="figures/barplot_phylum.png" ,units = 'in',width=9, height=6, res=1000)
wrap_plots(plot_phylum_relabundance, legend_phylum, nrow = 1, widths = c(2, 1))
#dev.off()

#######################################################################

#png(filename="figures/combined_barplot.png" ,units = 'in',width=9, height=6, res=1000)
wrap_plots(plot_phylum_relabundance, legend_phylum,plots[[1]], legend, nrow = 2, widths = c(2, 1))
#dev.off()

#######################################################################
# Getting top taxa on a species level - active from samsa2 
tse_active_species <- agglomerateByRank(tse_active, rank ="Species")

top_taxa_active <- getTopFeatures(tse_active_species, top = 19, assay.type = "relabundance")

# Renaming the "Phylum" rank to keep only top taxa and the rest to "Other"
species_renamed <- lapply(rowData(tse_active_species)$Species, function(x){
  if (x %in% top_taxa_active) {x} else {"Other"}
})
rowData(tse_active_species)$Species <- as.character(species_renamed)


plots_species <- plotAbundance(tse_active_species, rank = "Species", 
                              assay.type = "relabundance", features = "Date")

# Modify the legend of the first plot to be smaller
plots_species[[1]] <- plots_species[[1]] +
  theme(
    legend.key.size = unit(0.3, 'cm'),
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 10))

# Modify the legend of the second plot to be smaller
plots_species[[2]] <- plots_species[[2]] +
  theme(
    legend.key.height = unit(0.3, 'cm'),
    legend.key.width = unit(0.3, 'cm'),
    legend.text = element_text(size = 6),
    legend.title = element_text(size = 8),
    legend.direction = "vertical")

# Combine legends
legend_species <- wrap_plots(
  as_ggplot(get_legend(plots_species[[1]])),
  as_ggplot(get_legend(plots_species[[2]])),
  ncol = 1)

# Removes legends from the plots
plots_species[[1]] <- plots_species[[1]] + theme(legend.position = "none")
plots_species[[2]] <- plots_species[[2]] +
  theme(legend.position = "none", axis.title.x=element_blank())

# Combine plots
plot_species_relabundance <- wrap_plots(plots_species[[2]], plots_species[[1]], ncol = 1, heights = c(2, 10))
# Combine the plot with the legend

#png(filename="figures/barplot_species.png" ,units = 'in',width=9, height=6, res=1000)
#wrap_plots(plot_species_relabundance, legend_species, nrow = 1, widths = c(2, 1))
#dev.off()

#######################################################################
#Function 
# Agglomerate to phylum level
tse_pathway_cat <- agglomerateByPrevalence(tse_pathway, rank = "Kingdom")

# Add clr-transformation on samples
tse_pathway_cat <- transformAssay(
  tse_pathway_cat, assay.type = "counts", method = "relabundance", pseudocount = 1)
tse_pathway_cat <- transformAssay(tse_pathway_cat, assay.type = "relabundance", method = "clr")

# Add scale features (taxa)
tse_pathway_cat <- transformAssay(
  tse_pathway_cat, assay.type = "clr", MARGIN = "features", method = "standardize",
  name = "clr_z")

# Gets the assay table
mat <- assay(tse_pathway_cat, "clr_z")

# Create a heatmap with customized legend position
#png(filename="figures/heatmap_func_cat_2.png" ,units = 'in',width=9, height=6, res=1000)
#tiff(filename="figures/heatmap_func_cat_2.tiff" ,units = 'in',width=9, height=6, res=1000)
#Heatmap(mat, cluster_columns = FALSE,    cluster_rows = FALSE,
#        name = "value",
#        col = colorRampPalette(brewer.pal(11, "RdYlBu"))(256), # Diverging palette
#        heatmap_legend_param = list(title = "clr_z", 
#                                    legend_position = c(0.85, 0.95), 
#                                    direction = "vertical"), 
#        row_names_gp = gpar(fontsize = 10),
#        column_names_gp = gpar(fontsize = 10)
#        )
#dev.off()

#tiff(filename="figures/heatmap_func_cat_1.tiff" ,units = 'in',width=9, height=6, res=1000)
#png(filename="figures/heatmap_func_cat_1.png" ,units = 'in',width=9, height=6, res=1000)
#pheatmap(mat, cluster_rows = FALSE,fontsize_row = 8 ,fontsize_col = 9,
#         cluster_cols = FALSE,annotation_legend = TRUE, heatmap_legend_param = list(title = "clr_z"))
#dev.off()

#######################################################################
################################################################################
# Getting core functional category 
core_cat = as.data.frame(getPrevalence(
  tse_pathway, rank = "Kingdom",detection = 1, sort = TRUE, assay.type = "counts",
  as.relative = FALSE) %>% head(16))
# There are 16 core SEED Subsystem core functional category 


# Subset by feature (core funationa categories)
selected <- rowData(tse_pathway)$Kingdom %in% c("Motility and Chemotaxis", "Stress Response",
                                                "Membrane Transport","DNA Metabolism","Cell Wall and Capsule",
                                                "RNA Metabolism","Nucleosides and Nucleotides","Protein Metabolism",
                                                "Virulence, Disease and Defense","Clustering-based subsystems","Cofactors, Vitamins, Prosthetic Groups, Pigments",
                                                "Amino Acids and Derivatives","Fatty Acids, Lipids, and Isoprenoids",
                                                "Respiration","Carbohydrates") &
  !is.na(rowData(tse_pathway)$Kingdom)

tse_core_cat <- tse_pathway[selected, ]


#Now, to create bar plot 
tse_core_cat <- agglomerateByRank(tse_core_cat, rank ="Kingdom")


###
top_taxa_bacteria_phylum <- getTopFeatures(tse_phylum_bacteria, top = 10, assay.type = "relabundance")

# Renaming the "Phylum" rank to keep only top taxa and the rest to "Other"
phylum_renamed <- lapply(rowData(tse_phylum_bacteria)$Phylum, function(x){
  if (x %in% top_taxa_bacteria_phylum) {x} else {"Other"}
})
rowData(tse_phylum_bacteria)$Phylum <- as.character(phylum_renamed)
####
#tse_phylum_bacteria_sub <- agglomerateByRank(tse_phylum_bacteria, rank= "Phylum")


plots_cat <- plotAbundance(tse_core_cat, rank = "Kingdom", 
                              assay.type = "relabundance", features = "Date")

# Modify the legend of the first plot to be smaller
plots_cat[[1]] <- plots_cat[[1]] +
  theme(
    legend.key.size = unit(0.3, 'cm'),
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 10))

# Modify the legend of the second plot to be smaller
plots_cat[[2]] <- plots_cat[[2]] +
  theme(
    legend.key.height = unit(0.3, 'cm'),
    legend.key.width = unit(0.3, 'cm'),
    legend.text = element_text(size = 6),
    legend.title = element_text(size = 8),
    legend.direction = "vertical")

# Combine legends
legend_cat <- wrap_plots(
  as_ggplot(get_legend(plots_cat[[1]])),
  as_ggplot(get_legend(plots_cat[[2]])),
  ncol = 1)

# Removes legends from the plots
plots_cat[[1]] <- plots_cat[[1]] + theme(legend.position = "none")
plots_cat[[2]] <- plots_cat[[2]] +
  theme(legend.position = "none", axis.title.x=element_blank())

# Combine plots
plot_cat_relabundance <- wrap_plots(plots_cat[[2]], plots_cat[[1]], ncol = 1, heights = c(2, 10))
# Combine the plot with the legend

#png(filename="figures/barplot_phylum.png" ,units = 'in',width=9, height=6, res=1000)
wrap_plots(plot_cat_relabundance, legend_cat, nrow = 1, widths = c(2, 1))
#dev.off()

#######################################################################






################################################################################
# barplot of pathway funtional category known and unknown 



tse_fun_cat <- agglomerateByRank(tse_pathway, rank ="Kingdom")

fun_cat_known <- !(rowData(tse_pathway)$Kingdom %in% c(""))
selected_known = unique(as.data.frame(rowData(tse_pathway)$Kingdom)[fun_cat_known, ])

fun_cat_unknown <- (rowData(tse_pathway)$Kingdom %in% c(""))
selected_unknown = tse_pathway[fun_cat_unknown, ]
  
# Renaming the "known" rank to keep only top taxa and the rest to "Other"
class_renamed <- lapply(rowData(tse_class_bacteria)$Class, function(x){
  if (x %in% top_taxa_bacteria_class) {x} else {"Other"}
})
rowData(tse_class_bacteria)$Class <- as.character(class_renamed)

plots <- plotAbundance(tse_class_bacteria, rank = "Class", 
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

