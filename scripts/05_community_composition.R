#Composition barplot for core community based on class level

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
#wrap_plots(plot_phylum_relabundance, legend_phylum,plots[[1]], legend, nrow = 2, widths = c(2, 1))
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

# Agglomerate to phylum level
tse_pathway_des <- agglomerateByPrevalence(tse_pathway_des, rank = "Kingdom")

# Add clr-transformation on samples
tse_pathway_des <- transformAssay(
  tse_pathway_des, assay.type = "counts", method = "relabundance", pseudocount = 1)
tse_pathway_des <- transformAssay(tse_pathway_des, assay.type = "relabundance", method = "clr")

# Add scale features (taxa)
tse_pathway_des <- transformAssay(
  tse_pathway_des, assay.type = "clr", MARGIN = "features", method = "standardize",
  name = "clr_z")

# Gets the assay table
mat <- assay(tse_pathway_des, "clr_z")

# Creates the heatmap
#png(filename="figures/heatmap_des.png" ,units = 'in',width=9, height=6, res=1000)
#pheatmap(mat)
#dev.off()

