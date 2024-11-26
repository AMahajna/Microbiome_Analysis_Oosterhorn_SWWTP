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
#Creating alternative experiment for esach taxonomic level
#ranks = c('Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus','Species')
#setTaxonomyRanks(ranks)
#######################################################################
# Stress Response 
selected_stress <- rowData(tse_pathway)$Kingdom %in% c("Stress Response") &
  !is.na(rowData(tse_pathway)$Kingdom)
tse_stress <- tse_pathway[selected_stress, ]
################################################################################
# Core functional categories 
tse_func_cat <- agglomerateByRank(tse_pathway, rank ="Kingdom")

top_func_cat <- getTopFeatures(tse_func_cat, top = 10, assay.type = "relabundance")


# Renaming the "Phylum" rank to keep only top taxa and the rest to "Other"
func_cat_renamed <- lapply(rowData(tse_func_cat)$Kingdom, function(x){
  if (x %in% top_func_cat) {x} else {"Other"}
})
rowData(tse_func_cat)$Kingdom <- as.character(func_cat_renamed)

#######################################################################
tse_func_cat <- agglomerateByRank(tse_func_cat, rank ="Kingdom")


plots_func_cat <- plotAbundance(tse_func_cat, rank = "Kingdom", 
                       assay.type = "relabundance", features = "Date")

# Modify the legend of the first plot to be smaller
plots_func_cat[[1]] <- plots_func_cat[[1]] +
  theme(
    legend.key.size = unit(0.3, 'cm'),
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 10))

# Modify the legend of the second plot to be smaller
plots_func_cat[[2]] <- plots_func_cat[[2]] +
  theme(
    legend.key.height = unit(0.3, 'cm'),
    legend.key.width = unit(0.3, 'cm'),
    legend.text = element_text(size = 6),
    legend.title = element_text(size = 8),
    legend.direction = "vertical")

# Combine legends
legend_func_cat <- wrap_plots(
  as_ggplot(get_legend(plots_func_cat[[1]])),
  as_ggplot(get_legend(plots_func_cat[[2]])),
  ncol = 1)

# Removes legends from the plots
plots_func_cat[[1]] <- plots_func_cat[[1]] + theme(legend.position = "none")
plots_func_cat[[2]] <- plots_func_cat[[2]] +
  theme(legend.position = "none", axis.title.x=element_blank())

# Combine plots
plot_fun_cat <- wrap_plots(plots_func_cat[[2]] , plots_func_cat[[1]], ncol = 1, heights = c(2, 10))
# Combine the plot with the legend

barplot_func_cat = wrap_plots(plot_fun_cat, legend_func_cat, nrow = 1, widths = c(2, 1))

#png(filename="figures/barplot_func_cat.png" ,units = 'in',width=9, height=6, res=1000)
#print(barplot_func_cat)
#dev.off()

################################################################################
# stress related 

tse_stress_Phylum <- agglomerateByRank(tse_stress, rank ="Phylum")

#Delete top_taxa_bacteria_phylum <- getTopFeatures(tse_phylum_bacteria, top = 10, assay.type = "relabundance")

# Renaming the "Phylum" rank to keep only top taxa and the rest to "Other"
#phylum_renamed <- lapply(rowData(tse_phylum_bacteria)$Phylum, function(x){
#  if (x %in% top_taxa_bacteria_phylum) {x} else {"Other"}
#})
#rowData(tse_phylum_bacteria)$Phylum <- as.character(phylum_renamed)

#tse_phylum_bacteria_sub <- agglomerateByRank(tse_phylum_bacteria, rank= "Phylum")


plots_stress_phylum <- plotAbundance(tse_stress_Phylum, rank = "Phylum", 
                              assay.type = "relabundance", features = "Date")

# Modify the legend of the first plot to be smaller
plots_stress_phylum[[1]] <- plots_stress_phylum[[1]] +
  theme(
    legend.key.size = unit(0.3, 'cm'),
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 10))



# Modify the legend of the second plot to be smaller
plots_stress_phylum[[2]] <- plots_stress_phylum[[2]] +
  theme(
    legend.key.height = unit(0.3, 'cm'),
    legend.key.width = unit(0.3, 'cm'),
    legend.text = element_text(size = 6),
    legend.title = element_text(size = 8),
    legend.direction = "vertical")

# Combine legends
legend_stress_phylum <- wrap_plots(
  as_ggplot(get_legend(plots_stress_phylum[[1]])),
  as_ggplot(get_legend(plots_stress_phylum[[2]])),
  ncol = 1)

# Removes legends from the plots
plots_stress_phylum[[1]] <- plots_stress_phylum[[1]] + theme(legend.position = "none")
plots_stress_phylum[[2]] <- plots_stress_phylum[[2]] +
  theme(legend.position = "none", axis.title.x=element_blank())

# Combine plots
plot_stress_phylum_relabundance <- wrap_plots(plots_stress_phylum[[2]], plots_stress_phylum[[1]], ncol = 1, heights = c(2, 10))
# Combine the plot with the legend

barplot_stress_response = wrap_plots(plot_stress_phylum_relabundance, legend_stress_phylum, nrow = 1, widths = c(2, 1))


#png(filename="figures/barplot_stress_response.png" ,units = 'in',width=9, height=6, res=1000)
#print(barplot_stress_response)
#dev.off()

#######################################################################

#png(filename="figures/combined_barplot_func_cat_stress.png" ,units = 'in',width=9, height=6, res=1000)
#wrap_plots(plot_fun_cat, legend_func_cat,plots_stress_phylum[[1]], legend_stress_phylum, nrow = 2, widths = c(2, 1))
#dev.off()
