#Composition barplot for core community based on class level

# Load the global variable 
mae <- readRDS("mae.rds")
tse = mae[[1]]
tse_bacteria = mae[[2]]
tse_active = mae[[3]]
tse_pathway = mae[[4]]
tse_enzyme = mae[[5]]
tse_gene = mae[[6]]

# Getting top taxa on a Phylum level
tse_class_bacteria <- agglomerateByRank(tse_bacteria, rank ="Class")

core_class = getPrevalence(
  tse_class_bacteria, rank = "Class",detection = 1, sort = TRUE, assay.type = "counts",
  as.relative = FALSE) %>% head(12)


core_bacteria_class <-  names(core_class)


core_bacteria_class <- setdiff(core_bacteria_class, "c__") 
# Renaming the "Phylum" rank to keep only top taxa and the rest to "Other"
class_renamed <- lapply(rowData(tse_class_bacteria)$Class, function(x){
  if (x %in% core_bacteria_class) {x} else {"Other"}
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
#wrap_plots(plot_phylum_relabundance, legend_phylum, nrow = 1, widths = c(2, 1))
#dev.off()

#######################################################################

png(filename="figures/combined_barplot.png" ,units = 'in',width=9, height=6, res=1000)
tiff(filename="figures/combined_barplot.tiff" ,units = 'in',width=9, height=6, res=1000)

print(wrap_plots(plot_phylum_relabundance, legend_phylum,plots[[1]], legend, nrow = 2, widths = c(2, 1)))
dev.off()
