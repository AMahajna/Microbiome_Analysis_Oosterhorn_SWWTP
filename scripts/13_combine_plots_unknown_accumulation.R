# Load the global variable 
mae <- readRDS("mae.rds")
tse = mae[[1]]
tse_bacteria = mae[[2]]
tse_active = mae[[3]]
tse_pathway = mae[[4]]
tse_enzyme = mae[[5]]
tse_gene = mae[[6]]

# Getting top taxa on a Phylum level
tse_unknown <- agglomerateByRank(tse_pathway, rank ="Order")
unknown <- c("NO HIERARCHY")

# Renaming the "Phylum" rank to keep only top taxa and the rest to "Other"
unknown_renamed <- lapply(rowData(tse_unknown)$Order, function(x){
  if (x %in% unknown) {x} else {"Other"}
})
rowData(tse_unknown)$Order <- as.character(unknown_renamed)

# Visualizing the composition barplot, with samples order by "Bacteroidetes"
plots <- plotAbundance(tse_unknown, rank = "Order", 
                       assay.type = "relabundance", features = "Date")

# Modify the legend of the first plot to be smaller
plots[[1]] <- plots[[1]] +  # Change legend title
  theme(
    legend.key.size = unit(0.3, 'cm'),
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 0, , color = "white")
  )

# Modify the legend of the second plot to be smaller
plots[[2]] <- plots[[2]] +
  theme(
    legend.key.height = unit(0.3, 'cm'),
    legend.key.width = unit(0.3, 'cm'),
    legend.text = element_text(size = 6),
    legend.title = element_text(size = 0),
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
unknown_barplot = wrap_plots(plot_class_relabundance, legend, nrow = 1, widths = c(2, 1))
#dev.off()

################################################################################
# Getting top taxa on a Phylum level
tse <- agglomerateByRank(tse, rank ="Phylum")
top_taxa <- getTop(tse, top = 10, assay.type = "relabundance")

# Renaming the "Phylum" rank to keep only top taxa and the rest to "Other"
phylum_renamed <- lapply(rowData(tse)$Phylum, function(x){
  if (x %in% top_taxa) {x} else {"Other"}
})
rowData(tse)$Phylum_sub <- as.character(phylum_renamed)
# Agglomerate the data based on specified taxa
tse_sub <- agglomerateByVariable(tse, by = "rows", f = "Phylum_sub")

# Visualizing the composition barplot, with samples order by "Bacteroidetes"
plotAbundance(
  tse_sub, assay.type = "relabundance",
  order.row.by = "abund", order.col.by = "Bacteroidetes")

################################################################################
# Change the names of the last four columns
new_names <- c("Enzyme Accumulation", "SEED Subsystem Accumulation", "Description Accumulation", "SEED Subsystem Functional\nCategory Accumulation")
colnames(species_function_curve)[(ncol(species_function_curve)-3):ncol(species_function_curve)] <- new_names

# Reshape the dataframe to long format
df_long <- melt(
  species_function_curve,
  id.vars = "Species_Accumulation", 
  variable.name = "Y_Variable", 
  value.name = "Y_Value"
)

# Create the scatter plot with LOESS regression lines
curves = ggplot(df_long, aes(x = Species_Accumulation, y = Y_Value, color = Y_Variable)) +
  geom_point() +                           # Scatter plot
  geom_smooth(method = "loess", se = FALSE) + # LOESS regression lines
  labs(
    x = "Species Accumulation",
    y = "Accumulation Values",
    color = "Legend"
  ) +
  theme_minimal()

################################################################################

# Combine plots with cowplot
combined_plot <- plot_grid(
   curves, unknown_barplot,
  labels = c("A", "B"), # Labels for the plots
  nrow = 2,             # Correct argument: nrow (not nrows)
  label_size = 14       # Size of the labels
)

png("figures/combined_plot_unknowns_accumulation.png", width = 8, height = 6, units = "in", res = 1000)
print(combined_plot)
dev.off()

