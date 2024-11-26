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


# Stress Response 
selected_stress <- rowData(tse_pathway)$Kingdom %in% c("Stress Response") &
  !is.na(rowData(tse_pathway)$Kingdom)
tse_stress <- tse_pathway[selected_stress, ]


#######################################################################

tse_stress =  agglomerateByRank(tse_stress, rank = "Order", update.tree = TRUE)



#relative abundance for the top-20 phylum
top_20_log_10 = plotAbundanceDensity(tse_stress, layout = "jitter", 
                                     assay.type = "relabundance",
                                     n = 20, point_size=2, point_shape=19, colour_by="Season",
                                     point_alpha=0.5)+
  scale_x_log10(label=scales::percent)

#relative abundance for the top-20 phylum
top_20 = plotAbundanceDensity(tse_stress, layout = "jitter", 
                              assay.type = "relabundance",
                              n = 20, point_size=2, point_shape=19, colour_by="Season",
                              point_alpha=0.5)


#relative abundance for the top-5 species
#png(filename="figures/abundance_density_plot_active_5.png" ,units = 'in',width=9, height=6, res=1000)
top_5 = plotAbundanceDensity(tse_stress, layout = "density", 
                             assay.type = "relabundance",
                             n = 5, colour_by="Season", 
                             point_alpha=1/10) 
#dev.off()

combined_plot_density <- plot_grid(
  top_20_log_10, 
  top_5, 
  ncol = 1, 
  labels = c("A", "B"),  # Labels for the plots
  label_size = 12        # Customize label size
)

png(filename="figures/combined_abundance_density_plot_top20log10_top5_stress_genes.png" ,units = 'in',width=9, height=6, res=1000)
print(combined_plot_density)
dev.off()
