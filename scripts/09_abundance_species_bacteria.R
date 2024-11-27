################################################################################
##Abundance
# Load the global variable 
mae <- readRDS("mae.rds")
tse = mae[[1]]
tse_bacteria = mae[[2]]
tse_active = mae[[3]]
tse_pathway = mae[[4]]
tse_enzyme = mae[[5]]
tse_gene = mae[[6]]

#relative abundance for the top-20 phylum
top_20_log_10 = plotAbundanceDensity(tse_bacteria, layout = "jitter", 
                                     assay.type = "relabundance",
                                     n = 20, point_size=2, point_shape=19, colour_by="Season",
                                     point_alpha=0.5)+
  scale_x_log10(label=scales::percent)

#relative abundance for the top-20 phylum
top_20 = plotAbundanceDensity(tse_bacteria, layout = "jitter", 
                              assay.type = "relabundance",
                              n = 20, point_size=2, point_shape=19, colour_by="Season",
                              point_alpha=0.5)

#png(filename="figures/abundance_density_plot_active_top20.png" ,units = 'in',width=9, height=6, res=1000)
#print(top_20)
#dev.off()
#png(filename="figures/abundance_density_plot_active_top20_log10.png" ,units = 'in',width=9, height=6, res=1000)
#print(top_20_log_10)
#dev.off()

#relative abundance for the top-5 species
#png(filename="figures/abundance_density_plot_active_5.png" ,units = 'in',width=9, height=6, res=1000)
top_5 = plotAbundanceDensity(tse_bacteria, layout = "density", 
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

png(filename="figures/combined_abundance_species_bacteria_top20log10_top5.png" ,units = 'in',width=9, height=6, res=1000)
print(combined_plot_density)
dev.off()
################################################################################
