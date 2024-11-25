################################################################################
##Abundance: SAMSA2 Data for active bacterial community .

# Load the global variable 
mae <- readRDS("mae.rds")
tse = mae[[1]]
tse_bacteria = mae[[2]]
tse_active = mae[[3]]
tse_pathway = mae[[4]]
tse_enzyme = mae[[5]]
tse_gene = mae[[6]]


#relative abundance for the top-20 phylum
top_20_log_10 = plotAbundanceDensity(tse_active, layout = "jitter", 
                     assay.type = "relabundance",
                     n = 20, point_size=2, point_shape=19, colour_by="Season",
                     point_alpha=0.5)+
  scale_x_log10(label=scales::percent)

#relative abundance for the top-20 phylum
top_20 = plotAbundanceDensity(tse_active, layout = "jitter", 
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
top_5 = plotAbundanceDensity(tse_active, layout = "density", 
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

#png(filename="figures/combined_abundance_density_plot_top20log10_top5.png" ,units = 'in',width=9, height=6, res=1000)
#print(combined_plot_density)
#dev.off()
################################################################################
##Core community 
core_active_species = getPrevalence(
  tse_active, detection = 1, sort = TRUE, assay.type = "counts",
  as.relative = FALSE) %>% head(19)

cat("There are 19 core active species in the microbiome")
print(core_active_species)

df_core_active_species <- data.frame(
  Species = names(core_active_species),
  Value = as.numeric(core_active_species)
)

#write_xlsx(df_core_active_species, "output_data/core_active_species.xlsx")


################################################################################
################################################################################
################################################################################

##Core community from tse (generated using Kraken2 includes all Kingdoms)

################################################################################

core_phylum = getPrevalence(
  tse, rank = "Phylum",detection = 1, sort = TRUE, assay.type = "counts",
  as.relative = FALSE) %>% head(8)

core_species =getPrevalence(
  tse, rank = "Species",detection = 1, sort = TRUE, assay.type = "counts",
  as.relative = FALSE) %>% head(14)

core_phylum_bacteria = getPrevalence(
  tse_bacteria, rank = "Phylum",detection = 1, sort = TRUE, assay.type = "counts",
  as.relative = FALSE) %>% head(5)

core_species_bacteria = getPrevalence(
  tse_bacteria, rank = "Species",detection = 1, sort = TRUE, assay.type = "counts",
  as.relative = FALSE) %>% head(13)

core_class_bacteria = getPrevalence(
  tse_bacteria, rank = "Class",detection = 1, sort = TRUE, assay.type = "counts",
  as.relative = FALSE) %>% head(12)

mapTaxonomy(tse_bacteria, taxa = "s__flagellatus")
#Homonyms in Taxonomy 

################################################################################
##Prevalence  
#Create tse_phylum 
tse_phylum <- agglomerateByRank(tse, rank = "Phylum", update.tree = TRUE)
altExp(tse, "Phylum") <- tse_phylum

#Prevalence of Phylum in total community labeled by kingdom 
rowData(altExp(tse,"Phylum"))$prevalence <- 
  getPrevalence(altExp(tse,"Phylum"), detection = 1/100, 
                sort = FALSE,
                assay.type = "counts", as.relative = TRUE)

prevalence_plot = plotRowData(altExp(tse,"Phylum"), "prevalence", point_size=5,
                              colour_by = "Kingdom")
#png(filename="figures/prevalence.png" ,units = 'in',width=9, height=6, res=1000)
print(prevalence_plot)
#dev.off()

#Core community is made of 8 Phyla (5 bacteria and 3 Eukaryota)
sum(getPrevalence(altExp(tse,"Phylum"), detection = 1/100, sort = FALSE,
                  assay.type = "counts", as.relative = TRUE, prevalence = 1) == 1)

################################################################################
#Prevalence taxonomic tree

altExp(tse,"Phylum") <- agglomerateByRank(tse, "Phylum")

rowData(altExp(tse,"Phylum"))$prevalence <- 
  getPrevalence(altExp(tse,"Phylum"), detection = 1/100, 
                sort = FALSE,
                assay.type = "counts", as.relative = TRUE)


top_phyla <- getTopFeatures(
  altExp(tse,"Phylum"),
  method="sum",
  top=5L,
  assay.type="counts")
top_phyla_mean <- getTopFeatures(
  altExp(tse,"Phylum"),
  method="mean",
  top=5L,
  assay.type="counts")

#Number of ranks can also be decreased
x <- unsplitByRanks(tse, ranks = taxonomyRanks(tse)[1:7])
x <- addHierarchyTree(x)

prevalence_tree = plotRowTree(
  x[rowData(x)$Phylum %in% top_phyla,],
  edge_colour_by = "Phylum",
  tip_colour_by = "prevalence",
  node_colour_by = "prevalence")

#png(filename="figures/prevalence_tree.png" ,units = 'in',width=9, height=6, res=1000)
#print(prevalence_tree)
#dev.off()

################################################################################
